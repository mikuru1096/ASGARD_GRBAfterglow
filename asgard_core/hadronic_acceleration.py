from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping

import numpy as np

from src import constants

try:
    import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module
except ImportError:
    hadronic_fortran_module = None


ELEMENTARY_CHARGE_ESU = constants.para_e

ERG_PER_GEV = constants.para_gev2erg
ELECTRON_MASS_GEV = constants.para_m_e_gev
PROTON_MASS_GEV = constants.para_m_p_gev
NEUTRON_MASS_GEV = constants.para_m_n_gev
PION_CHARGED_MASS_GEV = constants.para_m_pi_charged_gev
MUON_MASS_GEV = constants.para_m_mu_gev

GEV_TO_GRAM = constants.para_gev2erg / (constants.para_c * constants.para_c)
_HAS_FORTRAN_ACCELERATION = hadronic_fortran_module is not None and hasattr(
    hadronic_fortran_module, "fs_hadronic_acceleration_shell"
)
ACCELERATION_BACKEND = "fortran_acceleration" if _HAS_FORTRAN_ACCELERATION else "unavailable"


@dataclass(frozen=True)
class SpeciesProperties:
    name: str
    mass_gev: float
    charge_number: int

    @property
    def mass_g(self) -> float:
        return self.mass_gev * GEV_TO_GRAM

    @property
    def abs_charge_esu(self) -> float:
        return abs(float(self.charge_number)) * ELEMENTARY_CHARGE_ESU


SPECIES_DB: dict[str, SpeciesProperties] = {
    "proton": SpeciesProperties("proton", PROTON_MASS_GEV, +1),
    "neutron": SpeciesProperties("neutron", NEUTRON_MASS_GEV, 0),
    "pion_plus": SpeciesProperties("pion_plus", PION_CHARGED_MASS_GEV, +1),
    "pion_minus": SpeciesProperties("pion_minus", PION_CHARGED_MASS_GEV, -1),
    "muon_plus": SpeciesProperties("muon_plus", MUON_MASS_GEV, +1),
    "muon_minus": SpeciesProperties("muon_minus", MUON_MASS_GEV, -1),
}


@dataclass(frozen=True)
class InjectionConfig:
    species: str
    luminosity_erg_s: float
    spectral_index: float
    gamma_min: float
    gamma_max: float
    gamma_cut: float | None = None


@dataclass(frozen=True)
class MaxEnergyEstimate:
    species: str
    gamma_max: float
    gamma_dyn: float
    gamma_syn: float
    gamma_ext: float | None


ExternalPhotonCoolingRate = Callable[[str, np.ndarray, Mapping[str, float]], np.ndarray]


def species_properties(species: str) -> SpeciesProperties:
    prop = SPECIES_DB.get(species)
    if prop is None:
        raise ValueError(f"Unsupported species: {species}")
    return prop


def proton_acceleration_timescale_s(
    gamma: np.ndarray,
    b_field_g: float,
    eta_acc: float,
) -> np.ndarray:
    return acceleration_timescale_s("proton", gamma, b_field_g, eta_acc)


def acceleration_timescale_s(
    species: str,
    gamma: np.ndarray,
    b_field_g: float,
    eta_acc: float,
) -> np.ndarray:
    """
    t_acc(gamma) = eta_acc * r_L / c
                 = eta_acc * gamma * m * c / (|q| * B)
    """
    gamma_arr = _as_positive_1d(gamma, "gamma")
    if not _HAS_FORTRAN_ACCELERATION:
        raise RuntimeError("Acceleration core must be provided by the Fortran backend.")
    t_acc, _t_syn, _q_inj, _gmax, _gdyn, _gsyn, _gext, _has_gext = hadronic_fortran_module.fs_hadronic_acceleration_shell(
        species,
        gamma_arr,
        float(b_field_g),
        float(eta_acc),
        1.0,
        2.0,
        float(gamma_arr[0]),
        float(gamma_arr[-1]),
        1.0,
        False,
        1.0,
        1.0,
        np.array([1.0, 2.0], dtype=float),
        np.array([1.0, 2.0], dtype=float),
        False,
    )
    return np.asarray(t_acc, dtype=float)


def external_photon_cooling_timescale_s(
    species: str,
    gamma: np.ndarray,
    cooling_rate: ExternalPhotonCoolingRate | None = None,
    context: Mapping[str, float] | None = None,
) -> np.ndarray:
    """
    External extension point for future photon-field cooling.
    cooling_rate(species, gamma, context) returns A_ext(gamma) = -dgamma/dt (>0),
    then t_ext(gamma) = gamma / A_ext(gamma).
    """
    _ = species_properties(species)
    gamma_arr = _as_positive_1d(gamma, "gamma")
    if cooling_rate is None:
        return np.full_like(gamma_arr, np.inf)
    ctx = {} if context is None else context
    rate = np.asarray(cooling_rate(species, gamma_arr, ctx), dtype=float)
    if rate.shape != gamma_arr.shape:
        raise ValueError("cooling_rate output must have the same shape as gamma.")
    if np.any(rate <= 0.0):
        raise ValueError("cooling_rate must be strictly positive when provided.")
    return gamma_arr / rate


def species_injection_operator(
    gamma: np.ndarray,
    config: InjectionConfig,
) -> np.ndarray:
    """
    Q(gamma) = Q0 * gamma^{-p} * exp(-gamma/gamma_cut)   (optional cutoff)
    L_inj = ∫_{gamma_min}^{gamma_max} Q(gamma) * gamma * m * c^2 dgamma
    Q0 = L_inj / (m*c^2 * ∫ shape(gamma)*gamma dgamma)
    """
    gamma_arr = _as_positive_1d(gamma, "gamma")
    if not _HAS_FORTRAN_ACCELERATION:
        raise RuntimeError("Acceleration core must be provided by the Fortran backend.")
    _t_acc, _t_syn, q_inj, _gmax, _gdyn, _gsyn, _gext, _has_gext = hadronic_fortran_module.fs_hadronic_acceleration_shell(
        config.species,
        gamma_arr,
        1.0,
        1.0,
        float(config.luminosity_erg_s),
        float(config.spectral_index),
        float(config.gamma_min),
        float(config.gamma_max),
        1.0 if config.gamma_cut is None else float(config.gamma_cut),
        config.gamma_cut is not None,
        1.0,
        1.0,
        np.array([1.0, 2.0], dtype=float),
        np.array([1.0, 2.0], dtype=float),
        False,
    )
    return np.asarray(q_inj, dtype=float)


def build_species_injection_operators(
    gamma: np.ndarray,
    configs: list[InjectionConfig],
) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for cfg in configs:
        out[cfg.species] = species_injection_operator(gamma, cfg)
    return out


def estimate_max_gamma(
    species: str,
    b_field_g: float,
    radius_cm: float,
    gamma_bulk: float,
    eta_acc: float,
    cooling_rate: ExternalPhotonCoolingRate | None = None,
    context: Mapping[str, float] | None = None,
) -> MaxEnergyEstimate:
    """
    gamma_dyn from t_acc = t_dyn, where t_dyn = R / (Gamma * c):
      gamma_dyn = |q| * B * t_dyn / (eta_acc * m * c)

    gamma_syn from t_acc = t_syn:
      gamma_syn = sqrt( 6*pi*|q|*m^2 / (eta_acc*sigma_T*m_e^2*B) )

    If external cooling is provided, gamma_ext from t_acc = t_ext
    by solving f(gamma) = t_acc(gamma) - t_ext(gamma) on a log grid.
    """
    if not _HAS_FORTRAN_ACCELERATION:
        raise RuntimeError("Acceleration core must be provided by the Fortran backend.")
    if cooling_rate is None:
        gamma_scan = np.array([1.0, 2.0], dtype=float)
        ext_rate = np.array([1.0, 2.0], dtype=float)
        has_external = False
    else:
        gamma_scan = np.logspace(0.0, 14.0, 2048)
        ext_rate = np.asarray(cooling_rate(species, gamma_scan, {} if context is None else context), dtype=float)
        has_external = True
    _t_acc, _t_syn, _q_inj, gamma_max, gamma_dyn, gamma_syn, gamma_ext, has_gamma_ext = hadronic_fortran_module.fs_hadronic_acceleration_shell(
        species,
        np.array([1.0, 2.0], dtype=float),
        float(b_field_g),
        float(eta_acc),
        1.0,
        2.0,
        1.0,
        2.0,
        1.0,
        False,
        float(radius_cm),
        float(gamma_bulk),
        gamma_scan,
        ext_rate,
        has_external,
    )
    return MaxEnergyEstimate(
        species=species,
        gamma_max=float(gamma_max),
        gamma_dyn=float(gamma_dyn),
        gamma_syn=float(gamma_syn),
        gamma_ext=float(gamma_ext) if bool(has_gamma_ext) else None,
    )


def _as_positive_1d(values: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be 1d.")
    if arr.size < 2:
        raise ValueError(f"{name} must contain at least two points.")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be > 0.")
    return arr
