from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import constants

try:
    import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module
except ImportError:
    hadronic_fortran_module = None


PROTON_MASS_GEV = constants.para_m_p_gev
PI0_MASS_GEV = constants.para_m_pi0_gev
MBARN_TO_CM2 = 1.0e-27
_HAS_FORTRAN_PP_DELTA = hadronic_fortran_module is not None and hasattr(
    hadronic_fortran_module, "fs_hadronic_pp_delta_shell"
)
PP_DELTA_BACKEND = "fortran_pp_delta" if _HAS_FORTRAN_PP_DELTA else "unavailable"


@dataclass(frozen=True)
class HadronicPPOutput:
    gamma_energy_gev: np.ndarray
    gamma_rate_per_gev: np.ndarray
    neutrino_energy_gev: np.ndarray
    neutrino_rate_per_gev: np.ndarray
    pair_energy_gev: np.ndarray
    pair_rate_per_gev: np.ndarray
    proton_energy_gev: np.ndarray
    proton_loss_rate: np.ndarray


def solve_pp_delta(
    proton_energy_gev: np.ndarray,
    proton_density_per_gev: np.ndarray,
    target_proton_density_cm3: float,
    gamma_energy_gev: np.ndarray,
    neutrino_energy_gev: np.ndarray,
    pair_energy_gev: np.ndarray,
    *,
    kappa_inelastic: float = 0.5,
    pion_energy_fraction: float = 0.17,
    neutral_pion_fraction: float = 1.0 / 3.0,
) -> HadronicPPOutput:
    """Minimal pp microphysics backend in delta-functional approximation.

    Implemented relations:
    1) Inelastic cross section (Kelner et al. 2006, Eq. 79 form):
       sigma_inel(Tp) = (34.3 + 1.88 L + 0.25 L^2) * [1 - (Tp_th/Tp)^4]^2 mb, L = ln(Tp/Tp_th).
    2) Collision frequency:
       t_pp^{-1}(Ep) = n_H * c * sigma_inel(Ep).
    3) Proton loss term (continuous approximation):
       dN_p/dt|loss = -kappa_inelastic * t_pp^{-1} * N_p.
    4) Secondary source terms (delta mapping E_s = x_s * E_p):
       Q_s(E_s) = (m_s / x_s) * t_pp^{-1}(E_p) * N_p(E_p), E_p = E_s / x_s.

    Validity:
    - Physical threshold enforced by sigma_inel(Tp < Tp_th) = 0.
    - Delta approximation is intended as a minimal backend for broad hadronic sweeps;
      accuracy degrades near threshold and for detailed spectral features.
    """

    e_p = _as_strictly_increasing(proton_energy_gev, "proton_energy_gev")
    n_p = _as_matching(proton_density_per_gev, e_p, "proton_density_per_gev")
    e_gamma = _as_strictly_increasing(gamma_energy_gev, "gamma_energy_gev")
    e_nu = _as_strictly_increasing(neutrino_energy_gev, "neutrino_energy_gev")
    e_pair = _as_strictly_increasing(pair_energy_gev, "pair_energy_gev")
    if target_proton_density_cm3 < 0.0:
        raise ValueError("target_proton_density_cm3 must be non-negative.")
    if not (0.0 < kappa_inelastic <= 1.0):
        raise ValueError("kappa_inelastic must be in (0, 1].")
    if not (0.0 < pion_energy_fraction < 1.0):
        raise ValueError("pion_energy_fraction must be in (0, 1).")
    if not (0.0 <= neutral_pion_fraction <= 1.0):
        raise ValueError("neutral_pion_fraction must be in [0, 1].")

    if not _HAS_FORTRAN_PP_DELTA:
        raise RuntimeError("pp delta core must be provided by the Fortran backend.")

    gamma_rate, neutrino_rate, pair_rate, proton_loss_rate = hadronic_fortran_module.fs_hadronic_pp_delta_shell(
        e_p,
        n_p,
        float(target_proton_density_cm3),
        e_gamma,
        e_nu,
        e_pair,
        float(kappa_inelastic),
        float(pion_energy_fraction),
        float(neutral_pion_fraction),
    )
    gamma_rate = np.asarray(gamma_rate, dtype=float)
    neutrino_rate = np.asarray(neutrino_rate, dtype=float)
    pair_rate = np.asarray(pair_rate, dtype=float)
    proton_loss_rate = np.asarray(proton_loss_rate, dtype=float)

    return HadronicPPOutput(
        gamma_energy_gev=e_gamma,
        gamma_rate_per_gev=gamma_rate,
        neutrino_energy_gev=e_nu,
        neutrino_rate_per_gev=neutrino_rate,
        pair_energy_gev=e_pair,
        pair_rate_per_gev=pair_rate,
        proton_energy_gev=e_p,
        proton_loss_rate=proton_loss_rate,
    )


def _pp_threshold_kinetic_energy_gev() -> float:
    return 2.0 * PI0_MASS_GEV + PI0_MASS_GEV * PI0_MASS_GEV / (2.0 * PROTON_MASS_GEV)


def _as_strictly_increasing(values: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1d array.")
    if arr.size < 2:
        raise ValueError(f"{name} must contain at least two points.")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive.")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing.")
    return arr


def _as_matching(values: np.ndarray, grid: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != grid.shape:
        raise ValueError(f"{name} must match grid shape.")
    return arr
