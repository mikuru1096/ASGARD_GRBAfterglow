from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import constants

try:
    import src.Hadronic.FS_hadronic_1d as hadronic_fortran_module
except ImportError:
    hadronic_fortran_module = None


_HAS_FORTRAN_PG_OPERATOR = hadronic_fortran_module is not None and hasattr(
    hadronic_fortran_module, "fs_hadronic_pgamma_operator_shell"
)
_HAS_FORTRAN_DECAY_OPERATOR = hadronic_fortran_module is not None and hasattr(
    hadronic_fortran_module, "fs_hadronic_decay_operator_shell"
)

HUMMER2010_OPERATOR_BACKEND = "fortran_operator" if _HAS_FORTRAN_PG_OPERATOR else "unavailable"
HUMMER2010_DECAY_BACKEND = "fortran_decay" if _HAS_FORTRAN_DECAY_OPERATOR else "unavailable"


PROTON_MASS_GEV = constants.para_m_p_gev
PI_PLUS_MASS_GEV = constants.para_m_pi_charged_gev
MUON_MASS_GEV = constants.para_m_mu_gev
CHARGED_PION_DECAY_S = 2.6033e-8
MUON_DECAY_S = 2.1969811e-6
GEV_TO_ERG = constants.para_gev2erg


@dataclass(frozen=True)
class Hummer2010Output:
    gamma_energy_gev: np.ndarray
    gamma_rate_per_gev: np.ndarray
    neutrino_energy_gev: np.ndarray
    neutrino_rate_per_gev: np.ndarray
    process_energy_gev: np.ndarray
    process_rate_per_gev: np.ndarray
    pion_plus_source_rate_per_gev: np.ndarray
    pion_minus_source_rate_per_gev: np.ndarray
    muon_plus_right_source_rate_per_gev: np.ndarray
    muon_plus_left_source_rate_per_gev: np.ndarray
    muon_minus_left_source_rate_per_gev: np.ndarray
    muon_minus_right_source_rate_per_gev: np.ndarray
    prompt_pion_neutrino_rate_per_gev: np.ndarray
    muon_neutrino_rate_per_gev: np.ndarray
    muon_electron_rate_per_gev: np.ndarray
    hadron_energy_gev: np.ndarray
    proton_reinjection_rate_per_gev: np.ndarray
    neutron_reinjection_rate_per_gev: np.ndarray
    proton_loss_rate: np.ndarray
    neutron_loss_rate: np.ndarray
    photon_loss_energy_gev: np.ndarray
    photon_loss_rate: np.ndarray


def solve_hummer2010_pgamma(
    proton_energy_gev: np.ndarray,
    proton_density_per_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    photon_density_per_gev: np.ndarray,
    gamma_energy_gev: np.ndarray,
    neutrino_energy_gev: np.ndarray,
    process_energy_gev: np.ndarray,
    neutron_density_per_gev: np.ndarray | None = None,
) -> Hummer2010Output:
    if not _HAS_FORTRAN_PG_OPERATOR:
        raise RuntimeError("Hummer 2010 photopion core must be provided by the Fortran backend.")
    if not _HAS_FORTRAN_DECAY_OPERATOR:
        raise RuntimeError("Hummer 2010 decay core must be provided by the Fortran backend.")

    e_p = _as_strictly_increasing_1d(proton_energy_gev, "proton_energy_gev")
    n_p = _as_matching_nonnegative_1d(proton_density_per_gev, e_p, "proton_density_per_gev")
    n_n = np.zeros_like(n_p) if neutron_density_per_gev is None else _as_matching_nonnegative_1d(
        neutron_density_per_gev, e_p, "neutron_density_per_gev"
    )
    e_gam_t = _as_strictly_increasing_1d(photon_energy_gev, "photon_energy_gev")
    n_gam_t = _as_matching_nonnegative_1d(photon_density_per_gev, e_gam_t, "photon_density_per_gev")
    e_gamma = _as_strictly_increasing_1d(gamma_energy_gev, "gamma_energy_gev")
    e_nu = _as_strictly_increasing_1d(neutrino_energy_gev, "neutrino_energy_gev")
    e_proc = _as_strictly_increasing_1d(process_energy_gev, "process_energy_gev")

    (
        qpi0_rate,
        qpip_rate,
        qpim_rate,
        qpro_reinj_rate,
        qntr_reinj_rate,
        afpg_pro,
        afpg_ntr,
        afpg_rad,
    ) = hadronic_fortran_module.fs_hadronic_pgamma_operator_shell(
        e_p,
        n_p,
        e_gam_t,
        n_gam_t,
        n_n,
    )
    qpi0_rate = np.asarray(qpi0_rate, dtype=float)
    qpip_rate = np.asarray(qpip_rate, dtype=float)
    qpim_rate = np.asarray(qpim_rate, dtype=float)

    (
        gamma_rate,
        process_rate,
        mu_plus_r_rate,
        mu_plus_l_rate,
        mu_minus_l_rate,
        mu_minus_r_rate,
        prompt_pion_neutrino_rate,
        muon_neutrino_rate,
        muon_electron_rate,
        neutrino_rate,
    ) = hadronic_fortran_module.fs_hadronic_decay_operator_shell(
        e_p,
        qpi0_rate,
        qpip_rate,
        qpim_rate,
        e_gamma,
        e_nu,
        e_proc,
    )

    return Hummer2010Output(
        gamma_energy_gev=e_gamma,
        gamma_rate_per_gev=np.asarray(gamma_rate, dtype=float),
        neutrino_energy_gev=e_nu,
        neutrino_rate_per_gev=np.asarray(neutrino_rate, dtype=float),
        process_energy_gev=e_proc,
        process_rate_per_gev=np.asarray(process_rate, dtype=float),
        pion_plus_source_rate_per_gev=qpip_rate,
        pion_minus_source_rate_per_gev=qpim_rate,
        muon_plus_right_source_rate_per_gev=np.asarray(mu_plus_r_rate, dtype=float),
        muon_plus_left_source_rate_per_gev=np.asarray(mu_plus_l_rate, dtype=float),
        muon_minus_left_source_rate_per_gev=np.asarray(mu_minus_l_rate, dtype=float),
        muon_minus_right_source_rate_per_gev=np.asarray(mu_minus_r_rate, dtype=float),
        prompt_pion_neutrino_rate_per_gev=np.asarray(prompt_pion_neutrino_rate, dtype=float),
        muon_neutrino_rate_per_gev=np.asarray(muon_neutrino_rate, dtype=float),
        muon_electron_rate_per_gev=np.asarray(muon_electron_rate, dtype=float),
        hadron_energy_gev=e_p,
        proton_reinjection_rate_per_gev=np.asarray(qpro_reinj_rate, dtype=float),
        neutron_reinjection_rate_per_gev=np.asarray(qntr_reinj_rate, dtype=float),
        proton_loss_rate=np.asarray(afpg_pro, dtype=float),
        neutron_loss_rate=np.asarray(afpg_ntr, dtype=float),
        photon_loss_energy_gev=e_gam_t,
        photon_loss_rate=np.asarray(afpg_rad, dtype=float),
    )


def _as_strictly_increasing_1d(values: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1d array.")
    if arr.size < 2:
        raise ValueError(f"{name} must contain at least two points.")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be positive.")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing.")
    return arr


def _as_matching_nonnegative_1d(values: np.ndarray, reference: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1 or arr.shape != reference.shape:
        raise ValueError(f"{name} must be a 1d array with the same shape as its energy grid.")
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must be non-negative.")
    return arr
