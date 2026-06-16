from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import constants
from asgard_core.hadronic_validation import as_matching, as_strictly_increasing_grid

import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module


AM3_C_CGS = constants.para_c
AM3_SIGMA_T_CGS = constants.para_sigmat
AM3_MASS_ELECTRON_GEV = constants.para_m_e_gev
AM3_MASS_PROTON_GEV = constants.para_m_p_gev
AM3_MASS_PION_CHARGED_GEV = constants.para_m_pi_charged_gev
AM3_MASS_MUON_GEV = constants.para_m_mu_gev
HADRONIC_IC_BACKEND = "fortran_hadronic_ic"


@dataclass(frozen=True)
class HadronicInverseComptonKernel:
    """AM3-compatible delta-kernel indices for hadronic IC emission."""

    dln_energy: float
    ind_min_energy_pho_hadgrid: int
    delta_e_p: np.ndarray
    jmax_p: np.ndarray
    delta_e_pi: np.ndarray
    jmax_pi: np.ndarray
    delta_e_mu: np.ndarray
    jmax_mu: np.ndarray


@dataclass(frozen=True)
class HadronicInverseComptonOutput:
    """Hadronic IC emissivity output (proton/pion/muon channels)."""

    photon_energy_gev: np.ndarray
    epsilon_p_ic: np.ndarray
    epsilon_pi_ic: np.ndarray
    epsilon_mu_ic: np.ndarray
    coeff_p_cgs: float
    coeff_pi_cgs: float
    coeff_mu_cgs: float
    kernel: HadronicInverseComptonKernel


def initialize_hadronic_inverse_compton_kernel(
    hadron_energy_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    ind_min_energy_pho_hadgrid: int,
) -> HadronicInverseComptonKernel:
    """Initialize AM3 HadIC small-step kernels.

    This follows AM3 `HadIC::initialize_kernels` exactly:
    gamma = E_had / m,
    deltaE = int(log(gamma^2)/dlnE),
    jmax1 = int(log(0.5*m/gamma)/dlnE) + ind_min_energy_pho_hadgrid,
    jmax = min(nbins_pho, jmax1 + deltaE).
    """

    e_had = as_strictly_increasing_grid(hadron_energy_gev, "hadron_energy_gev", require_finite=True)
    e_ph = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev", require_finite=True)
    zeros_ph = np.zeros_like(e_ph)
    zeros_had = np.zeros_like(e_had)
    (
        _eps_p,
        _eps_pi,
        _eps_mu,
        _coeff_p,
        _coeff_pi,
        _coeff_mu,
        dln_energy,
        delta_e_p,
        jmax_p,
        delta_e_pi,
        jmax_pi,
        delta_e_mu,
        jmax_mu,
    ) = hadronic_fortran_module.fs_hadronic_hadronic_ic_shell(
        e_had,
        e_ph,
        zeros_ph,
        zeros_had,
        zeros_had,
        zeros_had,
        zeros_had,
        zeros_had,
        zeros_had,
        zeros_had,
        int(ind_min_energy_pho_hadgrid),
    )

    return HadronicInverseComptonKernel(
        dln_energy=float(dln_energy),
        ind_min_energy_pho_hadgrid=int(ind_min_energy_pho_hadgrid),
        delta_e_p=np.asarray(delta_e_p, dtype=np.int64),
        jmax_p=np.asarray(jmax_p, dtype=np.int64),
        delta_e_pi=np.asarray(delta_e_pi, dtype=np.int64),
        jmax_pi=np.asarray(jmax_pi, dtype=np.int64),
        delta_e_mu=np.asarray(delta_e_mu, dtype=np.int64),
        jmax_mu=np.asarray(jmax_mu, dtype=np.int64),
    )


def solve_hadronic_inverse_compton(
    hadron_energy_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    photons_on_had_grid_per_gev: np.ndarray,
    protons_per_gev: np.ndarray,
    pion_plus_per_gev: np.ndarray,
    pion_minus_per_gev: np.ndarray,
    muon_minus_left_per_gev: np.ndarray,
    muon_minus_right_per_gev: np.ndarray,
    muon_plus_left_per_gev: np.ndarray,
    muon_plus_right_per_gev: np.ndarray,
    ind_min_energy_pho_hadgrid: int = 0,
    kernel: HadronicInverseComptonKernel | None = None,
) -> HadronicInverseComptonOutput:
    """Compute hadronic-grid IC emission in AM3 delta approximation.

    Input/Output units follow AM3's HadIC path:
    - energies in GeV
    - distributions as number density per GeV [cm^-3 GeV^-1]
    - output epsilon_* arrays are AM3-style emissivity terms with
      coeff = c * sigma_T * (m_e / m)^2 and one dlnE integration factor.
    """

    e_had = as_strictly_increasing_grid(hadron_energy_gev, "hadron_energy_gev", require_finite=True)
    e_ph = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev", require_finite=True)
    photons = as_matching(photons_on_had_grid_per_gev, e_ph, "photons_on_had_grid_per_gev", require_finite=True)
    protons = as_matching(protons_per_gev, e_had, "protons_per_gev", require_finite=True)
    pion_plus = as_matching(pion_plus_per_gev, e_had, "pion_plus_per_gev", require_finite=True)
    pion_minus = as_matching(pion_minus_per_gev, e_had, "pion_minus_per_gev", require_finite=True)
    muon_minus_left = as_matching(muon_minus_left_per_gev, e_had, "muon_minus_left_per_gev", require_finite=True)
    muon_minus_right = as_matching(muon_minus_right_per_gev, e_had, "muon_minus_right_per_gev", require_finite=True)
    muon_plus_left = as_matching(muon_plus_left_per_gev, e_had, "muon_plus_left_per_gev", require_finite=True)
    muon_plus_right = as_matching(muon_plus_right_per_gev, e_had, "muon_plus_right_per_gev", require_finite=True)

    (
        epsilon_p,
        epsilon_pi,
        epsilon_mu,
        coeff_p,
        coeff_pi,
        coeff_mu,
        dln_energy,
        delta_e_p,
        jmax_p,
        delta_e_pi,
        jmax_pi,
        delta_e_mu,
        jmax_mu,
    ) = hadronic_fortran_module.fs_hadronic_hadronic_ic_shell(
        e_had,
        e_ph,
        photons,
        protons,
        pion_plus,
        pion_minus,
        muon_minus_left,
        muon_minus_right,
        muon_plus_left,
        muon_plus_right,
        int(ind_min_energy_pho_hadgrid),
    )
    if kernel is None:
        kernel = HadronicInverseComptonKernel(
            dln_energy=float(dln_energy),
            ind_min_energy_pho_hadgrid=int(ind_min_energy_pho_hadgrid),
            delta_e_p=np.asarray(delta_e_p, dtype=np.int64),
            jmax_p=np.asarray(jmax_p, dtype=np.int64),
            delta_e_pi=np.asarray(delta_e_pi, dtype=np.int64),
            jmax_pi=np.asarray(jmax_pi, dtype=np.int64),
            delta_e_mu=np.asarray(delta_e_mu, dtype=np.int64),
            jmax_mu=np.asarray(jmax_mu, dtype=np.int64),
        )

    coeff_p_target = _am3_ic_coeff(AM3_MASS_PROTON_GEV)
    coeff_pi_target = _am3_ic_coeff(AM3_MASS_PION_CHARGED_GEV)
    coeff_mu_target = _am3_ic_coeff(AM3_MASS_MUON_GEV)
    epsilon_p = np.asarray(epsilon_p, dtype=float) * (coeff_p_target / float(coeff_p))
    epsilon_pi = np.asarray(epsilon_pi, dtype=float) * (coeff_pi_target / float(coeff_pi))
    epsilon_mu = np.asarray(epsilon_mu, dtype=float) * (coeff_mu_target / float(coeff_mu))

    return HadronicInverseComptonOutput(
        photon_energy_gev=e_ph,
        epsilon_p_ic=epsilon_p,
        epsilon_pi_ic=epsilon_pi,
        epsilon_mu_ic=epsilon_mu,
        coeff_p_cgs=float(coeff_p_target),
        coeff_pi_cgs=float(coeff_pi_target),
        coeff_mu_cgs=float(coeff_mu_target),
        kernel=kernel,
    )


def _am3_ic_coeff(mass_gev: float) -> float:
    mass_ratio = mass_gev / AM3_MASS_ELECTRON_GEV
    return AM3_C_CGS * AM3_SIGMA_T_CGS / mass_ratio / mass_ratio
