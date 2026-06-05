from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import constants

try:
    import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module
except ImportError:
    hadronic_fortran_module = None


AM3_C_CGS = constants.para_c
AM3_SIGMA_T_CGS = constants.para_sigmat
AM3_MASS_ELECTRON_GEV = constants.para_m_e_gev
AM3_MASS_PROTON_GEV = constants.para_m_p_gev
AM3_MASS_PION_CHARGED_GEV = constants.para_m_pi_charged_gev
AM3_MASS_MUON_GEV = constants.para_m_mu_gev
_HAS_FORTRAN_HADRONIC_IC = hadronic_fortran_module is not None and hasattr(
    hadronic_fortran_module, "fs_hadronic_hadronic_ic_shell"
)
HADRONIC_IC_BACKEND = "fortran_hadronic_ic" if _HAS_FORTRAN_HADRONIC_IC else "unavailable"


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

    if not _HAS_FORTRAN_HADRONIC_IC:
        raise RuntimeError("Hadronic IC core must be provided by the Fortran backend.")
    e_had = _as_strictly_increasing(hadron_energy_gev, "hadron_energy_gev")
    e_ph = _as_strictly_increasing(photon_energy_gev, "photon_energy_gev")
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

    e_had = _as_strictly_increasing(hadron_energy_gev, "hadron_energy_gev")
    e_ph = _as_strictly_increasing(photon_energy_gev, "photon_energy_gev")
    photons = _as_matching_1d(photons_on_had_grid_per_gev, e_ph, "photons_on_had_grid_per_gev")
    protons = _as_matching_1d(protons_per_gev, e_had, "protons_per_gev")
    pion_plus = _as_matching_1d(pion_plus_per_gev, e_had, "pion_plus_per_gev")
    pion_minus = _as_matching_1d(pion_minus_per_gev, e_had, "pion_minus_per_gev")
    muon_minus_left = _as_matching_1d(muon_minus_left_per_gev, e_had, "muon_minus_left_per_gev")
    muon_minus_right = _as_matching_1d(muon_minus_right_per_gev, e_had, "muon_minus_right_per_gev")
    muon_plus_left = _as_matching_1d(muon_plus_left_per_gev, e_had, "muon_plus_left_per_gev")
    muon_plus_right = _as_matching_1d(muon_plus_right_per_gev, e_had, "muon_plus_right_per_gev")

    if not _HAS_FORTRAN_HADRONIC_IC:
        raise RuntimeError("Hadronic IC core must be provided by the Fortran backend.")

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


def _build_species_kernel(
    hadron_energy_gev: np.ndarray,
    dln_energy: float,
    nbins_had: int,
    nbins_pho: int,
    ind_min_energy_pho_hadgrid: int,
    mass_gev: float,
) -> tuple[np.ndarray, np.ndarray]:
    delta_e = np.zeros(nbins_had, dtype=np.int64)
    jmax = np.ones(nbins_had, dtype=np.int64) * int(nbins_pho)
    for i in range(nbins_had):
        gamma = float(hadron_energy_gev[i] / mass_gev)
        d_e = int(np.log(gamma * gamma) / dln_energy)
        delta_e[i] = d_e
        jmax1 = int(np.log(0.5 * mass_gev / gamma) / dln_energy) + int(ind_min_energy_pho_hadgrid)
        candidate = jmax1 + d_e
        jmax[i] = int(nbins_pho) if candidate > nbins_pho else int(candidate)
    return delta_e, jmax


def _compute_emissivity_channel(
    photons_on_had_grid_per_gev: np.ndarray,
    hadron_density_per_gev: np.ndarray,
    delta_e: np.ndarray,
    jmax: np.ndarray,
    dln_energy: float,
    coeff_cgs: float,
) -> np.ndarray:
    nbins_pho = int(photons_on_had_grid_per_gev.size)
    nbins_had = int(hadron_density_per_gev.size)
    out = np.zeros(nbins_pho, dtype=float)
    for j in range(nbins_pho):
        z = 0.0
        for i in range(nbins_had):
            d_e = int(delta_e[i])
            if j < d_e or j > int(jmax[i]):
                continue
            src_idx = j - d_e
            if src_idx < 0 or src_idx >= nbins_pho:
                raise ValueError("Kernel maps to an out-of-grid photon source index.")
            z += float(photons_on_had_grid_per_gev[src_idx]) * float(hadron_density_per_gev[i])
        out[j] = z * dln_energy * coeff_cgs
    return out


def _am3_ic_coeff(mass_gev: float) -> float:
    mass_ratio = mass_gev / AM3_MASS_ELECTRON_GEV
    return AM3_C_CGS * AM3_SIGMA_T_CGS / mass_ratio / mass_ratio


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
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    return arr


def _as_matching_1d(values: np.ndarray, grid: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != grid.shape:
        raise ValueError(f"{name} must match grid shape.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    return arr


def _log_spacing(grid: np.ndarray, name: str) -> float:
    dln = np.diff(np.log(grid))
    if not np.allclose(dln, dln[0], rtol=1.0e-6, atol=1.0e-12):
        raise ValueError(f"{name} must be logarithmically uniform.")
    return float(dln[0])
