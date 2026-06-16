from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import constants
from asgard_core.hadronic_validation import as_matching, as_strictly_increasing_grid, log_spacing

import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module

AM3_C_CGS = constants.para_c
AM3_SIGMA_T_CGS = constants.para_sigmat
AM3_MASS_ELECTRON_GEV = constants.para_m_e_gev
AM3_MASS_PION_CHARGED_GEV = constants.para_m_pi_charged_gev
AM3_MASS_MUON_GEV = constants.para_m_mu_gev
AM3_B_CRIT_GAUSS = 4.41e13
AM3_ERG_TO_GEV = 624.0
AM3_PI = constants.pi
SECONDARY_RADIATION_BACKEND = "fortran_secondary_radiation"


@dataclass(frozen=True)
class SecondarySpeciesDistribution:
    pion_plus_per_gev: np.ndarray
    pion_minus_per_gev: np.ndarray
    muon_minus_left_per_gev: np.ndarray
    muon_minus_right_per_gev: np.ndarray
    muon_plus_left_per_gev: np.ndarray
    muon_plus_right_per_gev: np.ndarray


@dataclass(frozen=True)
class SecondaryTargetPhotonField:
    photon_energy_gev: np.ndarray
    photons_on_had_grid_per_gev: np.ndarray
    ind_min_energy_pho_hadgrid: int = 0


@dataclass(frozen=True)
class SecondarySynchrotronKernel:
    dln_energy: float
    kernel_pion: np.ndarray
    kernel_muon: np.ndarray
    magnetic_field_g: float


@dataclass(frozen=True)
class SecondaryInverseComptonKernel:
    dln_energy: float
    ind_min_energy_pho_hadgrid: int
    delta_e_pi: np.ndarray
    jmax_pi: np.ndarray
    delta_e_mu: np.ndarray
    jmax_mu: np.ndarray


@dataclass(frozen=True)
class SecondaryRadiationSpectrum:
    photon_energy_gev: np.ndarray
    pion_synch_rate_per_gev: np.ndarray
    muon_synch_rate_per_gev: np.ndarray
    pion_ic_rate_per_gev: np.ndarray
    muon_ic_rate_per_gev: np.ndarray
    synch_kernel: SecondarySynchrotronKernel
    ic_kernel: SecondaryInverseComptonKernel


def initialize_secondary_synchrotron_kernel(
    hadron_energy_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    magnetic_field_g: float,
) -> SecondarySynchrotronKernel:
    e_had = as_strictly_increasing_grid(hadron_energy_gev, "hadron_energy_gev", require_finite=True)
    e_ph = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev", require_finite=True)
    dln_had = log_spacing(e_had, "hadron_energy_gev")
    b = float(magnetic_field_g)
    if b <= 0.0:
        raise ValueError("magnetic_field_g must be > 0.")

    kernel_pi = np.zeros((e_ph.size, e_had.size), dtype=float)
    kernel_mu = np.zeros((e_ph.size, e_had.size), dtype=float)
    for i in range(e_ph.size):
        for j in range(e_had.size):
            kernel_pi[i, j] = _am3_syn_kernel_ultrarel(e_ph[i], e_had[j], AM3_MASS_PION_CHARGED_GEV, b)
            kernel_mu[i, j] = _am3_syn_kernel_ultrarel(e_ph[i], e_had[j], AM3_MASS_MUON_GEV, b)

    return SecondarySynchrotronKernel(
        dln_energy=float(dln_had),
        kernel_pion=kernel_pi,
        kernel_muon=kernel_mu,
        magnetic_field_g=b,
    )


def initialize_secondary_inverse_compton_kernel(
    hadron_energy_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    ind_min_energy_pho_hadgrid: int = 0,
) -> SecondaryInverseComptonKernel:
    e_had = as_strictly_increasing_grid(hadron_energy_gev, "hadron_energy_gev", require_finite=True)
    e_ph = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev", require_finite=True)
    dln_had = log_spacing(e_had, "hadron_energy_gev")
    dln_ph = log_spacing(e_ph, "photon_energy_gev")
    if not np.isclose(dln_had, dln_ph, rtol=1.0e-10, atol=0.0):
        raise ValueError("hadron_energy_gev and photon_energy_gev must share one logarithmic spacing dlnE.")

    delta_e_pi, jmax_pi = _build_ic_species_kernel(
        e_had,
        dln_had,
        e_ph.size,
        int(ind_min_energy_pho_hadgrid),
        AM3_MASS_PION_CHARGED_GEV,
    )
    delta_e_mu, jmax_mu = _build_ic_species_kernel(
        e_had,
        dln_had,
        e_ph.size,
        int(ind_min_energy_pho_hadgrid),
        AM3_MASS_MUON_GEV,
    )
    return SecondaryInverseComptonKernel(
        dln_energy=float(dln_had),
        ind_min_energy_pho_hadgrid=int(ind_min_energy_pho_hadgrid),
        delta_e_pi=delta_e_pi,
        jmax_pi=jmax_pi,
        delta_e_mu=delta_e_mu,
        jmax_mu=jmax_mu,
    )


def pion_synchrotron_rate(
    hadron_energy_gev: np.ndarray,
    species: SecondarySpeciesDistribution,
    photon_energy_gev: np.ndarray,
    magnetic_field_g: float,
    kernel: SecondarySynchrotronKernel | None = None,
) -> np.ndarray:
    e_had = as_strictly_increasing_grid(hadron_energy_gev, "hadron_energy_gev", require_finite=True)
    e_ph = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev", require_finite=True)
    pion_plus = as_matching(species.pion_plus_per_gev, e_had, "pion_plus_per_gev", require_finite=True)
    pion_minus = as_matching(species.pion_minus_per_gev, e_had, "pion_minus_per_gev", require_finite=True)
    if kernel is None:
        kernel = initialize_secondary_synchrotron_kernel(e_had, e_ph, magnetic_field_g)
    return kernel.dln_energy * (kernel.kernel_pion @ (pion_plus + pion_minus))


def muon_synchrotron_rate(
    hadron_energy_gev: np.ndarray,
    species: SecondarySpeciesDistribution,
    photon_energy_gev: np.ndarray,
    magnetic_field_g: float,
    kernel: SecondarySynchrotronKernel | None = None,
) -> np.ndarray:
    e_had = as_strictly_increasing_grid(hadron_energy_gev, "hadron_energy_gev", require_finite=True)
    e_ph = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev", require_finite=True)
    mu_ml = as_matching(species.muon_minus_left_per_gev, e_had, "muon_minus_left_per_gev", require_finite=True)
    mu_mr = as_matching(species.muon_minus_right_per_gev, e_had, "muon_minus_right_per_gev", require_finite=True)
    mu_pl = as_matching(species.muon_plus_left_per_gev, e_had, "muon_plus_left_per_gev", require_finite=True)
    mu_pr = as_matching(species.muon_plus_right_per_gev, e_had, "muon_plus_right_per_gev", require_finite=True)
    if kernel is None:
        kernel = initialize_secondary_synchrotron_kernel(e_had, e_ph, magnetic_field_g)
    return kernel.dln_energy * (kernel.kernel_muon @ (mu_ml + mu_mr + mu_pl + mu_pr))


def pion_inverse_compton_rate(
    hadron_energy_gev: np.ndarray,
    species: SecondarySpeciesDistribution,
    target: SecondaryTargetPhotonField,
    kernel: SecondaryInverseComptonKernel | None = None,
) -> np.ndarray:
    e_had = as_strictly_increasing_grid(hadron_energy_gev, "hadron_energy_gev", require_finite=True)
    e_ph = as_strictly_increasing_grid(target.photon_energy_gev, "photon_energy_gev", require_finite=True)
    photons = as_matching(target.photons_on_had_grid_per_gev, e_ph, "photons_on_had_grid_per_gev", require_finite=True)
    pion_plus = as_matching(species.pion_plus_per_gev, e_had, "pion_plus_per_gev", require_finite=True)
    pion_minus = as_matching(species.pion_minus_per_gev, e_had, "pion_minus_per_gev", require_finite=True)
    if kernel is None:
        kernel = initialize_secondary_inverse_compton_kernel(
            e_had,
            e_ph,
            int(target.ind_min_energy_pho_hadgrid),
        )
    coeff_pi = _am3_ic_coeff(AM3_MASS_PION_CHARGED_GEV)
    return _compute_ic_channel(photons, pion_plus + pion_minus, kernel.delta_e_pi, kernel.jmax_pi, kernel.dln_energy, coeff_pi)


def muon_inverse_compton_rate(
    hadron_energy_gev: np.ndarray,
    species: SecondarySpeciesDistribution,
    target: SecondaryTargetPhotonField,
    kernel: SecondaryInverseComptonKernel | None = None,
) -> np.ndarray:
    e_had = as_strictly_increasing_grid(hadron_energy_gev, "hadron_energy_gev", require_finite=True)
    e_ph = as_strictly_increasing_grid(target.photon_energy_gev, "photon_energy_gev", require_finite=True)
    photons = as_matching(target.photons_on_had_grid_per_gev, e_ph, "photons_on_had_grid_per_gev", require_finite=True)
    mu_ml = as_matching(species.muon_minus_left_per_gev, e_had, "muon_minus_left_per_gev", require_finite=True)
    mu_mr = as_matching(species.muon_minus_right_per_gev, e_had, "muon_minus_right_per_gev", require_finite=True)
    mu_pl = as_matching(species.muon_plus_left_per_gev, e_had, "muon_plus_left_per_gev", require_finite=True)
    mu_pr = as_matching(species.muon_plus_right_per_gev, e_had, "muon_plus_right_per_gev", require_finite=True)
    if kernel is None:
        kernel = initialize_secondary_inverse_compton_kernel(
            e_had,
            e_ph,
            int(target.ind_min_energy_pho_hadgrid),
        )
    coeff_mu = _am3_ic_coeff(AM3_MASS_MUON_GEV)
    return _compute_ic_channel(photons, mu_ml + mu_mr + mu_pl + mu_pr, kernel.delta_e_mu, kernel.jmax_mu, kernel.dln_energy, coeff_mu)


def solve_secondary_radiation_spectrum(
    hadron_energy_gev: np.ndarray,
    species: SecondarySpeciesDistribution,
    target: SecondaryTargetPhotonField,
    magnetic_field_g: float,
    synch_kernel: SecondarySynchrotronKernel | None = None,
    ic_kernel: SecondaryInverseComptonKernel | None = None,
) -> SecondaryRadiationSpectrum:
    e_had = as_strictly_increasing_grid(hadron_energy_gev, "hadron_energy_gev", require_finite=True)
    e_ph = as_strictly_increasing_grid(target.photon_energy_gev, "photon_energy_gev", require_finite=True)
    (
        pion_syn,
        muon_syn,
        pion_ic,
        muon_ic,
    ) = hadronic_fortran_module.fs_hadronic_secondary_radiation_shell(
        e_had,
        e_ph,
        np.asarray(species.pion_plus_per_gev, dtype=float),
        np.asarray(species.pion_minus_per_gev, dtype=float),
        np.asarray(species.muon_minus_left_per_gev, dtype=float),
        np.asarray(species.muon_minus_right_per_gev, dtype=float),
        np.asarray(species.muon_plus_left_per_gev, dtype=float),
        np.asarray(species.muon_plus_right_per_gev, dtype=float),
        np.asarray(target.photons_on_had_grid_per_gev, dtype=float),
        int(target.ind_min_energy_pho_hadgrid),
        float(magnetic_field_g),
    )
    local_synch_kernel = synch_kernel if synch_kernel is not None else initialize_secondary_synchrotron_kernel(e_had, e_ph, magnetic_field_g)
    local_ic_kernel = ic_kernel if ic_kernel is not None else initialize_secondary_inverse_compton_kernel(e_had, e_ph, int(target.ind_min_energy_pho_hadgrid))
    return SecondaryRadiationSpectrum(
        photon_energy_gev=e_ph,
        pion_synch_rate_per_gev=np.asarray(pion_syn, dtype=float),
        muon_synch_rate_per_gev=np.asarray(muon_syn, dtype=float),
        pion_ic_rate_per_gev=np.asarray(pion_ic, dtype=float),
        muon_ic_rate_per_gev=np.asarray(muon_ic, dtype=float),
        synch_kernel=local_synch_kernel,
        ic_kernel=local_ic_kernel,
    )


def _am3_syn_kernel_ultrarel(photon_energy_gev: float, particle_energy_gev: float, particle_mass_gev: float, magnetic_field_g: float) -> float:
    norm = (
        (3.0 * 1.732 / AM3_PI)
        * AM3_SIGMA_T_CGS
        * AM3_C_CGS
        * magnetic_field_g
        * AM3_B_CRIT_GAUSS
        / (8.0 * AM3_PI)
        * AM3_ERG_TO_GEV
        / particle_mass_gev
    )
    b_dimless = magnetic_field_g / AM3_B_CRIT_GAUSS
    mass_ratio = AM3_MASS_ELECTRON_GEV / particle_mass_gev
    xbar = photon_energy_gev * particle_mass_gev / (
        3.0 * particle_energy_gev * particle_energy_gev * b_dimless * mass_ratio * mass_ratio
    )
    two_xbar = 2.0 * xbar
    if two_xbar < 1.0e-2:
        return norm * 1.80842 * np.power(xbar, 0.33333) * np.power(2.0, -0.6666667)
    if two_xbar < 1.0:
        y = np.log10(two_xbar)
        poly = (
            -0.35775237
            - 0.83695385 * y
            - 1.1449608 * y * y
            - 0.68137283 * y * y * y
            - 0.22754737 * y * y * y * y
            - 0.031967334 * y * y * y * y * y
        )
        return norm * np.power(10.0, poly) / 2.0
    if two_xbar < 10.0:
        y = np.log10(two_xbar)
        poly = (
            -0.35842494
            - 0.79652041 * y
            - 1.6113032 * y * y
            + 0.26055213 * y * y * y
            - 1.6979017 * y * y * y * y
            + 0.032955035 * y * y * y * y * y
        )
        return norm * np.power(10.0, poly) / 2.0
    if two_xbar < 100.0:
        return norm * AM3_PI / 4.0 * np.exp(-two_xbar) * (1.0 - 99.0 / 162.0 / two_xbar)
    return 0.0


def _build_ic_species_kernel(
    hadron_energy_gev: np.ndarray,
    dln_energy: float,
    nbins_pho: int,
    ind_min_energy_pho_hadgrid: int,
    mass_gev: float,
) -> tuple[np.ndarray, np.ndarray]:
    delta_e = np.zeros(hadron_energy_gev.size, dtype=np.int64)
    jmax = np.ones(hadron_energy_gev.size, dtype=np.int64) * int(nbins_pho)
    for i in range(hadron_energy_gev.size):
        gamma = float(hadron_energy_gev[i] / mass_gev)
        d_e = int(np.log(gamma * gamma) / dln_energy)
        delta_e[i] = d_e
        jmax1 = int(np.log(0.5 * mass_gev / gamma) / dln_energy) + int(ind_min_energy_pho_hadgrid)
        candidate = jmax1 + d_e
        jmax[i] = int(nbins_pho) if candidate > nbins_pho else int(candidate)
    return delta_e, jmax


def _compute_ic_channel(
    photons_on_had_grid_per_gev: np.ndarray,
    hadron_density_per_gev: np.ndarray,
    delta_e: np.ndarray,
    jmax: np.ndarray,
    dln_energy: float,
    coeff_cgs: float,
) -> np.ndarray:
    out = np.zeros(photons_on_had_grid_per_gev.size, dtype=float)
    for j in range(photons_on_had_grid_per_gev.size):
        z = 0.0
        for i in range(hadron_density_per_gev.size):
            d_e = int(delta_e[i])
            if j < d_e or j > int(jmax[i]):
                continue
            src_idx = j - d_e
            if src_idx < 0 or src_idx >= photons_on_had_grid_per_gev.size:
                raise ValueError("IC kernel maps to an out-of-grid photon source index.")
            z += float(photons_on_had_grid_per_gev[src_idx]) * float(hadron_density_per_gev[i])
        out[j] = z * dln_energy * coeff_cgs
    return out


def _am3_ic_coeff(mass_gev: float) -> float:
    mass_ratio = mass_gev / AM3_MASS_ELECTRON_GEV
    return AM3_C_CGS * AM3_SIGMA_T_CGS / mass_ratio / mass_ratio
