from __future__ import annotations

from dataclasses import dataclass
from math import pi

import numpy as np

from prompt.eats import EATSGeometry, EATSNumerics, project_branch_flux
from prompt.internal_shock import BranchHistory, InternalShockSolution
from src import Radiation, constants
from src.Electron.electron_radiation import electron_radiation_kernel as electron_radiation_module


@dataclass(frozen=True)
class InternalShockMicrophysics:
    epsilon_e: float = 0.1
    epsilon_b: float = 0.01
    electron_index_p: float = 2.3
    accelerated_electron_fraction: float = 1.0
    acceleration_efficiency: float = 1.0


@dataclass(frozen=True)
class RadiationNumerics:
    electron_gamma_min: float = 1.0
    electron_gamma_max: float = 1.0e8
    num_electron_gamma: int = 161
    photon_frequency_min_hz: float = 1.0e8
    photon_frequency_max_hz: float = 1.0e28
    num_photon_frequency: int = 181
    index_syn_integr: int = 2
    num_threads: int = 1


@dataclass(frozen=True)
class BranchRadiation:
    gamma_e: np.ndarray
    seed_frequency_hz: np.ndarray
    d_n_dgamma: np.ndarray
    sync_luminosity: np.ndarray
    sync_seed: np.ndarray
    ssc_luminosity: np.ndarray
    ssc_seed: np.ndarray
    gamma_gamma_absorption: np.ndarray
    gamma_m: np.ndarray
    gamma_c: np.ndarray
    gamma_max: np.ndarray
    compactness: np.ndarray
    efficiency: np.ndarray


@dataclass(frozen=True)
class PromptObservedFlux:
    observer_frequency_hz: np.ndarray
    observer_time_s: np.ndarray
    fs_sync: np.ndarray
    fs_ssc: np.ndarray
    rs_sync: np.ndarray
    rs_ssc: np.ndarray
    total: np.ndarray
    fs_radiation: BranchRadiation
    rs_radiation: BranchRadiation


def compute_prompt_observed_flux(
    solution: InternalShockSolution,
    observer_frequency_hz: np.ndarray,
    observer_time_s: np.ndarray,
    *,
    microphysics: InternalShockMicrophysics | None = None,
    radiation_numerics: RadiationNumerics | None = None,
    eats_numerics: EATSNumerics | None = None,
) -> PromptObservedFlux:
    microphysics = InternalShockMicrophysics() if microphysics is None else microphysics
    radiation_numerics = RadiationNumerics() if radiation_numerics is None else radiation_numerics
    eats_numerics = EATSNumerics(num_threads=radiation_numerics.num_threads) if eats_numerics is None else eats_numerics
    _validate_microphysics(microphysics)
    _validate_radiation_numerics(radiation_numerics)
    observer_frequency_hz = np.asarray(observer_frequency_hz, dtype=float)
    observer_time_s = np.asarray(observer_time_s, dtype=float)
    if np.any(observer_frequency_hz <= 0.0) or np.any(observer_time_s < 0.0):
        raise ValueError("observer frequencies must be positive and observer times non-negative.")

    fs_radiation = compute_branch_radiation(solution.fs, microphysics, radiation_numerics)
    rs_radiation = compute_branch_radiation(solution.rs, microphysics, radiation_numerics)
    geometry = EATSGeometry(solution.redshift, solution.opening_angle_rad, solution.viewing_angle_rad)
    distance_prefactor = (1.0 + solution.redshift) / (4.0 * pi * solution.luminosity_distance_cm**2)
    fs_sync_source = fs_radiation.sync_luminosity * fs_radiation.gamma_gamma_absorption * distance_prefactor
    fs_ssc_source = fs_radiation.ssc_luminosity * fs_radiation.gamma_gamma_absorption * distance_prefactor
    rs_sync_source = rs_radiation.sync_luminosity * rs_radiation.gamma_gamma_absorption * distance_prefactor
    rs_ssc_source = rs_radiation.ssc_luminosity * rs_radiation.gamma_gamma_absorption * distance_prefactor

    fs_sync = _project(solution.fs, fs_sync_source, fs_radiation.seed_frequency_hz, observer_frequency_hz, observer_time_s, geometry, eats_numerics)
    fs_ssc = _project(solution.fs, fs_ssc_source, fs_radiation.seed_frequency_hz, observer_frequency_hz, observer_time_s, geometry, eats_numerics)
    rs_sync = _project(solution.rs, rs_sync_source, rs_radiation.seed_frequency_hz, observer_frequency_hz, observer_time_s, geometry, eats_numerics)
    rs_ssc = _project(solution.rs, rs_ssc_source, rs_radiation.seed_frequency_hz, observer_frequency_hz, observer_time_s, geometry, eats_numerics)
    return PromptObservedFlux(
        observer_frequency_hz=observer_frequency_hz,
        observer_time_s=observer_time_s,
        fs_sync=fs_sync,
        fs_ssc=fs_ssc,
        rs_sync=rs_sync,
        rs_ssc=rs_ssc,
        total=fs_sync + fs_ssc + rs_sync + rs_ssc,
        fs_radiation=fs_radiation,
        rs_radiation=rs_radiation,
    )


def compute_branch_radiation(
    branch: BranchHistory,
    microphysics: InternalShockMicrophysics,
    numerics: RadiationNumerics,
) -> BranchRadiation:
    gamma_e = np.logspace(np.log10(numerics.electron_gamma_min), np.log10(numerics.electron_gamma_max), numerics.num_electron_gamma)
    seed_frequency = np.logspace(
        np.log10(numerics.photon_frequency_min_hz),
        np.log10(numerics.photon_frequency_max_hz),
        numerics.num_photon_frequency,
    )
    num_nu = seed_frequency.size
    num_r = branch.radius_cm.size
    d_n_dgamma = np.zeros((gamma_e.size, num_r), dtype=float)
    gamma_m = np.zeros(num_r, dtype=float)
    gamma_c = np.zeros(num_r, dtype=float)
    gamma_max = np.zeros(num_r, dtype=float)
    compactness = np.zeros(num_r, dtype=float)
    efficiency = np.zeros(num_r, dtype=float)
    sync_luminosity = np.zeros((num_nu, num_r), dtype=float)
    sync_seed = np.zeros((num_nu, num_r), dtype=float)
    if branch.valid_shock:
        for i in range(num_r):
            if branch.electron_luminosity_comoving_erg_s[i] <= 0.0 or branch.total_b_g[i] <= 0.0:
                continue
            electron_distribution = _electron_distribution(branch, i, gamma_e, microphysics)
            d_n_dgamma[:, i] = electron_distribution
            gamma_m[i] = _gamma_m(branch, i, microphysics)
            gamma_c[i] = _gamma_c(branch, i)
            gamma_max[i] = _gamma_max(branch, i, microphysics)
            compactness[i] = _compactness(branch, i)
            efficiency[i] = _efficiency(branch, i, microphysics)
            p_syn, seed_syn = electron_radiation_module.get_syn_selected(
                numerics.index_syn_integr,
                float(branch.radius_cm[i]),
                float(branch.total_b_g[i]),
                numerics.num_threads,
                gamma_e,
                electron_distribution,
                seed_frequency,
            )
            sync_luminosity[:, i] = np.asarray(p_syn, dtype=float)
            sync_seed[:, i] = np.asarray(seed_syn, dtype=float)
    ssc_luminosity, ssc_seed = Radiation.ssc_spec(
        branch.radius_cm,
        gamma_e,
        d_n_dgamma,
        seed_frequency,
        sync_seed,
        numerics.num_threads,
    )
    tau_extra = np.zeros((num_nu, num_r), dtype=float)
    absorption = Radiation.annihilation(
        branch.gamma,
        branch.radius_cm,
        seed_frequency,
        sync_seed,
        ssc_seed,
        tau_extra,
        numerics.num_threads,
    )
    if not branch.valid_shock:
        absorption = np.ones_like(absorption)
    return BranchRadiation(
        gamma_e=gamma_e,
        seed_frequency_hz=seed_frequency,
        d_n_dgamma=d_n_dgamma,
        sync_luminosity=np.asarray(sync_luminosity, dtype=float),
        sync_seed=np.asarray(sync_seed, dtype=float),
        ssc_luminosity=np.asarray(ssc_luminosity, dtype=float),
        ssc_seed=np.asarray(ssc_seed, dtype=float),
        gamma_gamma_absorption=np.asarray(absorption, dtype=float),
        gamma_m=gamma_m,
        gamma_c=gamma_c,
        gamma_max=gamma_max,
        compactness=compactness,
        efficiency=efficiency,
    )


def _project(
    branch: BranchHistory,
    source_flux: np.ndarray,
    seed_frequency_hz: np.ndarray,
    observer_frequency_hz: np.ndarray,
    observer_time_s: np.ndarray,
    geometry: EATSGeometry,
    numerics: EATSNumerics,
) -> np.ndarray:
    return project_branch_flux(
        characteristic_time_s=branch.characteristic_time_s,
        gamma=branch.gamma,
        radius_cm=branch.radius_cm,
        source_flux=source_flux,
        seed_frequency_hz=seed_frequency_hz,
        observer_frequency_hz=observer_frequency_hz,
        observer_time_s=observer_time_s,
        geometry=geometry,
        numerics=numerics,
    )


def _electron_distribution(
    branch: BranchHistory,
    index: int,
    gamma_e: np.ndarray,
    microphysics: InternalShockMicrophysics,
) -> np.ndarray:
    gamma_m = _gamma_m(branch, index, microphysics)
    gamma_c = _gamma_c(branch, index)
    gamma_max = _gamma_max(branch, index, microphysics)
    spectrum = np.zeros_like(gamma_e)
    active = (gamma_e >= min(gamma_m, gamma_c)) & (gamma_e <= gamma_max)
    if np.any(active):
        if gamma_c < gamma_m:
            cooled = active & (gamma_e < gamma_m)
            high = active & (gamma_e >= gamma_m)
            spectrum[cooled] = gamma_e[cooled] ** -2.0
            spectrum[high] = gamma_m ** (microphysics.electron_index_p - 1.0) * gamma_e[high] ** (-(microphysics.electron_index_p + 1.0))
        else:
            low = active & (gamma_e < gamma_c)
            high = active & (gamma_e >= gamma_c)
            spectrum[low] = gamma_e[low] ** (-microphysics.electron_index_p)
            spectrum[high] = gamma_c * gamma_e[high] ** (-(microphysics.electron_index_p + 1.0))
        t_age_comoving = _branch_age_comoving(branch, index)
        target_energy = branch.electron_luminosity_comoving_erg_s[index] * t_age_comoving
        energy_integral = np.trapezoid(spectrum * (gamma_e - 1.0) * constants.para_m_energy, gamma_e)
        if energy_integral > 0.0 and target_energy > 0.0:
            spectrum = spectrum * (target_energy / energy_integral)
    return spectrum


def _gamma_m(branch: BranchHistory, index: int, microphysics: InternalShockMicrophysics) -> float:
    dmdt_comoving = branch.thermal_luminosity_comoving_erg_s[index] / (
        branch.jump.specific_internal_energy * constants.para_c**2
    )
    electron_rate = microphysics.accelerated_electron_fraction * dmdt_comoving / constants.para_m_p
    electron_power = branch.electron_luminosity_comoving_erg_s[index]
    mean_energy = electron_power / electron_rate
    return 1.0 + (microphysics.electron_index_p - 2.0) / (microphysics.electron_index_p - 1.0) * mean_energy / constants.para_m_energy


def _gamma_c(branch: BranchHistory, index: int) -> float:
    t_age_comoving = _branch_age_comoving(branch, index)
    return 6.0 * pi * constants.para_m_e * constants.para_c / (
        constants.para_sigmat * branch.total_b_g[index] ** 2 * t_age_comoving
    )


def _gamma_max(branch: BranchHistory, index: int, microphysics: InternalShockMicrophysics) -> float:
    return np.sqrt(6.0 * pi * constants.para_e / (constants.para_sigmat * branch.total_b_g[index] * microphysics.acceleration_efficiency))


def _branch_age_comoving(branch: BranchHistory, index: int) -> float:
    shell_index = index + 0.5
    return branch.jump.crossing_time_lab_s * shell_index / (branch.radius_cm.size * branch.gamma[index])


def _compactness(branch: BranchHistory, index: int) -> float:
    return (
        constants.para_sigmat
        * branch.thermal_luminosity_comoving_erg_s[index]
        / (4.0 * pi * branch.radius_cm[index] * constants.para_m_e * constants.para_c**3)
    )


def _efficiency(branch: BranchHistory, index: int, microphysics: InternalShockMicrophysics) -> float:
    return microphysics.epsilon_e * branch.jump.specific_internal_energy / (
        1.0 + branch.jump.specific_internal_energy + branch.jump.magnetic_fraction
    )


def _validate_microphysics(microphysics: InternalShockMicrophysics) -> None:
    if not (microphysics.electron_index_p > 2.0):
        raise ValueError("electron_index_p must exceed 2 for finite injected electron energy.")
    if microphysics.epsilon_e < 0.0 or microphysics.epsilon_b < 0.0:
        raise ValueError("epsilon_e and epsilon_b must be non-negative.")
    if microphysics.accelerated_electron_fraction <= 0.0 or microphysics.acceleration_efficiency <= 0.0:
        raise ValueError("electron fraction and acceleration efficiency must be positive.")


def _validate_radiation_numerics(numerics: RadiationNumerics) -> None:
    if numerics.num_electron_gamma < 8 or numerics.num_photon_frequency < 8:
        raise ValueError("radiation grids are too small for synchrotron/SSC quadrature.")
    if numerics.electron_gamma_min < 1.0 or numerics.electron_gamma_max <= numerics.electron_gamma_min:
        raise ValueError("electron gamma grid must be ordered and start at gamma >= 1.")
    if numerics.photon_frequency_min_hz <= 0.0 or numerics.photon_frequency_max_hz <= numerics.photon_frequency_min_hz:
        raise ValueError("photon frequency grid must be positive and ordered.")
