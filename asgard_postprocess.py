from __future__ import annotations

from typing import Optional

import numpy as np

from asgard_models import FitConfig, PhysicalSolution, SimulationSetup
from asgard_observables import (
    FITTING_BANDS,
    FITTING_FREQUENCIES_HZ,
    ZEROPOINT_FLUX,
    build_multiband_observer_frequencies,
    combine_multiband_flux,
)
from cal_chi2_light_curve import cal_chi2_light_curve as cal_chi2_lc
from src import Interpolation


def interpolate_observed_flux(
    setup: SimulationSetup,
    characteristic_time_s: np.ndarray,
    gamma: np.ndarray,
    radius_cm: np.ndarray,
    absorbed_spectral_flux: np.ndarray,
    frequencies_hz: np.ndarray,
    config: FitConfig,
) -> np.ndarray:
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    order = np.argsort(frequencies_hz)
    sorted_frequencies = frequencies_hz[order]
    flux_sorted = Interpolation.sed_interpolation(
        setup.boundary,
        characteristic_time_s,
        gamma,
        radius_cm,
        absorbed_spectral_flux,
        setup.seed_frequency_hz,
        sorted_frequencies,
        setup.observer_time_s,
        config.num_theta,
        config.num_phi,
        config.num_threads,
    )
    if np.array_equal(order, np.arange(order.shape[0])):
        return flux_sorted

    flux_matrix = np.empty_like(flux_sorted)
    flux_matrix[order] = flux_sorted
    return flux_matrix


def compute_observed_flux_matrix(
    setup: SimulationSetup,
    physical: PhysicalSolution,
    frequencies_hz: np.ndarray,
    config: FitConfig,
) -> np.ndarray:
    return interpolate_observed_flux(
        setup,
        physical.characteristic_time_s,
        physical.gamma,
        physical.radius_cm,
        physical.absorbed_spectral_flux,
        frequencies_hz,
        config,
    )


def compute_band_fluxes(
    setup: SimulationSetup,
    physical: PhysicalSolution,
    config: FitConfig,
) -> np.ndarray:
    num_xrt, band_frequencies = build_multiband_observer_frequencies()
    band_flux_matrix = compute_observed_flux_matrix(
        setup,
        physical,
        band_frequencies,
        config,
    )
    return combine_multiband_flux(band_flux_matrix, band_frequencies, num_xrt)


def compute_light_curve_redchi(
    bands_flux: np.ndarray,
    t_obs_s: np.ndarray,
    config: FitConfig,
) -> float:
    redchi = cal_chi2_lc(
        list(FITTING_BANDS),
        bands_flux[: len(FITTING_BANDS)],
        t_obs_s,
        FITTING_FREQUENCIES_HZ,
        config.rv,
        config.ebv,
        ZEROPOINT_FLUX,
        config.z,
        config.f_sys,
    )
    return float(redchi)


def build_spectrum_frequency_grid(config: FitConfig) -> np.ndarray:
    spec = config.spectrum_output
    return np.logspace(np.log10(spec.nu_min_hz), np.log10(spec.nu_max_hz), spec.num_nu_obs)


def compute_spectrum_flux(
    setup: SimulationSetup,
    physical: PhysicalSolution,
    config: FitConfig,
) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    if not config.spectrum_output.enabled:
        return None, None

    spectrum_freq_hz = build_spectrum_frequency_grid(config)
    spectrum_fnu = compute_observed_flux_matrix(
        setup,
        physical,
        spectrum_freq_hz,
        config,
    )
    return spectrum_freq_hz, spectrum_fnu
