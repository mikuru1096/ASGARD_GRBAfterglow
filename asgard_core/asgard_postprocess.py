from __future__ import annotations

from typing import Optional

import numpy as np

from asgard_core.asgard_config import FitConfig, PhysicalSolution, SimulationSetup
from asgard_core.asgard_observables import (
    FITTING_BANDS,
    FITTING_FREQUENCIES_HZ,
    ZEROPOINT_FLUX,
    build_multiband_observer_frequencies,
    combine_multiband_flux,
)
from asgard_core.support.chi2_light_curve import cal_chi2_light_curve as cal_chi2_lc
from asgard_core.support.chi2_spectrum import cal_chi2_spectrum as cal_chi2_spec
from asgard_core.asgard_paths import DATA_SPECTRUM_DIR
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
        config.lyman_ar,
    )
    return float(redchi)


def build_spectrum_dataset_names() -> tuple[str, ...]:
    data_dir = DATA_SPECTRUM_DIR
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory does not exist: {data_dir}")
    return tuple(path.stem for path in sorted(data_dir.iterdir()) if path.is_file())


def compute_spectrum_redchi(
    spectrum_fnu: np.ndarray,
    spectrum_freq_hz: np.ndarray,
    spectrum_names: Optional[tuple[str, ...] | list[str]] = None,
) -> float:
    spectrum_freq_hz = np.asarray(spectrum_freq_hz, dtype=float)
    model_curves = np.asarray(spectrum_fnu, dtype=float)
    if model_curves.ndim == 1:
        model_curves = model_curves.reshape(1, -1)
    elif model_curves.ndim != 2:
        raise ValueError("spectrum_fnu must be a 1D or 2D array.")

    if model_curves.shape[1] != spectrum_freq_hz.size:
        if model_curves.shape[0] == spectrum_freq_hz.size:
            model_curves = model_curves.T
        else:
            raise ValueError("Could not align spectrum_fnu with spectrum_freq_hz.")

    names = build_spectrum_dataset_names() if spectrum_names is None else tuple(spectrum_names)
    if model_curves.shape[0] != len(names):
        raise ValueError("Number of model spectra must match number of spectrum dataset names.")

    redchi = cal_chi2_spec(list(names), model_curves, spectrum_freq_hz)
    return float(redchi)


def select_spectrum_time_index(t_obs_s: np.ndarray, time_s: Optional[float]) -> int:
    t_obs_s = np.asarray(t_obs_s, dtype=float)
    if t_obs_s.ndim != 1 or t_obs_s.size == 0:
        raise ValueError("t_obs_s must be a non-empty 1D array.")
    if time_s is None:
        return int(t_obs_s.size - 1)
    return int(np.abs(t_obs_s - float(time_s)).argmin())


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
