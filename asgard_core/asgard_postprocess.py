from __future__ import annotations

from pathlib import Path

import numpy as np
from extinction import fitzpatrick99 as f99
from scipy import interpolate

from asgard_core.asgard_config import RuntimeConfig, SimulationSetup
from src import Interpolation, constants

PACKAGE_ROOT = Path(__file__).resolve().parent
DATA_LIGHT_CURVE_DIR = PACKAGE_ROOT / "data_light_curve"
DATA_SPECTRUM_DIR = PACKAGE_ROOT / "data_spectrum"
FITTING_BANDS = ("xrt", "optr", "optz", "opti", "optg", "9GHz", "5.5GHz", "3GHz")
OUTPUT_BANDS = FITTING_BANDS + ("1GeV", "1TeV")
POINT_SOURCE_BANDS = FITTING_BANDS[1:]
FITTING_FREQUENCIES_HZ = np.array([2.7e17, 4.63e14, 3.39e14, 4.01e14, 6.42e14, 9e9, 5.5e9, 3e9], dtype=float)
_NUM_XRT_BINS = 8
_MULTIBAND_OBSERVER_FREQUENCIES_HZ = np.concatenate(
    (
        np.logspace(np.log10(0.5 * constants.para_kev2hz), np.log10(10.0 * constants.para_kev2hz), _NUM_XRT_BINS),
        FITTING_FREQUENCIES_HZ[1:],
        np.array([1e24, 1e27], dtype=float),
    ),
    axis=0,
)
ZEROPOINT_FLUX = np.array([0.0, -48.6, -48.6, -48.6, -48.6, 0.0, 0.0, 0.0], dtype=float)


def build_multiband_observer_frequencies() -> tuple[int, np.ndarray]:
    return _NUM_XRT_BINS, _MULTIBAND_OBSERVER_FREQUENCIES_HZ.copy()


def combine_multiband_flux(flux_matrix: np.ndarray, frequencies_hz: np.ndarray, num_xrt: int) -> np.ndarray:
    num_point_source_bands = len(POINT_SOURCE_BANDS)
    xrt_energy_flux = np.trapezoid(flux_matrix[:num_xrt].T, frequencies_hz[:num_xrt], axis=1).reshape(1, -1)
    optical_and_radio = flux_matrix[num_xrt : num_xrt + num_point_source_bands] * 1e29
    gev_index = num_xrt + num_point_source_bands
    tev_index = gev_index + 1
    gev_energy_flux = (flux_matrix[gev_index] * frequencies_hz[gev_index]).reshape(1, -1)
    tev_energy_flux = (flux_matrix[tev_index] * frequencies_hz[tev_index]).reshape(1, -1)
    return np.vstack([xrt_energy_flux, optical_and_radio, gev_energy_flux, tev_energy_flux])


def interpolate_observed_flux(
    setup: SimulationSetup,
    characteristic_time_s: np.ndarray,
    gamma: np.ndarray,
    radius_cm: np.ndarray,
    absorbed_spectral_flux: np.ndarray,
    frequencies_hz: np.ndarray,
    config: RuntimeConfig,
) -> np.ndarray:
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    order = np.argsort(frequencies_hz)
    sorted_frequencies = frequencies_hz[order]
    geometry_kernel = str(config.geometry_kernel).lower()
    if geometry_kernel == "sed_adaptive_theta":
        interpolate = Interpolation.sed_adaptive_theta
        interpolate_args = (
            float(config.projection_adaptive_rtol),
            int(config.projection_adaptive_max_depth),
        )
    else:
        interpolate = Interpolation.sed_interpolation
        interpolate_args = ()
    flux_sorted = interpolate(
        setup.boundary,
        characteristic_time_s,
        gamma,
        radius_cm,
        absorbed_spectral_flux,
        setup.seed_frequency_hz,
        sorted_frequencies,
        setup.observer_time_s,
        *interpolate_args,
        config.eats_num_theta,
        config.eats_num_phi,
        config.num_threads,
    )
    if np.array_equal(order, np.arange(order.shape[0])):
        return flux_sorted

    flux_matrix = np.empty_like(flux_sorted)
    flux_matrix[order] = flux_sorted
    return flux_matrix


def compute_light_curve_redchi(
    bands_flux: np.ndarray,
    t_obs_s: np.ndarray,
    config: RuntimeConfig,
) -> float:
    redchi = _cal_chi2_light_curve(
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
    spectrum_names: tuple[str, ...] | list[str] | None = None,
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

    redchi = _cal_chi2_spectrum(list(names), model_curves, spectrum_freq_hz)
    return float(redchi)


def select_spectrum_time_index(t_obs_s: np.ndarray, time_s: float | None) -> int:
    t_obs_s = np.asarray(t_obs_s, dtype=float)
    if t_obs_s.ndim != 1 or t_obs_s.size == 0:
        raise ValueError("t_obs_s must be a non-empty 1D array.")
    if time_s is None:
        return int(t_obs_s.size - 1)
    return int(np.abs(t_obs_s - float(time_s)).argmin())


def build_spectrum_frequency_grid(config: RuntimeConfig) -> np.ndarray:
    spec = config.spectrum_output
    return np.logspace(np.log10(spec.nu_min_hz), np.log10(spec.nu_max_hz), spec.num_nu_obs)


def _cal_chi2_light_curve(
    bands_fit,
    model_curves,
    model_serial,
    frequency,
    Rv,
    Ebv,
    zeropointflux,
    z,
    f_sys,
    lyman_ar=0.0,
):
    data_dir = DATA_LIGHT_CURVE_DIR
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory does not exist: {data_dir}")

    chi2_total = 0.0
    model_interpolators = _build_model_interpolators(model_curves, model_serial)

    for data_file in data_dir.glob("*"):
        if not data_file.is_file():
            continue
        name = data_file.stem
        if name not in bands_fit:
            continue
        table = _load_observation_table(data_file)
        range_data, flux_data, flux_err = _parse_observation_data(table)
        band_idx = bands_fit.index(name)
        if name.startswith("opt"):
            flux_data, flux_err = _opt_extinction(
                flux_data,
                flux_err,
                frequency[band_idx],
                Rv,
                Ebv,
                zeropointflux[band_idx],
                redshift=z,
                lyman_ar=lyman_ar,
            )
        range_data = range_data * 86400
        _validate_model_range(range_data, model_serial)
        fit_flux = model_interpolators[band_idx](range_data)
        if np.any(np.isnan(fit_flux)):
            raise ValueError("Some data points are beyond the scope of the model.")
        sigma = np.where(fit_flux > flux_data, flux_err, -flux_err) if table.shape[1] == 6 else flux_err
        variance = (flux_data * 0.1) ** 2 if f_sys <= 0 else (flux_data * f_sys) ** 2 + sigma**2
        chi2_total += np.sum((fit_flux - flux_data) ** 2 / variance)
    return chi2_total


def _cal_chi2_spectrum(name_fit, model_curves, model_serial):
    data_dir = DATA_SPECTRUM_DIR
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory does not exist: {data_dir}")

    chi2_total = 0.0
    model_interpolators = _build_model_interpolators(model_curves, model_serial)
    for data_file in data_dir.glob("*"):
        if not data_file.is_file():
            continue
        name = data_file.stem
        if name not in name_fit:
            continue
        table = _load_observation_table(data_file)
        range_data, flux_data, flux_err = _parse_spectrum_data(table)
        _validate_model_range(range_data, model_serial)
        band_idx = name_fit.index(name)
        fit_flux = model_interpolators[band_idx](range_data)
        if np.any(np.isnan(fit_flux)):
            raise ValueError("Some data points are beyond the scope of the model.")
        sigma = np.where(fit_flux > flux_data, flux_err[0], flux_err[1]) if table.shape[1] == 6 else flux_err
        chi2_total += np.sum(((fit_flux - flux_data) / sigma) ** 2) / len(range_data)
    return chi2_total


def _load_observation_table(data_file: Path) -> np.ndarray:
    table = np.loadtxt(data_file)
    if table.ndim == 1:
        table = table.reshape(1, -1)
    return table


def _build_model_interpolators(model_curves: np.ndarray, model_serial: np.ndarray) -> list[interpolate.interp1d]:
    return [
        interpolate.interp1d(
            model_serial,
            model_curves[index, :],
            kind="linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        for index in range(len(model_curves))
    ]


def _validate_model_range(sample_points: np.ndarray, model_serial: np.ndarray) -> None:
    if np.min(sample_points) < model_serial[0] or np.max(sample_points) > model_serial[-1]:
        raise ValueError("The model curve cannot fully cover the data range.")


def _parse_observation_data(table):
    n_cols = table.shape[1] if table.ndim > 1 else table.size
    if n_cols == 6:
        range_data, flux_data = table[:, 0], table[:, 3]
        flux_err_up, flux_err_down = table[:, 4], -table[:, 5]
        flux_err = np.where(flux_data > table[:, 3], flux_err_up, flux_err_down)
    elif n_cols == 4:
        range_data = 0.5 * (table[:, 0] + table[:, 1]) if table[0, 0] < table[0, 1] else table[:, 0]
        flux_data, flux_err = table[:, 2], table[:, 3]
    elif n_cols == 3:
        range_data, flux_data, flux_err = table[:, 0], table[:, 1], table[:, 2]
    elif n_cols == 2:
        range_data, flux_data = table[:, 0], table[:, 1]
        flux_err = flux_data * 0.1
    else:
        raise ValueError(f"The observation data should be 2 to 6 columns. Currently, there are {n_cols} columns.")
    return range_data, flux_data, flux_err


def _parse_spectrum_data(table):
    n_cols = table.shape[1] if table.ndim > 1 else table.size
    if n_cols == 6:
        range_data = table[:, 0]
        flux_data = table[:, 3]
        flux_err = np.vstack((table[:, 4], np.abs(table[:, 5])))
    elif n_cols == 4:
        range_data = table[:, 0]
        flux_data = table[:, 2]
        flux_err = table[:, 3]
    elif n_cols == 3:
        range_data = table[:, 0]
        flux_data = table[:, 1]
        flux_err = table[:, 2]
    elif n_cols == 2:
        range_data = table[:, 0]
        flux_data = table[:, 1]
        flux_err = flux_data * 0.1
    else:
        raise ValueError(f"The observation data should be 2 to 6 columns. Currently, there are {n_cols} columns.")
    return range_data, flux_data, flux_err


def _frequency_to_wavelength_micron(frequency):
    return np.array([2.997e10 / frequency * 1e4], dtype=float)


def _flux_from_mag(mag_data, mag_err, zeropointflux):
    flux_data_deredden = 10 ** (0.4 * (zeropointflux - mag_data))
    flux_data_err = 0.4 * np.log(10.0) * flux_data_deredden * mag_err
    return flux_data_deredden, flux_data_err


def _opt_extinction(
    mag_data,
    mag_err,
    frequency,
    Rv,
    Ebv,
    zeropointflux,
    redshift=None,
    lyman_ar=0.0,
):
    Av = Rv * Ebv
    wave = np.array([2.997e10 / frequency * 1e8], dtype=float)
    mag_data_deredden = mag_data - f99(wave, Av, Rv)
    if redshift is not None and lyman_ar != 0.0:
        wave_in_mu_m = _frequency_to_wavelength_micron(frequency)
        if np.any((wave_in_mu_m > 0.6) & (wave_in_mu_m < 0.68)):
            mag_data_deredden = mag_data_deredden - lyman_ar
    return _flux_from_mag(mag_data_deredden, mag_err, zeropointflux)
