import warnings

import numpy as np

from asgard_core.asgard_paths import DATA_LIGHT_CURVE_DIR

from .chi2_utils import build_model_interpolators, load_observation_table, validate_model_range
from .extinc import opt_extinction


def cal_chi2_light_curve(bands_fit, model_curves, model_serial, frequency, Rv, Ebv, zeropointflux, z, f_sys, lyman_ar=0.0):
    data_dir = DATA_LIGHT_CURVE_DIR
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory does not exist: {data_dir}")

    chi2_total = 0.0
    model_interpolators = build_model_interpolators(model_curves, model_serial)

    for data_file in data_dir.glob("*"):
        if not data_file.is_file():
            continue
        name = data_file.stem
        if name not in bands_fit:
            continue
        try:
            table = load_observation_table(data_file)
            range_data, flux_data, flux_err = _parse_observation_data(table, name)
            band_idx = bands_fit.index(name)
            if name.startswith("opt"):
                flux_data, flux_err = opt_extinction(
                    flux_data,
                    flux_err,
                    frequency[band_idx],
                    Rv,
                    Ebv,
                    zeropointflux[band_idx],
                    redshift=z,
                    lyman_ar=lyman_ar,
                )
            range_data = _convert_time_units(range_data, name)
            validate_model_range(range_data, model_serial)
            fit_flux = model_interpolators[band_idx](range_data)
            if np.any(np.isnan(fit_flux)):
                raise ValueError("Some data points are beyond the scope of the model.")
            sigma = _get_uncertainties(flux_data, flux_err, fit_flux, table.shape[1])
            variance = _calculate_variance(flux_data, sigma, f_sys)
            chi2_total += np.sum((fit_flux - flux_data) ** 2 / variance)
        except Exception as exc:
            warnings.warn(f"An error occurred while processing the file {data_file.name}: {str(exc)}")
            continue
    return chi2_total


def _parse_observation_data(table, name):
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


def _convert_time_units(range_data, name):
    return range_data * 86400


def _get_uncertainties(flux_data, flux_err, fit_flux, n_cols):
    return np.where(fit_flux > flux_data, flux_err, -flux_err) if n_cols == 6 else flux_err


def _calculate_variance(flux_data, sigma, f_sys):
    return (flux_data * 0.1) ** 2 if f_sys <= 0 else (flux_data * f_sys) ** 2 + sigma**2
