import warnings

import numpy as np

from asgard_core.asgard_paths import DATA_SPECTRUM_DIR

from .chi2_utils import build_model_interpolators, load_observation_table, validate_model_range


def cal_chi2_spectrum(name_fit, model_curves, model_serial):
    data_dir = DATA_SPECTRUM_DIR
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory does not exist: {data_dir}")

    chi2_total = 0.0
    model_interpolators = build_model_interpolators(model_curves, model_serial)
    for data_file in data_dir.glob("*"):
        if not data_file.is_file():
            continue
        name = data_file.stem
        if name not in name_fit:
            continue
        try:
            table = load_observation_table(data_file)
            range_data, flux_data, flux_err = _parse_spectrum_data(table)
            validate_model_range(range_data, model_serial)
            band_idx = name_fit.index(name)
            fit_flux = model_interpolators[band_idx](range_data)
            if np.any(np.isnan(fit_flux)):
                raise ValueError("Some data points are beyond the scope of the model.")
            if table.shape[1] == 6:
                sigma = np.where(fit_flux > flux_data, flux_err[0], flux_err[1])
            else:
                sigma = flux_err
            chi2_total += np.sum(((fit_flux - flux_data) / sigma) ** 2) / len(range_data)
        except Exception as exc:
            warnings.warn(f"An error occurred while processing the file {data_file.name}: {exc}")
    return chi2_total


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
