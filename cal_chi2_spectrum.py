from pathlib import Path
import warnings

import numpy as np
from scipy import interpolate


def cal_chi2_spectrum(name_fit, model_curves, model_serial):
    data_dir = Path(__file__).parent / "data_spectrum"
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory does not exist: {data_dir}")

    chi2_total = 0.0
    model_interpolators = []
    for i in range(len(name_fit)):
        interp = interpolate.interp1d(model_serial, model_curves[i, :], kind="linear", bounds_error=False, fill_value=np.nan)
        model_interpolators.append(interp)

    for data_file in data_dir.glob("*"):
        if not data_file.is_file():
            continue

        name = data_file.stem
        if name not in name_fit:
            continue

        try:
            table = np.loadtxt(data_file)
            if table.ndim == 1:
                table = table.reshape(1, -1)

            range_data, flux_data, flux_err = _parse_spectrum_data(table)
            _validate_model_range(range_data, model_serial)

            band_idx = name_fit.index(name)
            fit_flux = model_interpolators[band_idx](range_data)
            if np.any(np.isnan(fit_flux)):
                raise ValueError("Some data points are beyond the scope of the model.")

            sigma = _get_spectrum_uncertainties(flux_data, flux_err, fit_flux, table.shape[1])
            num_data = len(range_data)
            temp_chi2 = np.sum(((fit_flux - flux_data) / sigma) ** 2) / num_data
            chi2_total += temp_chi2
        except Exception as exc:
            warnings.warn(f"An error occurred while processing the file {data_file.name}: {exc}")
            continue

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


def _validate_model_range(range_data, model_serial):
    if np.min(range_data) < model_serial[0] or np.max(range_data) > model_serial[-1]:
        raise ValueError("The model curve cannot fully cover the data range.")


def _get_spectrum_uncertainties(flux_data, flux_err, fit_flux, n_cols):
    if n_cols == 6:
        flux_err_up = flux_err[0]
        flux_err_down = flux_err[1]
        return np.where(fit_flux > flux_data, flux_err_up, flux_err_down)
    return flux_err
