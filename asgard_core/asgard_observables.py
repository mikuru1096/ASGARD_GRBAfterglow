from __future__ import annotations

import numpy as np

from src import constants


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
