from __future__ import annotations

from functools import lru_cache

import numpy as np

from src import constants
import numpy as np


FITTING_BANDS = ("xrt", "optr", "optz", "opti", "optg", "9GHz", "5.5GHz", "3GHz")
OUTPUT_BANDS = FITTING_BANDS + ("1GeV", "1TeV")
POINT_SOURCE_BANDS = FITTING_BANDS[1:]
FITTING_FREQUENCIES_HZ = np.array([2.7e17, 4.63e14, 3.39e14, 4.01e14, 6.42e14, 9e9, 5.5e9, 3e9], dtype=float)
POINT_SOURCE_FREQUENCIES_HZ = FITTING_FREQUENCIES_HZ[1:]
HIGH_ENERGY_FREQUENCIES_HZ = np.array([1e24, 1e27], dtype=float)
ZEROPOINT_FLUX = np.array([0.0, -48.6, -48.6, -48.6, -48.6, 0.0, 0.0, 0.0], dtype=float)
LIGHT_CURVE_PLOT_SPECS = (
    ("xrt", 1.0, "k", "$0.5-10$ keV"),
    ("optr", 0.3e-16, "#FFD700", "r $\\times$ 0.3e-16"),
    ("optz", 0.1e-16, "#FF8C00", "z $\\times$ 0.1e-16"),
    ("opti", 0.03e-16, "#FF4500", "i $\\times$ 0.03e-16"),
    ("optg", 0.01e-16, "#FF0000", "g $\\times$ 0.01e-16"),
    ("9GHz", 1.0e-19, "#FF1493", "9GHz $\\times$ 1e-19"),
    ("5.5GHz", 0.3e-19, "#8A2BE2", "5.5GHz $\\times$ 0.3e-19"),
    ("3GHz", 0.1e-19, "#4B0082", "3GHz $\\times$ 0.1e-19"),
    ("1GeV", 1.0e4, "#9400D3", "1e24 $\\times$ 1e4"),
    ("1TeV", 1.0e5, "#FF00FF", "1e27 $\\times$ 1e5"),
)


@lru_cache(maxsize=1)
def _cached_multiband_observer_frequencies() -> tuple[int, np.ndarray]:
    num_xrt = 8
    xrt_low = 0.5 * constants.para_kev2hz
    xrt_high = 10.0 * constants.para_kev2hz
    xrt = np.logspace(np.log10(xrt_low), np.log10(xrt_high), num_xrt)
    return num_xrt, np.concatenate((xrt, POINT_SOURCE_FREQUENCIES_HZ, HIGH_ENERGY_FREQUENCIES_HZ), axis=0)


def build_multiband_observer_frequencies() -> tuple[int, np.ndarray]:
    num_xrt, frequencies_hz = _cached_multiband_observer_frequencies()
    return num_xrt, frequencies_hz.copy()


def combine_multiband_flux(flux_matrix: np.ndarray, frequencies_hz: np.ndarray, num_xrt: int) -> np.ndarray:
    num_point_source_bands = len(POINT_SOURCE_BANDS)
    xrt_energy_flux = np.trapezoid(flux_matrix[:num_xrt].T, frequencies_hz[:num_xrt], axis=1).reshape(1, -1)
    optical_and_radio = flux_matrix[num_xrt : num_xrt + num_point_source_bands] * 1e29
    gev_index = num_xrt + num_point_source_bands
    tev_index = gev_index + 1
    gev_energy_flux = (flux_matrix[gev_index] * frequencies_hz[gev_index]).reshape(1, -1)
    tev_energy_flux = (flux_matrix[tev_index] * frequencies_hz[tev_index]).reshape(1, -1)
    return np.vstack([xrt_energy_flux, optical_and_radio, gev_energy_flux, tev_energy_flux])
