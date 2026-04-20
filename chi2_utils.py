from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import interpolate


def load_observation_table(data_file: Path) -> np.ndarray:
    table = np.loadtxt(data_file)
    if table.ndim == 1:
        table = table.reshape(1, -1)
    return table


def build_model_interpolators(model_curves: np.ndarray, model_serial: np.ndarray) -> list[interpolate.interp1d]:
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


def validate_model_range(sample_points: np.ndarray, model_serial: np.ndarray) -> None:
    if np.min(sample_points) < model_serial[0] or np.max(sample_points) > model_serial[-1]:
        raise ValueError("The model curve cannot fully cover the data range.")
