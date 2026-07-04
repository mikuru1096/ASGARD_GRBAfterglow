from pathlib import Path

import numpy as np

from src import constants

_EBL_DIR = Path(__file__).resolve().parent / "EBL"


def cal_ebl(z, v_obs, model="Dominguez11.txt"):
    """Return EBL attenuation exp(-tau) for observed frequencies."""
    table = np.loadtxt(_EBL_DIR / model)
    v_obs = np.asarray(v_obs, dtype=float)

    redshift = table[0, 1:]
    energy_hz = table[1:, 0] * constants.para_tev2hz
    tau_grid = table[1:, 1:]

    if z < 0.0 or z > redshift[-1]:
        raise ValueError(f"EBL model {model} supports 0 <= z <= {redshift[-1]}")
    if np.any(v_obs < 0.0) or np.any(v_obs > energy_hz[-1]):
        raise ValueError(f"EBL model {model} supports 0 <= nu <= {energy_hz[-1]} Hz")

    redshift = np.concatenate(([0.0], redshift))
    tau_grid = np.column_stack((np.zeros(tau_grid.shape[0]), tau_grid))
    tau_z = np.array([np.interp(z, redshift, tau_grid[i, :]) for i in range(tau_grid.shape[0])])
    tau_obs = np.interp(v_obs, energy_hz, tau_z, left=0.0)
    return np.exp(-tau_obs)
