from pathlib import Path

import numpy as np

from src import constants

_EBL_DIR = Path(__file__).resolve().parent / "EBL"


def cal_ebl(z, v_obs, model="Dominguez11.txt"):
    """Return EBL attenuation exp(-tau) for observed frequencies."""
    table = np.loadtxt(_EBL_DIR / model)
    v_obs = np.asarray(v_obs, dtype=float)

    redshifts = table[0, 1:]
    energies_hz = table[1:, 0] * constants.para_tev2hz
    tau_values = table[1:, 1:]

    if z <= redshifts[0]:
        tau_at_z = tau_values[:, 0]
    elif z >= redshifts[-1]:
        tau_at_z = tau_values[:, -1]
    else:
        tau_at_z = np.array([np.interp(z, redshifts, tau_values[i, :]) for i in range(tau_values.shape[0])])

    tau_obs = np.interp(v_obs, energies_hz, tau_at_z, left=tau_at_z[0], right=tau_at_z[-1])
    absorption = np.exp(-tau_obs)
    return np.where(v_obs < energies_hz[0], 1.0, np.where(v_obs > energies_hz[-1], 1.0e-30, absorption))
