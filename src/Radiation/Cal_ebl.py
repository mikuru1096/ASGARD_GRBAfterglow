from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d


PARA_TEV2HZ = 2.418e26


def cal_ebl(z, v_obs, model="Dominguez11.txt"):
    file_path = Path(__file__).parent / "EBL" / model
    table = np.loadtxt(file_path)

    redshifts = table[0, 1:]
    energies_hz = table[1:, 0] * PARA_TEV2HZ
    tau_values = table[1:, 1:]

    if z <= redshifts[0]:
        tau_at_z = tau_values[:, 0]
    elif z >= redshifts[-1]:
        tau_at_z = tau_values[:, -1]
    else:
        tau_at_z = np.array([np.interp(z, redshifts, tau_values[i, :]) for i in range(tau_values.shape[0])])

    tau_interp = interp1d(energies_hz, tau_at_z, kind="linear", bounds_error=False, fill_value=(tau_at_z[0], tau_at_z[-1]))
    tau_obs = tau_interp(v_obs)
    absorption = np.exp(-tau_obs)

    absorption[v_obs < energies_hz[0]] = 1.0
    absorption[v_obs > energies_hz[-1]] = 1.0e-30
    return absorption
