from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import Fitter, ISM, Model, ObsData, Observer, ParamDef, Radiation, Scale, Setups, TophatJet


def main() -> None:
    times_s = np.logspace(2.0, 4.0, 8)
    frequencies_hz = np.full(times_s.shape, 4.63e14, dtype=float)
    model = Model(
        medium=ISM(0.1),
        jet=TophatJet(E_iso=1.0e53, lf=1.0e2, theta_j=1.0e-1),
        observer=Observer(z=0.4, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.5, xi_N=1.0e-1, ssc=True, kn=False),
        setups=Setups(num_threads=4, num_r=180, num_theta=180, num_tobs=64),
    )
    truth = model.flux_density_grid(times_s, np.array([4.63e14], dtype=float))[0]

    data = ObsData()
    data.add_flux_density(times_s, frequencies_hz, truth, 0.05 * np.maximum(truth, 1.0e-99))

    fitter = Fitter(
        model=model,
        data=data,
        params=[
            ParamDef("log10_Eiso", "jet.E_iso", 52.0, 54.0, Scale.LOG10),
            ParamDef("p", "fwd_rad.p", 2.1, 2.9, Scale.LINEAR),
        ],
    )
    values = {"log10_Eiso": 53.0, "p": 2.5}
    flux = fitter.flux_density_grid(values, times_s, np.array([4.63e14], dtype=float))
    assert flux.shape == (1, times_s.shape[0])
    loglike = fitter.loglike(values)
    assert np.isfinite(loglike)

    print("PASS: ASGARD fitter check succeeded.")


if __name__ == "__main__":
    main()
