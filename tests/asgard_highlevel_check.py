from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import Fitter, ISM, Magnetar, Model, ObsData, Observer, ParamDef, Radiation, Scale, TophatJet, units


def main() -> None:
    times_s = np.logspace(2.0, 3.5, 5)
    freq_opt = np.full(times_s.shape, 4.63e14, dtype=float)

    model = Model(
        medium=ISM(n0=0.1),
        jet=TophatJet(E_iso=1.0e53, Gamma0=1.0e2, theta_c=1.0e-1, duration=10.0, magnetar=Magnetar(L0=1.0e48, t0=1.0e4, q=2.0)),
        observer=Observer(z=0.4, theta_obs=0.0, lumi_dist=7.0e27),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.5, xi_N=1.0e-1, ssc=False, kn=False),
        resolutions=(6, 12, 48),
    )
    truth = model.flux_density_grid(times_s, np.array([4.63e14], dtype=float))[0]
    band_flux = model.flux(1.0e3, 1.0e17, 1.0e19, num_points=32)

    data = ObsData()
    data.add_flux_density(t=times_s, nu=freq_opt, f_nu=truth, err=0.05 * np.maximum(truth, 1.0e-99))
    data.add_flux(t=1.0e3, nu_min=1.0e17, nu_max=1.0e19, value=band_flux, err=max(0.05 * band_flux, 1.0e-99), num_points=32)

    cfg = {
        "medium_type": "ism",
        "n0": 0.1,
        "jet_type": "tophat",
        "E_iso": 1.0e53,
        "Gamma0": 1.0e2,
        "theta_c": 1.0e-1,
        "z": 0.4,
        "lumi_dist": 7.0e27,
        "eps_e": 1.0e-1,
        "eps_B": 1.0e-3,
        "p": 2.5,
        "f_e": 1.0e-1,
        "ssc": False,
        "resolutions": (6, 12, 48),
    }
    fitter = Fitter(data=data, cfg=cfg)
    loglike = fitter.loglike({"E_iso": 53.0})
    assert np.isfinite(loglike)

    result = fitter.fit(
        [
            ParamDef("E_iso", 53.0, 53.0, Scale.log),
            ParamDef("p", 2.5, 2.5, Scale.fixed),
        ],
        nsteps=16,
        nburn=4,
        resolution=(6, 12, 48),
    )
    assert "E_iso" in result.best_params
    assert "p" in result.best_params
    assert result.top_k_params is not None
    assert result.top_k_params[0]["E_iso"] == result.best_params["E_iso"]

    screened = ObsData.logscale_screen(np.logspace(2.0, 6.0, 100), data_density=3.0)
    assert screened.ndim == 1
    assert screened.size > 0
    assert np.isclose(units.day, 86400.0)
    assert np.isclose(units.GHz, 1.0e9)
    assert np.isclose(units.mJy, 1.0e-3)

    print("PASS: ASGARD high-level check succeeded.")


if __name__ == "__main__":
    main()
