from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet


def _model(coupling: str | None) -> Model:
    setups = Setups(
        electron_solver="fullhide_1d",
        num_gam_e=18,
        num_nu=18,
        num_r=14,
        num_theta=8,
        num_tobs=14,
    )
    if coupling is not None:
        setups.electron_photon_coupling = coupling
    return Model(
        TophatJet(0.1, 1.0e52, 100.0),
        ISM(1.0),
        Observer(1.0e26, 0.05, 0.0),
        Radiation(0.1, 1.0e-3, 2.3, ssc=True),
        setups=setups,
    )


def main() -> None:
    times = np.array([1.0e4, 1.0e5], dtype=float)
    freqs = np.array([1.0e14, 1.0e18], dtype=float)
    default_flux = _model(None).flux_density(times, freqs).total
    explicit_flux = _model("separated").flux_density(times, freqs).total
    np.testing.assert_allclose(explicit_flux, default_flux, rtol=1.0e-14, atol=0.0)
    print("electron_photon_separated_regression_smoke: ok")


if __name__ == "__main__":
    main()
