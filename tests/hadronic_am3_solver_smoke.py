from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet


def _build_model(*, include_pg: bool = False, include_neutrino: bool = False) -> Model:
    return Model(
        TophatJet(0.1, 1.0e52, 300.0),
        ISM(1.0),
        Observer(1.0e26, 0.1, 0.0),
        Radiation(
            0.1,
            1.0e-3,
            2.3,
            epsilon_p=0.2,
            ssc=True,
            proton_synch=True,
            pg=include_pg,
            neutrino=include_neutrino,
        ),
        setups=Setups(
            electron_solver="fullhide_1d",
            hadronic_enabled=True,
            hadronic_solver="am3_1d",
            num_gam_e=24,
            num_gam_p=40,
            num_nu=40,
            num_nu_nu=24,
            num_r=24,
            num_theta=16,
            num_tobs=24,
        ),
    )


def main() -> None:
    details = _build_model().details(1.0e3, 1.0e6).fwd
    assert details.gamma_p is not None
    assert details.dN_dgamma_p is not None
    assert details.hadronic_gamma is not None
    assert np.all(np.isfinite(details.hadronic_gamma))
    assert details.neutrino_luminosity is None


if __name__ == "__main__":
    main()
