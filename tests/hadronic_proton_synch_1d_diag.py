from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

from src import constants


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet


def _model(epsilon_b: float, epsilon_p: float) -> Model:
    return Model(
        TophatJet(0.1, 1.0e52, 300.0),
        ISM(1.0),
        Observer(1.0e26, 0.1, 0.0),
        Radiation(0.1, epsilon_b, 2.3, epsilon_p=epsilon_p, proton_synch=True),
        setups=Setups(
            electron_solver="fullhide_1d",
            num_gam_e=24,
            num_gam_p=40,
            num_nu=48,
            num_r=28,
            num_theta=16,
            num_tobs=28,
            hadronic_enabled=True,
        ),
    )


def _peak_frequency(track) -> float:
    emission = np.asarray(track.hadronic_gamma, dtype=float)
    idx = np.unravel_index(np.argmax(emission), emission.shape)
    spectrum = np.logspace(
        np.log10(1.0e-8 * constants.para_ev2hz),
        np.log10(1.0e4 * constants.para_tev2hz),
        emission.shape[0],
    )
    return float(spectrum[idx[0]])


def main() -> None:
    weak = _model(1.0e-3, 0.05).details(1.0e3, 1.0e6).fwd
    strong = _model(1.0e-3, 0.2).details(1.0e3, 1.0e6).fwd
    assert float(np.sum(strong.hadronic_gamma)) > float(np.sum(weak.hadronic_gamma))

    low_b = _model(1.0e-4, 0.2).details(1.0e3, 1.0e6).fwd
    high_b = _model(1.0e-2, 0.2).details(1.0e3, 1.0e6).fwd
    assert _peak_frequency(high_b) >= _peak_frequency(low_b)


if __name__ == "__main__":
    main()
