from __future__ import annotations

from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import VegasAfterglow.api as va_api
from VegasAfterglow import ISM, Model, Observer, Radiation, Setups, TophatJet


def _build_model() -> Model:
    return Model(
        medium=ISM(0.1),
        jet=TophatJet(E_iso=1.0e53, Gamma0=100.0, theta_c=0.1),
        observer=Observer(lumi_dist=1.0e28, z=0.4, theta_obs=0.0, phi_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.5, xi_N=0.1, ssc=True, kn=False),
        setups=Setups(
            num_threads=8,
            num_gam_e=121,
            num_nu=121,
            num_r=140,
            num_theta=120,
            num_phi=1,
            num_tobs=64,
            observer_time_min_s=1.0e2,
            observer_time_max_s=1.0e6,
        ),
    )


def main() -> None:
    time_s = np.array([1.0e4], dtype=float)
    frequencies_hz = np.logspace(9.0, 21.0, 128)

    model = _build_model()
    with patch.object(va_api, "solve_model_state_from_setup", wraps=va_api.solve_model_state_from_setup) as solve_mock:
        adaptive = model.flux_density_grid(time_s, frequencies_hz)
        solve_count = solve_mock.call_count

    reference_model = _build_model()
    reference = reference_model._compute_raw(time_s, frequencies_hz)

    rel_total = np.max(np.abs(adaptive.total - reference.total) / np.maximum(np.abs(reference.total), 1.0e-99))
    rel_sync = np.max(np.abs(adaptive.fwd.sync - reference.fwd.sync) / np.maximum(np.abs(reference.fwd.sync), 1.0e-99))
    rel_ssc = np.max(np.abs(adaptive.fwd.ssc - reference.fwd.ssc) / np.maximum(np.abs(reference.fwd.ssc), 1.0e-99))

    assert solve_count == 1, solve_count
    assert rel_total <= 1.0e-2, rel_total
    assert rel_sync <= 1.0e-2, rel_sync
    assert rel_ssc <= 1.0e-2, rel_ssc
    print("PASS: adaptive spectrum check succeeded.")


if __name__ == "__main__":
    main()
