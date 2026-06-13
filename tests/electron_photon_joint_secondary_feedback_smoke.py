from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet
from ASGARD.api_model import _build_fit_config_for_patch, _solve_patch_state


def _joint_secondary_model() -> Model:
    return Model(
        TophatJet(0.1, 3.0e52, 140.0),
        ISM(3.0),
        Observer(1.0e26, 0.05, 0.0),
        Radiation(
            0.08,
            3.0e-3,
            2.25,
            epsilon_p=0.25,
            ssc=True,
            kn=True,
            proton_synch=True,
            pg=True,
            pp=True,
            bethe_heitler=True,
            pair_production=True,
            pgamma_scheme="hummer_2010_response",
        ),
        setups=Setups(
            electron_solver="fullhide_1d",
            electron_photon_coupling="joint",
            num_gam_e=18,
            num_gam_p=22,
            num_nu=20,
            num_nu_nu=12,
            num_r=14,
            num_theta=8,
            num_tobs=14,
            hadronic_enabled=True,
            hadronic_solver="am3_1d",
            pgamma_scheme="hummer_2010_response",
            pair_cascade_iterations=2,
        ),
    )


def main() -> None:
    model = _joint_secondary_model()
    times = np.array([3.0e4, 3.0e5], dtype=float)
    frequencies = np.array([1.0e14, 1.0e18], dtype=float)
    config = _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )
    state = _solve_patch_state(model, config, times, frequencies)
    assert state.adapter_reports["electron"].grid_semantics == "log-gamma-1d-joint-cooling"
    assert state.hadronic is not None
    assert state.hadronic.secondary_electron_source_r is not None
    assert np.all(np.isfinite(state.hadronic.secondary_electron_source_r))
    assert float(np.max(state.hadronic.secondary_electron_source_r)) > 0.0
    assert state.hadronic.l_had_pair_production is not None
    assert np.all(np.isfinite(state.hadronic.l_had_pair_production))
    assert float(np.max(state.hadronic.l_had_pair_production)) >= 0.0
    assert np.all(np.isfinite(state.observer.tau_pair))
    assert np.all(state.observer.tau_pair >= 0.0)
    assert np.all(np.isfinite(state.photon_field.hadronic_target_seed))
    assert np.all(np.isfinite(state.electron.d_n_gam_e))
    print("electron_photon_joint_secondary_feedback_smoke: ok")


if __name__ == "__main__":
    main()
