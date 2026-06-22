from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core import Model, Observer, UniformMedium, top_hat_jet
from asgard_core.api_model import _build_fit_config_for_patch, _solve_patch_state
from tests.public_api_builders import hadronic, numerics, radiation, solver_options, top_hat_model


def _joint_secondary_model() -> Model:
    return top_hat_model(
        jet=top_hat_jet(
            energy_iso_erg=3.0e52,
            initial_lorentz_factor=140.0,
            opening_angle_rad=0.1,
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        medium=UniformMedium(number_density_cm3=3.0),
        observer=Observer(z=0.05, viewing_angle_rad=0.0, viewing_azimuth_rad=0.0, luminosity_distance_cm=1.0e26),
        fwd_rad=radiation(
            epsilon_e=0.08,
            epsilon_B=3.0e-3,
            p=2.25,
            proton_energy_fraction=0.25,
            include_ssc=True,
            include_kn_correction=True,
            proton_synch=True,
            include_pgamma=True,
            pp=True,
            bethe_heitler=True,
            pair_production=True,
            pgamma_scheme="hummer_2010_response",
        ),
        numerics=numerics(
            num_electron_gamma=18,
            num_photon_frequency=20,
            num_radius=14,
            eats_num_theta=8,
            num_observer_time=14,
        ),
        solver_options=solver_options(
            electron_solver="fullhide_1d",
            electron_photon_coupling="joint",
            ssc_cooling_mode="numeric_ic_kn",
        ),
        hadronic=hadronic(
            enabled=True,
            solver="am3_1d",
            num_proton_gamma=22,
            num_neutrino_frequency=12,
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
