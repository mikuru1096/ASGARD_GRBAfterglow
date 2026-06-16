from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core import Model, Observer, top_hat_jet, UniformMedium
from asgard_core.api_model import _build_fit_config_for_patch, _solve_patch_state
from tests.public_api_builders import hadronic, numerics, radiation, solver_options, top_hat_model


def _joint_bh_model() -> Model:
    return top_hat_model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=120.0,
            opening_angle_rad=0.1,
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        medium=UniformMedium(number_density_cm3=1.0),
        observer=Observer(z=0.05, viewing_angle_rad=0.0, viewing_azimuth_rad=0.0, luminosity_distance_cm=1.0e26),
        fwd_rad=radiation(
            proton_energy_fraction=0.1,
            include_ssc=True,
            include_kn_correction=True,
            proton_synch=True,
            bethe_heitler=True,
            pgamma_scheme="hummer_2010_response",
        ),
        numerics=numerics(
            num_electron_gamma=20,
            num_photon_frequency=24,
            num_radius=18,
            num_theta=10,
            num_observer_time=18,
        ),
        solver_options=solver_options(
            electron_solver="fullhide_1d",
            electron_photon_coupling="joint",
            ssc_cooling_mode="numeric_ic_kn",
        ),
        hadronic=hadronic(
            enabled=True,
            solver="am3_1d",
            num_proton_gamma=24,
            num_neutrino_frequency=16,
            pgamma_scheme="hummer_2010_response",
        ),
    )


def _assert_shell_profile_smooth(name: str, values: np.ndarray, limit: float) -> None:
    profile = np.asarray(values, dtype=float)
    positive = profile[np.isfinite(profile) & (profile > 0.0)]
    assert positive.size >= 3, name
    adjacent_ratio = np.maximum(positive[1:] / positive[:-1], positive[:-1] / positive[1:])
    assert float(np.max(adjacent_ratio)) < limit, name


def main() -> None:
    model = _joint_bh_model()
    times = np.array([1.0e4, 1.0e5], dtype=float)
    frequencies = np.array([1.0e14, 1.0e18], dtype=float)
    result = model.flux_density(times, frequencies)
    assert np.all(np.isfinite(result.total))
    details = model.details(1.0e4, 1.0e14)
    assert details.fwd.dN_dgamma_e_bh is not None
    assert details.fwd.tau_bh is not None
    assert np.all(np.isfinite(details.fwd.dN_dgamma_e_bh))
    assert np.all(np.isfinite(details.fwd.tau_bh))
    assert float(np.max(details.fwd.tau_bh)) >= 0.0

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
    assert state.hadronic.tau_bh is not None
    np.testing.assert_allclose(state.observer.tau_extra, state.hadronic.tau_bh, rtol=1.0e-14, atol=0.0)
    assert np.all(np.isfinite(state.photon_field.hadronic_target_seed))
    synch_part_after_sink = state.photon_field.hadronic_target_seed - state.photon_field.hadronic_forward_ssc_seed
    positive_excess = synch_part_after_sink - state.electron.seed_syn
    subtraction_scale = np.maximum.reduce(
        [
            state.photon_field.hadronic_target_seed,
            state.photon_field.hadronic_forward_ssc_seed,
            state.electron.seed_syn,
            np.full_like(state.electron.seed_syn, 1.0e-300),
        ]
    )
    resolved_excess = positive_excess > 1.0e-12 * subtraction_scale
    assert not np.any(resolved_excess)
    sink_mask = (state.electron.seed_syn > 0.0) & (state.hadronic.tau_bh > 0.0)
    if np.any(sink_mask) and float(np.max(state.hadronic.tau_bh[sink_mask])) > 1.0e-10:
        assert float(np.max(state.electron.seed_syn[sink_mask] - synch_part_after_sink[sink_mask])) > 0.0
    _assert_shell_profile_smooth(
        "joint electron spectrum",
        np.trapezoid(state.electron.d_n_gam_e, state.electron.gam_e, axis=0),
        200.0,
    )
    _assert_shell_profile_smooth(
        "joint BH pair spectrum",
        np.trapezoid(state.hadronic.d_n_gam_e_bh, state.electron.gam_e, axis=0),
        200.0,
    )
    _assert_shell_profile_smooth(
        "joint photon target spectrum",
        np.trapezoid(state.photon_field.hadronic_target_seed, state.photon_field.seed_frequency_hz, axis=0),
        200.0,
    )
    print("electron_photon_joint_smoke: ok")


if __name__ == "__main__":
    main()
