from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core import Model, Observer, top_hat_jet
from tests.public_api_builders import numerics, observer_grid, radiation, reverse_shock, solver_options, top_hat_model


def _model(theta_obs: float = 0.0, gamma0: float = 80.0) -> Model:
    model = top_hat_model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=gamma0,
            opening_angle_rad=0.1,
            shell_duration_s=20.0,
            magnetar=None,
            spreading=False,
        ),
        observer=Observer(z=0.1, viewing_angle_rad=theta_obs, viewing_azimuth_rad=0.0, luminosity_distance_cm=1.0e26),
        fwd_rad=radiation(epsilon_e=0.1, epsilon_B=1.0e-4, p=2.3, accelerated_electron_fraction=0.2),
        rvs_rad=radiation(epsilon_e=0.15, epsilon_B=3.0e-3, p=2.35, accelerated_electron_fraction=1.0),
        numerics=numerics(
            num_radius=48,
            num_theta=8,
            num_observer_time=24,
            num_electron_gamma=41,
            num_photon_frequency=31,
        ),
        observer_grid=observer_grid(time_min_s=1.0e2, time_max_s=1.0e5),
        solver_options=solver_options(electron_solver="fullhide_1d", ssc_cooling_mode="none"),
        reverse_shock=reverse_shock(enabled=True, shell_duration_s=20.0),
    )
    model.setups.jump_r_cm = (3.0e16,)
    model.setups.jump_factor = (6.0,)
    model.setups.jump_width_log10 = (0.10,)
    return model


def _assert_offaxis_adaptive_eats_resolution() -> None:
    model = _model(theta_obs=0.05, gamma0=20.0)
    model.setups.num_theta = 4
    model.setups.num_phi = 1
    model._ensure_direct_adaptive_eats_resolution()
    gamma = model.jet.initial_lorentz_factor
    theta_j = model.jet.opening_angle_rad
    expected_theta = int(np.ceil(theta_j * model.setups.patch_sampling_beaming_resolution * gamma / model.setups.patch_sampling_beaming_factor)) + 1
    expected_phi = int(np.ceil(np.pi * model.setups.patch_sampling_beaming_resolution * gamma * np.sin(theta_j) / model.setups.patch_sampling_beaming_factor)) + 1
    assert model.setups.num_theta >= expected_theta
    assert model.setups.num_phi >= max(expected_phi, 2)


def main() -> None:
    times = np.logspace(2.0, 5.0, 8)
    freqs = np.array([1.0e10], dtype=float)
    model = _model()
    baseline = model.flux_density_grid(times, freqs)
    adaptive = model.flux_density_grid_adaptive(times, freqs)

    assert baseline.total.shape == (freqs.size, times.size)
    assert adaptive.total.shape[0] == freqs.size
    assert adaptive.total.shape[1] == adaptive.time_s.size
    assert adaptive.time_s.size > times.size
    assert np.all(np.diff(adaptive.time_s) > 0.0)
    assert np.isclose(adaptive.time_s[0], times[0])
    assert np.isclose(adaptive.time_s[-1], times[-1])
    assert np.all(np.isfinite(adaptive.total))
    for time_s in times:
        assert np.any(np.isclose(adaptive.time_s, time_s, rtol=1.0e-13, atol=0.0))
    assert np.any(adaptive.rev.sync > 0.0)
    details = model.details(times[0], times[-1])
    assert details.rev is not None
    assert details.rev.secondary_rs_event_active is not None
    assert np.any(details.rev.secondary_rs_event_active)
    assert details.rev.secondary_rs_start_radius is not None
    assert details.rev.secondary_rs_shock_end_radius is not None
    assert details.rev.secondary_rs_B is not None
    active = details.rev.secondary_rs_event_active
    assert np.all(details.rev.secondary_rs_shock_end_radius[active] >= details.rev.secondary_rs_start_radius[active])
    latest_end = float(np.max(details.rev.secondary_rs_shock_end_radius[active]))
    cooling_tail = (details.rev.radius > latest_end) & (details.rev.secondary_rs_B > 0.0)
    assert np.any(cooling_tail)
    tail_axis_times = details.rev.t_obs[cooling_tail]
    covered_tail = tail_axis_times[(tail_axis_times >= adaptive.time_s[0]) & (tail_axis_times <= adaptive.time_s[-1])]
    assert covered_tail.size > 0
    for time_s in covered_tail[: min(3, covered_tail.size)]:
        assert np.any(np.isclose(adaptive.time_s, time_s, rtol=1.0e-13, atol=0.0))
    _assert_offaxis_adaptive_eats_resolution()
    print("adaptive-observer-time-grid-smoke-ok")


if __name__ == "__main__":
    main()
