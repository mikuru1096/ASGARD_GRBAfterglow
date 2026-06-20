from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.api_observe import run_fit
from asgard_core.asgard_config import RuntimeConfig, ReverseShockConfig
from asgard_core.asgard_setup import build_simulation_setup
from asgard_core.asgard_runtime import solve_dynamics
from asgard_core.asgard_state import make_query_setup, solve_state_from_setup


def _config(index_y: int, *, sigma: float | None = None) -> RuntimeConfig:
    reverse_kwargs = {}
    if sigma is not None:
        reverse_kwargs["sigma"] = sigma
    return RuntimeConfig(
        index_y=index_y,
        index_syn_integr=2,
        num_threads=1,
        num_gam_e=61,
        num_nu=41,
        num_r=40,
        num_theta=24,
        num_tobs=24,
        reverse=True,
        reverse_shock=ReverseShockConfig(
            enabled=True,
            delta_t_s=10.0,
            epsilon_e=0.1,
            epsilon_b=1.0e-2,
            p=2.4,
            f_e=1.0,
            **reverse_kwargs,
        ),
    )


def _run(index_y: int) -> None:
    config = _config(index_y)
    result = run_fit(config)
    assert np.all(np.isfinite(result.bands_flux))

    setup = make_query_setup(config, np.logspace(2.0, 5.0, 6), np.array([1.0e9, 1.0e14]))
    state = solve_state_from_setup(config, setup)
    rs = state.dynamics.reverse_shock
    assert rs is not None
    assert rs.causality is not None
    assert rs.causality.global_reverse_shock_allowed
    assert rs.causality.pressure_balance_condition_seen
    assert rs.causality.local_fast_condition_seen
    assert rs.causality.reverse_shock_started
    assert rs.causality.criteria_agree
    assert rs.causality.reference_crossing_radius_cm > 0.0
    assert rs.causality.pressure_balance_start_radius_cm > 0.0
    assert rs.causality.pressure_balance_start_ratio >= 1.0
    assert rs.causality.local_start_radius_cm > 0.0
    assert rs.causality.actual_start_pressure_ratio >= 1.0
    assert 0.0 < rs.causality.actual_start_contact_fraction <= 1.0
    assert np.all(np.isfinite(rs.magnetic_field_g))
    assert np.all(np.isfinite(rs.internal_energy_erg))
    assert np.all(np.isfinite(rs.comoving_volume_cm3))
    assert np.all(np.isfinite(rs.gamma34))
    assert np.isfinite(rs.ordered_magnetic_cross_g)
    assert np.all(rs.internal_energy_erg > 0.0)
    assert np.all(rs.comoving_volume_cm3 > 0.0)
    assert np.all(rs.magnetic_field_g >= 0.0)
    assert np.all(rs.gamma34 >= 1.0)
    post = np.asarray(state.dynamics.radius, dtype=float) >= float(rs.r_cross)
    if np.count_nonzero(post) > 2:
        thermal_energy = rs.internal_energy_erg[post] / rs.swept_mass_g[post]
        assert np.all(np.diff(thermal_energy) <= 0.0)


def _run_sigma_zero_baseline() -> None:
    baseline = run_fit(_config(0))
    sigma_zero = run_fit(_config(0, sigma=0.0))
    assert np.allclose(sigma_zero.bands_flux, baseline.bands_flux, rtol=0.0, atol=0.0)


def _run_magnetized_interface() -> None:
    config = _config(0, sigma=1.0e-3)
    setup = make_query_setup(config, np.logspace(2.0, 5.0, 6), np.array([1.0e9, 1.0e14]))
    state = solve_state_from_setup(config, setup)
    rs = state.dynamics.reverse_shock
    assert rs is not None
    assert rs.causality is not None
    assert rs.causality.global_reverse_shock_allowed
    assert isinstance(rs.causality.global_reverse_shock_allowed, bool)
    assert isinstance(rs.causality.pressure_balance_condition_seen, bool)
    assert isinstance(rs.causality.local_fast_condition_seen, bool)
    assert isinstance(rs.causality.reverse_shock_started, bool)
    assert np.isfinite(rs.causality.contact_radius_cm)
    assert rs.causality.reference_crossing_radius_cm > 0.0
    assert rs.causality.pressure_balance_condition_seen
    assert rs.causality.pressure_balance_start_radius_cm > 0.0
    assert rs.causality.pressure_balance_start_ratio >= 1.0
    assert np.isfinite(rs.ordered_magnetic_cross_g)
    assert rs.ordered_magnetic_cross_g > 0.0
    assert np.all(np.isfinite(rs.magnetic_field_g))
    assert np.all(np.isfinite(rs.internal_energy_erg))
    assert np.all(np.isfinite(rs.comoving_volume_cm3))
    assert np.all(np.isfinite(rs.gamma34))
    injected = np.asarray(rs.swept_mass_g, dtype=float) > float(rs.swept_mass_g[0]) * (1.0 + 1.0e-6)
    assert np.any(injected)
    assert np.all(rs.magnetic_field_g[injected] > 0.0)
    assert np.all(rs.internal_energy_erg > 0.0)
    assert np.all(rs.comoving_volume_cm3 > 0.0)


def _run_magnetized_delayed_reverse_shock() -> None:
    config = _config(0, sigma=1.0)
    setup = make_query_setup(config, np.logspace(2.0, 5.0, 6), np.array([1.0e9, 1.0e14]))
    state = solve_state_from_setup(config, setup)
    rs = state.dynamics.reverse_shock
    assert rs is not None
    assert rs.causality is not None
    assert rs.causality.global_reverse_shock_allowed
    assert rs.causality.pressure_balance_condition_seen
    assert rs.causality.pressure_balance_start_ratio >= 1.0
    assert rs.causality.local_fast_condition_seen
    assert rs.causality.reverse_shock_started
    assert rs.causality.criteria_agree
    assert rs.causality.actual_start_pressure_ratio >= 1.0
    assert 0.0 < rs.causality.actual_start_contact_fraction <= 1.0
    assert rs.causality.fast_wave_crossing_radius_cm > 0.0
    assert rs.t_cross > 0.0
    assert rs.ordered_magnetic_cross_g > 0.0
    assert np.max(rs.magnetic_field_g) > 0.0
    radius = np.asarray(state.dynamics.radius, dtype=float)
    gamma = np.asarray(state.dynamics.r_gamma, dtype=float)
    swept_reverse = np.asarray(rs.swept_mass_g, dtype=float)
    assert np.all(np.isfinite(radius))
    assert np.all(np.isfinite(gamma))
    assert np.all(np.isfinite(swept_reverse))
    assert np.all(np.diff(radius) > 0.0)
    assert np.all(gamma > 1.0)
    assert np.all(np.diff(swept_reverse) >= 0.0)
    assert np.max(swept_reverse) <= rs.m3_cross_g * (1.0 + 1.0e-8)


def _run_dense_radius_grid() -> None:
    config = _config(0)
    config.num_r = 220
    setup = build_simulation_setup(config)
    dynamics = solve_dynamics(setup.boundary, config)
    radius = np.asarray(dynamics.radius, dtype=float)
    assert np.all(np.diff(radius) > 0.0)
    assert radius[1] - radius[0] > np.spacing(radius[0])


def _run_requested_frequency_seed_bounds() -> None:
    config = _config(2)
    config.eta_0 = 300.0
    config.z = 0.1
    requested = np.array([1.0e8, 1.0e30])
    default_setup = build_simulation_setup(config)
    query_setup = make_query_setup(config, np.array([1.0e3, 1.0e5]), requested)

    beta0 = np.sqrt(1.0 - config.eta_0**-2)
    min_doppler_den = config.eta_0 * (1.0 - beta0)
    max_doppler_den = config.eta_0 * (1.0 + beta0)
    redshift_factor = 1.0 + config.z
    projected_min = query_setup.seed_frequency_hz[0] / (min_doppler_den * redshift_factor)
    projected_max = query_setup.seed_frequency_hz[-1] / (max_doppler_den * redshift_factor)

    assert query_setup.seed_frequency_hz[0] < default_setup.seed_frequency_hz[0]
    assert projected_min < requested[0]
    assert projected_max > requested[-1]


def main() -> None:
    for index_y in (0, 2):
        _run(index_y)
    _run_sigma_zero_baseline()
    _run_magnetized_interface()
    _run_magnetized_delayed_reverse_shock()
    _run_dense_radius_grid()
    _run_requested_frequency_seed_bounds()
    print("reverse-shock-smoke-ok")


if __name__ == "__main__":
    main()
