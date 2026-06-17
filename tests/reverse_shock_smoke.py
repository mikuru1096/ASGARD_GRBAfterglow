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
    config = _config(0, sigma=1.0e-2)
    setup = make_query_setup(config, np.logspace(2.0, 5.0, 6), np.array([1.0e9, 1.0e14]))
    state = solve_state_from_setup(config, setup)
    rs = state.dynamics.reverse_shock
    assert rs is not None
    assert np.isfinite(rs.ordered_magnetic_cross_g)
    assert rs.ordered_magnetic_cross_g > 0.0
    assert np.all(np.isfinite(rs.magnetic_field_g))
    assert np.all(np.isfinite(rs.internal_energy_erg))
    assert np.all(np.isfinite(rs.comoving_volume_cm3))
    assert np.all(np.isfinite(rs.gamma34))
    assert np.all(rs.magnetic_field_g > 0.0)
    assert np.all(rs.internal_energy_erg > 0.0)
    assert np.all(rs.comoving_volume_cm3 > 0.0)


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
    _run_dense_radius_grid()
    _run_requested_frequency_seed_bounds()
    print("reverse-shock-smoke-ok")


if __name__ == "__main__":
    main()
