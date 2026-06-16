from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_setup import build_simulation_setup
from asgard_core.asgard_runtime import solve_dynamics


def _event_config(num_r: int) -> RuntimeConfig:
    return RuntimeConfig(
        num_r=num_r,
        num_tobs=num_r,
        t_obs_min_log10=1.0,
        t_obs_max_log10=6.0,
        eta_0=180.0,
        e_iso=3.0e52,
        d_ne=0.8,
        z=0.2,
        initial_radius_cm=1.0e14,
        e_inj_t1=80.0,
        e_inj_t2=1.5e3,
        l_inj_0=2.0e47,
        q_inj=-0.4,
        density_profile_radius_cm=(1.0e15, 3.0e15, 9.0e15, 2.7e16),
        density_profile_n_cm3=(0.8, 12.0, 1.6, 4.0),
    )


def _solve(config: RuntimeConfig):
    setup = build_simulation_setup(config)
    return solve_dynamics(setup.boundary, config)


def _assert_event_is_inside_step(axis: np.ndarray, event_value: float) -> None:
    inside = (axis[:-1] < event_value) & (event_value < axis[1:])
    assert np.any(inside)


def main() -> None:
    coarse = _solve(_event_config(36))
    z = _event_config(36).z
    source_time = coarse.r_tobs / (1.0 + z)

    for event_time in (80.0 / (1.0 + z), 1.5e3 / (1.0 + z)):
        _assert_event_is_inside_step(source_time, event_time)
    for event_radius in (3.0e15, 9.0e15):
        _assert_event_is_inside_step(coarse.radius, event_radius)

    assert np.all(np.isfinite(coarse.r_gamma))
    assert np.all(np.isfinite(coarse.radius))
    assert np.all(np.isfinite(coarse.swept_mass_g))
    assert np.all(np.diff(coarse.radius) > 0.0)
    assert np.all(np.diff(coarse.swept_mass_g) >= 0.0)
    assert np.all(coarse.r_gamma > 1.0)
    assert coarse.swept_mass_g[-1] > coarse.swept_mass_g[0]
    print("dynamics-event-split-smoke-ok")


if __name__ == "__main__":
    main()
