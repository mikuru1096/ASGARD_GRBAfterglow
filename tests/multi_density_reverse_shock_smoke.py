from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.asgard_config import FitConfig, ReverseShockConfig
from asgard_core.asgard_physics_utils import ambient_density
from asgard_core.asgard_state import make_query_setup, solve_state_from_setup


def _base_config() -> FitConfig:
    return FitConfig(
        index_y=0,
        index_syn_integr=2,
        num_threads=1,
        num_gam_e=81,
        num_nu=48,
        num_r=96,
        num_theta=16,
        num_tobs=64,
        t_obs_min_log10=2.0,
        t_obs_max_log10=6.0,
        eta_0=120.0,
        e_iso=1.0e53,
        d_ne=1.0,
        a_star=-1.0,
        electron_solver="fullhide_1d",
        reverse=True,
        include_forward_ssc=False,
        reverse_shock=ReverseShockConfig(
            enabled=True,
            delta_t_s=30.0,
            epsilon_e=0.12,
            epsilon_b=3.0e-3,
            p=2.35,
            f_e=1.0,
            include_ssc=False,
            include_cross_zone_ic=False,
        ),
    )


def _run_single_bump_equivalence() -> None:
    radius = np.logspace(15.0, 18.0, 64)
    scalar = _base_config()
    scalar.r_tr = 3.0e16
    scalar.f_jump = 6.0
    scalar.f_wide = 0.08
    arrays = _base_config()
    arrays.jump_r_cm = (scalar.r_tr,)
    arrays.jump_factor = (scalar.f_jump,)
    arrays.jump_width_log10 = (scalar.f_wide,)
    assert np.allclose(ambient_density(radius, scalar), ambient_density(radius, arrays), rtol=0.0, atol=0.0)


def _run_multi_bump_reverse() -> None:
    config = _base_config()
    config.jump_r_cm = (2.5e16, 7.0e16)
    config.jump_factor = (5.0, 8.0)
    config.jump_width_log10 = (0.10, 0.12)
    times = np.logspace(2.0, 6.0, 18)
    freqs = np.array([1.0e10, 1.0e14], dtype=float)
    setup = make_query_setup(config, times, freqs)
    state = solve_state_from_setup(config, setup)
    assert state.reverse_emission is not None
    assert state.reverse_emission.secondary_rs is not None
    secondary = state.reverse_emission.secondary_rs
    assert np.all(np.isfinite(state.reverse_emission.l_syn_spec))
    assert np.any(secondary.luminosity_syn > 0.0)
    _assert_secondary_injection_inside_candidate_windows(state.dynamics.radius, secondary.dissipated_energy_density, config)
    _assert_secondary_event_diagnostics(secondary, config)
    for values in (
        secondary.gamma_contact,
        secondary.pressure_3,
        secondary.gamma_43,
        secondary.dissipated_energy_density,
        secondary.magnetic_field_g,
        secondary.nu_m,
        secondary.nu_c,
        secondary.nu_a,
    ):
        assert np.all(np.isfinite(values))
    reservoir = secondary.magnetic_field_g > 0.0
    injection = secondary.dissipated_energy_density > 0.0
    assert np.any(reservoir)
    assert np.all(secondary.gamma_contact[injection] > 1.0)
    assert np.all(secondary.gamma_43[injection] >= 1.0)
    assert np.all(secondary.pressure_3[injection] > 0.0)
    assert np.all(secondary.dissipated_energy_density[injection] > 0.0)
    assert np.isclose(
        np.sum(secondary.electron_injected_energy_erg),
        config.reverse_shock.epsilon_e * np.sum(secondary.dissipated_energy_erg),
        rtol=1.0e-13,
        atol=0.0,
    )
    _assert_log_smooth(secondary.gamma_contact[injection], 0.5)
    _assert_log_smooth(secondary.pressure_3[injection], 3.0)
    _assert_log_smooth(secondary.gamma_43[injection], 0.5)
    _assert_log_smooth(secondary.magnetic_field_g[reservoir], 3.0)
    assert state.components.rev_sync is not None
    _assert_log_smooth(state.components.fwd_sync[0], 2.0)
    _assert_log_smooth(state.components.rev_sync[0], 5.0)


def _run_single_bump_secondary() -> None:
    config = _base_config()
    config.jump_r_cm = (3.0e16,)
    config.jump_factor = (6.0,)
    config.jump_width_log10 = (0.10,)
    times = np.logspace(2.0, 6.0, 14)
    freqs = np.array([1.0e10], dtype=float)
    state = solve_state_from_setup(config, make_query_setup(config, times, freqs))
    assert state.reverse_emission is not None
    assert state.reverse_emission.secondary_rs is not None
    assert np.any(state.reverse_emission.secondary_rs.luminosity_syn > 0.0)
    _assert_secondary_event_diagnostics(state.reverse_emission.secondary_rs, config)
    _assert_secondary_injection_inside_candidate_windows(
        state.dynamics.radius,
        state.reverse_emission.secondary_rs.dissipated_energy_density,
        config,
    )


def _run_disabled_branch_rejections() -> None:
    cases = [
        ("non-synch electron cooling", lambda cfg: setattr(cfg, "index_y", 1)),
        ("forward SSC", lambda cfg: setattr(cfg, "include_forward_ssc", True)),
        ("wind medium", lambda cfg: setattr(cfg, "a_star", 0.1)),
        ("2D electron solver", lambda cfg: setattr(cfg, "electron_solver", "fullhide_2d")),
        ("structured observer path", lambda cfg: setattr(cfg, "geometry_kernel", "chi_eats_2d")),
        ("forward hadronic", lambda cfg: setattr(cfg.hadronic, "enabled", True)),
        ("reverse hadronic", lambda cfg: setattr(cfg.hadronic, "reverse_enabled", True)),
        ("RS SSC", lambda cfg: setattr(cfg.reverse_shock, "include_ssc", True)),
        ("cross-zone IC", lambda cfg: setattr(cfg.reverse_shock, "include_cross_zone_ic", True)),
    ]
    for label, mutate in cases:
        config = _base_config()
        config.jump_r_cm = (3.0e16,)
        config.jump_factor = (6.0,)
        config.jump_width_log10 = (0.10,)
        mutate(config)
        try:
            solve_state_from_setup(config, make_query_setup(config, np.logspace(2.0, 3.0, 3), np.array([1.0e10])))
        except (NotImplementedError, ValueError):
            continue
        raise AssertionError(f"multi-density reverse shock accepted {label}")


def _assert_log_smooth(values: np.ndarray, max_jump: float) -> None:
    positive = np.asarray(values, dtype=float)
    positive = positive[np.isfinite(positive) & (positive > 0.0)]
    if positive.size < 3:
        return
    jumps = np.abs(np.diff(np.log10(positive)))
    assert np.max(jumps) < max_jump


def _assert_secondary_injection_inside_candidate_windows(radius: np.ndarray, injection: np.ndarray, config: FitConfig) -> None:
    active = np.asarray(injection, dtype=float) > 0.0
    if not np.any(active):
        return
    allowed = np.zeros_like(active, dtype=bool)
    for radius_j, width_j in zip(config.jump_r_cm, config.jump_width_log10):
        width_cm = float(width_j) * float(radius_j)
        allowed |= (radius >= float(radius_j) - 4.0 * width_cm) & (radius < float(radius_j))
    assert not np.any(active & ~allowed)


def _assert_secondary_event_diagnostics(secondary, config: FitConfig) -> None:
    assert secondary.event_active.shape == (len(config.jump_r_cm),)
    assert secondary.start_radius_cm.shape == secondary.event_active.shape
    assert secondary.shock_end_radius_cm.shape == secondary.event_active.shape
    assert np.any(secondary.event_active)
    active = secondary.event_active
    assert np.all(secondary.shock_end_radius_cm[active] >= secondary.start_radius_cm[active])
    assert np.all(secondary.shock_end_tobs_axis_s[active] >= secondary.start_tobs_axis_s[active])
    for i_jump, is_active in enumerate(active):
        if not is_active:
            continue
        radius_j = float(config.jump_r_cm[i_jump])
        width_cm = float(config.jump_width_log10[i_jump]) * radius_j
        assert radius_j - 4.0 * width_cm <= secondary.start_radius_cm[i_jump] < radius_j
        assert radius_j - 4.0 * width_cm <= secondary.shock_end_radius_cm[i_jump] < radius_j


def main() -> None:
    _run_single_bump_equivalence()
    _run_disabled_branch_rejections()
    _run_single_bump_secondary()
    _run_multi_bump_reverse()
    print("multi-density-reverse-shock-smoke-ok")


if __name__ == "__main__":
    main()
