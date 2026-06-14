from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD.api_model import _make_details
from asgard_core.asgard_config import FitConfig, ReverseShockConfig
from asgard_core.asgard_numpy import trapezoid
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


def _run_zero_amplitude_bump_baseline() -> None:
    base = _base_config()
    explicit_zero = _base_config()
    explicit_zero.jump_r_cm = (3.0e16,)
    explicit_zero.jump_factor = (1.0,)
    explicit_zero.jump_width_log10 = (0.10,)
    times = np.logspace(2.0, 5.0, 10)
    freqs = np.array([1.0e10, 1.0e14], dtype=float)
    base_state = solve_state_from_setup(base, make_query_setup(base, times, freqs))
    zero_state = solve_state_from_setup(explicit_zero, make_query_setup(explicit_zero, times, freqs))
    assert zero_state.reverse_emission is not None
    assert zero_state.reverse_emission.secondary_rs is None
    assert zero_state.dynamics.reverse_shock is not None
    assert zero_state.dynamics.reverse_shock.secondary_branch_swept_mass_g is not None
    assert np.all(zero_state.dynamics.reverse_shock.secondary_branch_swept_mass_g == 0.0)
    assert np.all(zero_state.dynamics.reverse_shock.secondary_branch_internal_energy_erg == 0.0)
    assert np.all(zero_state.dynamics.reverse_shock.secondary_branch_comoving_volume_cm3 == 0.0)
    assert np.allclose(zero_state.dynamics.radius, base_state.dynamics.radius, rtol=0.0, atol=0.0)
    assert np.allclose(zero_state.dynamics.r_gamma, base_state.dynamics.r_gamma, rtol=0.0, atol=0.0)
    assert np.allclose(zero_state.components.total, base_state.components.total, rtol=0.0, atol=0.0)


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
    assert state.dynamics.reverse_shock is not None
    dyn_branch_m = state.dynamics.reverse_shock.secondary_branch_swept_mass_g
    dyn_branch_u = state.dynamics.reverse_shock.secondary_branch_internal_energy_erg
    dyn_branch_v = state.dynamics.reverse_shock.secondary_branch_comoving_volume_cm3
    assert dyn_branch_m is not None and dyn_branch_u is not None and dyn_branch_v is not None
    assert dyn_branch_m.shape == (len(config.jump_r_cm), state.dynamics.radius.size)
    assert np.all(np.isfinite(dyn_branch_m))
    assert np.all(np.isfinite(dyn_branch_u))
    assert np.all(np.isfinite(dyn_branch_v))
    assert np.any(dyn_branch_m > 0.0)
    assert np.any(dyn_branch_u > 0.0)
    assert np.any(dyn_branch_v > 0.0)
    assert np.all(np.diff(dyn_branch_m, axis=1) >= -1.0e-20 * np.maximum(1.0, dyn_branch_m[:, 1:]))
    details = _make_details(state.components, patches=[{"phi": 0.0, "theta": 0.0, "weight": 1.0}], state=state)
    assert details.rev is not None
    assert details.rev.secondary_rs_gamma_e is secondary.gam_e
    assert details.rev.secondary_rs_dN_dgamma_e is secondary.d_n_gam_e
    assert details.rev.secondary_rs_branch_swept_mass_g is secondary.branch_swept_mass_g
    assert details.rev.secondary_rs_branch_internal_energy_erg is secondary.branch_internal_energy_erg
    assert np.all(np.isfinite(state.reverse_emission.l_syn_spec))
    assert np.any(secondary.luminosity_syn > 0.0)
    _assert_secondary_injection_inside_candidate_windows(state.dynamics.radius, secondary.dissipated_energy_density, config)
    _assert_secondary_event_diagnostics(secondary, config)
    for values in (
        secondary.gamma_contact,
        secondary.pressure_3,
        secondary.gamma_43,
        secondary.beta_rs,
        secondary.dissipated_energy_density,
        secondary.swept_mass_g,
        secondary.internal_energy_erg,
        secondary.comoving_volume_cm3,
        secondary.pressure_total,
        secondary.enthalpy_density_total,
        secondary.branch_swept_mass_g,
        secondary.branch_internal_energy_erg,
        secondary.branch_comoving_volume_cm3,
        secondary.branch_magnetic_field_g,
        secondary.magnetic_field_g,
        secondary.gam_e,
        secondary.d_n_gam_e,
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
    beta_contact = np.sqrt(1.0 - secondary.gamma_contact[injection] ** -2)
    beta_upstream = np.sqrt(1.0 - state.dynamics.r_gamma[injection] ** -2)
    assert np.all(secondary.beta_rs[injection] > beta_contact)
    assert np.all(secondary.beta_rs[injection] < beta_upstream)
    assert np.all(secondary.dissipated_energy_density[injection] > 0.0)
    assert np.all(np.diff(secondary.swept_mass_g) >= 0.0)
    assert secondary.branch_swept_mass_g.shape == (len(config.jump_r_cm), secondary.swept_mass_g.size)
    assert np.allclose(np.sum(secondary.branch_swept_mass_g, axis=0), secondary.swept_mass_g, rtol=1.0e-13, atol=0.0)
    assert np.allclose(
        np.sum(secondary.branch_internal_energy_erg, axis=0),
        secondary.internal_energy_erg,
        rtol=1.0e-13,
        atol=0.0,
    )
    assert np.allclose(
        np.sum(secondary.branch_comoving_volume_cm3, axis=0),
        secondary.comoving_volume_cm3,
        rtol=1.0e-13,
        atol=0.0,
    )
    last_injection = int(np.flatnonzero(injection)[-1])
    if last_injection + 1 < secondary.magnetic_field_g.shape[0]:
        assert secondary.swept_mass_g[last_injection + 1] > 0.0
        assert secondary.internal_energy_erg[last_injection + 1] > 0.0
        assert secondary.comoving_volume_cm3[last_injection + 1] > 0.0
        assert secondary.pressure_total[last_injection + 1] > 0.0
        assert secondary.enthalpy_density_total[last_injection + 1] > 0.0
        assert secondary.magnetic_field_g[last_injection + 1] > 0.0
        assert float(trapezoid(secondary.d_n_gam_e[:, last_injection + 1], secondary.gam_e)) > 0.0
        assert np.any(secondary.luminosity_syn[:, last_injection + 1] > 0.0)
    assert np.all(secondary.pressure_total[reservoir] > 0.0)
    assert np.all(secondary.enthalpy_density_total[reservoir] > secondary.pressure_total[reservoir])
    assert np.allclose(
        secondary.pressure_total[reservoir],
        secondary.internal_energy_erg[reservoir] / (3.0 * secondary.comoving_volume_cm3[reservoir]),
        rtol=1.0e-13,
        atol=0.0,
    )
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


def _run_single_bump_event_grid_convergence() -> None:
    starts = []
    ends = []
    for num_r in (72, 96, 128):
        config = _base_config()
        config.num_r = num_r
        config.num_gam_e = 61
        config.num_nu = 32
        config.num_theta = 12
        config.num_tobs = 48
        config.jump_r_cm = (3.0e16,)
        config.jump_factor = (6.0,)
        config.jump_width_log10 = (0.10,)
        state = solve_state_from_setup(
            config,
            make_query_setup(config, np.logspace(2.0, 6.0, 10), np.array([1.0e10])),
        )
        assert state.reverse_emission is not None
        assert state.reverse_emission.secondary_rs is not None
        starts.append(float(state.reverse_emission.secondary_rs.start_radius_cm[0]))
        ends.append(float(state.reverse_emission.secondary_rs.shock_end_radius_cm[0]))
    assert np.allclose(starts, starts[-1], rtol=2.0e-3, atol=0.0)
    assert np.allclose(ends, ends[-1], rtol=2.0e-3, atol=0.0)


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
    _run_zero_amplitude_bump_baseline()
    _run_disabled_branch_rejections()
    _run_single_bump_secondary()
    _run_single_bump_event_grid_convergence()
    _run_multi_bump_reverse()
    print("multi-density-reverse-shock-smoke-ok")


if __name__ == "__main__":
    main()
