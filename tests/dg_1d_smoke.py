from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.asgard_config import RuntimeConfig, ReverseShockConfig
from asgard_core.asgard_state import solve_state


def _config(solver: str) -> RuntimeConfig:
    return RuntimeConfig(
        electron_solver=solver,
        index_y=0,
        include_forward_ssc=False,
        num_threads=1,
        num_gam_e=121,
        num_nu=31,
        num_r=24,
        num_theta=12,
        num_tobs=8,
    )


def _jump_config(solver: str) -> RuntimeConfig:
    return RuntimeConfig(
        electron_solver=solver,
        index_y=0,
        include_forward_ssc=False,
        num_threads=1,
        num_gam_e=121,
        num_nu=31,
        num_r=32,
        num_theta=12,
        num_tobs=8,
        jump_r_cm=(3.0e16, 1.2e17),
        jump_factor=(8.0, 0.35),
        jump_width_log10=(0.025, 0.035),
    )


def _reverse_config(solver: str) -> RuntimeConfig:
    return RuntimeConfig(
        electron_solver=solver,
        index_y=0,
        include_forward_ssc=False,
        num_threads=1,
        num_gam_e=121,
        num_nu=25,
        num_r=32,
        num_theta=10,
        num_tobs=8,
        reverse=True,
        reverse_shock=ReverseShockConfig(
            enabled=True,
            delta_t_s=10.0,
            epsilon_e=0.1,
            epsilon_b=1.0e-2,
            p=2.4,
            f_e=1.0,
        ),
    )


def _electron_count_x(state) -> float:
    gam = np.asarray(state.electron.gam_e, dtype=float)
    dndg = np.asarray(state.electron.d_n_gam_e[:, -1], dtype=float)
    x = np.log10(gam)
    dnx = dndg * gam * np.log(10.0)
    return float(np.trapezoid(dnx, x))


def _reverse_count_x(state) -> float:
    rs = state.dynamics.reverse_shock
    assert rs is not None
    gam = np.asarray(rs.gam_e, dtype=float)
    dndg = np.asarray(rs.d_n_gam_e[:, -1], dtype=float)
    x = np.log10(gam)
    dnx = dndg * gam * np.log(10.0)
    return float(np.trapezoid(dnx, x))


def _assert_no_grid_oscillation(gamma: np.ndarray, dndg: np.ndarray) -> None:
    gamma = np.asarray(gamma, dtype=float)
    dndg = np.asarray(dndg, dtype=float)
    tail = gamma * gamma * dndg
    finite_positive = np.isfinite(tail) & (tail > 0.0)
    if not np.any(finite_positive):
        raise AssertionError("empty positive electron spectrum")
    threshold = float(np.nanmax(tail[finite_positive])) * 1.0e-9
    active = finite_positive & (tail >= threshold)
    if int(np.count_nonzero(active)) < 9:
        return
    active_idx = np.flatnonzero(active)
    y = np.log10(tail[active_idx])
    if np.any(np.diff(active_idx) > 1):
        raise AssertionError("DG electron spectrum has a zero gap inside the active support")
    slope_sign = np.sign(np.diff(y))
    turning_points = np.count_nonzero(slope_sign[1:] * slope_sign[:-1] < 0.0)
    if turning_points > 2:
        raise AssertionError(f"DG electron spectrum has grid-scale sawtooth turns: {turning_points}")
    curvature = np.abs(y[2:] - 2.0 * y[1:-1] + y[:-2])
    if curvature.size and float(np.max(curvature)) > 1.2:
        raise AssertionError(f"DG electron spectrum has grid-scale ringing: {float(np.max(curvature)):.3f} dex")
    hole = (~active[1:-1]) & active[:-2] & active[2:]
    if np.any(hole):
        raise AssertionError("DG electron spectrum has an isolated zero hole")


def main() -> None:
    times = np.array([1.0e4, 1.0e5], dtype=float)
    reference = solve_state(_config("fullhide_1d"), times)
    dg_state = solve_state(_config("dg_1d"), times)

    report = dg_state.adapter_reports["electron"]
    assert report.solver == "dg_1d"
    assert report.grid_semantics == "log-gamma-1d-dg"
    assert report.status == "ok"

    electron = dg_state.electron
    assert electron.gam_e.shape == (121,)
    assert electron.d_n_gam_e.shape == (121, 24)
    assert electron.l_syn_spec.shape == (31, 24)
    assert electron.seed_syn.shape == (31, 24)
    assert np.all(np.isfinite(electron.d_n_gam_e))
    assert np.all(np.isfinite(electron.l_syn_spec))
    assert np.min(electron.d_n_gam_e) >= 0.0
    assert np.max(electron.l_syn_spec) > 0.0
    _assert_no_grid_oscillation(electron.gam_e, electron.d_n_gam_e[:, -1])

    count_ratio = _electron_count_x(dg_state) / _electron_count_x(reference)
    assert 0.8 < count_ratio < 1.3

    luminosity_ratio = np.max(electron.l_syn_spec) / np.max(reference.electron.l_syn_spec)
    assert 0.2 < float(luminosity_ratio) < 5.0

    jump_reference = solve_state(_jump_config("fullhide_1d"), times)
    jump_dg_state = solve_state(_jump_config("dg_1d"), times)
    assert np.all(np.isfinite(jump_dg_state.electron.d_n_gam_e))
    assert np.all(np.isfinite(jump_dg_state.electron.l_syn_spec))
    assert np.min(jump_dg_state.electron.d_n_gam_e) >= 0.0
    _assert_no_grid_oscillation(jump_dg_state.electron.gam_e, jump_dg_state.electron.d_n_gam_e[:, -1])

    jump_count_ratio = _electron_count_x(jump_dg_state) / _electron_count_x(jump_reference)
    assert 0.9 < jump_count_ratio < 1.1

    jump_luminosity_ratio = np.max(jump_dg_state.electron.l_syn_spec) / np.max(jump_reference.electron.l_syn_spec)
    assert 0.7 < float(jump_luminosity_ratio) < 1.35

    reverse_reference = solve_state(_reverse_config("fullhide_1d"), np.array([1.0e2, 1.0e4], dtype=float))
    reverse_dg = solve_state(_reverse_config("dg_1d"), np.array([1.0e2, 1.0e4], dtype=float))
    reverse_rs = reverse_dg.dynamics.reverse_shock
    assert reverse_rs is not None
    assert np.all(np.isfinite(reverse_rs.d_n_gam_e))
    assert np.min(reverse_rs.d_n_gam_e) >= 0.0
    assert np.max(reverse_rs.d_n_gam_e) > 0.0
    _assert_no_grid_oscillation(reverse_rs.gam_e, reverse_rs.d_n_gam_e[:, -1])
    reverse_count_ratio = _reverse_count_x(reverse_dg) / _reverse_count_x(reverse_reference)
    assert 0.5 < float(reverse_count_ratio) < 2.0


if __name__ == "__main__":
    main()
