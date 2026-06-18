from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_state import solve_state


def _config(solver: str) -> RuntimeConfig:
    return RuntimeConfig(
        electron_solver=solver,
        index_y=0,
        include_forward_ssc=False,
        num_threads=1,
        num_gam_e=41,
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
        num_gam_e=41,
        num_nu=31,
        num_r=32,
        num_theta=12,
        num_tobs=8,
        jump_r_cm=(3.0e16, 1.2e17),
        jump_factor=(8.0, 0.35),
        jump_width_log10=(0.025, 0.035),
    )


def _electron_count_x(state) -> float:
    gam = np.asarray(state.electron.gam_e, dtype=float)
    dndg = np.asarray(state.electron.d_n_gam_e[:, -1], dtype=float)
    x = np.log10(gam)
    dnx = dndg * gam * np.log(10.0)
    return float(np.trapezoid(dnx, x))


def main() -> None:
    times = np.array([1.0e4, 1.0e5], dtype=float)
    reference = solve_state(_config("fullhide_1d"), times)
    dg_state = solve_state(_config("dg_1d"), times)

    report = dg_state.adapter_reports["electron"]
    assert report.solver == "dg_1d"
    assert report.grid_semantics == "log-gamma-1d-dg"
    assert report.status == "ok"

    electron = dg_state.electron
    assert electron.gam_e.shape == (41,)
    assert electron.d_n_gam_e.shape == (41, 24)
    assert electron.l_syn_spec.shape == (31, 24)
    assert electron.seed_syn.shape == (31, 24)
    assert np.all(np.isfinite(electron.d_n_gam_e))
    assert np.all(np.isfinite(electron.l_syn_spec))
    assert np.min(electron.d_n_gam_e) >= 0.0
    assert np.max(electron.l_syn_spec) > 0.0

    count_ratio = _electron_count_x(dg_state) / _electron_count_x(reference)
    assert 0.8 < count_ratio < 1.3

    luminosity_ratio = np.max(electron.l_syn_spec) / np.max(reference.electron.l_syn_spec)
    assert 0.2 < float(luminosity_ratio) < 5.0

    jump_reference = solve_state(_jump_config("fullhide_1d"), times)
    jump_dg_state = solve_state(_jump_config("dg_1d"), times)
    assert np.all(np.isfinite(jump_dg_state.electron.d_n_gam_e))
    assert np.all(np.isfinite(jump_dg_state.electron.l_syn_spec))
    assert np.min(jump_dg_state.electron.d_n_gam_e) >= 0.0

    jump_count_ratio = _electron_count_x(jump_dg_state) / _electron_count_x(jump_reference)
    assert 0.9 < jump_count_ratio < 1.1

    jump_luminosity_ratio = np.max(jump_dg_state.electron.l_syn_spec) / np.max(jump_reference.electron.l_syn_spec)
    assert 0.7 < float(jump_luminosity_ratio) < 1.35


if __name__ == "__main__":
    main()
