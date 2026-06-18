from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.asgard_config import RuntimeConfig, ReverseShockConfig
from asgard_core.asgard_state import solve_state


def _config(solver: str, *, medium: str) -> RuntimeConfig:
    wind = medium == "wind"
    return RuntimeConfig(
        electron_solver=solver,
        a_star=3.0e-3 if wind else -1.0,
        index_y=0,
        include_forward_ssc=False,
        num_threads=1,
        num_gam_e=121,
        num_nu=25,
        num_r=34,
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


def _assert_reverse_state(state) -> None:
    rs = state.dynamics.reverse_shock
    assert rs is not None
    assert np.all(np.isfinite(rs.d_n_gam_e))
    assert np.min(rs.d_n_gam_e) >= 0.0
    assert np.max(rs.d_n_gam_e) > 0.0
    assert state.reverse_emission is not None
    assert np.all(np.isfinite(state.reverse_emission.l_syn_spec))
    assert np.max(state.reverse_emission.l_syn_spec) > 0.0


def _count_ratio(dg_state, ref_state) -> float:
    dg = dg_state.dynamics.reverse_shock
    ref = ref_state.dynamics.reverse_shock
    assert dg is not None and ref is not None
    x_dg = np.log10(np.asarray(dg.gam_e, dtype=float))
    x_ref = np.log10(np.asarray(ref.gam_e, dtype=float))
    n_dg = np.asarray(dg.d_n_gam_e[:, -1], dtype=float) * dg.gam_e * np.log(10.0)
    n_ref = np.asarray(ref.d_n_gam_e[:, -1], dtype=float) * ref.gam_e * np.log(10.0)
    return float(np.trapezoid(n_dg, x_dg) / np.trapezoid(n_ref, x_ref))


def _run_case(medium: str) -> None:
    times = np.array([1.0e2, 1.0e4], dtype=float)
    reference = solve_state(_config("fullhide_1d", medium=medium), times)
    dg_state = solve_state(_config("dg_1d", medium=medium), times)
    _assert_reverse_state(reference)
    _assert_reverse_state(dg_state)
    ratio = _count_ratio(dg_state, reference)
    assert 0.5 < ratio < 2.0


def main() -> None:
    _run_case("ism")
    _run_case("wind")
    print("reverse-shared-solver-smoke-ok")


if __name__ == "__main__":
    main()
