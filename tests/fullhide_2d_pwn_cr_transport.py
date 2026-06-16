from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_setup import build_boundary
from asgard_core.asgard_state import solve_state
from tests._bench_common import run_case


def _config(**overrides) -> RuntimeConfig:
    return RuntimeConfig(**(dict(
        electron_solver="fullhide_2d",
        num_gam_e=10,
        num_chi=5,
        num_nu=10,
        num_r=10,
        num_theta=8,
        num_tobs=5,
        include_forward_ssc=False,
        electron_adaptive_substeps=False,
    ) | overrides))


def case_boundary_layout() -> dict[str, object]:
    legacy = build_boundary(_config(), 1.0e28)
    pwn = build_boundary(_config(fullhide2d_transport_model="pwn_cr_v1"), 1.0e28)
    free = build_boundary(
        _config(fullhide2d_transport_model="pwn_cr_v1", fullhide2d_escape_mode="free_outer"),
        1.0e28,
    )
    assert pwn.shape == (legacy.size + 3,)
    assert np.array_equal(pwn[:-3], legacy)
    assert pwn[-3] == 1.0
    assert pwn[-2] == 0.0
    assert pwn[-1] == 0.0
    assert free[-1] == 1.0
    return {"legacy_len": int(legacy.size), "pwn_len": int(pwn.size), "r0": float(pwn[26])}


def _electron_count(config: RuntimeConfig) -> tuple[float, bool]:
    state = solve_state(config, np.array([1.0e2, 1.2e2]))
    dnde = np.asarray(state.electron.d_n_gam_e, dtype=float)
    gamma = np.asarray(state.electron.gam_e, dtype=float)
    total = float(np.trapezoid(dnde[:, -1], gamma))
    return total, bool(np.isfinite(dnde).all())


def case_pwn_closed_smoke() -> dict[str, object]:
    total, finite = _electron_count(_config(fullhide2d_transport_model="pwn_cr_v1"))
    assert finite
    assert np.isfinite(total)
    return {"electron_count_last": total}


def case_free_outer_not_larger() -> dict[str, object]:
    closed_total, closed_finite = _electron_count(_config(fullhide2d_transport_model="pwn_cr_v1"))
    free_total, free_finite = _electron_count(
        _config(fullhide2d_transport_model="pwn_cr_v1", fullhide2d_escape_mode="free_outer")
    )
    assert closed_finite and free_finite
    assert free_total <= closed_total * (1.0 + 1.0e-12)
    return {"closed": closed_total, "free_outer": free_total}


def case_stochastic_zero_equivalence() -> dict[str, object]:
    implicit_zero = _config(fullhide2d_transport_model="pwn_cr_v1")
    explicit_zero = _config(fullhide2d_transport_model="pwn_cr_v1", fullhide2d_stochastic_accel_norm=0.0)
    state_a = solve_state(implicit_zero, np.array([1.0e2, 1.2e2]))
    state_b = solve_state(explicit_zero, np.array([1.0e2, 1.2e2]))
    a = np.asarray(state_a.electron.d_n_gam_e, dtype=float)
    b = np.asarray(state_b.electron.d_n_gam_e, dtype=float)
    assert np.array_equal(a, b)
    return {"shape": list(a.shape)}


def _max_adjacent_dex(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    lhs = arr[:-1]
    rhs = arr[1:]
    mask = np.isfinite(lhs) & np.isfinite(rhs) & (lhs > 0.0) & (rhs > 0.0)
    if not np.any(mask):
        return 0.0
    return float(np.max(np.abs(np.log10(rhs[mask]) - np.log10(lhs[mask]))))


def case_pwn_spectral_smoothness() -> dict[str, object]:
    config = _config(fullhide2d_transport_model="pwn_cr_v1", num_r=16, num_tobs=8)
    state = solve_state(config, np.array([1.0e2, 1.2e2, 1.5e2]))
    jumps = {
        "nu_m": _max_adjacent_dex(np.asarray(state.electron.nu_m, dtype=float)),
        "nu_c": _max_adjacent_dex(np.asarray(state.electron.nu_c, dtype=float)),
        "nu_a": _max_adjacent_dex(np.asarray(state.electron.nu_a, dtype=float)),
    }
    assert np.isfinite(np.asarray(list(jumps.values()), dtype=float)).all()
    assert max(jumps.values()) < 4.0
    return jumps


def case_stochastic_nonzero_smoke() -> dict[str, object]:
    config = _config(fullhide2d_transport_model="pwn_cr_v1", fullhide2d_stochastic_accel_norm=1.0e-3)
    total, finite = _electron_count(config)
    assert finite
    assert np.isfinite(total)
    return {"electron_count_last": total}


def main() -> None:
    results = [
        run_case("[1/6] pwn_cr:boundary_layout", case_boundary_layout),
        run_case("[2/6] pwn_cr:closed_smoke", case_pwn_closed_smoke),
        run_case("[3/6] pwn_cr:free_outer_not_larger", case_free_outer_not_larger),
        run_case("[4/6] pwn_cr:stochastic_zero_equivalence", case_stochastic_zero_equivalence),
        run_case("[5/6] pwn_cr:spectral_smoothness", case_pwn_spectral_smoothness),
        run_case("[6/6] pwn_cr:stochastic_nonzero_smoke", case_stochastic_nonzero_smoke),
    ]
    if any(item["status"] == "FAIL" for item in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
