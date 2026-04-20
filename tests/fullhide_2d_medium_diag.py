from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_config import FitConfig
from ASGARD.api import _build_model_from_fit_config, _make_details
from asgard_state import solve_state


def _max_adjacent_dex_jump(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    lhs = arr[:-1]
    rhs = arr[1:]
    mask = np.isfinite(lhs) & np.isfinite(rhs) & (lhs > 0.0) & (rhs > 0.0)
    if not np.any(mask):
        return float("nan")
    return float(np.max(np.abs(np.log10(rhs[mask]) - np.log10(lhs[mask]))))


def _first_nonfinite_index(values: np.ndarray) -> int | None:
    arr = np.asarray(values, dtype=float)
    bad = np.where(~np.isfinite(arr))[0]
    if bad.size == 0:
        return None
    return int(bad[0])


def _first_bad_shell_index(values: np.ndarray) -> int | None:
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 0:
        return None if np.isfinite(arr) else 0
    finite_mask = np.isfinite(arr)
    if arr.ndim > 1:
        axis = tuple(range(arr.ndim - 1))
        finite_mask = np.all(finite_mask, axis=axis)
    bad = np.where(~finite_mask)[0]
    if bad.size == 0:
        return None
    return int(bad[0])


def _first_nonfinite_multi_index(values: np.ndarray) -> list[int] | None:
    arr = np.asarray(values, dtype=float)
    bad = np.argwhere(~np.isfinite(arr))
    if bad.size == 0:
        return None
    return [int(x) for x in bad[0]]


def _jump_info(values: np.ndarray) -> dict[str, object] | None:
    arr = np.asarray(values, dtype=float)
    lhs = arr[:-1]
    rhs = arr[1:]
    mask = np.isfinite(lhs) & np.isfinite(rhs) & (lhs > 0.0) & (rhs > 0.0)
    if not np.any(mask):
        return None
    idx = np.where(mask)[0]
    jumps = np.abs(np.log10(rhs[mask]) - np.log10(lhs[mask]))
    j = int(np.argmax(jumps))
    return {
        "pair": [int(idx[j]), int(idx[j] + 1)],
        "dex": float(jumps[j]),
        "lhs": float(lhs[mask][j]),
        "rhs": float(rhs[mask][j]),
    }


def main() -> None:
    config = FitConfig(
        electron_solver="fullhide_2d",
        num_gam_e=16,
        num_chi=8,
        num_nu=24,
        num_r=32,
        num_theta=24,
        num_tobs=12,
        plot_lc=False,
        show_plots=False,
    )
    model = _build_model_from_fit_config(config)
    times_s = model.default_times()
    state = solve_state(config, times_s)
    details = _make_details(state.components, patches=[{"phi": 0.0, "theta": 0.0, "weight": 1.0}], state=state)
    gamma_last = float(np.asarray(state.dynamics.r_gamma, dtype=float)[-1])
    radius_last = float(np.asarray(state.dynamics.radius, dtype=float)[-1])
    d_ne = float(getattr(config, "d_ne", 1.0))
    eps_b = float(getattr(config, "epsilon_b", 0.0))
    if float(getattr(config, "a_star", 0.0)) > 0.0:
        d_ne = float(config.a_star) * 3.0e35 / max(radius_last * radius_last, 1.0)
    db_min = 0.39 * np.sqrt(max(eps_b * d_ne * gamma_last * max(gamma_last - 1.0, 0.0), 0.0))
    gam_e_max_max = np.inf if db_min <= 0.0 else 3.0 * 9.1094e-28 * 2.9979e10 / np.sqrt(8.0 * db_min * (4.8032e-10 ** 3))

    summary = {
        "r_tobs_first_nonfinite": _first_nonfinite_index(np.asarray(state.dynamics.r_tobs, dtype=float)),
        "r_gamma_first_nonfinite": _first_nonfinite_index(np.asarray(state.dynamics.r_gamma, dtype=float)),
        "radius_first_nonfinite": _first_nonfinite_index(np.asarray(state.dynamics.radius, dtype=float)),
        "r_gamma_tail": np.asarray(state.dynamics.r_gamma, dtype=float)[-8:].tolist(),
        "radius_tail": np.asarray(state.dynamics.radius, dtype=float)[-8:].tolist(),
        "gamma_last": gamma_last,
        "radius_last": radius_last,
        "db_min_estimate": db_min,
        "gam_e_max_max_estimate": gam_e_max_max,
        "state_dnde_shape": list(np.asarray(state.electron.d_n_gam_e, dtype=float).shape),
        "state_dnde_chi_shape": list(np.asarray(state.electron.d_n_gam_e_chi, dtype=float).shape),
        "nu_m_max_adj_dex": _max_adjacent_dex_jump(np.asarray(state.electron.nu_m, dtype=float)),
        "nu_c_max_adj_dex": _max_adjacent_dex_jump(np.asarray(state.electron.nu_c, dtype=float)),
        "nu_a_max_adj_dex": _max_adjacent_dex_jump(np.asarray(state.electron.nu_a, dtype=float)),
        "nu_m_jump": _jump_info(np.asarray(state.electron.nu_m, dtype=float)),
        "nu_c_jump": _jump_info(np.asarray(state.electron.nu_c, dtype=float)),
        "nu_a_jump": _jump_info(np.asarray(state.electron.nu_a, dtype=float)),
        "nu_m_tail": np.asarray(state.electron.nu_m, dtype=float)[-6:].tolist(),
        "nu_c": np.asarray(state.electron.nu_c, dtype=float).tolist(),
        "nu_c_first_nonfinite": _first_nonfinite_index(np.asarray(state.electron.nu_c, dtype=float)),
        "nu_c_tail": np.asarray(state.electron.nu_c, dtype=float)[-4:].tolist(),
        "nu_a_tail": np.asarray(state.electron.nu_a, dtype=float)[-6:].tolist(),
        "nu_c_finite_mask_tail": np.isfinite(np.asarray(state.electron.nu_c, dtype=float))[-8:].astype(int).tolist(),
        "state_nu_c_tail": np.asarray(state.electron.nu_c, dtype=float)[-8:].tolist(),
        "details_nu_c_tail": np.asarray(details.fwd.nu_c, dtype=float)[-8:].tolist(),
        "state_nu_m_tail": np.asarray(state.electron.nu_m, dtype=float)[-8:].tolist(),
        "state_nu_a_tail": np.asarray(state.electron.nu_a, dtype=float)[-8:].tolist(),
        "state_dnde_finite": bool(np.isfinite(np.asarray(state.electron.d_n_gam_e, dtype=float)).all()),
        "state_dnde_chi_finite": bool(np.isfinite(np.asarray(state.electron.d_n_gam_e_chi, dtype=float)).all()),
        "state_dnde_first_bad_shell": _first_bad_shell_index(np.asarray(state.electron.d_n_gam_e, dtype=float)),
        "state_dnde_chi_first_bad_shell": _first_bad_shell_index(np.asarray(state.electron.d_n_gam_e_chi, dtype=float)),
        "state_dnde_shell_finite_tail": np.all(np.isfinite(np.asarray(state.electron.d_n_gam_e, dtype=float)), axis=0)[-8:].astype(int).tolist(),
        "state_dnde_chi_shell_finite_tail": np.all(np.isfinite(np.asarray(state.electron.d_n_gam_e_chi, dtype=float)), axis=(0, 1))[-8:].astype(int).tolist(),
        "state_dnde_first_nonfinite_index": _first_nonfinite_multi_index(np.asarray(state.electron.d_n_gam_e, dtype=float)),
        "state_dnde_chi_first_nonfinite_index": _first_nonfinite_multi_index(np.asarray(state.electron.d_n_gam_e_chi, dtype=float)),
    }
    print(summary)


if __name__ == "__main__":
    main()
