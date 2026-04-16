from __future__ import annotations

import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_component_backend import extract_physical_solution_from_state, solve_model_state_from_setup
from asgard_postprocess import compute_band_fluxes
from asgard_presets import build_baseline_config
from asgard_runtime import solve_dynamics, solve_electron
from asgard_setup import build_simulation_setup
from asgard_paths import asgard_doc_path


OUTPUT_JSON = asgard_doc_path("order_convergence.json")
ELECTRON_ORDER_TARGET = 2.0
RADIATION_ORDER_TARGET = 2.0


@dataclass
class OrderResult:
    name: str
    category: str
    target_order: float
    hypothesis: str
    refinement_ratio: float
    coarse_error: float
    medium_error: float
    measured_order: float
    passed: bool
    extra: dict


def _json_safe(value):
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def _resolution_triplet(solver: str) -> list[int]:
    return [61, 101, 161]


def _log_interp(x_new: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return np.interp(np.log10(x_new), np.log10(x), y)


def _log_interp_positive(x_new: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(x_new, dtype=float)
    mask = np.isfinite(y) & (y > 0.0)
    if np.count_nonzero(mask) < 2:
        return out
    x_src = x[mask]
    y_src = y[mask]
    lo = max(float(np.min(x_new)), float(x_src[0]))
    hi = min(float(np.max(x_new)), float(x_src[-1]))
    if not (hi > lo):
        return out
    tgt = (x_new >= lo) & (x_new <= hi)
    if not np.any(tgt):
        return out
    out[tgt] = 10.0 ** np.interp(np.log10(x_new[tgt]), np.log10(x_src), np.log10(y_src))
    return out


def _max_rel(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    scale = np.maximum(np.abs(b), 1.0e-99)
    return float(np.max(np.abs(a - b) / scale))


def _estimate_order(err_coarse: float, err_medium: float, ratio: float = 2.0) -> float:
    if err_coarse == 0.0 and err_medium == 0.0:
        return float("inf")
    if err_medium == 0.0 and err_coarse > 0.0:
        return float("inf")
    if err_coarse <= 0.0 or err_medium <= 0.0:
        return float("nan")
    return float(np.log(err_coarse / err_medium) / np.log(ratio))


def _log_grid_ratio_from_gamma(gamma_coarse: np.ndarray, gamma_fine: np.ndarray) -> float:
    x_coarse = np.log10(np.asarray(gamma_coarse, dtype=float))
    x_fine = np.log10(np.asarray(gamma_fine, dtype=float))
    h_coarse = float(np.mean(np.diff(x_coarse)))
    h_fine = float(np.mean(np.diff(x_fine)))
    return h_coarse / h_fine


def _direct_case(solver: str, n: int) -> dict[str, np.ndarray | int]:
    cfg = build_baseline_config(
        electron_solver=solver,
        num_gam_e=n,
        num_nu=n,
        num_r=96,
        include_forward_ssc=False,
        index_y=0,
    )
    setup = build_simulation_setup(cfg)
    dyn = solve_dynamics(setup.boundary, cfg)
    ele = solve_electron(setup.boundary, dyn, setup.seed_frequency_hz, cfg)
    return {
        "gam_e": np.asarray(ele.gam_e, dtype=float),
        "d_n_gam_e": np.asarray(ele.d_n_gam_e, dtype=float),
        "num_r": int(cfg.num_r),
    }


def _band_case(solver: str, n: int) -> dict[str, np.ndarray]:
    cfg = build_baseline_config(
        electron_solver=solver,
        num_gam_e=n,
        num_nu=n,
        num_r=96,
        num_tobs=32,
        include_forward_ssc=True,
    )
    setup = build_simulation_setup(cfg)
    state = solve_model_state_from_setup(cfg, setup)
    physical = extract_physical_solution_from_state(state)
    return {
        "bands_flux": np.asarray(compute_band_fluxes(setup, physical, cfg), dtype=float),
    }


def _dynamics_order() -> list[OrderResult]:
    num_r_values = [64, 128, 256]
    sols: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for num_r in num_r_values:
        cfg = build_baseline_config(num_r=num_r, include_forward_ssc=False, index_y=0)
        setup = build_simulation_setup(cfg)
        dyn = solve_dynamics(setup.boundary, cfg)
        sols[num_r] = (
            np.asarray(dyn.r_tobs, dtype=float),
            np.asarray(dyn.r_gamma, dtype=float),
            np.asarray(dyn.radius, dtype=float),
        )

    t_min = max(float(sols[64][0][1]), float(sols[128][0][1]), float(sols[256][0][1]))
    t_max = min(float(sols[64][0][-1]), float(sols[128][0][-1]), float(sols[256][0][-1]))
    t_common = np.logspace(np.log10(t_min), np.log10(t_max), 48)
    gamma_64 = _log_interp(t_common, sols[64][0], sols[64][1])
    gamma_128 = _log_interp(t_common, sols[128][0], sols[128][1])
    gamma_256 = _log_interp(t_common, sols[256][0], sols[256][1])
    radius_64 = _log_interp(t_common, sols[64][0], sols[64][2])
    radius_128 = _log_interp(t_common, sols[128][0], sols[128][2])
    radius_256 = _log_interp(t_common, sols[256][0], sols[256][2])

    err_gamma_c = _max_rel(gamma_64, gamma_128)
    err_gamma_m = _max_rel(gamma_128, gamma_256)
    err_radius_c = _max_rel(radius_64, radius_128)
    err_radius_m = _max_rel(radius_128, radius_256)
    order_gamma = _estimate_order(err_gamma_c, err_gamma_m)
    order_radius = _estimate_order(err_radius_c, err_radius_m)
    return [
        OrderResult(
            name="dynamics-forward-gamma",
            category="dynamics",
            target_order=0.0,
            hypothesis="动力学 Gamma 的粗阶数至少应为正。",
            refinement_ratio=2.0,
            coarse_error=err_gamma_c,
            medium_error=err_gamma_m,
            measured_order=order_gamma,
            passed=np.isfinite(order_gamma) and order_gamma > 0.0,
            extra={"num_r": num_r_values},
        ),
        OrderResult(
            name="dynamics-forward-radius",
            category="dynamics",
            target_order=0.0,
            hypothesis="动力学半径 R(t_obs) 的粗阶数至少应为正。",
            refinement_ratio=2.0,
            coarse_error=err_radius_c,
            medium_error=err_radius_m,
            measured_order=order_radius,
            passed=np.isfinite(order_radius) and order_radius > 0.0,
            extra={"num_r": num_r_values},
        ),
    ]


def _solver_spectrum_order(solver: str) -> list[OrderResult]:
    num_values = _resolution_triplet(solver)
    shell_fracs = [0.3, 0.5, 0.8]
    direct = {n: _direct_case(solver, n) for n in num_values}
    ratio = _log_grid_ratio_from_gamma(direct[num_values[0]]["gam_e"], direct[num_values[1]]["gam_e"])

    def _whole_shape_error(test_case: dict, ref_case: dict) -> float:
        gamma_ref = ref_case["gam_e"]
        num_r_ref = ref_case["num_r"]
        errs: list[float] = []
        for frac in shell_fracs:
            shell_ref = min(max(int(round(frac * (num_r_ref - 1))), 1), num_r_ref - 1)
            shell_test = min(max(int(round(frac * (test_case["num_r"] - 1))), 1), test_case["num_r"] - 1)
            y_ref = ref_case["d_n_gam_e"][:, shell_ref]
            y_test = _log_interp_positive(gamma_ref, test_case["gam_e"], test_case["d_n_gam_e"][:, shell_test])
            peak = float(np.max(y_ref))
            if not np.isfinite(peak) or peak <= 0.0:
                continue
            mask = y_ref > 1.0e-8 * peak
            if np.count_nonzero(mask) < 4:
                continue
            xlog = np.log10(gamma_ref[mask])
            numer = np.trapezoid(np.abs(y_test[mask] - y_ref[mask]), xlog)
            denom = np.trapezoid(np.abs(y_ref[mask]), xlog)
            if denom > 0.0:
                errs.append(float(numer / denom))
        return float(max(errs)) if errs else float("nan")

    order_value = _estimate_order(
        _whole_shape_error(direct[num_values[0]], direct[num_values[1]]),
        _whole_shape_error(direct[num_values[1]], direct[num_values[2]]),
        ratio=ratio,
    )
    return [
        OrderResult(
            name=f"{solver}-electron-spectrum",
            category="electron",
            target_order=ELECTRON_ORDER_TARGET,
            hypothesis="整条电子谱的粗阶数检查，只保留整体谱形误差。",
            refinement_ratio=ratio,
            coarse_error=float("nan"),
            medium_error=float("nan"),
            measured_order=order_value,
            passed=np.isfinite(order_value) and order_value > ELECTRON_ORDER_TARGET,
            extra={"num_gam_e_num_nu": num_values, "shell_fractions": shell_fracs},
        )
    ]


def _solver_radiation_order(solver: str) -> list[OrderResult]:
    num_values = _resolution_triplet(solver)
    bands = {n: _band_case(solver, n) for n in num_values}
    ratio = 2.0
    order_value = _estimate_order(
        _max_rel(bands[num_values[0]]["bands_flux"], bands[num_values[1]]["bands_flux"]),
        _max_rel(bands[num_values[1]]["bands_flux"], bands[num_values[2]]["bands_flux"]),
        ratio=ratio,
    )
    return [
        OrderResult(
            name=f"{solver}-radiation-bands-flux",
            category="radiation",
            target_order=RADIATION_ORDER_TARGET,
            hypothesis="观测端只保留多波段总通量的粗阶数检查。",
            refinement_ratio=ratio,
            coarse_error=float("nan"),
            medium_error=float("nan"),
            measured_order=order_value,
            passed=np.isfinite(order_value) and order_value > RADIATION_ORDER_TARGET,
            extra={"num_gam_e_num_nu": num_values},
        )
    ]


def _build_summary(results: list[OrderResult]) -> dict[str, object]:
    summary: dict[str, object] = {}
    for category in ("electron", "radiation", "dynamics"):
        cat_results = [item for item in results if item.category == category]
        if not cat_results:
            continue
        finite_orders = [item.measured_order for item in cat_results if np.isfinite(item.measured_order)]
        summary[category] = {
            "count": len(cat_results),
            "failed": int(sum(not item.passed for item in cat_results)),
            "min_measured_order": None if not finite_orders else float(min(finite_orders)),
        }
    return summary


def main() -> None:
    checks: list[OrderResult] = []
    checks.extend(_dynamics_order())
    for solver in ("fullhide", "slc1", "charint"):
        checks.extend(_solver_spectrum_order(solver))
        checks.extend(_solver_radiation_order(solver))

    payload = {
        "targets": {
            "electron": ELECTRON_ORDER_TARGET,
            "radiation": RADIATION_ORDER_TARGET,
        },
        "summary": _build_summary(checks),
        "results": [_json_safe(asdict(item)) for item in checks],
        "total": len(checks),
        "failed": int(sum(not item.passed for item in checks)),
    }
    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if payload["failed"] != 0:
        raise SystemExit(1)
    print("PASS: order convergence check succeeded.")


if __name__ == "__main__":
    main()
