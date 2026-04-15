from __future__ import annotations

import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_presets import build_baseline_config
from asgard_runtime import solve_dynamics, solve_electron
from asgard_setup import build_simulation_setup
from asgard_component_backend import solve_model_state_from_setup, extract_physical_solution_from_state
from asgard_setup import build_simulation_setup
from asgard_postprocess import compute_band_fluxes


OUTPUT_JSON = ROOT / "output" / "vegasafterglow_doc" / "order_convergence.json"
CACHE_DIR = ROOT / "output" / "vegasafterglow_doc" / "order_convergence_cache"
ORDER_CACHE_VERSION = 9


@dataclass
class OrderResult:
    name: str
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


def _cache_path(kind: str, name: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / f"{kind}_{name}.npz"


def _load_npz_cache(path: Path) -> dict | None:
    if not path.exists():
        return None
    data = np.load(path, allow_pickle=True)
    version = int(np.asarray(data["version"]).item())
    if version != ORDER_CACHE_VERSION:
        return None
    return {key: data[key] for key in data.files if key != "version"}


def _save_npz_cache(path: Path, **arrays) -> None:
    np.savez_compressed(path, version=np.array(ORDER_CACHE_VERSION), **arrays)


def _observed_case(solver: str, n: int) -> dict[str, np.ndarray]:
    path = _cache_path("observed", f"{solver}_{n}")
    cached = _load_npz_cache(path)
    if cached is not None:
        return {
            "bands_flux": np.asarray(cached["bands_flux"], dtype=float),
            "nu_a": np.asarray(cached["nu_a"], dtype=float),
        }

    cfg = build_baseline_config(
        electron_solver=solver,
        num_gam_e=n,
        num_nu=n,
        num_r=160,
        num_tobs=64,
        include_forward_ssc=True,
    )
    setup = build_simulation_setup(cfg)
    state = solve_model_state_from_setup(cfg, setup)
    physical = extract_physical_solution_from_state(state)
    bands_flux = compute_band_fluxes(setup, physical, cfg)
    payload = {
        "bands_flux": np.asarray(bands_flux, dtype=float),
        "nu_a": np.asarray(physical.nu_a, dtype=float),
    }
    _save_npz_cache(path, **payload)
    return payload


def _direct_case(solver: str, n: int) -> dict[str, np.ndarray | int]:
    path = _cache_path("direct", f"{solver}_{n}")
    cached = _load_npz_cache(path)
    if cached is not None:
        return {
            "gam_e": np.asarray(cached["gam_e"], dtype=float),
            "d_n_gam_e": np.asarray(cached["d_n_gam_e"], dtype=float),
            "num_r": int(np.asarray(cached["num_r"]).item()),
        }

    cfg = build_baseline_config(
        electron_solver=solver,
        num_gam_e=n,
        num_nu=n,
        num_r=160,
        include_forward_ssc=False,
        index_y=0,
    )
    setup = build_simulation_setup(cfg)
    dyn = solve_dynamics(setup.boundary, cfg)
    ele = solve_electron(setup.boundary, dyn, setup.seed_frequency_hz, cfg)
    payload = {
        "gam_e": np.asarray(ele.gam_e, dtype=float),
        "d_n_gam_e": np.asarray(ele.d_n_gam_e, dtype=float),
        "num_r": np.array(int(cfg.num_r)),
    }
    _save_npz_cache(path, **payload)
    return {
        "gam_e": payload["gam_e"],
        "d_n_gam_e": payload["d_n_gam_e"],
        "num_r": int(payload["num_r"].item()),
    }


def _log_interp(x_new: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return np.interp(np.log10(x_new), np.log10(x), y)


def _max_rel(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    scale = np.maximum(np.abs(b), 1.0e-99)
    return float(np.max(np.abs(a - b) / scale))


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


def _electron_support_metrics(gamma: np.ndarray, spectrum: np.ndarray) -> tuple[float, float, float, float]:
    gamma = np.asarray(gamma, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    peak_idx = int(np.argmax(spectrum))
    peak_gamma = float(gamma[peak_idx])
    peak_value = float(spectrum[peak_idx])
    if not np.isfinite(peak_value) or peak_value <= 0.0:
        return float("nan"), float("nan"), float("nan"), float("nan")
    threshold = 1.0e-8 * peak_value
    mask = spectrum > threshold
    if np.count_nonzero(mask) < 4:
        return peak_gamma, float("nan"), float("nan"), peak_value
    support_idx = np.where(mask)[0]
    g_lo = _electron_support_edge(gamma, spectrum, threshold, support_idx[0], side="low")
    g_hi = _electron_support_edge(gamma, spectrum, threshold, support_idx[-1], side="high")
    return peak_gamma, g_lo, g_hi, peak_value


def _electron_support_edge(
    gamma: np.ndarray,
    spectrum: np.ndarray,
    threshold: float,
    idx_inside: int,
    side: str,
) -> float:
    gamma = np.asarray(gamma, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)

    if side == "low":
        idx0 = max(idx_inside - 1, 0)
        idx1 = idx_inside
    else:
        idx0 = idx_inside
        idx1 = min(idx_inside + 1, gamma.size - 1)

    y0 = float(spectrum[idx0])
    y1 = float(spectrum[idx1])
    x0 = float(np.log10(gamma[idx0]))
    x1 = float(np.log10(gamma[idx1]))

    if idx0 == idx1 or y0 <= 0.0 or y1 <= 0.0 or x1 == x0:
        return float(gamma[idx_inside])

    if (y0 - threshold) * (y1 - threshold) > 0.0:
        return float(gamma[idx_inside])

    logy0 = float(np.log10(y0))
    logy1 = float(np.log10(y1))
    logt = float(np.log10(threshold))
    if logy1 == logy0:
        return float(gamma[idx_inside])

    frac = (logt - logy0) / (logy1 - logy0)
    frac = float(np.clip(frac, 0.0, 1.0))
    x_edge = x0 + frac * (x1 - x0)
    return float(10.0**x_edge)


def _aligned_shape_error(
    gamma_ref: np.ndarray,
    y_ref: np.ndarray,
    gamma_test: np.ndarray,
    y_test: np.ndarray,
) -> float:
    peak_ref, g_lo_ref, g_hi_ref, peak_ref_val = _electron_support_metrics(gamma_ref, y_ref)
    peak_test, _, _, peak_test_val = _electron_support_metrics(gamma_test, y_test)
    if not np.isfinite(peak_ref) or not np.isfinite(peak_test):
        return float("nan")
    if not np.isfinite(peak_ref_val) or not np.isfinite(peak_test_val):
        return float("nan")
    if peak_ref_val <= 0.0 or peak_test_val <= 0.0:
        return float("nan")

    shift = peak_test / peak_ref
    gamma_query = np.asarray(gamma_ref, dtype=float) * shift
    y_shifted = _log_interp_positive(gamma_query, gamma_test, y_test)
    y_ref_n = np.asarray(y_ref, dtype=float) / peak_ref_val
    y_shifted_n = np.asarray(y_shifted, dtype=float) / peak_test_val
    mask = (gamma_ref >= g_lo_ref) & (gamma_ref <= g_hi_ref) & (y_ref_n > 1.0e-8)
    if np.count_nonzero(mask) < 4:
        return float("nan")
    xlog = np.log10(gamma_ref[mask])
    numer = np.trapezoid(np.abs(y_shifted_n[mask] - y_ref_n[mask]), xlog)
    denom = np.trapezoid(np.abs(y_ref_n[mask]), xlog)
    if denom <= 0.0:
        return float("nan")
    return float(numer / denom)


def _dynamics_order() -> list[OrderResult]:
    num_r_values = [80, 160, 320]
    sols: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for num_r in num_r_values:
        path = _cache_path("dynamics", f"{num_r}")
        cached = _load_npz_cache(path)
        if cached is None:
            cfg = build_baseline_config(num_r=num_r, include_forward_ssc=False, index_y=0)
            setup = build_simulation_setup(cfg)
            dyn = solve_dynamics(setup.boundary, cfg)
            payload = {
                "r_tobs": np.asarray(dyn.r_tobs, dtype=float),
                "r_gamma": np.asarray(dyn.r_gamma, dtype=float),
                "radius": np.asarray(dyn.radius, dtype=float),
            }
            _save_npz_cache(path, **payload)
            cached = payload
        sols[num_r] = (
            np.asarray(cached["r_tobs"], dtype=float),
            np.asarray(cached["r_gamma"], dtype=float),
            np.asarray(cached["radius"], dtype=float),
        )

    t_min = max(float(sols[80][0][1]), float(sols[160][0][1]), float(sols[320][0][1]))
    t_max = min(float(sols[80][0][-1]), float(sols[160][0][-1]), float(sols[320][0][-1]))
    t_common = np.logspace(np.log10(t_min), np.log10(t_max), 64)

    gamma_80 = _log_interp(t_common, sols[80][0], sols[80][1])
    gamma_160 = _log_interp(t_common, sols[160][0], sols[160][1])
    gamma_320 = _log_interp(t_common, sols[320][0], sols[320][1])
    radius_80 = _log_interp(t_common, sols[80][0], sols[80][2])
    radius_160 = _log_interp(t_common, sols[160][0], sols[160][2])
    radius_320 = _log_interp(t_common, sols[320][0], sols[320][2])

    err_gamma_c = _max_rel(gamma_80, gamma_160)
    err_gamma_m = _max_rel(gamma_160, gamma_320)
    err_radius_c = _max_rel(radius_80, radius_160)
    err_radius_m = _max_rel(radius_160, radius_320)
    order_gamma = _estimate_order(err_gamma_c, err_gamma_m)
    order_radius = _estimate_order(err_radius_c, err_radius_m)

    return [
        OrderResult(
            name="dynamics-forward-gamma",
            hypothesis="动力学观测系 Gamma 的有效收敛阶应为正，若显著低于预期，说明步长控制、分段物理或后处理在降阶。",
            refinement_ratio=2.0,
            coarse_error=err_gamma_c,
            medium_error=err_gamma_m,
            measured_order=order_gamma,
            passed=np.isfinite(order_gamma) and order_gamma > 0.0,
            extra={"num_r": num_r_values},
        ),
        OrderResult(
            name="dynamics-forward-radius",
            hypothesis="动力学观测系半径 R(t_obs) 的有效收敛阶应为正，若偏低，需要继续拆分动力学和观测时间映射误差。",
            refinement_ratio=2.0,
            coarse_error=err_radius_c,
            medium_error=err_radius_m,
            measured_order=order_radius,
            passed=np.isfinite(order_radius) and order_radius > 0.0,
            extra={"num_r": num_r_values},
        ),
    ]


def _solver_observed_order(solver: str) -> list[OrderResult]:
    num_values = [20, 32, 40] if solver == "slc1_mmg2" else [61, 121, 241]
    results = {}
    for n in num_values:
        results[n] = _observed_case(solver, n)

    flux_c = _max_rel(results[num_values[0]]["bands_flux"], results[num_values[1]]["bands_flux"])
    flux_m = _max_rel(results[num_values[1]]["bands_flux"], results[num_values[2]]["bands_flux"])
    nua_c = _max_rel(results[num_values[0]]["nu_a"], results[num_values[1]]["nu_a"])
    nua_m = _max_rel(results[num_values[1]]["nu_a"], results[num_values[2]]["nu_a"])

    order_flux = _estimate_order(flux_c, flux_m)
    order_nua = _estimate_order(nua_c, nua_m)
    hypothesis = {
        "fullhide": "fullhide 是默认一阶隐式电子求解器，观测量的有效收敛阶至少应为正。",
        "t2g1": "t2g1 目标是高于 fullhide 的时间精度，观测量收敛阶应高于 fullhide。",
        "slc1": "slc1 是实验性的守恒半拉格朗日电子求解器，若方向正确，观测量收敛阶至少应为正且速度应优于 fullhide。",
        "slc1_mmg2": "slc1_mmg2 是固定网格数的 moving-mesh 半拉格朗日电子求解器，目标是在 20-40 网格下把特征电子能量位置和观测量收敛阶抬到二阶附近。",
        "weno5": "weno5 若空间离散真正生效，观测量的有效收敛阶不应劣于 fullhide。",
    }[solver]

    return [
        OrderResult(
            name=f"{solver}-bands-flux",
            hypothesis=hypothesis,
            refinement_ratio=2.0,
            coarse_error=flux_c,
            medium_error=flux_m,
            measured_order=order_flux,
            passed=np.isfinite(order_flux) and order_flux > 0.0,
            extra={"num_gam_e_num_nu": num_values},
        ),
        OrderResult(
            name=f"{solver}-nu_a",
            hypothesis="nu_a 是 SSA 派生量；只要链路稳定，它的有效收敛阶应为正。",
            refinement_ratio=2.0,
            coarse_error=nua_c,
            medium_error=nua_m,
            measured_order=order_nua,
            passed=np.isfinite(order_nua) and order_nua > 0.0,
            extra={"num_gam_e_num_nu": num_values},
        ),
    ]


def _solver_spectrum_order(solver: str) -> list[OrderResult]:
    num_values = [20, 32, 40] if solver == "slc1_mmg2" else [61, 121, 241]
    shell_fracs = [0.3, 0.5, 0.8]
    direct = {}
    for n in num_values:
        direct[n] = _direct_case(solver, n)

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
            y_ref = np.asarray(y_ref, dtype=float)
            peak = float(np.max(y_ref))
            if not np.isfinite(peak) or peak <= 0.0:
                continue
            mask = y_ref > 1.0e-8 * peak
            if np.count_nonzero(mask) < 4:
                continue
            xlog = np.log10(gamma_ref[mask])
            numer = np.trapezoid(np.abs(y_test[mask] - y_ref[mask]), xlog)
            denom = np.trapezoid(np.abs(y_ref[mask]), xlog)
            if denom <= 0.0:
                continue
            errs.append(float(numer / denom))
        return float(max(errs)) if errs else float("nan")

    def _support_metric_error(test_case: dict, ref_case: dict, metric: str) -> float:
        num_r_ref = ref_case["num_r"]
        errs: list[float] = []
        for frac in shell_fracs:
            shell_ref = min(max(int(round(frac * (num_r_ref - 1))), 1), num_r_ref - 1)
            shell_test = min(max(int(round(frac * (test_case["num_r"] - 1))), 1), test_case["num_r"] - 1)
            peak_ref, g_lo_ref, g_hi_ref, _ = _electron_support_metrics(ref_case["gam_e"], ref_case["d_n_gam_e"][:, shell_ref])
            peak_test, g_lo_test, g_hi_test, _ = _electron_support_metrics(test_case["gam_e"], test_case["d_n_gam_e"][:, shell_test])
            value_ref = {"peak_gamma": peak_ref, "g_lo": g_lo_ref, "g_hi": g_hi_ref}[metric]
            value_test = {"peak_gamma": peak_test, "g_lo": g_lo_test, "g_hi": g_hi_test}[metric]
            if not np.isfinite(value_ref) or not np.isfinite(value_test):
                continue
            if value_ref <= 0.0:
                continue
            errs.append(float(abs(value_test - value_ref) / value_ref))
        return float(max(errs)) if errs else float("nan")

    def _aligned_shape_metric(test_case: dict, ref_case: dict) -> float:
        num_r_ref = ref_case["num_r"]
        errs: list[float] = []
        for frac in shell_fracs:
            shell_ref = min(max(int(round(frac * (num_r_ref - 1))), 1), num_r_ref - 1)
            shell_test = min(max(int(round(frac * (test_case["num_r"] - 1))), 1), test_case["num_r"] - 1)
            err = _aligned_shape_error(
                ref_case["gam_e"],
                ref_case["d_n_gam_e"][:, shell_ref],
                test_case["gam_e"],
                test_case["d_n_gam_e"][:, shell_test],
            )
            if np.isfinite(err):
                errs.append(float(err))
        return float(max(errs)) if errs else float("nan")

    spec_c = _whole_shape_error(direct[num_values[0]], direct[num_values[1]])
    spec_m = _whole_shape_error(direct[num_values[1]], direct[num_values[2]])
    order_spec = _estimate_order(spec_c, spec_m, ratio=ratio)

    peak_c = _support_metric_error(direct[num_values[0]], direct[num_values[1]], "peak_gamma")
    peak_m = _support_metric_error(direct[num_values[1]], direct[num_values[2]], "peak_gamma")
    order_peak = _estimate_order(peak_c, peak_m, ratio=ratio)

    glo_c = _support_metric_error(direct[num_values[0]], direct[num_values[1]], "g_lo")
    glo_m = _support_metric_error(direct[num_values[1]], direct[num_values[2]], "g_lo")
    order_glo = _estimate_order(glo_c, glo_m, ratio=ratio)

    ghi_c = _support_metric_error(direct[num_values[0]], direct[num_values[1]], "g_hi")
    ghi_m = _support_metric_error(direct[num_values[1]], direct[num_values[2]], "g_hi")
    order_ghi = _estimate_order(ghi_c, ghi_m, ratio=ratio)

    aligned_c = _aligned_shape_metric(direct[num_values[0]], direct[num_values[1]])
    aligned_m = _aligned_shape_metric(direct[num_values[1]], direct[num_values[2]])
    order_aligned = _estimate_order(aligned_c, aligned_m, ratio=ratio)

    return [
        OrderResult(
            name=f"{solver}-electron-spectrum",
            hypothesis="整条电子谱的混合有效阶数，包含峰位漂移、支撑区漂移和谱形误差。",
            refinement_ratio=ratio,
            coarse_error=spec_c,
            medium_error=spec_m,
            measured_order=order_spec,
            passed=np.isfinite(order_spec) and order_spec > 0.0,
            extra={"num_gam_e_num_nu": num_values, "shell_fractions": shell_fracs},
        ),
        OrderResult(
            name=f"{solver}-electron-peak-gamma",
            hypothesis="谱峰位置误差的有效阶数；若该项低，主问题在峰位锁网格。",
            refinement_ratio=ratio,
            coarse_error=peak_c,
            medium_error=peak_m,
            measured_order=order_peak,
            passed=np.isfinite(order_peak) and order_peak > 0.0,
            extra={"num_gam_e_num_nu": num_values, "shell_fractions": shell_fracs},
        ),
        OrderResult(
            name=f"{solver}-electron-support-low",
            hypothesis="低能支撑边界 g_lo 误差的有效阶数；若该项低，主问题在低能端边界和注入/冷却切换。",
            refinement_ratio=ratio,
            coarse_error=glo_c,
            medium_error=glo_m,
            measured_order=order_glo,
            passed=np.isfinite(order_glo) and order_glo > 0.0,
            extra={"num_gam_e_num_nu": num_values, "shell_fractions": shell_fracs},
        ),
        OrderResult(
            name=f"{solver}-electron-support-high",
            hypothesis="高能支撑边界 g_hi 误差的有效阶数；若该项低，主问题在高能截止和高能耗散。",
            refinement_ratio=ratio,
            coarse_error=ghi_c,
            medium_error=ghi_m,
            measured_order=order_ghi,
            passed=np.isfinite(order_ghi) and order_ghi > 0.0,
            extra={"num_gam_e_num_nu": num_values, "shell_fractions": shell_fracs},
        ),
        OrderResult(
            name=f"{solver}-electron-shape-aligned",
            hypothesis="对齐峰位后的谱形误差有效阶数；若该项仍低，主问题在谱形离散本身。",
            refinement_ratio=ratio,
            coarse_error=aligned_c,
            medium_error=aligned_m,
            measured_order=order_aligned,
            passed=np.isfinite(order_aligned) and order_aligned > 0.0,
            extra={"num_gam_e_num_nu": num_values, "shell_fractions": shell_fracs},
        ),
    ]


def main() -> None:
    checks: list[OrderResult] = []
    checks.extend(_solver_observed_order("slc1_mmg2"))
    checks.extend(_solver_spectrum_order("slc1_mmg2"))

    payload = {
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
