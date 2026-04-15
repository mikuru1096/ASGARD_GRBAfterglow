from __future__ import annotations

import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_component_backend import (
    extract_physical_solution_from_state,
    observe_flux_grid_from_state,
    solve_model_state_from_setup,
)
from asgard_postprocess import compute_band_fluxes
from asgard_presets import build_baseline_config
from asgard_runtime import solve_dynamics, solve_electron
from asgard_setup import build_simulation_setup


OUTPUT_JSON = ROOT / "output" / "asgard_doc" / "order_convergence.json"
CACHE_DIR = ROOT / "output" / "asgard_doc" / "order_convergence_cache"
ORDER_CACHE_VERSION = 11
ELECTRON_ORDER_TARGET = 2.0
RADIATION_ORDER_TARGET = 2.0
OBSERVER_FREQUENCIES_HZ = np.logspace(9.0, 26.0, 48)


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


def _resolution_triplet(solver: str) -> list[int]:
    return [20, 32, 40] if solver == "slc1_mmg2" else [61, 121, 241]


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


def _matrix_spectrum_error(
    freq_ref: np.ndarray,
    mat_ref: np.ndarray,
    freq_test: np.ndarray,
    mat_test: np.ndarray,
) -> float:
    freq_ref = np.asarray(freq_ref, dtype=float)
    freq_test = np.asarray(freq_test, dtype=float)
    mat_ref = np.asarray(mat_ref, dtype=float)
    mat_test = np.asarray(mat_test, dtype=float)
    if mat_ref.shape[1] != mat_test.shape[1]:
        raise ValueError("Radiation matrices must share the same shell/time dimension.")
    interp = np.zeros_like(mat_ref)
    for idx in range(mat_ref.shape[1]):
        interp[:, idx] = _log_interp_positive(freq_ref, freq_test, mat_test[:, idx])
    return _max_rel(interp, mat_ref)


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
    frac = float(np.clip((logt - logy0) / (logy1 - logy0), 0.0, 1.0))
    return float(10.0 ** (x0 + frac * (x1 - x0)))


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
    gamma_query = np.asarray(gamma_ref, dtype=float) * (peak_test / peak_ref)
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


def _radiation_case(solver: str, n: int) -> dict[str, np.ndarray]:
    path = _cache_path("radiation", f"{solver}_{n}")
    cached = _load_npz_cache(path)
    if cached is not None:
        return {
            "seed_frequency_hz": np.asarray(cached["seed_frequency_hz"], dtype=float),
            "fwd_sync": np.asarray(cached["fwd_sync"], dtype=float),
            "fwd_ssc": np.asarray(cached["fwd_ssc"], dtype=float),
            "observer_total": np.asarray(cached["observer_total"], dtype=float),
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
    observed = observe_flux_grid_from_state(state, setup.observer_time_s, OBSERVER_FREQUENCIES_HZ)
    payload = {
        "seed_frequency_hz": np.asarray(setup.seed_frequency_hz, dtype=float),
        "fwd_sync": np.asarray(state.component_spectra.fwd_sync, dtype=float),
        "fwd_ssc": np.asarray(state.component_spectra.fwd_ssc, dtype=float),
        "observer_total": np.asarray(observed.components["total"], dtype=float),
    }
    _save_npz_cache(path, **payload)
    return payload


def _band_case(solver: str, n: int) -> dict[str, np.ndarray]:
    path = _cache_path("bands", f"{solver}_{n}")
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
    payload = {
        "bands_flux": np.asarray(compute_band_fluxes(setup, physical, cfg), dtype=float),
        "nu_a": np.asarray(physical.nu_a, dtype=float),
    }
    _save_npz_cache(path, **payload)
    return payload


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
            category="dynamics",
            target_order=0.0,
            hypothesis="动力学观测系 Gamma 的有效收敛阶至少应为正。",
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
            hypothesis="动力学观测系半径 R(t_obs) 的有效收敛阶至少应为正。",
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
            if np.isfinite(value_ref) and np.isfinite(value_test) and value_ref > 0.0:
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

    metrics = {
        "electron-spectrum": _estimate_order(
            _whole_shape_error(direct[num_values[0]], direct[num_values[1]]),
            _whole_shape_error(direct[num_values[1]], direct[num_values[2]]),
            ratio=ratio,
        ),
        "electron-peak-gamma": _estimate_order(
            _support_metric_error(direct[num_values[0]], direct[num_values[1]], "peak_gamma"),
            _support_metric_error(direct[num_values[1]], direct[num_values[2]], "peak_gamma"),
            ratio=ratio,
        ),
        "electron-support-low": _estimate_order(
            _support_metric_error(direct[num_values[0]], direct[num_values[1]], "g_lo"),
            _support_metric_error(direct[num_values[1]], direct[num_values[2]], "g_lo"),
            ratio=ratio,
        ),
        "electron-support-high": _estimate_order(
            _support_metric_error(direct[num_values[0]], direct[num_values[1]], "g_hi"),
            _support_metric_error(direct[num_values[1]], direct[num_values[2]], "g_hi"),
            ratio=ratio,
        ),
        "electron-shape-aligned": _estimate_order(
            _aligned_shape_metric(direct[num_values[0]], direct[num_values[1]]),
            _aligned_shape_metric(direct[num_values[1]], direct[num_values[2]]),
            ratio=ratio,
        ),
    }
    hypotheses = {
        "electron-spectrum": "整条电子谱的混合有效阶数，包含峰位漂移、支撑区漂移和谱形误差。",
        "electron-peak-gamma": "谱峰位置误差的有效阶数，直接反映电子谱主峰跟踪精度。",
        "electron-support-low": "低能支撑边界 g_lo 误差的有效阶数。",
        "electron-support-high": "高能支撑边界 g_hi 误差的有效阶数。",
        "electron-shape-aligned": "对齐峰位后的谱形误差有效阶数，用于排除纯平移误差。",
    }
    return [
        OrderResult(
            name=f"{solver}-{metric}",
            category="electron",
            target_order=ELECTRON_ORDER_TARGET,
            hypothesis=hypotheses[metric],
            refinement_ratio=ratio,
            coarse_error=float("nan"),
            medium_error=float("nan"),
            measured_order=order_value,
            passed=np.isfinite(order_value) and order_value > ELECTRON_ORDER_TARGET,
            extra={"num_gam_e_num_nu": num_values, "shell_fractions": shell_fracs},
        )
        for metric, order_value in metrics.items()
    ]


def _solver_radiation_order(solver: str) -> list[OrderResult]:
    num_values = _resolution_triplet(solver)
    cases = {n: _radiation_case(solver, n) for n in num_values}
    bands = {n: _band_case(solver, n) for n in num_values}
    ratio = 2.0
    metrics = {
        "radiation-sync-spectrum": _estimate_order(
            _matrix_spectrum_error(
                cases[num_values[1]]["seed_frequency_hz"],
                cases[num_values[1]]["fwd_sync"],
                cases[num_values[0]]["seed_frequency_hz"],
                cases[num_values[0]]["fwd_sync"],
            ),
            _matrix_spectrum_error(
                cases[num_values[2]]["seed_frequency_hz"],
                cases[num_values[2]]["fwd_sync"],
                cases[num_values[1]]["seed_frequency_hz"],
                cases[num_values[1]]["fwd_sync"],
            ),
            ratio=ratio,
        ),
        "radiation-ssc-spectrum": _estimate_order(
            _matrix_spectrum_error(
                cases[num_values[1]]["seed_frequency_hz"],
                cases[num_values[1]]["fwd_ssc"],
                cases[num_values[0]]["seed_frequency_hz"],
                cases[num_values[0]]["fwd_ssc"],
            ),
            _matrix_spectrum_error(
                cases[num_values[2]]["seed_frequency_hz"],
                cases[num_values[2]]["fwd_ssc"],
                cases[num_values[1]]["seed_frequency_hz"],
                cases[num_values[1]]["fwd_ssc"],
            ),
            ratio=ratio,
        ),
        "radiation-observer-total": _estimate_order(
            _max_rel(cases[num_values[0]]["observer_total"], cases[num_values[1]]["observer_total"]),
            _max_rel(cases[num_values[1]]["observer_total"], cases[num_values[2]]["observer_total"]),
            ratio=ratio,
        ),
        "radiation-bands-flux": _estimate_order(
            _max_rel(bands[num_values[0]]["bands_flux"], bands[num_values[1]]["bands_flux"]),
            _max_rel(bands[num_values[1]]["bands_flux"], bands[num_values[2]]["bands_flux"]),
            ratio=ratio,
        ),
        "radiation-nu_a": _estimate_order(
            _max_rel(bands[num_values[0]]["nu_a"], bands[num_values[1]]["nu_a"]),
            _max_rel(bands[num_values[1]]["nu_a"], bands[num_values[2]]["nu_a"]),
            ratio=ratio,
        ),
    }
    hypotheses = {
        "radiation-sync-spectrum": "前向 synchrotron 共动谱的有效阶数。",
        "radiation-ssc-spectrum": "前向 SSC 共动谱的有效阶数。",
        "radiation-observer-total": "observer-side 总观测矩阵的有效阶数，覆盖插值与 EATS 投影。",
        "radiation-bands-flux": "多波段光变聚合后的有效阶数。",
        "radiation-nu_a": "SSA 派生特征频率 nu_a 的有效阶数。",
    }
    return [
        OrderResult(
            name=f"{solver}-{metric}",
            category="radiation",
            target_order=RADIATION_ORDER_TARGET,
            hypothesis=hypotheses[metric],
            refinement_ratio=ratio,
            coarse_error=float("nan"),
            medium_error=float("nan"),
            measured_order=order_value,
            passed=np.isfinite(order_value) and order_value > RADIATION_ORDER_TARGET,
            extra={"num_gam_e_num_nu": num_values, "observer_num_nu": int(OBSERVER_FREQUENCIES_HZ.size)},
        )
        for metric, order_value in metrics.items()
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
    for solver in ("fullhide", "slc1", "slc1_mmg2"):
        checks.extend(_solver_spectrum_order(solver))
        checks.extend(_solver_radiation_order(solver))

    summary = _build_summary(checks)
    payload = {
        "targets": {
            "electron": ELECTRON_ORDER_TARGET,
            "radiation": RADIATION_ORDER_TARGET,
        },
        "summary": summary,
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
