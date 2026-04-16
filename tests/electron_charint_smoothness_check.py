from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_paths import asgard_doc_path
from asgard_presets import build_baseline_config
from asgard_runtime import solve_dynamics, solve_electron
from asgard_setup import build_simulation_setup


OUTPUT_JSON = asgard_doc_path("electron_charint_smoothness.json")
MAX_LOG_JUMP = 0.25
SAWTOOTH_TOL_DEX = 0.03


def _support_edge(gamma: np.ndarray, spectrum: np.ndarray, threshold: float, idx_inside: int) -> float:
    gamma = np.asarray(gamma, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    idx0 = max(idx_inside - 1, 0)
    idx1 = idx_inside
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


def _support_metrics(gamma: np.ndarray, spectrum: np.ndarray) -> tuple[float, float]:
    gamma = np.asarray(gamma, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    peak_idx = int(np.argmax(spectrum))
    peak_gamma = float(gamma[peak_idx])
    peak_value = float(spectrum[peak_idx])
    if not np.isfinite(peak_value) or peak_value <= 0.0:
        return float("nan"), float("nan")
    mask = spectrum > 1.0e-8 * peak_value
    if not np.any(mask):
        return peak_gamma, float("nan")
    idx = int(np.where(mask)[0][0])
    return peak_gamma, _support_edge(gamma, spectrum, 1.0e-8 * peak_value, idx)


def _metric_report(values: np.ndarray) -> dict[str, float | bool]:
    values = np.asarray(values, dtype=float)
    ok = np.isfinite(values) & (values > 0.0)
    report: dict[str, float | bool] = {
        "count": int(values.size),
        "positive_finite": bool(np.all(ok)),
        "min": float(np.min(values)) if values.size else float("nan"),
        "max": float(np.max(values)) if values.size else float("nan"),
        "max_adjacent_log_jump": float("nan"),
        "has_sawtooth": True,
    }
    if values.size < 3 or not np.all(ok):
        return report
    logv = np.log10(values)
    diffs = np.diff(logv)
    report["max_adjacent_log_jump"] = float(np.max(np.abs(diffs)))
    strong = np.abs(diffs) > SAWTOOTH_TOL_DEX
    alternating = (diffs[1:] * diffs[:-1] < 0.0) & strong[1:] & strong[:-1]
    report["has_sawtooth"] = bool(np.any(alternating))
    return report


def _passed(report: dict[str, float | bool]) -> bool:
    jump = report["max_adjacent_log_jump"]
    return bool(
        report["positive_finite"]
        and np.isfinite(jump)
        and float(jump) < MAX_LOG_JUMP
        and not report["has_sawtooth"]
    )


def main() -> None:
    payload: dict[str, object] = {
        "solver": "charint",
        "max_log_jump_dex": MAX_LOG_JUMP,
        "sawtooth_tol_dex": SAWTOOTH_TOL_DEX,
        "index_y": {},
    }
    failed = 0

    for index_y in (1, 2, 3):
        cfg = build_baseline_config(
            electron_solver="charint",
            num_gam_e=81,
            num_nu=81,
            num_r=160,
            include_forward_ssc=False,
            index_y=index_y,
        )
        setup = build_simulation_setup(cfg)
        dyn = solve_dynamics(setup.boundary, cfg)
        ele = solve_electron(setup.boundary, dyn, setup.seed_frequency_hz, cfg)

        shell_ids = np.arange(1, cfg.num_r - 1, dtype=int)
        peak_gamma = np.zeros(shell_ids.size, dtype=float)
        g_lo = np.zeros(shell_ids.size, dtype=float)
        for idx, shell in enumerate(shell_ids):
            peak_gamma[idx], g_lo[idx] = _support_metrics(ele.gam_e, ele.d_n_gam_e[:, shell])
        nu_a = np.asarray(ele.nu_a[shell_ids], dtype=float)

        peak_report = _metric_report(peak_gamma)
        g_lo_report = _metric_report(g_lo)
        nu_a_report = _metric_report(nu_a)
        passed = _passed(peak_report) and _passed(g_lo_report) and _passed(nu_a_report)
        if not passed:
            failed += 1

        payload["index_y"][str(index_y)] = {
            "passed": bool(passed),
            "peak_gamma": peak_report,
            "g_lo": g_lo_report,
            "nu_a": nu_a_report,
        }

    payload["failed"] = failed
    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if failed != 0:
        raise SystemExit(1)
    print("PASS: charint smoothness check succeeded.")


if __name__ == "__main__":
    main()
