from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import src.Electron.FS_electron_charint as charint_module


OUTPUT_JSON = ROOT / "output" / "asgard_doc" / "electron_charint_exact.json"
REL_TOL = 1.0e-12
ABS_TOL = 1.0e-12


def _max_rel(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = np.maximum(np.abs(b), 1.0e-99)
    return float(np.max(np.abs(a - b) / denom))


def _run_case(name: str, d_n_x_in: np.ndarray, d_f1: np.ndarray, ddr: float, a_u: float, b_u: float) -> dict[str, float | bool]:
    x_edge = np.linspace(0.0, 1.0, d_n_x_in.size + 1, dtype=float)
    d_n_x_out = np.asarray(
        charint_module.fs_electron_charint_affine_step_test(
            ddr,
            x_edge,
            float(a_u),
            float(b_u),
            np.asarray(d_f1, dtype=float),
            np.asarray(d_n_x_in, dtype=float),
        ),
        dtype=float,
    )
    expected = np.asarray(d_n_x_in, dtype=float) + float(ddr) * np.asarray(d_f1, dtype=float)
    dx = np.diff(x_edge)
    particle_out = float(np.sum(d_n_x_out * dx))
    particle_expected = float(np.sum(expected * dx))
    rel_err = _max_rel(d_n_x_out, expected)
    abs_err = float(np.max(np.abs(d_n_x_out - expected)))
    particle_err = 0.0 if particle_expected == 0.0 else abs(particle_out - particle_expected) / abs(particle_expected)
    passed = (
        np.all(np.isfinite(d_n_x_out))
        and np.all(d_n_x_out >= 0.0)
        and rel_err <= REL_TOL
        and abs_err <= ABS_TOL
        and particle_err <= REL_TOL
    )
    return {
        "name": name,
        "ddr": float(ddr),
        "a_u": float(a_u),
        "b_u": float(b_u),
        "rel_error": rel_err,
        "abs_error": abs_err,
        "particle_rel_error": particle_err,
        "passed": bool(passed),
    }


def main() -> None:
    num = 24
    ddr = 2.5e-1
    cases = [
        _run_case("transport-only-static", np.full(num, 3.25), np.zeros(num), ddr, 0.0, 0.0),
        _run_case("source-only-static", np.zeros(num), np.full(num, 1.75), ddr, 0.0, 0.0),
        _run_case("transport-plus-source-static", np.full(num, 2.0), np.full(num, 0.5), ddr, 0.0, 0.0),
    ]
    payload = {
        "target_rel_tol": REL_TOL,
        "target_abs_tol": ABS_TOL,
        "cases": cases,
        "failed": int(sum(not item["passed"] for item in cases)),
    }
    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if payload["failed"] != 0:
        raise SystemExit(1)
    print("PASS: charint affine exact check succeeded.")


if __name__ == "__main__":
    main()
