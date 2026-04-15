from __future__ import annotations

import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import src.Electron.FS_electron_fullhide as fullhide_module


OUTPUT_JSON = ROOT / "output" / "asgard_doc" / "electron_operator_scheme.json"

ecommon = fullhide_module.electron_common


@dataclass
class SchemeOrderResult:
    name: str
    refinement_ratio: float
    coarse_error: float
    medium_error: float
    measured_order: float
    passed: bool
    extra: dict


def _estimate_order(err_coarse: float, err_medium: float, ratio: float) -> float:
    if err_coarse == 0.0 and err_medium == 0.0:
        return float("inf")
    if err_medium == 0.0 and err_coarse > 0.0:
        return float("inf")
    if err_coarse <= 0.0 or err_medium <= 0.0:
        return float("nan")
    return float(np.log(err_coarse / err_medium) / np.log(ratio))


def _grid(num_gam: int, x_lo: float = 0.0, x_hi: float = 6.0) -> tuple[np.ndarray, np.ndarray]:
    x = np.linspace(x_lo, x_hi, num_gam)
    dx = float(x[1] - x[0])
    x_edge = np.empty(num_gam + 1, dtype=float)
    x_edge[1:-1] = 0.5 * (x[:-1] + x[1:])
    x_edge[0] = x[0] - 0.5 * dx
    x_edge[-1] = x[-1] + 0.5 * dx
    gamma = 10.0**x
    return x_edge, gamma


def _gaussian_cell_average(x_edge: np.ndarray, x0: float, sigma: float) -> np.ndarray:
    out = np.empty(x_edge.size - 1, dtype=float)
    root2 = math.sqrt(2.0)
    norm = sigma * math.sqrt(math.pi / 2.0)
    for i in range(out.size):
        xa = (x_edge[i] - x0) / (root2 * sigma)
        xb = (x_edge[i + 1] - x0) / (root2 * sigma)
        integ = norm * (math.erf(xb) - math.erf(xa))
        out[i] = integ / (x_edge[i + 1] - x_edge[i])
    return out


def _run_fullhide_operator(d_nx_init: np.ndarray, dx: float, lam: float, nsteps: int) -> np.ndarray:
    d_nx = np.asarray(d_nx_init, dtype=float).copy()
    r_loc = 1.0e30
    ddr = 1.0
    dEL_mean = np.full(d_nx.size - 1, lam * dx / ddr, dtype=float)
    zero_src = np.zeros_like(d_nx)
    for _ in range(nsteps):
        d_nx = np.asarray(ecommon.electron_fullhide_step(r_loc, ddr, dx, dEL_mean, zero_src, d_nx), dtype=float)
    return d_nx


def _run_t2g1_operator(d_nx_init: np.ndarray, dx: float, lam: float, nsteps: int) -> np.ndarray:
    d_nx = np.asarray(d_nx_init, dtype=float).copy()
    d_nx_prev = d_nx.copy()
    r_loc = 1.0e30
    ddr = 1.0
    dEL_mean = np.full(d_nx.size - 1, lam * dx / ddr, dtype=float)
    for istep in range(nsteps):
        up = -(ddr / dx) * (dEL_mean + 1.0 / r_loc / math.log(10.0))
        if istep < 2:
            principal, temp1 = ecommon.electron_prepare_implicit_coeffs(1.0, up)
            temp2 = d_nx / principal
        else:
            principal, temp1 = ecommon.electron_prepare_implicit_coeffs(1.5, up)
            temp2 = (2.0 * d_nx - 0.5 * d_nx_prev) / principal
        x = np.asarray(ecommon.electron_backward_sweep(temp1, temp2), dtype=float)
        d_nx_prev = d_nx
        d_nx = x
    return d_nx


def _scheme_error(num_gam: int, solver: str, lam: float, shift_total: float, x0: float, sigma: float) -> tuple[float, dict]:
    x_edge, gamma = _grid(num_gam)
    x_center = np.log10(gamma)
    dx = float(x_center[1] - x_center[0])
    d_nx_init = _gaussian_cell_average(x_edge, x0, sigma)
    nsteps = max(1, int(round(shift_total / (lam * dx))))
    x_target = x0 - lam * nsteps * dx
    d_nx_exact = _gaussian_cell_average(x_edge, x_target, sigma)

    if solver == "fullhide":
        d_nx_num = _run_fullhide_operator(d_nx_init, dx, lam, nsteps)
    elif solver == "t2g1":
        d_nx_num = _run_t2g1_operator(d_nx_init, dx, lam, nsteps)
    else:
        raise ValueError(solver)

    denom = float(np.trapezoid(np.abs(d_nx_exact), x_center))
    err = float(np.trapezoid(np.abs(d_nx_num - d_nx_exact), x_center) / denom)
    return err, {
        "nsteps": nsteps,
        "peak_init_gamma": float(gamma[int(np.argmax(d_nx_init))]),
        "peak_exact_gamma": float(gamma[int(np.argmax(d_nx_exact))]),
        "peak_num_gamma": float(gamma[int(np.argmax(d_nx_num))]),
    }


def main() -> None:
    num_values = [61, 121, 241]
    lam = 0.5
    shift_total = 0.8
    x0 = 3.2
    sigma = 0.08
    payload: dict[str, object] = {"results": [], "cases": {}}

    for solver in ("fullhide", "t2g1"):
        errs = {}
        extras = {}
        for n in num_values:
            err, extra = _scheme_error(n, solver, lam, shift_total, x0, sigma)
            errs[n] = err
            extras[n] = extra
        ratio = (num_values[1] - 1) / (num_values[0] - 1)
        order = _estimate_order(errs[61], errs[121], ratio)
        result = SchemeOrderResult(
            name=f"{solver}-operator-only",
            refinement_ratio=float(ratio),
            coarse_error=errs[61],
            medium_error=errs[121],
            measured_order=order,
            passed=bool(np.isfinite(order) and order > 0.8),
            extra={"num_gam_e": num_values, "lam": lam, "shift_total": shift_total, "peak_trace": extras},
        )
        payload["results"].append(asdict(result))
        payload["cases"][solver] = {"errors": errs, "peak_trace": extras}

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
