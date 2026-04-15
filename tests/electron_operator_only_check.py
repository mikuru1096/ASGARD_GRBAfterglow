from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_presets import build_baseline_config
from asgard_runtime import _minimum_electron_lorentz_factor, solve_dynamics, solve_electron
from asgard_setup import build_simulation_setup
import src.Electron.FS_electron_fullhide as fullhide_module


OUTPUT_JSON = ROOT / "output" / "asgard_doc" / "electron_operator_only.json"


ecommon = fullhide_module.electron_common


def _support_metrics(gamma: np.ndarray, spectrum: np.ndarray) -> dict[str, float]:
    peak_idx = int(np.argmax(spectrum))
    peak_gamma = float(gamma[peak_idx])
    peak_value = float(spectrum[peak_idx])
    mask = spectrum > 1.0e-8 * peak_value
    if np.any(mask):
        idx = np.where(mask)[0]
        g_lo = float(gamma[idx[0]])
        g_hi = float(gamma[idx[-1]])
    else:
        g_lo = float("nan")
        g_hi = float("nan")
    return {
        "peak_gamma": peak_gamma,
        "peak_value": peak_value,
        "g_lo": g_lo,
        "g_hi": g_hi,
        "peak_index": peak_idx,
    }


def _log_edges(gamma: np.ndarray) -> np.ndarray:
    x = np.log10(np.asarray(gamma, dtype=float))
    edges = np.empty(x.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = x[0] - 0.5 * (x[1] - x[0])
    edges[-1] = x[-1] + 0.5 * (x[-1] - x[-2])
    return edges


def _gaussian_cell_average(x_edges: np.ndarray, x0: float, sigma: float) -> np.ndarray:
    # cell-average of exp(-(x-x0)^2 / (2 sigma^2)) over each x-cell
    avg = np.empty(x_edges.size - 1, dtype=float)
    root2 = math.sqrt(2.0)
    norm = sigma * math.sqrt(math.pi / 2.0)
    for i in range(avg.size):
        xa = (x_edges[i] - x0) / (root2 * sigma)
        xb = (x_edges[i + 1] - x0) / (root2 * sigma)
        integ = norm * (math.erf(xb) - math.erf(xa))
        avg[i] = integ / (x_edges[i + 1] - x_edges[i])
    return avg


def _shell_context(solver: str, n_gam: int, shell_frac: float) -> dict[str, object]:
    cfg = build_baseline_config(
        electron_solver=solver,
        num_gam_e=n_gam,
        num_nu=n_gam,
        num_r=160,
        num_threads=1,
        include_forward_ssc=False,
        index_y=0,
    )
    setup = build_simulation_setup(cfg)
    dyn = solve_dynamics(setup.boundary, cfg)
    ele = solve_electron(setup.boundary, dyn, setup.seed_frequency_hz, cfg)

    shell = min(max(int(round(shell_frac * (cfg.num_r - 1))), 2), cfg.num_r - 1)
    gamma = np.asarray(ele.gam_e, dtype=float)
    d_x = float(np.log10(gamma[1] / gamma[0]))
    d_n_prev = np.asarray(ele.d_n_gam_e[:, shell - 1], dtype=float)
    dN_x_prev = d_n_prev * gamma * math.log(10.0)

    boundary = np.asarray(setup.boundary, dtype=float)
    r_prev = float(dyn.radius[shell - 1])
    r_cur = float(dyn.radius[shell])
    r_gamma_loc = 0.5 * float(dyn.r_gamma[shell - 1] + dyn.r_gamma[shell])

    a_star = float(boundary[11])
    d_ne_ism = float(boundary[10])
    r0 = float(boundary[-1])
    r_tr = float(boundary[20])
    f_jump = float(boundary[21])
    f_wide = float(boundary[22])

    d_ne_shell = float(ecommon.electron_external_density(a_star, d_ne_ism, r_prev, r0, r_tr, f_jump, f_wide, 1))
    db = 0.39 * math.sqrt(float(boundary[5]) * d_ne_shell * (r_gamma_loc * (r_gamma_loc - 1.0)))
    f_r = (1.35e-19) / math.sqrt(1.0 - 1.0 / r_gamma_loc**2) / r_gamma_loc * db**2 / math.pi
    gam_e_max = 3.0 * fullhide_module.constants.para_m_energy / math.sqrt(
        8.0 * db * fullhide_module.constants.para_e**3
    )
    dDR = 0.1 / (f_r * gam_e_max + 1.333 / (r_cur + r_prev))
    dDD = r_cur - r_prev
    L1 = max(100, min(1000, int(dDD / dDR)))
    dDR = dDD / L1

    gam_e_m = float(_minimum_electron_lorentz_factor(cfg, r_gamma_loc, gam_e_max))
    gam_e_c = 7.7e8 * (1.0 + float(boundary[7])) / r_gamma_loc / db**2 / float(dyn.r_tobs[shell])
    p_syn, seed_syn = fullhide_module.get_y.get_syn_selected(
        cfg.index_syn_integr,
        r_prev,
        db,
        1,
        gamma,
        d_n_prev,
        np.asarray(setup.seed_frequency_hz, dtype=float),
    )
    dEl = np.asarray(
        fullhide_module.get_y.get_forward_cooling(
            cfg.index_y,
            float(boundary[4]),
            float(boundary[5]),
            float(boundary[6]),
            db,
            gam_e_m,
            gam_e_c,
            gam_e_max,
            r_prev,
            r_gamma_loc,
            math.sqrt(1.0 - 1.0 / r_gamma_loc**2),
            d_ne_shell,
            1,
            gamma,
            np.asarray(setup.seed_frequency_hz, dtype=float),
            p_syn,
            seed_syn,
        ),
        dtype=float,
    )
    dEL_mean = np.asarray(ecommon.electron_loss_mean(dEl), dtype=float)

    return {
        "solver": solver,
        "gamma": gamma,
        "shell": shell,
        "d_x": d_x,
        "r_prev": r_prev,
        "dDR": dDR,
        "L1": L1,
        "dEL_mean": dEL_mean,
    }


def _run_fullhide_operator(ctx: dict[str, object], dN_x_init: np.ndarray) -> np.ndarray:
    dN_x = np.asarray(dN_x_init, dtype=float).copy()
    r_loc = float(ctx["r_prev"])
    for _ in range(int(ctx["L1"])):
        r_loc += float(ctx["dDR"])
        dN_x = np.asarray(
            ecommon.electron_fullhide_step(
                r_loc,
                float(ctx["dDR"]),
                float(ctx["d_x"]),
                np.asarray(ctx["dEL_mean"], dtype=float),
                np.zeros_like(dN_x),
                dN_x,
            ),
            dtype=float,
        )
    gamma = np.asarray(ctx["gamma"], dtype=float)
    return dN_x / gamma / math.log(10.0)


def _run_t2g1_operator(ctx: dict[str, object], dN_x_init: np.ndarray) -> np.ndarray:
    dN_x = np.asarray(dN_x_init, dtype=float).copy()
    dN_x_prev = dN_x.copy()
    r_loc = float(ctx["r_prev"])
    dEL_mean = np.asarray(ctx["dEL_mean"], dtype=float)
    up = -(float(ctx["dDR"]) / float(ctx["d_x"])) * (dEL_mean + 1.0 / r_loc / math.log(10.0))
    zero_src = np.zeros_like(dN_x)
    for l in range(int(ctx["L1"])):
        r_loc += float(ctx["dDR"])
        up = -(float(ctx["dDR"]) / float(ctx["d_x"])) * (dEL_mean + 1.0 / r_loc / math.log(10.0))
        if l < 2:
            principal, temp1 = ecommon.electron_prepare_implicit_coeffs(1.0, up)
            temp2 = dN_x / principal + float(ctx["dDR"]) * zero_src / principal
        else:
            principal, temp1 = ecommon.electron_prepare_implicit_coeffs(1.5, up)
            temp2 = (2.0 * dN_x - 0.5 * dN_x_prev) / principal
        x = np.asarray(ecommon.electron_backward_sweep(temp1, temp2), dtype=float)
        dN_x_prev = dN_x
        dN_x = x
    gamma = np.asarray(ctx["gamma"], dtype=float)
    return dN_x / gamma / math.log(10.0)


def main() -> None:
    shell_frac = 0.5
    n_values = [61, 121, 241]
    solvers = ["fullhide", "t2g1"]
    payload: dict[str, object] = {"shell_frac": shell_frac, "solvers": {}}

    for solver in solvers:
        solver_block = {}
        for n in n_values:
            ctx = _shell_context(solver, n, shell_frac)
            gamma = np.asarray(ctx["gamma"], dtype=float)
            x_edges = _log_edges(gamma)
            x0 = math.log10(6.055655981579369e4)
            sigma = 0.08
            dN_x0 = _gaussian_cell_average(x_edges, x0, sigma)
            if solver == "fullhide":
                d_n_out = _run_fullhide_operator(ctx, dN_x0)
            else:
                d_n_out = _run_t2g1_operator(ctx, dN_x0)
            solver_block[str(n)] = {
                "L1": int(ctx["L1"]),
                "initial": _support_metrics(gamma, dN_x0 / gamma / math.log(10.0)),
                "final": _support_metrics(gamma, d_n_out),
            }
        payload["solvers"][solver] = solver_block

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
