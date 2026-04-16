from __future__ import annotations

import json
import math
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
import src.Electron.FS_electron_fullhide as fullhide_module
from src import constants


OUTPUT_JSON = asgard_doc_path("electron_clipping.json")


ecommon = fullhide_module.electron_common
get_y = fullhide_module.get_y


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
    }


def _rel_diff(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = np.maximum(np.abs(b), 1.0e-99)
    return float(np.max(np.abs(a - b) / denom))


def _python_backward_no_clip(temp1: np.ndarray, temp2: np.ndarray) -> np.ndarray:
    x = np.empty_like(temp2, dtype=float)
    x[-1] = temp2[-1]
    for i in range(temp2.size - 2, -1, -1):
        x[i] = temp2[i] - temp1[i] * x[i + 1]
    return x


def _python_fullhide_step_no_clip(
    r_loc: float,
    ddr: float,
    d_x: float,
    dEL_mean: np.ndarray,
    dF1: np.ndarray,
    dN_x_in: np.ndarray,
) -> np.ndarray:
    temp3 = dEL_mean + 1.0 / r_loc / math.log(10.0)
    up = -(ddr / d_x) * temp3
    principal, temp1 = ecommon.electron_prepare_implicit_coeffs(1.0, up)
    temp2 = dN_x_in / principal + ddr * dF1 / principal
    return _python_backward_no_clip(np.asarray(temp1, dtype=float), np.asarray(temp2, dtype=float))


def _shell_context(n_gam: int, shell_frac: float) -> dict[str, object]:
    cfg = build_baseline_config(
        electron_solver="fullhide",
        num_gam_e=n_gam,
        num_nu=121,
        num_r=160,
        num_threads=1,
    )
    setup = build_simulation_setup(cfg)
    dyn = solve_dynamics(setup.boundary, cfg)
    ele = solve_electron(setup.boundary, dyn, setup.seed_frequency_hz, cfg)

    shell = min(max(int(round(shell_frac * (cfg.num_r - 1))), 2), cfg.num_r - 1)
    gamma = np.asarray(ele.gam_e, dtype=float)
    d_x = float(np.log10(gamma[1] / gamma[0]))
    d_n_prev = np.asarray(ele.d_n_gam_e[:, shell - 1], dtype=float)
    d_n_actual = np.asarray(ele.d_n_gam_e[:, shell], dtype=float)
    dN_x_prev = d_n_prev * gamma * math.log(10.0)

    r_prev = float(dyn.radius[shell - 1])
    r_cur = float(dyn.radius[shell])
    r_gamma_loc = 0.5 * float(dyn.r_gamma[shell - 1] + dyn.r_gamma[shell])
    beta_gam = math.sqrt(1.0 - 1.0 / r_gamma_loc**2)

    boundary = np.asarray(setup.boundary, dtype=float)
    epsilon_e = float(boundary[4])
    epsilon_b = float(boundary[5])
    p = float(boundary[6])
    z = float(boundary[7])
    d_ne_ism = float(boundary[10])
    a_star = float(boundary[11])
    f_e = float(boundary[15])
    r_tr = float(boundary[20])
    f_jump = float(boundary[21])
    f_wide = float(boundary[22])
    r0 = float(boundary[-1])

    d_ne_shell = float(ecommon.electron_external_density(a_star, d_ne_ism, r_prev, r0, r_tr, f_jump, f_wide, 1))
    db = 0.39 * math.sqrt(epsilon_b * d_ne_shell * (r_gamma_loc * (r_gamma_loc - 1.0)))
    gam_e_max = 3.0 * constants.para_m_energy / math.sqrt(8.0 * db * constants.para_e**3)
    temp_gam = epsilon_e / f_e * constants.para_m_p_div_m_e * (r_gamma_loc - 1.0)
    gam_e_m = float(ecommon.electron_gamma_m_exact(p, temp_gam, gam_e_max))
    gam_e_c = 7.7e8 * (1.0 + z) / r_gamma_loc / db**2 / float(dyn.r_tobs[shell])

    f_r = (1.35e-19) / beta_gam / r_gamma_loc * db**2 / math.pi
    dDR = 0.1 / (f_r * gam_e_max + 1.333 / (r_cur + r_prev))
    dDD = r_cur - r_prev
    L1 = max(100, min(1000, int(dDD / dDR)))
    dDR = dDD / L1

    p_syn, seed_syn = get_y.get_syn_selected(
        cfg.index_syn_integr,
        r_prev,
        db,
        cfg.num_threads,
        gamma,
        d_n_prev,
        np.asarray(setup.seed_frequency_hz, dtype=float),
    )
    dEl = get_y.get_forward_cooling(
        cfg.index_y,
        epsilon_e,
        epsilon_b,
        p,
        db,
        gam_e_m,
        gam_e_c,
        gam_e_max,
        r_prev,
        r_gamma_loc,
        beta_gam,
        d_ne_shell,
        cfg.num_threads,
        gamma,
        np.asarray(setup.seed_frequency_hz, dtype=float),
        p_syn,
        seed_syn,
    )
    dEL_mean_base = np.asarray(ecommon.electron_loss_mean(dEl), dtype=float)

    return {
        "cfg": cfg,
        "gamma": gamma,
        "d_x": d_x,
        "shell": shell,
        "r_prev": r_prev,
        "r_cur": r_cur,
        "r_gamma_loc": r_gamma_loc,
        "d_ne_shell": d_ne_shell,
        "a_star": a_star,
        "d_ne_ism": d_ne_ism,
        "r0": r0,
        "r_tr": r_tr,
        "f_jump": f_jump,
        "f_wide": f_wide,
        "epsilon_e": epsilon_e,
        "epsilon_b": epsilon_b,
        "f_e": f_e,
        "p": p,
        "dDR": dDR,
        "L1": L1,
        "dN_x_prev": dN_x_prev,
        "d_n_actual": d_n_actual,
        "dEL_mean_base": dEL_mean_base,
        "is_uniform_density": bool(a_star <= 0.0 and f_jump == 1.0),
    }


def _step_source(ctx: dict[str, object], r_loc: float, d_ne: float) -> tuple[np.ndarray, np.ndarray]:
    gamma = ctx["gamma"]
    epsilon_e = float(ctx["epsilon_e"])
    epsilon_b = float(ctx["epsilon_b"])
    f_e = float(ctx["f_e"])
    p = float(ctx["p"])
    r_gamma_loc = float(ctx["r_gamma_loc"])
    dDR = float(ctx["dDR"])
    db = 0.39 * math.sqrt(epsilon_b * d_ne * (r_gamma_loc * (r_gamma_loc - 1.0)))
    gam_e_max = 3.0 * constants.para_m_energy / math.sqrt(8.0 * db * constants.para_e**3)
    temp_gam = epsilon_e / f_e * constants.para_m_p_div_m_e * (r_gamma_loc - 1.0)
    gam_e_m = float(ecommon.electron_gamma_m_exact(p, temp_gam, gam_e_max))
    gam_e_m_p = (1.0 - p) / (gam_e_max ** (1.0 - p) - gam_e_m ** (1.0 - p))
    q = float(ecommon.electron_injection_prefactor(r_loc, dDR, d_ne, f_e, gam_e_m_p))
    dF1 = np.asarray(ecommon.electron_build_source_term(gamma, gam_e_m, gam_e_max, q, p), dtype=float)
    if float(ctx["d_ne_shell"]) > 0.0:
        dEL_mean = np.asarray(ctx["dEL_mean_base"], dtype=float) * (d_ne / float(ctx["d_ne_shell"]))
    else:
        dEL_mean = np.asarray(ctx["dEL_mean_base"], dtype=float)
    return dEL_mean, dF1


def _run_shell(ctx: dict[str, object], use_clip: bool) -> np.ndarray:
    gamma = np.asarray(ctx["gamma"], dtype=float)
    dN_x = np.asarray(ctx["dN_x_prev"], dtype=float).copy()
    r_loc = float(ctx["r_prev"])
    for _ in range(int(ctx["L1"])):
        r_loc += float(ctx["dDR"])
        if bool(ctx["is_uniform_density"]):
            d_ne = float(ctx["d_ne_shell"])
        else:
            d_ne = float(
                ecommon.electron_external_density(
                    float(ctx["a_star"]),
                    float(ctx["d_ne_ism"]),
                    r_loc,
                    float(ctx["r0"]),
                    float(ctx["r_tr"]),
                    float(ctx["f_jump"]),
                    float(ctx["f_wide"]),
                    1,
                )
            )
        dEL_mean, dF1 = _step_source(ctx, r_loc, d_ne)
        if use_clip:
            dN_x = np.asarray(
                ecommon.electron_fullhide_step(r_loc, float(ctx["dDR"]), float(ctx["d_x"]), dEL_mean, dF1, dN_x),
                dtype=float,
            )
        else:
            dN_x = _python_fullhide_step_no_clip(r_loc, float(ctx["dDR"]), float(ctx["d_x"]), dEL_mean, dF1, dN_x)
    return dN_x / gamma / math.log(10.0)


def main() -> None:
    shell_fracs = [0.3, 0.5]
    n_values = [61, 121, 241]
    payload: dict[str, object] = {"shells": []}

    for shell_frac in shell_fracs:
        shell_block = {"shell_frac": shell_frac, "grid": {}}
        for n in n_values:
            ctx = _shell_context(n, shell_frac)
            clipped = _run_shell(ctx, use_clip=True)
            unclipped = _run_shell(ctx, use_clip=False)
            actual = np.asarray(ctx["d_n_actual"], dtype=float)
            negatives = int(np.count_nonzero(unclipped < 0.0))

            shell_block["grid"][str(n)] = {
                "shell_index": int(ctx["shell"]),
                "L1": int(ctx["L1"]),
                "clipped_matches_solver": _rel_diff(clipped, actual),
                "unclipped_vs_clipped": _rel_diff(unclipped, clipped),
                "negative_cells_unclipped": negatives,
                "clipped": _support_metrics(np.asarray(ctx["gamma"], dtype=float), clipped),
                "unclipped": _support_metrics(np.asarray(ctx["gamma"], dtype=float), np.maximum(unclipped, 0.0)),
                "solver": _support_metrics(np.asarray(ctx["gamma"], dtype=float), actual),
            }
        payload["shells"].append(shell_block)

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
