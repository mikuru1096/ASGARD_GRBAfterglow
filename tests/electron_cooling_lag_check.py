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
from asgard_runtime import solve_dynamics, solve_electron
from asgard_setup import build_simulation_setup
import src.Electron.FS_electron_fullhide as fullhide_module
from src import constants


OUTPUT_JSON = ROOT / "output" / "vegasafterglow_doc" / "electron_cooling_lag.json"


get_y = fullhide_module.get_y
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
    }


def _rel_diff(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = np.maximum(np.abs(b), 1.0e-99)
    return float(np.max(np.abs(a - b) / denom))


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
    eta_0 = float(boundary[0])
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
        "config": cfg,
        "setup": setup,
        "gamma": gamma,
        "d_x": d_x,
        "shell": shell,
        "r_prev": r_prev,
        "r_cur": r_cur,
        "r_gamma_loc": r_gamma_loc,
        "beta_gam": beta_gam,
        "epsilon_e": epsilon_e,
        "epsilon_b": epsilon_b,
        "p": p,
        "z": z,
        "d_ne_ism": d_ne_ism,
        "a_star": a_star,
        "f_e": f_e,
        "r_tr": r_tr,
        "f_jump": f_jump,
        "f_wide": f_wide,
        "r0": r0,
        "d_ne_shell": d_ne_shell,
        "gam_e_c_shell": gam_e_c,
        "dDR": dDR,
        "L1": L1,
        "dN_x_prev": dN_x_prev,
        "d_n_actual": d_n_actual,
        "dEL_mean_base": dEL_mean_base,
        "v_seed": np.asarray(setup.seed_frequency_hz, dtype=float),
        "is_uniform_density": bool(a_star <= 0.0 and f_jump == 1.0),
        "index_y": cfg.index_y,
        "index_syn_integr": cfg.index_syn_integr,
    }


def _step_source_and_fields(ctx: dict[str, object], r_loc: float, d_ne: float) -> tuple[float, float, float, np.ndarray]:
    gamma = ctx["gamma"]
    p = float(ctx["p"])
    epsilon_e = float(ctx["epsilon_e"])
    epsilon_b = float(ctx["epsilon_b"])
    f_e = float(ctx["f_e"])
    r_gamma_loc = float(ctx["r_gamma_loc"])
    dDR = float(ctx["dDR"])

    db = 0.39 * math.sqrt(epsilon_b * d_ne * (r_gamma_loc * (r_gamma_loc - 1.0)))
    gam_e_max = 3.0 * constants.para_m_energy / math.sqrt(8.0 * db * constants.para_e**3)
    temp_gam = epsilon_e / f_e * constants.para_m_p_div_m_e * (r_gamma_loc - 1.0)
    gam_e_m = float(ecommon.electron_gamma_m_exact(p, temp_gam, gam_e_max))
    gam_e_m_p = (1.0 - p) / (gam_e_max ** (1.0 - p) - gam_e_m ** (1.0 - p))
    q = float(ecommon.electron_injection_prefactor(r_loc, dDR, d_ne, f_e, gam_e_m_p))
    dF1 = np.asarray(ecommon.electron_build_source_term(gamma, gam_e_m, gam_e_max, q, p), dtype=float)
    return db, gam_e_max, gam_e_m, dF1


def _run_shell_step(ctx: dict[str, object], refresh_cooling: bool) -> dict[str, object]:
    gamma = ctx["gamma"]
    dN_x = np.asarray(ctx["dN_x_prev"], dtype=float).copy()
    r_loc = float(ctx["r_prev"])
    d_x = float(ctx["d_x"])
    dDR = float(ctx["dDR"])
    d_ne_shell = float(ctx["d_ne_shell"])
    dEL_mean_base = np.asarray(ctx["dEL_mean_base"], dtype=float)
    d_ne_ism = float(ctx["d_ne_ism"])
    a_star = float(ctx["a_star"])
    r0 = float(ctx["r0"])
    r_tr = float(ctx["r_tr"])
    f_jump = float(ctx["f_jump"])
    f_wide = float(ctx["f_wide"])
    epsilon_e = float(ctx["epsilon_e"])
    epsilon_b = float(ctx["epsilon_b"])
    p = float(ctx["p"])
    gam_e_c_shell = float(ctx["gam_e_c_shell"])
    r_gamma_loc = float(ctx["r_gamma_loc"])
    beta_gam = float(ctx["beta_gam"])
    v_seed = np.asarray(ctx["v_seed"], dtype=float)
    index_y = int(ctx["index_y"])
    index_syn_integr = int(ctx["index_syn_integr"])
    is_uniform_density = bool(ctx["is_uniform_density"])
    num_threads = 1

    for _ in range(int(ctx["L1"])):
        r_loc = r_loc + dDR
        if is_uniform_density:
            d_ne = d_ne_shell
        else:
            d_ne = float(ecommon.electron_external_density(a_star, d_ne_ism, r_loc, r0, r_tr, f_jump, f_wide, 1))

        db, gam_e_max, gam_e_m, dF1 = _step_source_and_fields(ctx, r_loc, d_ne)

        if refresh_cooling:
            d_n_current = dN_x / gamma / math.log(10.0)
            p_syn, seed_syn = get_y.get_syn_selected(
                index_syn_integr,
                r_loc,
                db,
                num_threads,
                gamma,
                d_n_current,
                v_seed,
            )
            dEl = get_y.get_forward_cooling(
                index_y,
                epsilon_e,
                epsilon_b,
                p,
                db,
                gam_e_m,
                gam_e_c_shell,
                gam_e_max,
                r_loc,
                r_gamma_loc,
                beta_gam,
                d_ne,
                num_threads,
                gamma,
                v_seed,
                p_syn,
                seed_syn,
            )
            dEL_mean = np.asarray(ecommon.electron_loss_mean(dEl), dtype=float)
        else:
            if d_ne_shell > 0.0:
                dEL_mean = dEL_mean_base * (d_ne / d_ne_shell)
            else:
                dEL_mean = dEL_mean_base

        dN_x = np.asarray(ecommon.electron_fullhide_step(r_loc, dDR, d_x, dEL_mean, dF1, dN_x), dtype=float)

    d_n_out = dN_x / gamma / math.log(10.0)
    metrics = _support_metrics(gamma, d_n_out)
    metrics["spectrum"] = d_n_out
    return metrics


def main() -> None:
    shell_fracs = [0.3, 0.5]
    n_values = [61, 121, 241]
    payload: dict[str, object] = {"shells": []}

    for shell_frac in shell_fracs:
        shell_block = {"shell_frac": shell_frac, "grid": {}}
        for n in n_values:
            ctx = _shell_context(n, shell_frac)
            frozen = _run_shell_step(ctx, refresh_cooling=False)
            refreshed = _run_shell_step(ctx, refresh_cooling=True)
            actual = np.asarray(ctx["d_n_actual"], dtype=float)
            gamma = np.asarray(ctx["gamma"], dtype=float)
            shell_block["grid"][str(n)] = {
                "shell_index": int(ctx["shell"]),
                "L1": int(ctx["L1"]),
                "frozen_matches_solver": _rel_diff(frozen["spectrum"], actual),
                "refresh_vs_frozen": _rel_diff(refreshed["spectrum"], frozen["spectrum"]),
                "refresh_vs_solver": _rel_diff(refreshed["spectrum"], actual),
                "frozen": {
                    "peak_gamma": frozen["peak_gamma"],
                    "g_lo": frozen["g_lo"],
                    "g_hi": frozen["g_hi"],
                    "peak_to_gamma_grid_index": int(np.argmax(frozen["spectrum"])),
                },
                "refreshed": {
                    "peak_gamma": refreshed["peak_gamma"],
                    "g_lo": refreshed["g_lo"],
                    "g_hi": refreshed["g_hi"],
                    "peak_to_gamma_grid_index": int(np.argmax(refreshed["spectrum"])),
                },
                "solver": {
                    **_support_metrics(gamma, actual),
                    "peak_to_gamma_grid_index": int(np.argmax(actual)),
                },
            }
        payload["shells"].append(shell_block)

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
