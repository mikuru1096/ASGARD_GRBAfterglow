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
from asgard_runtime import _ambient_density, _minimum_electron_lorentz_factor, solve_dynamics, solve_electron
from asgard_setup import build_simulation_setup
from src import constants


OUTPUT_JSON = asgard_doc_path("electron_initial_projection.json")


def _log_cell_edges(gamma: np.ndarray) -> np.ndarray:
    x = np.log10(np.asarray(gamma, dtype=float))
    edges = np.empty(x.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = x[0] - 0.5 * (x[1] - x[0])
    edges[-1] = x[-1] + 0.5 * (x[-1] - x[-2])
    return edges


def _initial_cell_center_spectrum(
    gamma: np.ndarray,
    p: float,
    gamma_m: float,
    gamma_c: float,
    gamma_max: float,
) -> np.ndarray:
    out = np.zeros_like(gamma, dtype=float)
    if gamma_m > gamma_c:
        q1 = gamma_c
        for i, g in enumerate(gamma):
            if gamma_c > g or gamma_max < g:
                continue
            if gamma_m > g:
                out[i] = q1 * g ** (-2.0)
            else:
                out[i] = q1 * gamma_m ** (p - 1.0) * g ** (-(p + 1.0))
    else:
        q1 = gamma_m ** (p - 1.0)
        for i, g in enumerate(gamma):
            if gamma_m > g or gamma_max < g:
                continue
            if gamma_c > g:
                out[i] = q1 * g ** (-p)
            else:
                out[i] = q1 * gamma_c * g ** (-(p + 1.0))
    return out


def _integrate_powerlaw_dnx(a: float, slope: float, x_lo: float, x_hi: float) -> float:
    if x_hi <= x_lo:
        return 0.0
    expo = 1.0 - slope
    if abs(expo) < 1.0e-14:
        return a * math.log(10.0) * (x_hi - x_lo)
    return a * (10.0 ** (expo * x_hi) - 10.0 ** (expo * x_lo)) / expo


def _fill_dnx_interval(
    gamma: np.ndarray,
    x_edges: np.ndarray,
    a: float,
    slope: float,
    gamma_lo: float,
    gamma_hi: float,
    out_dnx: np.ndarray,
) -> None:
    if gamma_hi <= gamma_lo:
        return
    x_lo = math.log10(gamma_lo)
    x_hi = math.log10(gamma_hi)
    for i in range(gamma.size):
        cell_lo = max(x_edges[i], x_lo)
        cell_hi = min(x_edges[i + 1], x_hi)
        if cell_hi > cell_lo:
            dx_cell = x_edges[i + 1] - x_edges[i]
            out_dnx[i] += _integrate_powerlaw_dnx(a, slope, cell_lo, cell_hi) / dx_cell


def _initial_cell_integrated_spectrum(
    gamma: np.ndarray,
    p: float,
    gamma_m: float,
    gamma_c: float,
    gamma_max: float,
) -> np.ndarray:
    x_edges = _log_cell_edges(gamma)
    dnx = np.zeros_like(gamma, dtype=float)
    if gamma_m > gamma_c:
        q1 = gamma_c
        _fill_dnx_interval(gamma, x_edges, q1, 2.0, gamma_c, min(gamma_m, gamma_max), dnx)
        _fill_dnx_interval(
            gamma,
            x_edges,
            q1 * gamma_m ** (p - 1.0),
            p + 1.0,
            max(gamma_c, gamma_m),
            gamma_max,
            dnx,
        )
    else:
        q1 = gamma_m ** (p - 1.0)
        _fill_dnx_interval(gamma, x_edges, q1, p, gamma_m, min(gamma_c, gamma_max), dnx)
        _fill_dnx_interval(
            gamma,
            x_edges,
            q1 * gamma_c,
            p + 1.0,
            max(gamma_m, gamma_c),
            gamma_max,
            dnx,
        )
    return dnx / (gamma * math.log(10.0))


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


def _nearest_grid_metric(gamma: np.ndarray, target: float) -> dict[str, float]:
    idx = int(np.argmin(np.abs(np.log10(gamma) - math.log10(target))))
    nearest = float(gamma[idx])
    return {
        "nearest_gamma": nearest,
        "nearest_index": idx,
        "rel_error": float(abs(nearest - target) / target) if target > 0.0 else float("nan"),
        "log_distance": float(abs(math.log10(nearest) - math.log10(target))) if target > 0.0 else float("nan"),
    }


def main() -> None:
    num_values = [61, 121, 241]
    cfg = build_baseline_config(
        electron_solver="fullhide",
        num_gam_e=121,
        num_nu=121,
        num_r=160,
        include_forward_ssc=False,
        index_y=0,
    )
    setup = build_simulation_setup(cfg)
    dyn = solve_dynamics(setup.boundary, cfg)

    radius_ini = float(dyn.radius[0])
    gamma_bulk = float(dyn.r_gamma[0])
    d_ne = float(_ambient_density(radius_ini, cfg))
    db = 0.39 * math.sqrt(cfg.epsilon_b * d_ne * (gamma_bulk * (gamma_bulk - 1.0)))
    gam_e_max = 3.0 * constants.para_m_energy / math.sqrt(8.0 * db * constants.para_e**3)
    db_min = 0.39 * math.sqrt(
        cfg.epsilon_b * d_ne * (float(dyn.r_gamma[-1]) * (float(dyn.r_gamma[-1]) - 1.0))
    )
    gam_e_max_max = 3.0 * constants.para_m_energy / math.sqrt(8.0 * db_min * constants.para_e**3)
    gam_e_m = float(_minimum_electron_lorentz_factor(cfg, gamma_bulk, gam_e_max))
    gam_e_c = 7.7e8 / (1.0 + math.sqrt(cfg.epsilon_e / cfg.epsilon_b)) / gamma_bulk / db**2 / (
        float(dyn.r_tobs[0]) / 2.0
    )

    payload = {
        "gamma_m_theory": gam_e_m,
        "gamma_c_theory": gam_e_c,
        "gamma_max_theory": gam_e_max,
        "grid": {},
    }

    for n in num_values:
        ref_cfg = build_baseline_config(
            electron_solver="fullhide",
            num_gam_e=n,
            num_nu=n,
            num_r=160,
            include_forward_ssc=False,
            index_y=0,
        )
        ref_setup = build_simulation_setup(ref_cfg)
        ele = solve_electron(ref_setup.boundary, dyn, ref_setup.seed_frequency_hz, ref_cfg)
        gamma = np.asarray(ele.gam_e, dtype=float)
        cell_center = _initial_cell_center_spectrum(gamma, ref_cfg.p, gam_e_m, gam_e_c, gam_e_max)
        cell_integrated = _initial_cell_integrated_spectrum(gamma, ref_cfg.p, gam_e_m, gam_e_c, gam_e_max)
        payload["grid"][str(n)] = {
            "cell_center": {
                **_support_metrics(gamma, cell_center),
                "nearest_to_gamma_m": _nearest_grid_metric(gamma, gam_e_m),
            },
            "cell_integrated": {
                **_support_metrics(gamma, cell_integrated),
                "nearest_to_gamma_m": _nearest_grid_metric(gamma, gam_e_m),
            },
        }

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
