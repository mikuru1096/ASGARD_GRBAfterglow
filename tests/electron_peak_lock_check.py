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
from asgard_runtime import _ambient_density, _minimum_electron_lorentz_factor, solve_dynamics, solve_electron
from asgard_setup import build_simulation_setup
from src import constants


OUTPUT_JSON = asgard_doc_path("electron_peak_lock_diagnostic.json")


def _support_metrics(gamma: np.ndarray, spectrum: np.ndarray) -> dict[str, float]:
    peak_idx = int(np.argmax(spectrum))
    peak_gamma = float(gamma[peak_idx])
    peak_value = float(spectrum[peak_idx])
    mask = spectrum > 1.0e-8 * peak_value
    if np.any(mask):
        support_idx = np.where(mask)[0]
        g_lo = float(gamma[support_idx[0]])
        g_hi = float(gamma[support_idx[-1]])
    else:
        g_lo = float("nan")
        g_hi = float("nan")
    return {
        "peak_gamma": peak_gamma,
        "peak_value": peak_value,
        "g_lo": g_lo,
        "g_hi": g_hi,
    }


def _nearest_grid_metric(gamma_grid: np.ndarray, target: float) -> dict[str, float]:
    idx = int(np.argmin(np.abs(np.log10(gamma_grid) - np.log10(target))))
    nearest = float(gamma_grid[idx])
    return {
        "nearest_gamma": nearest,
        "nearest_index": idx,
        "rel_error": float(abs(nearest - target) / target) if target > 0.0 else float("nan"),
        "log_distance": float(abs(np.log10(nearest) - np.log10(target))) if target > 0.0 else float("nan"),
    }


def main() -> None:
    num_values = [61, 121, 241]
    shell_fracs = [0.3, 0.5, 0.8]
    solvers = ["fullhide", "slc1", "charint", "t2g1", "weno5"]

    base_cfg = build_baseline_config(
        electron_solver="fullhide",
        num_gam_e=121,
        num_nu=121,
        num_r=160,
        include_forward_ssc=False,
        index_y=0,
    )
    base_setup = build_simulation_setup(base_cfg)
    base_dyn = solve_dynamics(base_setup.boundary, base_cfg)

    payload: dict[str, dict] = {}
    for solver in solvers:
        solver_rows = []
        for n in num_values:
            cfg = build_baseline_config(
                electron_solver=solver,
                num_gam_e=n,
                num_nu=n,
                num_r=160,
                include_forward_ssc=False,
                index_y=0,
            )
            setup = build_simulation_setup(cfg)
            ele = solve_electron(base_setup.boundary, base_dyn, setup.seed_frequency_hz, cfg)
            solver_rows.append(
                {
                    "n": n,
                    "gam_e": np.asarray(ele.gam_e, dtype=float),
                    "d_n_gam_e": np.asarray(ele.d_n_gam_e, dtype=float),
                }
            )

        shell_table = []
        for frac in shell_fracs:
            shell = min(max(int(round(frac * (base_cfg.num_r - 1))), 1), base_cfg.num_r - 1)
            radius_loc = float(base_dyn.radius[shell - 1])
            gamma_bulk = 0.5 * float(base_dyn.r_gamma[shell - 1] + base_dyn.r_gamma[shell])
            d_ne = float(_ambient_density(radius_loc, base_cfg))
            db = 0.39 * np.sqrt(base_cfg.epsilon_b * d_ne * (gamma_bulk * (gamma_bulk - 1.0)))
            gam_e_max = 3.0 * constants.para_m_energy / np.sqrt(8.0 * db * constants.para_e**3)
            gam_e_m = float(_minimum_electron_lorentz_factor(base_cfg, gamma_bulk, gam_e_max))
            gam_e_c = 7.7e8 * (1.0 + base_cfg.z) / gamma_bulk / db**2 / float(base_dyn.r_tobs[shell])

            row = {
                "shell_frac": frac,
                "shell_index": shell,
                "radius_cm": radius_loc,
                "gamma_bulk": gamma_bulk,
                "gamma_m_theory": gam_e_m,
                "gamma_c_theory": gam_e_c,
                "grid": {},
            }
            for item in solver_rows:
                gamma = item["gam_e"]
                spectrum = item["d_n_gam_e"][:, shell]
                metrics = _support_metrics(gamma, spectrum)
                metrics["nearest_to_gamma_m"] = _nearest_grid_metric(gamma, gam_e_m)
                metrics["nearest_to_gamma_c"] = _nearest_grid_metric(gamma, gam_e_c)
                metrics["peak_to_gamma_m"] = float(metrics["peak_gamma"] / gam_e_m)
                metrics["peak_to_gamma_c"] = float(metrics["peak_gamma"] / gam_e_c)
                row["grid"][str(item["n"])] = metrics
            shell_table.append(row)
        payload[solver] = {"shells": shell_table}

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
