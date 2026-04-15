from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_presets import build_baseline_config
from asgard_runtime import solve_dynamics
from asgard_setup import build_simulation_setup
import src.Electron.FS_electron_fullhide as fullhide_module
import src.Electron.FS_electron_slc1 as slc1_module


OUTPUT_JSON = ROOT / "output" / "asgard_doc" / "mmg2_front_sharpness.json"
SHELL_FRACTIONS = (0.3, 0.5, 0.8)


def _front_metrics_from_gamma(gamma: np.ndarray, spectrum: np.ndarray) -> dict[str, float]:
    gamma = np.asarray(gamma, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    peak_idx = int(np.argmax(spectrum))
    peak_gamma = float(gamma[peak_idx])
    peak_value = float(spectrum[peak_idx])
    x = np.log10(np.maximum(gamma, 1.0))
    q = np.log10(1.0 + np.maximum(spectrum, 0.0))

    if peak_idx < spectrum.size - 1:
        dx = np.maximum(x[peak_idx + 1 :] - x[peak_idx:-1], 1.0e-30)
        drop = (q[peak_idx:-1] - q[peak_idx + 1 :]) / dx
        front_rel = int(np.argmax(drop))
        front_idx = peak_idx + front_rel
        max_drop = float(drop[front_rel])
        front_gamma = float(np.sqrt(gamma[front_idx] * gamma[front_idx + 1]))
    else:
        front_idx = peak_idx
        max_drop = 0.0
        front_gamma = peak_gamma

    active = np.where(spectrum >= 1.0e-6 * max(peak_value, 1.0e-99))[0]
    g_hi = float(gamma[active[-1]]) if active.size else peak_gamma
    return {
        "peak_gamma": peak_gamma,
        "peak_value": peak_value,
        "front_gamma": front_gamma,
        "front_drop_dlogn_dlogg": max_drop,
        "high_support_gamma": g_hi,
        "front_width_dex": float(max(np.log10(g_hi) - np.log10(peak_gamma), 0.0)),
    }


def _front_metrics_from_work_grid(x_edge_log10: np.ndarray, d_n_x: np.ndarray) -> dict[str, float]:
    x_edge = np.asarray(x_edge_log10, dtype=float)
    d_n_x = np.asarray(d_n_x, dtype=float)
    x_center = 0.5 * (x_edge[:-1] + x_edge[1:])
    gamma = np.power(10.0, x_center)
    return _front_metrics_from_gamma(gamma, d_n_x / (gamma * np.log(10.0)))


def main() -> None:
    config = build_baseline_config(
        num_threads=1,
        num_r=120,
        num_nu=121,
        num_tobs=96,
        epsilon_e=0.2,
        epsilon_b=1.0e-5,
        d_ne=10.0,
        a_star=-1.0,
    )
    setup = build_simulation_setup(config)
    dynamics = solve_dynamics(setup.boundary, config)
    v_seed = np.asfortranarray(setup.seed_frequency_hz)

    fullhide_cfg = build_baseline_config(**vars(config))
    fullhide_cfg.num_gam_e = 161
    gam_full, dng_full, *_ = fullhide_module.fs_electron_fullhide(
        setup.boundary,
        dynamics.r_tobs,
        dynamics.r_gamma,
        dynamics.radius,
        v_seed,
        fullhide_cfg.num_gam_e,
        fullhide_cfg.index_y,
        fullhide_cfg.index_syn_integr,
        fullhide_cfg.num_threads,
        0,
        fullhide_cfg.electron_substep_rtol,
        fullhide_cfg.electron_substep_min,
        fullhide_cfg.electron_substep_max,
    )
    gam_slc1, dng_slc1, *_ = slc1_module.fs_electron_slc1(
        setup.boundary,
        dynamics.r_tobs,
        dynamics.r_gamma,
        dynamics.radius,
        v_seed,
        32,
        config.index_y,
        config.index_syn_integr,
        config.num_threads,
    )
    (
        gam_mmg2,
        dng_mmg2,
        *_unused,
        work_x_edge_log10,
        work_d_n_x,
    ) = slc1_module.fs_electron_slc1_mmg2(
        setup.boundary,
        dynamics.r_tobs,
        dynamics.r_gamma,
        dynamics.radius,
        v_seed,
        32,
        config.index_y,
        config.index_syn_integr,
        config.num_threads,
    )

    payload: dict[str, object] = {"shells": []}
    for frac in SHELL_FRACTIONS:
        shell = min(max(int(round(frac * (config.num_r - 1))), 1), config.num_r - 1)
        row = {
            "shell_frac": float(frac),
            "shell_index": int(shell),
            "fullhide": _front_metrics_from_gamma(gam_full, dng_full[:, shell]),
            "slc1": _front_metrics_from_gamma(gam_slc1, dng_slc1[:, shell]),
            "slc1_mmg2_public": _front_metrics_from_gamma(gam_mmg2, dng_mmg2[:, shell]),
            "slc1_mmg2_work": _front_metrics_from_work_grid(work_x_edge_log10[:, shell], work_d_n_x[:, shell]),
        }
        payload["shells"].append(row)

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
