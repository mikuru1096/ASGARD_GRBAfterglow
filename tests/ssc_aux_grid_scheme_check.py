from __future__ import annotations

import json
import sys
from pathlib import Path
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_presets import build_baseline_config
from asgard_runtime import solve_dynamics
from asgard_setup import build_simulation_setup
from src import Radiation
import src.Electron.FS_electron_slc1 as slc1_module


OUTPUT_JSON = ROOT / "output" / "vegasafterglow_doc" / "ssc_aux_grid_scheme.json"
AUX_COUNTS = (40, 64, 121)


def _shell_support_bounds(work_x_edge_log10: np.ndarray, work_d_n_x: np.ndarray) -> tuple[float, float]:
    x_lo = np.inf
    x_hi = -np.inf
    for i_shell in range(work_d_n_x.shape[1]):
        q_shell = np.asarray(work_d_n_x[:, i_shell], dtype=float)
        peak = float(np.max(q_shell))
        if not np.isfinite(peak) or peak <= 0.0:
            continue
        active = np.where(q_shell >= 1.0e-12 * peak)[0]
        if active.size == 0:
            continue
        x_lo = min(x_lo, float(work_x_edge_log10[active[0], i_shell]))
        x_hi = max(x_hi, float(work_x_edge_log10[active[-1] + 1, i_shell]))
    if not np.isfinite(x_lo) or not np.isfinite(x_hi) or x_hi <= x_lo:
        x_lo = float(np.min(work_x_edge_log10))
        x_hi = float(np.max(work_x_edge_log10))
    return x_lo, x_hi


def _project_shell_conservative(x_old_edge: np.ndarray, q_old: np.ndarray, x_new_edge: np.ndarray) -> np.ndarray:
    n_old = q_old.shape[0]
    n_new = x_new_edge.shape[0] - 1
    q_new_int = np.zeros(n_new, dtype=float)
    i_old = 0
    i_new = 0
    while i_old < n_old and i_new < n_new:
        x_lo = max(x_old_edge[i_old], x_new_edge[i_new])
        x_hi = min(x_old_edge[i_old + 1], x_new_edge[i_new + 1])
        if x_hi > x_lo:
            q_new_int[i_new] += q_old[i_old] * (x_hi - x_lo)
        if x_old_edge[i_old + 1] <= x_new_edge[i_new + 1]:
            i_old += 1
        else:
            i_new += 1
    widths = np.diff(x_new_edge)
    return q_new_int / widths


def _build_aux_edges(work_x_edge_log10: np.ndarray, work_d_n_x: np.ndarray, num_aux: int) -> np.ndarray:
    x_lo, x_hi = _shell_support_bounds(work_x_edge_log10, work_d_n_x)
    margin = 0.08
    return np.linspace(x_lo - margin, x_hi + margin, num_aux + 1, dtype=float)


def _project_work_to_aux_conservative(
    work_x_edge_log10: np.ndarray,
    work_d_n_x: np.ndarray,
    num_aux: int,
) -> tuple[np.ndarray, np.ndarray]:
    x_aux_edge = _build_aux_edges(work_x_edge_log10, work_d_n_x, num_aux)
    d_n_x_aux = np.empty((num_aux, work_d_n_x.shape[1]), dtype=float)
    for i_shell in range(work_d_n_x.shape[1]):
        d_n_x_aux[:, i_shell] = _project_shell_conservative(
            np.asarray(work_x_edge_log10[:, i_shell], dtype=float),
            np.asarray(work_d_n_x[:, i_shell], dtype=float),
            x_aux_edge,
        )
    gam_aux = 10.0 ** (0.5 * (x_aux_edge[:-1] + x_aux_edge[1:]))
    d_n_gam_aux = d_n_x_aux / (gam_aux[:, None] * np.log(10.0))
    return gam_aux, np.asfortranarray(d_n_gam_aux)


def _project_work_to_aux_center_log(
    work_x_edge_log10: np.ndarray,
    work_d_n_x: np.ndarray,
    num_aux: int,
    renormalize: bool,
) -> tuple[np.ndarray, np.ndarray]:
    x_aux_edge = _build_aux_edges(work_x_edge_log10, work_d_n_x, num_aux)
    x_aux = 0.5 * (x_aux_edge[:-1] + x_aux_edge[1:])
    q_aux = np.zeros((num_aux, work_d_n_x.shape[1]), dtype=float)
    tiny = 1.0e-300
    for i_shell in range(work_d_n_x.shape[1]):
        x_old = 0.5 * (work_x_edge_log10[:-1, i_shell] + work_x_edge_log10[1:, i_shell])
        q_old = np.asarray(work_d_n_x[:, i_shell], dtype=float)
        y_old = np.log(np.maximum(q_old, tiny))
        y_new = np.interp(x_aux, x_old, y_old, left=np.log(tiny), right=np.log(tiny))
        q_new = np.exp(y_new)
        if renormalize:
            widths_old = np.diff(work_x_edge_log10[:, i_shell])
            widths_new = np.diff(x_aux_edge)
            integ_old = float(np.sum(q_old * widths_old))
            integ_new = float(np.sum(q_new * widths_new))
            if integ_new > 0.0:
                q_new *= integ_old / integ_new
        q_aux[:, i_shell] = q_new
    gam_aux = 10.0 ** x_aux
    d_n_gam_aux = q_aux / (gam_aux[:, None] * np.log(10.0))
    return np.asfortranarray(gam_aux), np.asfortranarray(d_n_gam_aux)


def _timed_ssc_public(
    radius: np.ndarray,
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
    v_seed: np.ndarray,
    seed_syn: np.ndarray,
) -> tuple[float, np.ndarray]:
    t0 = perf_counter()
    p_ssc, _ = Radiation.ssc_spec(radius, gam_e, d_n_gam_e, v_seed, seed_syn, 1)
    return perf_counter() - t0, np.asarray(p_ssc)


def _timed_ssc_nonuniform(
    radius: np.ndarray,
    work_x_edge_log10: np.ndarray,
    work_d_n_x: np.ndarray,
    v_seed: np.ndarray,
    seed_syn: np.ndarray,
) -> tuple[float, np.ndarray]:
    t0 = perf_counter()
    p_ssc, _ = Radiation.ssc_spec_nonuniform(radius, work_x_edge_log10, work_d_n_x, v_seed, seed_syn, 1)
    return perf_counter() - t0, np.asarray(p_ssc)


def _timed_ssc_aux(
    radius: np.ndarray,
    work_x_edge_log10: np.ndarray,
    work_d_n_x: np.ndarray,
    v_seed: np.ndarray,
    seed_syn: np.ndarray,
    num_aux: int,
    projection: str,
) -> tuple[float, np.ndarray]:
    t0 = perf_counter()
    if projection == "conservative":
        gam_aux, d_n_gam_aux = _project_work_to_aux_conservative(work_x_edge_log10, work_d_n_x, num_aux)
    elif projection == "center_log":
        gam_aux, d_n_gam_aux = _project_work_to_aux_center_log(
            work_x_edge_log10, work_d_n_x, num_aux, renormalize=False
        )
    elif projection == "center_log_renorm":
        gam_aux, d_n_gam_aux = _project_work_to_aux_center_log(
            work_x_edge_log10, work_d_n_x, num_aux, renormalize=True
        )
    else:
        raise ValueError(f"unknown projection: {projection}")
    p_ssc, _ = Radiation.ssc_spec(radius, np.asfortranarray(gam_aux), d_n_gam_aux, v_seed, seed_syn, 1)
    return perf_counter() - t0, np.asarray(p_ssc)


def _summary_against_reference(candidate: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    ref = np.asarray(reference, dtype=float)
    cand = np.asarray(candidate, dtype=float)
    mask = np.isfinite(ref) & np.isfinite(cand) & (np.abs(ref) > 0.0)
    rel = np.abs((cand[mask] - ref[mask]) / ref[mask])
    ref_threshold = 1.0e-6 * float(np.max(np.abs(ref)))
    core = np.isfinite(ref) & np.isfinite(cand) & (np.abs(ref) >= ref_threshold)
    rel_core = np.abs((cand[core] - ref[core]) / ref[core])
    return {
        "median_rel": float(np.median(rel)),
        "max_rel": float(np.max(rel)),
        "core_median_rel": float(np.median(rel_core)),
        "core_max_rel": float(np.max(rel_core)),
    }


def main() -> None:
    config = build_baseline_config(
        num_threads=1,
        num_gam_e=32,
        num_nu=121,
        num_r=120,
        num_tobs=96,
        epsilon_e=0.2,
        epsilon_b=1.0e-5,
        d_ne=10.0,
        a_star=-1.0,
        electron_solver="slc1_mmg2",
    )
    setup = build_simulation_setup(config)
    dynamics = solve_dynamics(setup.boundary, config)

    (
        gam_e,
        d_n_gam_e,
        p_syn,
        seed_syn,
        nu_m,
        nu_c,
        nu_a,
        work_x_edge_log10,
        work_d_n_x,
    ) = slc1_module.fs_electron_slc1_mmg2(
        setup.boundary,
        dynamics.r_tobs,
        dynamics.r_gamma,
        dynamics.radius,
        np.asfortranarray(setup.seed_frequency_hz),
        config.num_gam_e,
        config.index_y,
        config.index_syn_integr,
        config.num_threads,
    )

    radius = np.asfortranarray(dynamics.radius)
    v_seed = np.asfortranarray(setup.seed_frequency_hz)
    gam_e = np.asfortranarray(gam_e)
    d_n_gam_e = np.asfortranarray(d_n_gam_e)
    seed_syn = np.asfortranarray(seed_syn)
    work_x_edge_log10 = np.asfortranarray(work_x_edge_log10)
    work_d_n_x = np.asfortranarray(work_d_n_x)

    public_seconds, public_ssc = _timed_ssc_public(radius, gam_e, d_n_gam_e, v_seed, seed_syn)
    nonuniform_seconds, nonuniform_ssc = _timed_ssc_nonuniform(radius, work_x_edge_log10, work_d_n_x, v_seed, seed_syn)

    payload: dict[str, object] = {
        "config": {
            "num_gam_e_work": int(config.num_gam_e),
            "num_nu": int(config.num_nu),
            "num_r": int(config.num_r),
            "epsilon_e": float(config.epsilon_e),
            "epsilon_b": float(config.epsilon_b),
            "d_ne": float(config.d_ne),
            "a_star": float(config.a_star),
        },
        "reference": {
            "scheme": "direct_nonuniform_ssc",
            "seconds": float(nonuniform_seconds),
        },
        "public_uniform": {
            "seconds": float(public_seconds),
            **_summary_against_reference(public_ssc, nonuniform_ssc),
            "speedup_vs_nonuniform": float(nonuniform_seconds / public_seconds),
        },
        "aux_uniform": {},
    }

    for projection in ("conservative", "center_log", "center_log_renorm"):
        payload["aux_uniform"][projection] = {}
        for num_aux in AUX_COUNTS:
            aux_seconds, aux_ssc = _timed_ssc_aux(
                radius,
                work_x_edge_log10,
                work_d_n_x,
                v_seed,
                seed_syn,
                num_aux,
                projection,
            )
            payload["aux_uniform"][projection][str(num_aux)] = {
                "seconds": float(aux_seconds),
                **_summary_against_reference(aux_ssc, nonuniform_ssc),
                "speedup_vs_nonuniform": float(nonuniform_seconds / aux_seconds),
            }

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
