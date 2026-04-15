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
from asgard_component_backend import build_query_setup, observe_spectra_from_setup, solve_model_state_from_setup
from asgard_runtime import solve_dynamics, solve_electron
from asgard_ssc import compute_ssc_auxiliary_grid
from asgard_setup import build_simulation_setup
from src import Radiation
import src.Electron.FS_electron_slc1 as slc1_module


OUTPUT_JSON = ROOT / "output" / "asgard_doc" / "ssc_aux_grid_scheme.json"
AUX_COUNTS = (40, 64, 121)


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
) -> tuple[float, np.ndarray]:
    t0 = perf_counter()
    p_ssc, _ = compute_ssc_auxiliary_grid(
        radius,
        work_x_edge_log10,
        work_d_n_x,
        v_seed,
        seed_syn,
        1,
        num_auxiliary_gamma=num_aux,
    )
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


def _high_energy_smoothness(values: np.ndarray) -> dict[str, float]:
    arr = np.asarray(values, dtype=float)
    mask = np.isfinite(arr) & (arr > 0.0)
    if np.count_nonzero(mask) < 3:
        return {"max_curvature_dex": float("inf"), "rms_curvature_dex": float("inf")}
    logy = np.log10(arr[mask])
    curvature = np.diff(logy, n=2)
    return {
        "max_curvature_dex": float(np.max(np.abs(curvature))),
        "rms_curvature_dex": float(np.sqrt(np.mean(curvature * curvature))),
    }


def _observed_lightcurve_summary(config, solver_name: str, times_s: np.ndarray, frequencies_hz: np.ndarray) -> dict[str, np.ndarray]:
    query_config = build_baseline_config(**vars(config))
    query_config.electron_solver = solver_name
    if solver_name == "fullhide":
        query_config.num_gam_e = 161
    setup = build_query_setup(query_config, times_s, frequencies_hz)
    state = solve_model_state_from_setup(query_config, setup)
    observed = observe_spectra_from_setup(query_config, state.component_spectra, setup, frequencies_hz)
    return {key: np.asarray(value, dtype=float) for key, value in observed.items() if value is not None}


def _compute_reference_ssc(config, solver_name: str) -> np.ndarray:
    reference_config = build_baseline_config(**vars(config))
    reference_config.electron_solver = solver_name
    reference_config.num_gam_e = 161
    reference_config.num_nu = config.num_nu
    reference_config.num_r = config.num_r
    reference_config.num_tobs = config.num_tobs
    setup = build_simulation_setup(reference_config)
    dynamics = solve_dynamics(setup.boundary, reference_config)
    electron = solve_electron(setup.boundary, dynamics, setup.seed_frequency_hz, reference_config)
    ssc, _ = Radiation.ssc_spec(
        np.asfortranarray(dynamics.radius),
        np.asfortranarray(electron.gam_e),
        np.asfortranarray(electron.d_n_gam_e),
        np.asfortranarray(setup.seed_frequency_hz),
        np.asfortranarray(electron.seed_syn),
        1,
    )
    return np.asarray(ssc, dtype=float)


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
    fullhide_ref = _compute_reference_ssc(config, "fullhide")
    slc1_ref = _compute_reference_ssc(config, "slc1")
    high_energy_index = int(np.argmax(v_seed >= 2.4e26)) if np.any(v_seed >= 2.4e26) else v_seed.shape[0] - 1
    observer_time_s = np.logspace(2.0, 8.0, 80)
    observer_freq_hz = np.array([1.0e18, 2.4e26], dtype=float)
    observed_fullhide = _observed_lightcurve_summary(config, "fullhide", observer_time_s, observer_freq_hz)
    observed_slc1 = _observed_lightcurve_summary(config, "slc1", observer_time_s, observer_freq_hz)
    observed_mmg2 = _observed_lightcurve_summary(config, "slc1_mmg2", observer_time_s, observer_freq_hz)

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
        "validation_basis": "current_dirty_tree_state",
        "reference": {
            "physical": {
                "fullhide_highres": {"num_gam_e": 161, "num_nu": int(config.num_nu)},
                "slc1_highres": {"num_gam_e": 161, "num_nu": int(config.num_nu)},
            },
            "diagnostic_only": {
                "scheme": "direct_nonuniform_ssc",
                "seconds": float(nonuniform_seconds),
            },
        },
        "public_uniform": {
            "seconds": float(public_seconds),
            **_summary_against_reference(public_ssc, nonuniform_ssc),
            "speedup_vs_nonuniform": float(nonuniform_seconds / public_seconds),
            "vs_fullhide_highres": _summary_against_reference(public_ssc, fullhide_ref),
            "vs_slc1_highres": _summary_against_reference(public_ssc, slc1_ref),
        },
        "aux_uniform": {},
        "diagnostic_nonuniform": {
            "vs_fullhide_highres": _summary_against_reference(nonuniform_ssc, fullhide_ref),
            "vs_slc1_highres": _summary_against_reference(nonuniform_ssc, slc1_ref),
        },
        "observed_forward_ssc_lightcurve": {
            "times_s": {
                "min": float(observer_time_s[0]),
                "max": float(observer_time_s[-1]),
                "count": int(observer_time_s.size),
            },
            "frequencies_hz": {
                "xray": float(observer_freq_hz[0]),
                "tev": float(observer_freq_hz[1]),
            },
            "smoothness": {
                "fullhide": {
                    "xray": _high_energy_smoothness(observed_fullhide["fwd_ssc"][0]),
                    "tev": _high_energy_smoothness(observed_fullhide["fwd_ssc"][1]),
                },
                "slc1": {
                    "xray": _high_energy_smoothness(observed_slc1["fwd_ssc"][0]),
                    "tev": _high_energy_smoothness(observed_slc1["fwd_ssc"][1]),
                },
                "slc1_mmg2": {
                    "xray": _high_energy_smoothness(observed_mmg2["fwd_ssc"][0]),
                    "tev": _high_energy_smoothness(observed_mmg2["fwd_ssc"][1]),
                },
            },
            "slc1_mmg2_vs_fullhide": {
                "xray": _summary_against_reference(observed_mmg2["fwd_ssc"][0], observed_fullhide["fwd_ssc"][0]),
                "tev": _summary_against_reference(observed_mmg2["fwd_ssc"][1], observed_fullhide["fwd_ssc"][1]),
            },
            "slc1_mmg2_vs_slc1": {
                "xray": _summary_against_reference(observed_mmg2["fwd_ssc"][0], observed_slc1["fwd_ssc"][0]),
                "tev": _summary_against_reference(observed_mmg2["fwd_ssc"][1], observed_slc1["fwd_ssc"][1]),
            },
        },
    }

    for num_aux in AUX_COUNTS:
        aux_seconds, aux_ssc = _timed_ssc_aux(
            radius,
            work_x_edge_log10,
            work_d_n_x,
            v_seed,
            seed_syn,
            num_aux,
        )
        payload["aux_uniform"][str(num_aux)] = {
            "seconds": float(aux_seconds),
            "speedup_vs_nonuniform": float(nonuniform_seconds / aux_seconds),
            "vs_public_uniform": _summary_against_reference(aux_ssc, public_ssc),
            "vs_nonuniform": _summary_against_reference(aux_ssc, nonuniform_ssc),
            "vs_fullhide_highres": _summary_against_reference(aux_ssc, fullhide_ref),
            "vs_slc1_highres": _summary_against_reference(aux_ssc, slc1_ref),
            "tev_shell_smoothness": _high_energy_smoothness(aux_ssc[high_energy_index]),
        }

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
