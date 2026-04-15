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
from asgard_component_backend import solve_model_state_from_setup, extract_physical_solution_from_state
from asgard_setup import build_simulation_setup
from asgard_postprocess import compute_band_fluxes
import src.Electron.FS_electron_fullhide as fullhide_module
import src.Electron.FS_electron_slc1 as slc1_module
import src.Electron.FS_electron_t2g1 as t2g1_module
import src.Electron.FS_electron_weno5 as weno5_module


OUTPUT_JSON = ROOT / "output" / "asgard_doc" / "electron_solver_comparison.json"


def _run_direct_fullhide(boundary, r_tobs, r_gamma, radius, v_seed, config, adaptive: bool):
    start = perf_counter()
    gam_e, d_n_gam_e, p_syn, seed_syn, nu_m, nu_c, nu_a = fullhide_module.fs_electron_fullhide(
        boundary,
        r_tobs,
        r_gamma,
        radius,
        v_seed,
        config.num_gam_e,
        config.index_y,
        config.index_syn_integr,
        config.num_threads,
        1 if adaptive else 0,
        config.electron_substep_rtol,
        config.electron_substep_min,
        config.electron_substep_max,
    )
    elapsed = perf_counter() - start
    return {
        "seconds": elapsed,
        "gam_e": np.asarray(gam_e),
        "d_n_gam_e": np.asarray(d_n_gam_e),
        "p_syn": np.asarray(p_syn),
        "seed_syn": np.asarray(seed_syn),
        "nu_m": np.asarray(nu_m),
        "nu_c": np.asarray(nu_c),
        "nu_a": np.asarray(nu_a),
    }


def _run_direct_solver(name: str, boundary, r_tobs, r_gamma, radius, v_seed, config):
    if name == "fullhide":
        return _run_direct_fullhide(boundary, r_tobs, r_gamma, radius, v_seed, config, adaptive=False)
    if name == "fullhide_adaptive":
        return _run_direct_fullhide(boundary, r_tobs, r_gamma, radius, v_seed, config, adaptive=True)
    start = perf_counter()
    if name == "t2g1":
        gam_e, d_n_gam_e, p_syn, seed_syn, nu_m, nu_c, nu_a = t2g1_module.fs_electron_t2g1(
            boundary,
            r_tobs,
            r_gamma,
            radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
        )
        elapsed = perf_counter() - start
        return {
            "seconds": elapsed,
            "gam_e": np.asarray(gam_e),
            "d_n_gam_e": np.asarray(d_n_gam_e),
            "p_syn": np.asarray(p_syn),
            "seed_syn": np.asarray(seed_syn),
            "nu_m": np.asarray(nu_m),
            "nu_c": np.asarray(nu_c),
            "nu_a": np.asarray(nu_a),
        }
    if name == "slc1":
        gam_e, d_n_gam_e, p_syn, seed_syn, nu_m, nu_c, nu_a = slc1_module.fs_electron_slc1(
            boundary,
            r_tobs,
            r_gamma,
            radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
        )
        elapsed = perf_counter() - start
        return {
            "seconds": elapsed,
            "gam_e": np.asarray(gam_e),
            "d_n_gam_e": np.asarray(d_n_gam_e),
            "p_syn": np.asarray(p_syn),
            "seed_syn": np.asarray(seed_syn),
            "nu_m": np.asarray(nu_m),
            "nu_c": np.asarray(nu_c),
            "nu_a": np.asarray(nu_a),
        }
    if name == "slc1_mmg2":
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
            boundary,
            r_tobs,
            r_gamma,
            radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
        )
        elapsed = perf_counter() - start
        return {
            "seconds": elapsed,
            "gam_e": np.asarray(gam_e),
            "d_n_gam_e": np.asarray(d_n_gam_e),
            "p_syn": np.asarray(p_syn),
            "seed_syn": np.asarray(seed_syn),
            "nu_m": np.asarray(nu_m),
            "nu_c": np.asarray(nu_c),
            "nu_a": np.asarray(nu_a),
            "work_x_edge_log10": np.asarray(work_x_edge_log10),
            "work_d_n_x": np.asarray(work_d_n_x),
        }
    if name == "weno5":
        gam_e, d_n_gam_e, p_syn, seed_syn = weno5_module.fs_electron_weno5(
            boundary,
            r_tobs,
            r_gamma,
            radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.num_threads,
        )
        elapsed = perf_counter() - start
        return {
            "seconds": elapsed,
            "gam_e": np.asarray(gam_e),
            "d_n_gam_e": np.asarray(d_n_gam_e),
            "p_syn": np.asarray(p_syn),
            "seed_syn": np.asarray(seed_syn),
        }
    raise ValueError(name)


def _run_observed_solver(name: str, base_config):
    config = build_baseline_config(**vars(base_config))
    config.electron_solver = "fullhide" if name == "fullhide_adaptive" else name
    config.electron_adaptive_substeps = name == "fullhide_adaptive"
    start = perf_counter()
    setup = build_simulation_setup(config)
    state = solve_model_state_from_setup(config, setup)
    physical = extract_physical_solution_from_state(state)
    bands_flux = compute_band_fluxes(setup, physical, config)
    elapsed = perf_counter() - start
    return {
        "seconds": elapsed,
        "bands_flux": np.asarray(bands_flux),
        "nu_m": np.asarray(physical.nu_m),
        "nu_c": np.asarray(physical.nu_c),
        "nu_a": np.asarray(physical.nu_a),
    }


def _rel_diff(candidate: np.ndarray, baseline: np.ndarray) -> float:
    return float(np.nanmax(np.abs(candidate - baseline) / np.maximum(np.abs(baseline), 1.0e-99)))


def main() -> None:
    config = build_baseline_config(
        num_threads=1,
        num_r=120,
        num_nu=121,
        num_gam_e=121,
        num_tobs=64,
        electron_substep_rtol=2.0e-2,
        electron_substep_min=25,
        electron_substep_max=150,
    )
    setup = build_simulation_setup(config)
    dynamics = solve_dynamics(setup.boundary, config)
    boundary = setup.boundary
    v_seed = np.asfortranarray(setup.seed_frequency_hz)

    direct_fullhide = _run_direct_solver("fullhide", boundary, dynamics.r_tobs, dynamics.r_gamma, dynamics.radius, v_seed, config)
    direct_adaptive = _run_direct_solver("fullhide_adaptive", boundary, dynamics.r_tobs, dynamics.r_gamma, dynamics.radius, v_seed, config)
    direct_t2g1 = _run_direct_solver("t2g1", boundary, dynamics.r_tobs, dynamics.r_gamma, dynamics.radius, v_seed, config)
    direct_slc1 = _run_direct_solver("slc1", boundary, dynamics.r_tobs, dynamics.r_gamma, dynamics.radius, v_seed, config)
    direct_slc1_mmg2 = _run_direct_solver("slc1_mmg2", boundary, dynamics.r_tobs, dynamics.r_gamma, dynamics.radius, v_seed, config)
    direct_weno5 = _run_direct_solver("weno5", boundary, dynamics.r_tobs, dynamics.r_gamma, dynamics.radius, v_seed, config)

    obs_fullhide = _run_observed_solver("fullhide", config)
    obs_adaptive = _run_observed_solver("fullhide_adaptive", config)
    obs_t2g1 = _run_observed_solver("t2g1", config)
    obs_slc1 = _run_observed_solver("slc1", config)
    obs_slc1_mmg2 = _run_observed_solver("slc1_mmg2", config)
    obs_weno5 = _run_observed_solver("weno5", config)

    payload = {
        "baseline": "fullhide",
        "direct": {
            "fullhide": {"seconds": direct_fullhide["seconds"]},
            "fullhide_adaptive": {
                "seconds": direct_adaptive["seconds"],
                "rel_d_n_gam_e": _rel_diff(direct_adaptive["d_n_gam_e"], direct_fullhide["d_n_gam_e"]),
                "rel_p_syn": _rel_diff(direct_adaptive["p_syn"], direct_fullhide["p_syn"]),
                "rel_seed_syn": _rel_diff(direct_adaptive["seed_syn"], direct_fullhide["seed_syn"]),
                "rel_nu_m": _rel_diff(direct_adaptive["nu_m"], direct_fullhide["nu_m"]),
                "rel_nu_c": _rel_diff(direct_adaptive["nu_c"], direct_fullhide["nu_c"]),
                "rel_nu_a": _rel_diff(direct_adaptive["nu_a"], direct_fullhide["nu_a"]),
            },
            "t2g1": {
                "seconds": direct_t2g1["seconds"],
                "rel_d_n_gam_e": _rel_diff(direct_t2g1["d_n_gam_e"], direct_fullhide["d_n_gam_e"]),
                "rel_p_syn": _rel_diff(direct_t2g1["p_syn"], direct_fullhide["p_syn"]),
                "rel_seed_syn": _rel_diff(direct_t2g1["seed_syn"], direct_fullhide["seed_syn"]),
                "rel_nu_m": _rel_diff(direct_t2g1["nu_m"], direct_fullhide["nu_m"]),
                "rel_nu_c": _rel_diff(direct_t2g1["nu_c"], direct_fullhide["nu_c"]),
                "rel_nu_a": _rel_diff(direct_t2g1["nu_a"], direct_fullhide["nu_a"]),
            },
            "slc1": {
                "seconds": direct_slc1["seconds"],
                "rel_d_n_gam_e": _rel_diff(direct_slc1["d_n_gam_e"], direct_fullhide["d_n_gam_e"]),
                "rel_p_syn": _rel_diff(direct_slc1["p_syn"], direct_fullhide["p_syn"]),
                "rel_seed_syn": _rel_diff(direct_slc1["seed_syn"], direct_fullhide["seed_syn"]),
                "rel_nu_m": _rel_diff(direct_slc1["nu_m"], direct_fullhide["nu_m"]),
                "rel_nu_c": _rel_diff(direct_slc1["nu_c"], direct_fullhide["nu_c"]),
                "rel_nu_a": _rel_diff(direct_slc1["nu_a"], direct_fullhide["nu_a"]),
            },
            "slc1_mmg2": {
                "seconds": direct_slc1_mmg2["seconds"],
                "rel_d_n_gam_e": _rel_diff(direct_slc1_mmg2["d_n_gam_e"], direct_fullhide["d_n_gam_e"]),
                "rel_p_syn": _rel_diff(direct_slc1_mmg2["p_syn"], direct_fullhide["p_syn"]),
                "rel_seed_syn": _rel_diff(direct_slc1_mmg2["seed_syn"], direct_fullhide["seed_syn"]),
                "rel_nu_m": _rel_diff(direct_slc1_mmg2["nu_m"], direct_fullhide["nu_m"]),
                "rel_nu_c": _rel_diff(direct_slc1_mmg2["nu_c"], direct_fullhide["nu_c"]),
                "rel_nu_a": _rel_diff(direct_slc1_mmg2["nu_a"], direct_fullhide["nu_a"]),
            },
            "weno5": {
                "seconds": direct_weno5["seconds"],
                "rel_d_n_gam_e": _rel_diff(direct_weno5["d_n_gam_e"], direct_fullhide["d_n_gam_e"]),
                "rel_p_syn": _rel_diff(direct_weno5["p_syn"], direct_fullhide["p_syn"]),
                "rel_seed_syn": _rel_diff(direct_weno5["seed_syn"], direct_fullhide["seed_syn"]),
            },
        },
        "observed": {
            "fullhide": {"seconds": obs_fullhide["seconds"]},
            "fullhide_adaptive": {
                "seconds": obs_adaptive["seconds"],
                "rel_bands_flux": _rel_diff(obs_adaptive["bands_flux"], obs_fullhide["bands_flux"]),
                "rel_nu_m": _rel_diff(obs_adaptive["nu_m"], obs_fullhide["nu_m"]),
                "rel_nu_c": _rel_diff(obs_adaptive["nu_c"], obs_fullhide["nu_c"]),
                "rel_nu_a": _rel_diff(obs_adaptive["nu_a"], obs_fullhide["nu_a"]),
            },
            "t2g1": {
                "seconds": obs_t2g1["seconds"],
                "rel_bands_flux": _rel_diff(obs_t2g1["bands_flux"], obs_fullhide["bands_flux"]),
                "rel_nu_m": _rel_diff(obs_t2g1["nu_m"], obs_fullhide["nu_m"]),
                "rel_nu_c": _rel_diff(obs_t2g1["nu_c"], obs_fullhide["nu_c"]),
                "rel_nu_a": _rel_diff(obs_t2g1["nu_a"], obs_fullhide["nu_a"]),
            },
            "slc1": {
                "seconds": obs_slc1["seconds"],
                "rel_bands_flux": _rel_diff(obs_slc1["bands_flux"], obs_fullhide["bands_flux"]),
                "rel_nu_m": _rel_diff(obs_slc1["nu_m"], obs_fullhide["nu_m"]),
                "rel_nu_c": _rel_diff(obs_slc1["nu_c"], obs_fullhide["nu_c"]),
                "rel_nu_a": _rel_diff(obs_slc1["nu_a"], obs_fullhide["nu_a"]),
            },
            "slc1_mmg2": {
                "seconds": obs_slc1_mmg2["seconds"],
                "rel_bands_flux": _rel_diff(obs_slc1_mmg2["bands_flux"], obs_fullhide["bands_flux"]),
                "rel_nu_m": _rel_diff(obs_slc1_mmg2["nu_m"], obs_fullhide["nu_m"]),
                "rel_nu_c": _rel_diff(obs_slc1_mmg2["nu_c"], obs_fullhide["nu_c"]),
                "rel_nu_a": _rel_diff(obs_slc1_mmg2["nu_a"], obs_fullhide["nu_a"]),
            },
            "weno5": {
                "seconds": obs_weno5["seconds"],
                "rel_bands_flux": _rel_diff(obs_weno5["bands_flux"], obs_fullhide["bands_flux"]),
                "rel_nu_m": _rel_diff(obs_weno5["nu_m"], obs_fullhide["nu_m"]),
                "rel_nu_c": _rel_diff(obs_weno5["nu_c"], obs_fullhide["nu_c"]),
                "rel_nu_a": _rel_diff(obs_weno5["nu_a"], obs_fullhide["nu_a"]),
            },
        },
    }

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
