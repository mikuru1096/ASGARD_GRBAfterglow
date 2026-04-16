from __future__ import annotations

import json
import sys
from pathlib import Path
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_component_backend import extract_physical_solution_from_state, solve_model_state_from_setup
from asgard_postprocess import compute_band_fluxes
from asgard_presets import build_baseline_config
from asgard_runtime import solve_dynamics
from asgard_setup import build_simulation_setup
import src.Electron.FS_electron_charint as charint_module
import src.Electron.FS_electron_fullhide as fullhide_module
import src.Electron.FS_electron_slc1 as slc1_module


OUTPUT_JSON = ROOT / "output" / "asgard_doc" / "electron_solver_comparison.json"
SOLVERS = ("fullhide", "slc1", "charint")
NUM_GAM_BY_SOLVER = {
    "fullhide": 121,
    "slc1": 121,
    "charint": 121,
}


def _solver_config(base_config, solver: str):
    config = build_baseline_config(**vars(base_config))
    config.electron_solver = solver
    config.electron_adaptive_substeps = False
    config.num_gam_e = NUM_GAM_BY_SOLVER[solver]
    return config


def _run_direct_solver(name: str, boundary, r_tobs, r_gamma, radius, v_seed, config):
    start = perf_counter()
    if name == "fullhide":
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
            0,
            config.electron_substep_rtol,
            config.electron_substep_min,
            config.electron_substep_max,
        )
    elif name == "slc1":
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
    else:
        gam_e, d_n_gam_e, p_syn, seed_syn, nu_m, nu_c, nu_a = charint_module.fs_electron_charint(
            boundary,
            r_tobs,
            r_gamma,
            radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
            1 if config.electron_adaptive_substeps else 0,
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


def _run_observed_solver(name: str, base_config):
    config = _solver_config(base_config, name)
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
    v_seed = np.asfortranarray(setup.seed_frequency_hz)
    direct = {
        solver: _run_direct_solver(
            solver,
            setup.boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            _solver_config(config, solver),
        )
        for solver in SOLVERS
    }
    observed = {solver: _run_observed_solver(solver, config) for solver in SOLVERS}
    baseline_direct = direct["fullhide"]
    baseline_observed = observed["fullhide"]

    payload = {
        "baseline": "fullhide",
        "solvers": list(SOLVERS),
        "direct": {
            solver: {
                "seconds": direct[solver]["seconds"],
                **(
                    {}
                    if solver == "fullhide"
                    else {
                        "num_gam_e": NUM_GAM_BY_SOLVER[solver],
                        "rel_nu_m": _rel_diff(direct[solver]["nu_m"], baseline_direct["nu_m"]),
                        "rel_nu_c": _rel_diff(direct[solver]["nu_c"], baseline_direct["nu_c"]),
                        "rel_nu_a": _rel_diff(direct[solver]["nu_a"], baseline_direct["nu_a"]),
                    }
                ),
            }
            for solver in SOLVERS
        },
        "observed": {
            solver: {
                "seconds": observed[solver]["seconds"],
                **(
                    {}
                    if solver == "fullhide"
                    else {
                        "num_gam_e": NUM_GAM_BY_SOLVER[solver],
                        "rel_bands_flux": _rel_diff(observed[solver]["bands_flux"], baseline_observed["bands_flux"]),
                        "rel_nu_m": _rel_diff(observed[solver]["nu_m"], baseline_observed["nu_m"]),
                        "rel_nu_c": _rel_diff(observed[solver]["nu_c"], baseline_observed["nu_c"]),
                        "rel_nu_a": _rel_diff(observed[solver]["nu_a"], baseline_observed["nu_a"]),
                    }
                ),
            }
            for solver in SOLVERS
        },
    }

    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
