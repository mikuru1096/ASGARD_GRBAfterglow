from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet, Wind
from ASGARD.api_observe import _build_fit_config_for_patch, _solve_patch_state
from asgard_core.asgard_state import project_flux_grid


PRODUCTION_GRID = {
    "num_r": 300,
    "num_theta": 300,
    "num_phi": 1,
    "num_tobs": 72,
    "num_gam_e": 96,
    "num_nu": 128,
}


def _build_model(medium_name: str, threads: int) -> Model:
    if medium_name == "ism":
        medium = ISM(n_ism=1.0)
    elif medium_name == "wind":
        medium = Wind(A_star=0.1, n_ism=1.0)
    else:
        raise ValueError(f"Unsupported medium: {medium_name}")

    return Model(
        TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0),
        medium,
        Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.3, ssc=True, kn=True),
        setups=Setups(
            fwd_ssc=True,
            ssc_cooling=True,
            electron_solver="fullhide_1d",
            num_threads=threads,
            num_r=PRODUCTION_GRID["num_r"],
            num_theta=PRODUCTION_GRID["num_theta"],
            num_phi=PRODUCTION_GRID["num_phi"],
            num_tobs=PRODUCTION_GRID["num_tobs"],
            num_gam_e=PRODUCTION_GRID["num_gam_e"],
            num_nu=PRODUCTION_GRID["num_nu"],
            observer_time_min_s=1.0e2,
            observer_time_max_s=1.0e8,
            electron_adaptive_substeps=False,
        ),
    )


def _run_case(medium_name: str, threads: int) -> dict[str, float | int | str]:
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    model = _build_model(medium_name, threads)
    times_s = np.logspace(2.0, 8.0, PRODUCTION_GRID["num_tobs"])
    frequencies_hz = np.logspace(9.0, 22.0, PRODUCTION_GRID["num_nu"])
    config = _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=0.0,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )

    timings: dict[str, float] = {}
    start = perf_counter()
    state = _solve_patch_state(model, config, times_s, frequencies_hz, timings=timings)
    obs_state = project_flux_grid(state, times_s, frequencies_hz, timings=timings)
    flux = obs_state.components["total"]
    total_s = perf_counter() - start
    if not np.isfinite(flux).all() or not (flux >= 0.0).all():
        raise RuntimeError(f"{medium_name} threads={threads}: invalid flux")

    max_log_curv = 0.0
    for row in flux[[0, 32, 64, 96, 127], :]:
        good = row > 0.0
        if good.sum() > 4:
            max_log_curv = max(max_log_curv, float(np.max(np.abs(np.diff(np.log10(row[good]), 2)))))

    return {
        "medium": medium_name,
        "threads": threads,
        "total_s": total_s,
        "dynamics_s": sum(v for k, v in timings.items() if k.startswith("Dynamics.")),
        "electron_s": sum(v for k, v in timings.items() if k.startswith("Electron.")),
        "ssc_s": sum(v for k, v in timings.items() if k.startswith("Radiation.ssc_spec")),
        "annihilation_s": timings.get("Radiation.annihilation", 0.0),
        "interp_s": sum(v for k, v in timings.items() if k.startswith("Interpolation.sed_interpolation")),
        "max_log_curv": max_log_curv,
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Production-grid fullhide_1d runtime benchmark.")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--media", nargs="+", choices=("ism", "wind"), default=("ism", "wind"))
    parser.add_argument("--threads", nargs="+", type=int, default=(1, 8))
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    rows: list[dict[str, float | int | str]] = []
    for medium_name in args.media:
        for threads in args.threads:
            print(f"production {medium_name} threads={threads}", flush=True)
            row = _run_case(medium_name, threads)
            print(row, flush=True)
            rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(args.output)


if __name__ == "__main__":
    main()
