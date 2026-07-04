from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path
import sys
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core import Model, Observer, UniformMedium, WindMedium, top_hat_jet
from asgard_core.api_model import _direct_tophat_patch_config, _solve_patch_state
from asgard_core.asgard_paths import DOC_ROOT
from asgard_core.asgard_state import project_flux_grid
from tests.public_api_builders import hadronic, numerics, observer_grid, radiation, reverse_shock, solver_options


OUTPUT_DIR = DOC_ROOT / "runtime_benchmark"

GRID_MODES = {
    "standard": {"num_r": 80, "num_theta": 80, "num_tobs": 48, "num_gam_e": 81, "num_nu": 81, "num_chi": 8, "num_query_t": 32, "num_query_nu": 32},
}

STAGE_ORDER = [
    "dynamics",
    "electron_synch",
    "ssc",
    "annihilation",
    "sed_eats_interpolation",
    "python_middle_layer",
]

STAGE_LABELS = {
    "dynamics": "Dynamics",
    "electron_synch": "Electron + synch",
    "ssc": "SSC",
    "annihilation": "Gamma-gamma",
    "sed_eats_interpolation": "SED/EATS interp.",
    "python_middle_layer": "Python layer",
}

PALETTE = {
    1: "#0072B2",
    8: "#D55E00",
}


@dataclass(frozen=True)
class CaseSpec:
    medium: str
    dimension: str
    solver: str
    threads: int


def _build_model(case: CaseSpec, grid: dict[str, int]) -> Model:
    if case.medium == "ism":
        medium = UniformMedium(number_density_cm3=1.0)
    elif case.medium == "wind":
        medium = WindMedium(a_star=0.1, density_floor_cm3=1.0, density_cap_cm3=None)
    else:
        raise ValueError(f"Unsupported medium: {case.medium}")

    is_2d = case.solver == "fullhide_2d"
    return Model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=300.0,
            opening_angle_rad=0.1,
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        medium=medium,
        observer=Observer(
            z=0.1,
            viewing_angle_rad=0.0,
            viewing_azimuth_rad=0.0,
            luminosity_distance_cm=1.0e26,
        ),
        fwd_rad=radiation(
            epsilon_e=0.1,
            epsilon_B=1.0e-3,
            p=2.3,
            include_ssc=True,
            include_kn_correction=True,
        ),
        rvs_rad=None,
        numerics=numerics(
            num_radius=grid["num_r"],
            eats_num_theta=grid["num_theta"],
            eats_num_phi=1,
            num_observer_time=grid["num_tobs"],
            num_electron_gamma=grid["num_gam_e"],
            num_photon_frequency=grid["num_nu"],
            downstream_num_chi=grid["num_chi"] if is_2d else None,
            num_threads=case.threads,
            electron_adaptive_substeps=False,
            initial_radius_cm=1.0e14,
        ),
        observer_grid=observer_grid(time_min_s=1.0e2, time_max_s=1.0e7),
        solver_options=solver_options(
            electron_solver=case.solver,
            geometry_projection="chi_eats_2d" if is_2d else "sed_legacy",
            ssc_cooling_mode="numeric_ic_kn",
        ),
        reverse_shock=reverse_shock(),
        hadronic=hadronic(),
    )


def _query_grid(grid: dict[str, int]) -> tuple[np.ndarray, np.ndarray]:
    times_s = np.logspace(2.0, 7.0, grid["num_query_t"])
    frequencies_hz = np.logspace(9.0, 22.0, grid["num_query_nu"])
    return times_s, frequencies_hz


def _config_for_direct_model(model: Model):
    return _direct_tophat_patch_config(model)


def _time_high_level(case: CaseSpec, grid: dict[str, int], times_s: np.ndarray, frequencies_hz: np.ndarray) -> float:
    model = _build_model(case, grid)
    start = perf_counter()
    result = model.flux_density_grid(times_s, frequencies_hz)
    elapsed = perf_counter() - start
    if result.total.shape != (frequencies_hz.size, times_s.size):
        raise RuntimeError(f"Unexpected flux grid shape: {result.total.shape}")
    if not np.all(np.isfinite(result.total)):
        raise RuntimeError("Non-finite flux value in runtime benchmark output.")
    return elapsed


def _time_breakdown(case: CaseSpec, grid: dict[str, int], times_s: np.ndarray, frequencies_hz: np.ndarray) -> dict[str, float]:
    model = _build_model(case, grid)
    config = _config_for_direct_model(model)
    timings: dict[str, float] = {}

    start = perf_counter()
    state = _solve_patch_state(model, config, times_s, frequencies_hz, timings=timings)
    project_flux_grid(state, times_s, frequencies_hz, timings=timings)
    total = perf_counter() - start

    dynamics = sum(value for key, value in timings.items() if key.startswith("Dynamics."))
    electron_synch = sum(value for key, value in timings.items() if key.startswith("Electron."))
    ssc = sum(value for key, value in timings.items() if key.startswith("Radiation.ssc_spec"))
    annihilation = timings.get("Radiation.annihilation", 0.0)
    interpolation = sum(value for key, value in timings.items() if key.startswith("Interpolation.sed_interpolation"))
    measured = dynamics + electron_synch + ssc + annihilation + interpolation

    return {
        "total_low_level": total,
        "dynamics": dynamics,
        "electron_synch": electron_synch,
        "ssc": ssc,
        "annihilation": annihilation,
        "sed_eats_interpolation": interpolation,
        "python_middle_layer": total - measured,
    }


def _run_one(case: CaseSpec, grid: dict[str, int], repeats: int, warmup: int) -> list[dict[str, float | int | str]]:
    times_s, frequencies_hz = _query_grid(grid)
    for _ in range(warmup):
        _time_breakdown(case, grid, times_s, frequencies_hz)

    rows: list[dict[str, float | int | str]] = []
    for rep in range(repeats):
        high_level_total = _time_high_level(case, grid, times_s, frequencies_hz)
        breakdown = _time_breakdown(case, grid, times_s, frequencies_hz)
        rows.append(
            {
                "medium": case.medium,
                "dimension": case.dimension,
                "solver": case.solver,
                "threads": case.threads,
                "repeat": rep,
                "high_level_total": high_level_total,
                **breakdown,
            }
        )
    return rows


def _aggregate(rows: list[dict[str, float | int | str]]) -> list[dict[str, float | int | str]]:
    keys = {(row["medium"], row["dimension"], row["solver"], row["threads"]) for row in rows}
    out: list[dict[str, float | int | str]] = []
    for medium, dimension, solver, threads in sorted(keys):
        subset = [row for row in rows if (row["medium"], row["dimension"], row["solver"], row["threads"]) == (medium, dimension, solver, threads)]
        item: dict[str, float | int | str] = {
            "medium": medium,
            "dimension": dimension,
            "solver": solver,
            "threads": int(threads),
            "repeats": len(subset),
        }
        for name in ["high_level_total", "total_low_level", *STAGE_ORDER]:
            values = np.array([float(row[name]) for row in subset], dtype=float)
            item[f"{name}_median_s"] = float(np.median(values))
            item[f"{name}_min_s"] = float(np.min(values))
            item[f"{name}_max_s"] = float(np.max(values))
        out.append(item)
    return out


def _write_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _plot_breakdown(summary: list[dict[str, float | int | str]], path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14.0, 8.8), sharey=True, constrained_layout=True)
    case_pairs = [("ism", "1d"), ("ism", "2d"), ("wind", "1d"), ("wind", "2d")]
    x = np.arange(len(STAGE_ORDER), dtype=float)
    width = 0.36

    for ax, (medium, dimension) in zip(axes.flat, case_pairs):
        for offset, threads in [(-width / 2.0, 1), (width / 2.0, 8)]:
            row = next(
                item
                for item in summary
                if item["medium"] == medium and item["dimension"] == dimension and int(item["threads"]) == threads
            )
            values = np.array([float(row[f"{stage}_median_s"]) for stage in STAGE_ORDER], dtype=float)
            ax.bar(x + offset, values, width=width, color=PALETTE[threads], label=f"{threads} thread" if threads == 1 else "8 threads")
        ax.set_title(f"{medium.upper()} {dimension.upper()}")
        ax.set_xticks(x, [STAGE_LABELS[stage] for stage in STAGE_ORDER], rotation=28, ha="right")
        ax.set_yscale("log")
        ax.set_ylabel("median time [s]")
        ax.grid(axis="y", which="both", alpha=0.25)
        ax.legend(loc="upper right", fontsize=8)

    fig.suptitle("ASGARD fixed-grid runtime breakdown", fontsize=14)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def _plot_totals(summary: list[dict[str, float | int | str]], path: Path) -> None:
    labels = ["ISM 1D", "ISM 2D", "Wind 1D", "Wind 2D"]
    cases = [("ism", "1d"), ("ism", "2d"), ("wind", "1d"), ("wind", "2d")]
    x = np.arange(len(cases), dtype=float)
    width = 0.36

    fig, ax = plt.subplots(figsize=(9.8, 5.2), constrained_layout=True)
    serial_values = []
    parallel_values = []
    speedups = []
    for medium, dimension in cases:
        serial = next(item for item in summary if item["medium"] == medium and item["dimension"] == dimension and int(item["threads"]) == 1)
        parallel = next(item for item in summary if item["medium"] == medium and item["dimension"] == dimension and int(item["threads"]) == 8)
        serial_total = float(serial["high_level_total_median_s"])
        parallel_total = float(parallel["high_level_total_median_s"])
        serial_values.append(serial_total)
        parallel_values.append(parallel_total)
        speedups.append(serial_total / parallel_total)

    ax.bar(x - width / 2.0, serial_values, width=width, color=PALETTE[1], label="1 thread")
    ax.bar(x + width / 2.0, parallel_values, width=width, color=PALETTE[8], label="8 threads")
    ax.set_xticks(x, labels)
    ax.set_yscale("log")
    ax.set_ylabel("end-to-end median time [s]")
    ax.set_title("ASGARD fixed-grid end-to-end runtime")
    ax.grid(axis="y", which="both", alpha=0.25)
    ax.legend()

    for xi, value, speedup in zip(x, parallel_values, speedups):
        ax.text(xi + width / 2.0, value, f"x{speedup:.2f}", ha="center", va="bottom", fontsize=8)

    fig.savefig(path, dpi=220)
    plt.close(fig)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="ASGARD fixed-grid runtime breakdown benchmark.")
    parser.add_argument("--grid", choices=tuple(GRID_MODES), default="standard")
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    if args.repeats < 1:
        raise ValueError("--repeats must be positive.")
    if args.warmup < 0:
        raise ValueError("--warmup must be non-negative.")

    grid = GRID_MODES[args.grid]
    cases = [
        CaseSpec(medium=medium, dimension=dimension, solver=solver, threads=threads)
        for medium in ("ism", "wind")
        for dimension, solver in (("1d", "fullhide_1d"), ("2d", "fullhide_2d"))
        for threads in (1, 8)
    ]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, float | int | str]] = []
    print(f"grid={args.grid} repeats={args.repeats} warmup={args.warmup}", flush=True)
    for case in cases:
        os.environ["OMP_NUM_THREADS"] = str(case.threads)
        os.environ["OPENBLAS_NUM_THREADS"] = "1"
        print(f"{case.medium} {case.dimension} {case.solver} threads={case.threads}", flush=True)
        rows.extend(_run_one(case, grid, args.repeats, args.warmup))

    summary = _aggregate(rows)
    raw_csv = args.output_dir / "runtime_breakdown_raw.csv"
    summary_csv = args.output_dir / "runtime_breakdown_summary.csv"
    _write_csv(raw_csv, rows)
    _write_csv(summary_csv, summary)

    breakdown_png = args.output_dir / "runtime_breakdown_bars.png"
    totals_png = args.output_dir / "runtime_total_bars.png"
    _plot_breakdown(summary, breakdown_png)
    _plot_totals(summary, totals_png)

    print(raw_csv)
    print(summary_csv)
    print(breakdown_png)
    print(totals_png)


if __name__ == "__main__":
    main()
