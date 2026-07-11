"""Cold-solve, projection, and cached-query scaling benchmark."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from threading import Event, Thread
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core import Model, Observer, UniformMedium, WindMedium, top_hat_jet
from asgard_core.api_model import _direct_tophat_patch_config, _solve_patch_state
from asgard_core.asgard_state import project_flux
from scripts.benchmarks.benchmark_common import DATA_ROOT, FIGURE_ROOT, PALETTE, environment, plot_style, save_figure, summary, write_json
from tests.public_api_builders import hadronic, numerics, observer_grid, radiation, reverse_shock, solver_options


GRID = {
    "quick": dict(num_r=32, num_theta=24, num_tobs=24, num_gam_e=31, num_nu=31, num_chi=4, num_query_t=12, num_query_nu=12),
    "formal": dict(num_r=80, num_theta=80, num_tobs=48, num_gam_e=81, num_nu=81, num_chi=8, num_query_t=32, num_query_nu=32),
}
STAGES = ("cold_solve", "projection", "warm_query")
LABELS = {"cold_solve": "Cold solve", "projection": "Observer projection", "warm_query": "Warm cached query"}
COLORS = {1: PALETTE["black"], 2: PALETTE["sky"], 4: PALETTE["blue"], 8: PALETTE["vermillion"]}


def _model(medium_name: str, dimension: str, solver: str, threads: int, grid: dict[str, int]) -> Model:
    medium = UniformMedium(number_density_cm3=1.0) if medium_name == "ism" else WindMedium(a_star=0.1, density_floor_cm3=1.0, density_cap_cm3=None)
    two_d = dimension == "2d"
    return Model(
        jet=top_hat_jet(energy_iso_erg=1.0e52, initial_lorentz_factor=300.0, opening_angle_rad=0.1, shell_duration_s=None, magnetar=None, spreading=False),
        medium=medium,
        observer=Observer(z=0.1, viewing_angle_rad=0.0, viewing_azimuth_rad=0.0, luminosity_distance_cm=1.0e26),
        fwd_rad=radiation(epsilon_e=0.1, epsilon_B=1.0e-3, p=2.3, include_ssc=True),
        rvs_rad=None,
        numerics=numerics(num_radius=grid["num_r"], eats_num_theta=grid["num_theta"], eats_num_phi=1, num_observer_time=grid["num_tobs"], num_electron_gamma=grid["num_gam_e"], num_photon_frequency=grid["num_nu"], downstream_num_chi=grid["num_chi"] if two_d else None, num_threads=threads, electron_adaptive_substeps=False, initial_radius_cm=1.0e14),
        observer_grid=observer_grid(time_min_s=1.0e2, time_max_s=1.0e7),
        solver_options=solver_options(electron_solver=solver, geometry_projection="chi_eats_2d" if two_d else "sed_legacy", ssc_cooling_mode="numeric_ic_kn"),
        reverse_shock=reverse_shock(),
        hadronic=hadronic(),
    )


def _grid(grid: dict[str, int]) -> tuple[np.ndarray, np.ndarray]:
    return np.logspace(2.0, 7.0, grid["num_query_t"]), np.logspace(9.0, 22.0, grid["num_query_nu"])


def _rss() -> float:
    pages = int(Path("/proc/self/statm").read_text(encoding="ascii").split()[1])
    return pages * os.sysconf("SC_PAGE_SIZE") / 1024.0**2


def _rss_peak(stop: Event, samples: list[float]) -> None:
    while not stop.wait(0.01):
        samples.append(_rss())


def _sample(medium: str, dimension: str, solver: str, threads: int, grid: dict[str, int]) -> dict[str, float]:
    times, frequencies = _grid(grid)
    stop, rss = Event(), [_rss()]
    monitor = Thread(target=_rss_peak, args=(stop, rss))
    monitor.start()
    start = perf_counter()
    model = _model(medium, dimension, solver, threads, grid)
    cold_result = model.flux_density_grid(times, frequencies)
    cold = perf_counter() - start
    config = _direct_tophat_patch_config(model)
    state = _solve_patch_state(model, config, times, frequencies)
    start = perf_counter()
    result = project_flux(state, times, frequencies, mode="total_only")
    projection = perf_counter() - start
    start = perf_counter()
    cached = model.flux_density_grid(times, frequencies)
    warm = perf_counter() - start
    total = np.asarray(result.components["total"])
    if total.shape != np.asarray(cold_result.total).shape or not np.array_equal(cached.total, cold_result.total):
        raise RuntimeError("invalid runtime benchmark projection")
    stop.set()
    monitor.join()
    return {"cold_solve": cold, "projection": projection, "warm_query": warm, "peak_rss_mib": max(rss)}


def _run(medium: str, dimension: str, solver: str, threads: int, grid: dict[str, int], repeats: int, warmup: int) -> list[dict[str, object]]:
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    for _ in range(warmup):
        _sample(medium, dimension, solver, threads, grid)
    return [dict(medium=medium, dimension=dimension, solver=solver, threads=threads, repeat=repeat, **_sample(medium, dimension, solver, threads, grid)) for repeat in range(repeats)]


def _aggregate(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    keys = sorted({(row["medium"], row["dimension"], row["solver"], row["threads"]) for row in rows})
    out = []
    for medium, dimension, solver, threads in keys:
        subset = [row for row in rows if (row["medium"], row["dimension"], row["solver"], row["threads"]) == (medium, dimension, solver, threads)]
        item: dict[str, object] = dict(medium=medium, dimension=dimension, solver=solver, threads=threads, repeats=len(subset))
        for name in (*STAGES, "peak_rss_mib"):
            stats = summary([row[name] for row in subset])
            item.update({f"{name}_{key}": value for key, value in stats.items() if key != "count"})
        out.append(item)
    return out


def _csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _plot(rows: list[dict[str, object]], path: Path) -> None:
    cases = sorted({(str(row["medium"]), str(row["dimension"])) for row in rows})
    fig, axes = plt.subplots(1, len(cases), figsize=(3.25 * len(cases), 3.0), sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    for axis, (medium, dimension) in zip(axes, cases):
        subset = [row for row in rows if row["medium"] == medium and row["dimension"] == dimension]
        for stage, marker, linestyle in zip(STAGES, ("o", "s", "^"), ("-", "--", ":")):
            x = [int(row["threads"]) for row in subset]
            y = [row[f"{stage}_median"] for row in subset]
            axis.plot(x, y, color=PALETTE["black"], marker=marker, markerfacecolor="none", linestyle=linestyle, label=LABELS[stage])
            axis.scatter(x, y, c=[COLORS[value] for value in x], marker=marker, zorder=3)
        axis.set(xlabel="OpenMP threads", title=f"{medium.upper()} {dimension.upper()}", xticks=[row["threads"] for row in subset], yscale="log")
        axis.grid(alpha=0.18)
    axes[0].set_ylabel("Wall time [s]")
    axes[-1].legend()
    save_figure(fig, path)
    plt.close(fig)


def _args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("quick", "formal"), default="quick")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    return parser.parse_args()


def main() -> None:
    args = _args()
    formal = args.mode == "formal"
    repeats, warmup = (7, 1) if formal else (1, 1)
    threads = (1, 2, 4, 8) if formal else (1, 2)
    pairs = (("ism", "1d", "fullhide_1d"), ("ism", "2d", "fullhide_2d"), ("wind", "1d", "fullhide_1d"), ("wind", "2d", "fullhide_2d")) if formal else (("ism", "1d", "fullhide_1d"),)
    data_dir = args.output_dir or DATA_ROOT / "runtime_scaling"
    figure_dir = args.figure_dir or FIGURE_ROOT / "runtime_scaling"
    rows = [row for medium, dimension, solver in pairs for thread in threads for row in _run(medium, dimension, solver, thread, GRID[args.mode], repeats, warmup)]
    aggregate = _aggregate(rows)
    _csv(data_dir / "runtime_raw.csv", rows)
    _csv(data_dir / "runtime_summary.csv", aggregate)
    write_json(data_dir / "metadata.json", environment(args.mode, threads=max(threads), grid=GRID[args.mode], repeats=repeats) | {"warmup": warmup, "boundaries": {"cold_solve": "model construction through particle/radiation state solve", "projection": "first observer projection from solved state", "warm_query": "second projection from the same solved state"}})
    plt.rcParams.update(plot_style())
    _plot(aggregate, figure_dir / "runtime_scaling")


if __name__ == "__main__":
    main()
