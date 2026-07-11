from __future__ import annotations

import argparse
import csv
import json
import platform
import subprocess
import sys
from pathlib import Path
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core import Model, Observer, UniformMedium, top_hat_jet
from tests.public_api_builders import hadronic, numerics, observer_grid, radiation, reverse_shock, solver_options


MODES = {
    "quick": {"levels": (24, 36, 48), "reference": 72, "query": 24, "repeats": 1},
    "formal": {"levels": (48, 96, 192), "reference": 384, "query": 48, "repeats": 7},
}
BASE = {"radius": 48, "electron": 49, "photon": 41, "angle": 32}
DIMENSION_KEY = {
    "radius": "num_radius",
    "electron": "num_electron_gamma",
    "photon": "num_photon_frequency",
    "angle": "eats_num_theta",
}
M_E_C2 = 9.1093837015e-28 * 2.99792458e10**2
XI_E = 0.1
EPS_E = 0.1


def model(grid: dict[str, int], *, density_jump: bool = False) -> Model:
    item = Model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=100.0,
            opening_angle_rad=0.1,
            shell_duration_s=10.0,
            magnetar=None,
            spreading=False,
        ),
        medium=UniformMedium(number_density_cm3=1.0),
        observer=Observer(z=0.1, viewing_angle_rad=0.0, viewing_azimuth_rad=0.0, luminosity_distance_cm=1.0e26),
        fwd_rad=radiation(epsilon_e=EPS_E, accelerated_electron_fraction=XI_E),
        rvs_rad=radiation(epsilon_e=EPS_E, accelerated_electron_fraction=XI_E),
        numerics=numerics(
            num_radius=grid["radius"],
            num_electron_gamma=grid["electron"],
            num_photon_frequency=grid["photon"],
            eats_num_theta=grid["angle"],
            num_observer_time=max(24, grid["radius"] // 2),
            num_threads=1,
            initial_radius_cm=1.0e13,
        ),
        observer_grid=observer_grid(time_min_s=1.0, time_max_s=1.0e8),
        solver_options=solver_options(),
        reverse_shock=reverse_shock(enabled=True, shell_duration_s=10.0, upstream_sigma=0.1),
        hadronic=hadronic(),
    )
    if density_jump:
        item.setups.jump_r_cm = (1.0e15, 1.0e16, 1.0e17)
        item.setups.jump_factor = (100.0, 100.0, 100.0)
        item.setups.jump_width = (0.1, 0.1, 0.1)
    return item


def query(n: int) -> tuple[np.ndarray, np.ndarray]:
    return np.logspace(1.0, 7.0, n), np.array([1.0e9, 5.0e14, 1.0e18])


def solve(grid: dict[str, int], nquery: int, *, density_jump: bool = False) -> tuple[dict[str, np.ndarray], object, float]:
    times, frequencies = query(nquery)
    item = model(grid, density_jump=density_jump)
    start = perf_counter()
    result = item.flux_density_grid(times, frequencies)
    flux = {
        "total": np.asarray(result.total, dtype=float),
        "forward": np.asarray(result.fwd.sync + result.fwd.ssc, dtype=float),
        "reverse": np.asarray(result.rev.sync + result.rev.ssc, dtype=float),
    }
    details = item.details(float(times[0]), float(times[-1]))
    elapsed = perf_counter() - start
    if not all(np.all(np.isfinite(values)) for values in flux.values()):
        raise RuntimeError("non-finite end-to-end flux")
    return flux, details, elapsed


def relative_metrics(value: np.ndarray, reference: np.ndarray) -> tuple[float, float, float]:
    peak = np.max(np.abs(reference))
    active = np.abs(reference) >= peak * 1.0e-10
    if not np.any(active):
        raise RuntimeError("reference component has no cell above its fixed 1e-10 peak mask")
    error = np.abs(value[active] - reference[active]) / np.abs(reference[active])
    return float(np.median(error)), float(np.percentile(error, 95.0)), float(np.max(error))


def electron_number_error(track) -> float:
    gamma = np.asarray(track.gamma_e, dtype=float)
    dnde = np.asarray(track.dN_dgamma_e, dtype=float)
    number = np.trapezoid(dnde, gamma, axis=0)
    expected = XI_E * np.asarray(track.N_p, dtype=float)
    active = expected > 0.0
    return float(np.max(np.abs(number[active] / expected[active] - 1.0)))


def secondary_budget(details) -> dict[str, float]:
    rev = details.rev
    if rev is None or rev.secondary_rs_dissipated_energy_erg is None:
        raise RuntimeError("density-jump run returned no secondary-RS energy registers")
    dissipated = np.asarray(rev.secondary_rs_dissipated_energy_erg, dtype=float)
    injected = np.asarray(rev.secondary_rs_electron_injected_energy_erg, dtype=float)
    swept = np.asarray(rev.secondary_rs_swept_mass_g, dtype=float)
    internal = np.asarray(rev.secondary_rs_internal_energy_erg, dtype=float)
    active = dissipated > 0.0
    if not np.any(active):
        raise RuntimeError("density-jump run activated no secondary reverse shock")
    return {
        "electron_energy_fractional_error": float(np.max(np.abs(injected[active] / (EPS_E * dissipated[active]) - 1.0))),
        "min_swept_mass_increment_g": float(np.min(np.diff(swept[active]))) if np.count_nonzero(active) > 1 else 0.0,
        "min_internal_energy_erg": float(np.min(internal[active])),
        "active_shells": int(np.count_nonzero(active)),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(mode: str, output: Path) -> None:
    spec = MODES[mode]
    output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for dimension in DIMENSION_KEY:
        reference_grid = dict(BASE)
        reference_grid[dimension] = spec["reference"]
        reference, _, _ = solve(reference_grid, spec["query"])
        next_grid = dict(BASE)
        next_grid[dimension] = spec["levels"][-1]
        next_solution, _, _ = solve(next_grid, spec["query"])
        uncertainty = {name: relative_metrics(next_solution[name], reference[name])[1] for name in reference}
        for level in spec["levels"]:
            grid = dict(BASE)
            grid[dimension] = level
            samples = [solve(grid, spec["query"]) for _ in range(spec["repeats"])]
            for component in ("total", "forward", "reverse"):
                median, p95, maximum = relative_metrics(samples[-1][0][component], reference[component])
                rows.append({
                    "dimension": dimension,
                    "component": component,
                    "level": level,
                    "reference_level": spec["reference"],
                    "reference_uncertainty_p95": uncertainty[component],
                    "relative_error_median": median,
                    "relative_error_p95": p95,
                    "relative_error_max": maximum,
                    "wall_time_median_s": float(np.median([sample[2] for sample in samples])),
                    "electron_number_error_max": electron_number_error(samples[-1][1].fwd),
                })
    jump_grid = {key: spec["levels"][-1] for key in BASE}
    _, jump_details, jump_wall = solve(jump_grid, spec["query"], density_jump=True)
    budget = secondary_budget(jump_details)
    write_csv(output / "convergence.csv", rows)
    write_csv(output / "secondary_rs_budget.csv", [{"wall_time_s": jump_wall, **budget}])
    metadata = {
        "schema": "asgard.numerical_validation.v1",
        "mode": mode,
        "git_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, encoding="utf-8"
        ).strip(),
        "python": sys.version,
        "platform": platform.platform(),
        "threads": 1,
        "query_shape": [3, spec["query"]],
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate end-to-end ASGARD convergence and conservation evidence.")
    parser.add_argument("--mode", choices=MODES, default="quick")
    parser.add_argument("--output-dir", type=Path, default=Path("/tmp/asgard_benchmarks/numerical_validation"))
    args = parser.parse_args()
    run(args.mode, args.output_dir)


if __name__ == "__main__":
    main()
