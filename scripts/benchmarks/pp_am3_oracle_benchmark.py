from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from pathlib import Path
from threading import Event, Thread
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tests.public_api_builders import hadronic, numerics, radiation, top_hat_model
from scripts.benchmarks.benchmark_common import summary


MODELS = ("sibyll", "qgsjet", "geant4", "pythia8")
GRIDS = {"quick": (101,), "formal": (201, 401)}
REPEATS = {"quick": 1, "formal": 7}


def build(model_name: str | None, size: int):
    updates = {
        "enabled": True,
        "solver": "am3_1d",
        "num_proton_gamma": size,
        "num_neutrino_frequency": size,
        "pgamma_scheme": "hummer_2010_response",
    }
    if model_name is not None:
        updates["pp_gamma_model"] = model_name
    had = hadronic(**updates)
    return top_hat_model(
        fwd_rad=radiation(
            proton_energy_fraction=0.01,
            include_pgamma=False,
            pp=True,
            neutrino=True,
            pgamma_scheme="hummer_2010_response",
        ),
        hadronic=had,
        numerics=numerics(
            num_radius=8,
            num_electron_gamma=41,
            num_photon_frequency=size,
            num_observer_time=12,
            eats_num_theta=8,
            num_threads=1,
        ),
    )


def rss_mib() -> float:
    pages = int(Path("/proc/self/statm").read_text(encoding="ascii").split()[1])
    return pages * os.sysconf("SC_PAGE_SIZE") / 1024.0**2


def rss_peak(stop: Event, samples: list[float]) -> None:
    while not stop.wait(0.01):
        samples.append(rss_mib())


def solve(model_name: str | None, size: int) -> tuple[np.ndarray, np.ndarray, float, float]:
    item = build(model_name, size)
    stop, rss = Event(), [rss_mib()]
    monitor = Thread(target=rss_peak, args=(stop, rss))
    monitor.start()
    start = perf_counter()
    track = item.details(1.0e2, 1.0e5).fwd
    elapsed = perf_counter() - start
    stop.set()
    monitor.join()
    gamma = np.asarray(track.l_had_pg_gamma_spec, dtype=float)
    neutrino = np.asarray(track.neutrino_luminosity, dtype=float)
    if not np.all(np.isfinite(gamma)) or not np.all(np.isfinite(neutrino)):
        raise RuntimeError(f"non-finite pp output for {model_name or 'default'}")
    if np.any(gamma < 0.0) or np.any(neutrino < 0.0):
        raise RuntimeError(f"negative pp output for {model_name or 'default'}")
    return gamma, neutrino, elapsed, max(rss)


def read_oracle(path: Path) -> dict[tuple[str, int, int], float]:
    required = {"model", "grid_size", "flat_index", "qgamma"}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if set(reader.fieldnames or ()) != required:
            raise RuntimeError(f"AM3 oracle columns must be exactly {sorted(required)}")
        rows = list(reader)
    return {(row["model"], int(row["grid_size"]), int(row["flat_index"])): float(row["qgamma"]) for row in rows}


def oracle_error(model_name: str, size: int, values: np.ndarray, oracle: dict[tuple[str, int, int], float]) -> tuple[float, float, int]:
    expected = np.array([oracle[(model_name, size, index)] for index in range(values.size)], dtype=float)
    actual = values.ravel()
    active = (expected >= np.max(expected) * 1.0e-10) & (actual > 0.0)
    if not np.any(active):
        raise RuntimeError(f"AM3 oracle and ASGARD have no common positive significant domain for {model_name}")
    error = np.abs(actual[active] - expected[active]) / expected[active]
    return float(np.median(error)), float(np.max(error)), int(np.count_nonzero(active))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(mode: str, output: Path, oracle_path: Path | None, delta_path: Path | None) -> None:
    if mode == "formal" and (oracle_path is None or delta_path is None):
        raise RuntimeError("formal pp benchmark requires official --am3-oracle-csv and parent-version --delta-oracle-csv")
    am3_version = importlib.metadata.version("am3") if oracle_path is not None else None
    oracle = read_oracle(oracle_path) if oracle_path is not None else None
    delta_oracle = read_oracle(delta_path) if delta_path is not None else None
    output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    spectra: list[dict[str, object]] = []
    for size in GRIDS[mode]:
        delta_gamma, delta_nu, _, _ = solve("delta", size)
        delta_bitwise = ""
        if delta_oracle is not None:
            expected = np.array([delta_oracle[("delta", size, index)] for index in range(delta_gamma.size)])
            delta_bitwise = bool(np.array_equal(delta_gamma.ravel(), expected))
            if not delta_bitwise:
                raise RuntimeError("delta output is not bitwise identical to the parent-version formal pp oracle")
        for model_name in MODELS:
            samples = [solve(model_name, size) for _ in range(REPEATS[mode])]
            gamma, neutrino = samples[-1][:2]
            if not np.array_equal(neutrino, delta_nu):
                raise RuntimeError(f"{model_name} changed the pp neutrino chain")
            median_error = maximum_error = oracle_points = ""
            if oracle is not None:
                median_error, maximum_error, oracle_points = oracle_error(model_name, size, gamma, oracle)
            timing = summary([sample[2] for sample in samples])
            memory = summary([sample[3] for sample in samples])
            rows.append({
                "model": model_name,
                "grid_size": size,
                "grid_cells": size * size,
                "wall_time_median_s": timing["median"],
                "wall_time_iqr_s": timing["iqr"],
                "wall_time_p95_s": timing["p95"],
                "peak_rss_mib_median": memory["median"],
                "peak_rss_mib_iqr": memory["iqr"],
                "pp_gamma_sum": float(np.sum(gamma)),
                "pp_gamma_max": float(np.max(gamma)),
                "oracle_relative_error_median": median_error,
                "oracle_relative_error_max": maximum_error,
                "oracle_significant_points": oracle_points,
                "delta_parent_bitwise": delta_bitwise,
                "neutrino_delta_bitwise": True,
            })
            spectra.extend(
                {"model": model_name, "grid_size": size, "flat_index": index, "qgamma": value}
                for index, value in enumerate(gamma.ravel())
            )
    write_csv(output / "pp_models.csv", rows)
    write_csv(output / "pp_spectra.csv", spectra)
    metadata = {
        "schema": "asgard.pp_am3_oracle.v1",
        "mode": mode,
        "git_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, encoding="utf-8"
        ).strip(),
        "python": sys.version,
        "platform": platform.platform(),
        "threads": 1,
        "repeats": REPEATS[mode],
        "am3_oracle": str(oracle_path) if oracle_path else None,
        "delta_oracle": str(delta_path) if delta_path else None,
        "am3_python_version": am3_version,
        "pp_output_semantics": "p-gamma is disabled, so l_had_pg_gamma_spec contains pp gamma only",
        "oracle_mask": "oracle >= 1e-10 oracle peak and ASGARD > 0",
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark ASGARD pp gamma models against an official AM3 oracle table.")
    parser.add_argument("--mode", choices=GRIDS, default="quick")
    parser.add_argument("--output-dir", type=Path, default=Path("/tmp/asgard_benchmarks/pp_am3"))
    parser.add_argument("--am3-oracle-csv", type=Path)
    parser.add_argument("--delta-oracle-csv", type=Path)
    args = parser.parse_args()
    run(args.mode, args.output_dir, args.am3_oracle_csv, args.delta_oracle_csv)


if __name__ == "__main__":
    main()
