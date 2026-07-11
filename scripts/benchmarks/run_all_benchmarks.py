"""Run every ASGARD benchmark from raw calculation through paper figures."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

from benchmark_common import DATA_ROOT, FIGURE_ROOT, ROOT, environment, manifest, write_json


def command(script: str, *args: object) -> list[str]:
    return [sys.executable, str(ROOT / script), *(str(arg) for arg in args)]


def suites(mode: str) -> list[tuple[str, list[str]]]:
    formal = mode == "formal"
    pp_args: tuple[object, ...] = ()
    if formal:
        am3 = os.environ.get("ASGARD_AM3_ORACLE")
        delta = os.environ.get("ASGARD_DELTA_ORACLE")
        if am3 is None or delta is None:
            raise RuntimeError("formal pp benchmark requires ASGARD_AM3_ORACLE and ASGARD_DELTA_ORACLE")
        pp_args = ("--am3-oracle-csv", am3, "--delta-oracle-csv", delta)
    return [
        ("convergence", command("scripts/benchmarks/numerical_validation_benchmark.py", "--mode", mode, "--output-dir", DATA_ROOT / "convergence")),
        ("chi_eats", command("scripts/benchmarks/chi_eats_2d_benchmark.py", "--mode", mode, "--only", "convergence", "--output-dir", DATA_ROOT / "chi_eats", "--figure-dir", FIGURE_ROOT / "chi_eats")),
        ("runtime", command("scripts/benchmarks/runtime_breakdown_benchmark.py", "--mode", mode, "--output-dir", DATA_ROOT / "runtime", "--figure-dir", FIGURE_ROOT / "runtime")),
        ("pp_am3", command("scripts/benchmarks/pp_am3_oracle_benchmark.py", "--mode", mode, "--output-dir", DATA_ROOT / "pp_am3", *pp_args)),
        ("magnetized_rs_vegas_ism", command("scripts/benchmarks/magnetized_rs_vegas_lc_compare.py", "--mode", mode, "--medium", "ism", "--output", FIGURE_ROOT / "magnetized_rs_vegas" / "ism", "--data-dir", DATA_ROOT / "magnetized_rs_vegas" / "ism")),
        ("magnetized_rs_vegas_wind", command("scripts/benchmarks/magnetized_rs_vegas_lc_compare.py", "--mode", mode, "--medium", "wind", "--output", FIGURE_ROOT / "magnetized_rs_vegas" / "wind", "--data-dir", DATA_ROOT / "magnetized_rs_vegas" / "wind")),
        ("density_jump", command("scripts/benchmarks/triple_density_jump_rs_fs_tophat.py", "--mode", mode, "--output", FIGURE_ROOT / "density_jump" / "triple_density_jump", "--data-dir", DATA_ROOT / "density_jump")),
        ("skymap", command("tests/benchmark_skymap_centroid_motion.py", "--mode", mode, "--data-dir", DATA_ROOT / "skymap", "--figure-dir", FIGURE_ROOT / "skymap")),
        ("theta_decay", command("tests/benchmark_theta_j_multiples_magnetic_decay.py", "--mode", mode, "--data-dir", DATA_ROOT / "theta_decay", "--figure-dir", FIGURE_ROOT / "theta_decay")),
        ("angular_physical", command("tests/angular_sampling_physical_compare.py", "--mode", mode, "--output-dir", DATA_ROOT / "angular_physical", "--figure-dir", FIGURE_ROOT / "angular_physical")),
        ("angular_reference", command("tests/angular_sampling_full_reference_compare.py", "--mode", mode, "--output-dir", DATA_ROOT / "angular_reference", "--figure-dir", FIGURE_ROOT / "angular_reference")),
        ("magnetized_sigma_ism", command("tests/magnetized_rs_sigma_benchmark.py", "--medium", "ism", "--mode", mode, "--output-dir", FIGURE_ROOT / "magnetized_sigma" / "ism", "--data-dir", DATA_ROOT / "magnetized_sigma")),
        ("magnetized_sigma_wind", command("tests/magnetized_rs_sigma_benchmark.py", "--medium", "wind", "--mode", mode, "--output-dir", FIGURE_ROOT / "magnetized_sigma" / "wind", "--data-dir", DATA_ROOT / "magnetized_sigma")),
        ("cross_code", command("tests/generate_cross_code_benchmarks.py")),
        ("vegas_comparison", command("tests/vegas_afterglow_comparison.py", "--mode", mode, "--data-dir", DATA_ROOT / "vegas_comparison", "--figure-dir", FIGURE_ROOT / "vegas_comparison")),
        ("paper_figures", command("tests/generate_paper_figures.py")),
    ]


def run(mode: str) -> None:
    status = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=ROOT,
        check=True,
        text=True,
        encoding="utf-8",
        capture_output=True,
    ).stdout
    if mode == "formal" and status:
        raise RuntimeError("formal benchmarks require a clean tracked worktree")
    env = os.environ | {
        "ASGARD_BENCHMARK_MODE": mode,
        "ASGARD_BENCHMARK_DATA": str(DATA_ROOT / "cross_code"),
        "MPLBACKEND": "Agg",
    }
    repeats = 7 if mode == "formal" else 1
    system = environment(mode, threads=int(env.get("OMP_NUM_THREADS", "1")), grid="suite-specific", repeats=repeats)
    records = []
    for name, argv in suites(mode):
        print(f"[{name}] {' '.join(argv)}", flush=True)
        start = time.perf_counter()
        subprocess.run(argv, cwd=ROOT, env=env, check=True)
        records.append({"suite": name, "command": argv, "wall_seconds": time.perf_counter() - start})
    write_json(DATA_ROOT / "run_all.json", {"environment": system, "suites": records})
    write_json(DATA_ROOT / "manifest.json", {"environment": system, **manifest()})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("quick", "formal"), default="quick")
    args = parser.parse_args()
    run(args.mode)


if __name__ == "__main__":
    main()
