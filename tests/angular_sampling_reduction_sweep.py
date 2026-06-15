from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
import time

import matplotlib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _parse_grid(value: str) -> tuple[int, int]:
    left, right = value.lower().split("x", 1)
    return int(left), int(right)


def _make_model(theta_obs: float, *, sampling: str, patch_theta: int, patch_phi: int):
    from ASGARD import GaussianJet, ISM, Model, Observer, Radiation, Setups

    theta_c = 0.08
    return Model(
        GaussianJet(E_iso=1.0e52, Gamma0=120.0, theta_c=theta_c, theta_max=0.24),
        ISM(1.0),
        Observer(1.0e26, 0.1, theta_obs),
        Radiation(0.1, 1.0e-3, 2.3, ssc=False),
        setups=Setups(
            structured_backend="python_patch",
            patch_sampling=sampling,
            patch_projection="surface_element",
            patch_theta=patch_theta,
            patch_phi=patch_phi,
            num_gam_e=21,
            num_nu=17,
            num_r=16,
            num_theta=12,
            num_phi=12,
            num_tobs=12,
            num_threads=1,
            electron_solver="fullhide_1d",
            electron_adaptive_substeps=False,
        ),
    )


def _metrics(candidate: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(candidate) & np.isfinite(reference) & (reference > 0.0)
    rel = candidate[mask] / reference[mask] - 1.0
    abs_rel = np.abs(rel)
    return {
        "max_abs": float(np.max(abs_rel)),
        "p95_abs": float(np.percentile(abs_rel, 95.0)),
        "median_abs": float(np.median(abs_rel)),
    }


def _write_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _plot_summary(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summary = [row for row in rows if row["theta_obs_over_theta_c"] == "all"]
    patch_count = np.asarray([row["patch_count"] for row in summary], dtype=float)
    max_abs = np.asarray([row["max_abs"] for row in summary], dtype=float)
    p95_abs = np.asarray([row["p95_abs"] for row in summary], dtype=float)
    median_abs = np.asarray([row["median_abs"] for row in summary], dtype=float)
    reduction = np.asarray([row["cost_reduction"] for row in summary], dtype=float)

    order = np.argsort(patch_count)
    fig, ax = plt.subplots(figsize=(6.2, 4.0), constrained_layout=True)
    ax.plot(patch_count[order], max_abs[order], marker="o", color="#D55E00", label="max")
    ax.plot(patch_count[order], p95_abs[order], marker="s", color="#0072B2", label="p95")
    ax.plot(patch_count[order], median_abs[order], marker="^", color="#009E73", label="median")
    for x, y, r in zip(patch_count[order], p95_abs[order], reduction[order]):
        ax.annotate(f"{100.0 * r:.0f}% saved", (x, y), textcoords="offset points", xytext=(4, 5), fontsize=8)
    ax.set_xlabel("candidate patch count")
    ax.set_ylabel(r"$|F_{candidate}/F_{uniform300}-1|$")
    ax.set_yscale("log")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-grid", default="20x15")
    parser.add_argument("--candidate-sampling", default="dominant_region_ioka_v1")
    parser.add_argument("--candidate-grids", nargs="+", default=["8x8", "10x10", "12x12", "15x10"])
    parser.add_argument("--theta-c", type=float, default=0.08)
    parser.add_argument("--ratios", type=float, nargs="+", default=[0.0, 0.5, 1.0, 1.5, 2.0])
    parser.add_argument("--time-min", type=float, default=1.0e3)
    parser.add_argument("--time-max", type=float, default=1.0e6)
    parser.add_argument("--num-times", type=int, default=25)
    parser.add_argument("--freqs", type=float, nargs="+", default=[1.0e10, 1.0e14, 1.0e17])
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output" / "asgard_doc" / "angular_sampling_compare")
    args = parser.parse_args()

    reference_theta, reference_phi = _parse_grid(args.reference_grid)
    reference_count = reference_theta * reference_phi
    candidate_grids = [_parse_grid(item) for item in args.candidate_grids]
    times = np.logspace(np.log10(args.time_min), np.log10(args.time_max), int(args.num_times))
    freqs = np.asarray(args.freqs, dtype=float)
    ratios = np.asarray(args.ratios, dtype=float)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    reference_flux = []
    reference_start = time.perf_counter()
    for ratio in ratios:
        model = _make_model(
            float(ratio * args.theta_c),
            sampling="uniform",
            patch_theta=reference_theta,
            patch_phi=reference_phi,
        )
        reference_flux.append(np.asarray(model.flux_density_grid(times, freqs).total, dtype=float))
    reference_elapsed = time.perf_counter() - reference_start
    reference_flux = np.asarray(reference_flux, dtype=float)

    rows: list[dict[str, float | int | str]] = []
    for patch_theta, patch_phi in candidate_grids:
        patch_count = patch_theta * patch_phi
        candidate_flux = []
        candidate_start = time.perf_counter()
        for ratio in ratios:
            model = _make_model(
                float(ratio * args.theta_c),
                sampling=args.candidate_sampling,
                patch_theta=patch_theta,
                patch_phi=patch_phi,
            )
            candidate_flux.append(np.asarray(model.flux_density_grid(times, freqs).total, dtype=float))
        elapsed = time.perf_counter() - candidate_start
        candidate_flux = np.asarray(candidate_flux, dtype=float)
        common = {
            "sampling": args.candidate_sampling,
            "grid": f"{patch_theta}x{patch_phi}",
            "patch_count": patch_count,
            "cost_reduction": 1.0 - patch_count / reference_count,
            "candidate_elapsed_s": elapsed,
            "reference_elapsed_s": reference_elapsed,
        }
        for i_ratio, ratio in enumerate(ratios):
            row = {"theta_obs_over_theta_c": float(ratio), **common, **_metrics(candidate_flux[i_ratio], reference_flux[i_ratio])}
            rows.append(row)
            print(json.dumps(row, ensure_ascii=False))
        row = {"theta_obs_over_theta_c": "all", **common, **_metrics(candidate_flux, reference_flux)}
        rows.append(row)
        print(json.dumps(row, ensure_ascii=False))

    safe_sampling = args.candidate_sampling.replace("_", "-")
    csv_path = args.output_dir / f"reduction_sweep_{safe_sampling}_ref{reference_theta}x{reference_phi}.csv"
    png_path = args.output_dir / f"reduction_sweep_{safe_sampling}_ref{reference_theta}x{reference_phi}.png"
    _write_csv(csv_path, rows)
    _plot_summary(png_path, rows)
    print(json.dumps({"csv": str(csv_path), "png": str(png_path)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
