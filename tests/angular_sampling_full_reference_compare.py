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

from historical_angular_patch import make_gaussian_model, patch_flux_grid
from scripts.benchmarks.benchmark_common import PALETTE, plot_style, save_figure


def _make_model(
    theta_obs: float,
    *,
    gamma0: float,
    theta_c: float,
    theta_max: float,
    backend: str,
    sampling: str,
    patch_theta: int,
    patch_phi: int,
    beaming_factor: float = 3.0,
    beaming_resolution: float = 8.0,
):
    del backend, sampling, patch_theta, patch_phi
    return make_gaussian_model(
        theta_obs,
        gamma0=gamma0,
        theta_c=theta_c,
        theta_max=theta_max,
        beaming_factor=beaming_factor,
        beaming_resolution=beaming_resolution,
    )


def _python_surface_axisym_reference(
    theta_obs: float,
    times: np.ndarray,
    freqs: np.ndarray,
    *,
    gamma0: float,
    theta_c: float,
    theta_max: float,
    patch_theta: int,
    patch_phi: int,
) -> np.ndarray:
    reference_model = make_gaussian_model(theta_obs, gamma0=gamma0, theta_c=theta_c, theta_max=theta_max)
    flux, _grid = patch_flux_grid(
        reference_model,
        times,
        freqs,
        sampling="uniform",
        patch_theta=patch_theta,
        patch_phi=patch_phi,
    )
    return flux


def _metrics(candidate: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(candidate) & np.isfinite(reference) & (reference > 0.0)
    if not np.any(mask):
        raise AssertionError("no positive finite reference flux samples")
    rel = candidate[mask] / reference[mask] - 1.0
    abs_rel = np.abs(rel)
    return {
        "max_abs": float(np.max(abs_rel)),
        "p95_abs": float(np.percentile(abs_rel, 95.0)),
        "median_abs": float(np.median(abs_rel)),
    }


def _positive_ylim(values: list[np.ndarray]) -> tuple[float, float]:
    positive = np.concatenate([np.ravel(np.asarray(value, dtype=float)) for value in values])
    positive = positive[np.isfinite(positive) & (positive > 0.0)]
    if positive.size == 0:
        raise AssertionError("plot has no positive finite flux values")
    return float(np.min(positive)) / 2.0, float(np.max(positive)) * 2.0


def _write_metrics_csv(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=["theta_obs_over_theta_c", "max_abs", "p95_abs", "median_abs"],
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_lightcurve_csv(
    path: Path,
    ratios: np.ndarray,
    times: np.ndarray,
    freqs: np.ndarray,
    reference_flux: np.ndarray,
    candidate_flux: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["theta_obs_over_theta_c", "time_s", "frequency_hz", "reference", "candidate", "relative"])
        for i_ratio, ratio in enumerate(ratios):
            for i_freq, freq in enumerate(freqs):
                for i_time, time_s in enumerate(times):
                    reference = float(reference_flux[i_ratio, i_freq, i_time])
                    candidate = float(candidate_flux[i_ratio, i_freq, i_time])
                    writer.writerow([float(ratio), float(time_s), float(freq), reference, candidate, candidate / reference - 1.0])


def _plot_lightcurves(
    path: Path,
    ratios: np.ndarray,
    times: np.ndarray,
    freqs: np.ndarray,
    reference_flux: np.ndarray,
    candidate_flux: np.ndarray,
    candidate_label: str,
) -> None:
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(plot_style())
    fig, axes = plt.subplots(
        len(ratios),
        len(freqs),
        figsize=(3.6 * len(freqs), 2.3 * len(ratios)),
        sharex=True,
        constrained_layout=True,
    )
    axes = np.asarray(axes).reshape(len(ratios), len(freqs))
    for i_ratio, ratio in enumerate(ratios):
        for i_freq, freq in enumerate(freqs):
            ax = axes[i_ratio, i_freq]
            reference = reference_flux[i_ratio, i_freq]
            candidate = candidate_flux[i_ratio, i_freq]
            ax.loglog(times, reference, color=PALETTE["blue"], marker="o", markevery=4, label="uniform reference")
            ax.loglog(times, candidate, color=PALETTE["vermillion"], ls="--", marker="s", markevery=4, label=candidate_label)
            ax.set_ylim(*_positive_ylim([reference, candidate]))
            if i_ratio == 0:
                ax.set_title(f"{freq:.0e} Hz")
            if i_freq == 0:
                ax.set_ylabel(rf"$\theta_{{obs}}/\theta_c={ratio:g}$" + "\n" + r"$F_\nu$")
            if i_ratio == len(ratios) - 1:
                ax.set_xlabel("observer time [s]")
            if i_ratio == 0 and i_freq == len(freqs) - 1:
                ax.legend(frameon=False)
    outputs = save_figure(fig, path)
    plt.close(fig)
    return outputs


def _plot_relative(
    path: Path,
    ratios: np.ndarray,
    times: np.ndarray,
    freqs: np.ndarray,
    reference_flux: np.ndarray,
    candidate_flux: np.ndarray,
) -> None:
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(plot_style())
    fig, axes = plt.subplots(
        len(ratios),
        len(freqs),
        figsize=(3.6 * len(freqs), 2.0 * len(ratios)),
        sharex=True,
        constrained_layout=True,
    )
    axes = np.asarray(axes).reshape(len(ratios), len(freqs))
    for i_ratio, ratio in enumerate(ratios):
        for i_freq, freq in enumerate(freqs):
            ax = axes[i_ratio, i_freq]
            rel = candidate_flux[i_ratio, i_freq] / reference_flux[i_ratio, i_freq] - 1.0
            ax.semilogx(times, rel, color=PALETTE["green"], marker="o", markevery=4)
            ax.axhline(0.0, color="0.4", lw=0.8, ls=":")
            pad = max(float(np.max(np.abs(rel))) * 1.15, 1.0e-3)
            ax.set_ylim(-pad, pad)
            if i_ratio == 0:
                ax.set_title(f"{freq:.0e} Hz")
            if i_freq == 0:
                ax.set_ylabel(rf"$\theta_{{obs}}/\theta_c={ratio:g}$" + "\n" + r"$F_{ad}/F_{ref}-1$")
            if i_ratio == len(ratios) - 1:
                ax.set_xlabel("observer time [s]")
    outputs = save_figure(fig, path)
    plt.close(fig)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("quick", "formal"), default="formal")
    parser.add_argument("--theta-c", type=float, default=0.08)
    parser.add_argument("--theta-max", type=float, default=0.24)
    parser.add_argument("--gamma0", type=float, default=120.0)
    parser.add_argument("--beaming-factor", type=float, default=3.0)
    parser.add_argument("--beaming-resolution", type=float, default=8.0)
    parser.add_argument("--ratios", type=float, nargs="+")
    parser.add_argument("--time-min", type=float, default=1.0e3)
    parser.add_argument("--time-max", type=float, default=1.0e6)
    parser.add_argument("--num-times", type=int)
    parser.add_argument("--freqs", type=float, nargs="+")
    parser.add_argument("--candidate-sampling", default="dominant_region_ioka_time_v1")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    args = parser.parse_args()

    formal = args.mode == "formal"
    ref_theta, ref_phi = ((300, 48) if formal else (30, 12))
    cand_theta, cand_phi = ((20, 15) if formal else (8, 6))
    args.num_times = args.num_times or (25 if formal else 5)
    args.ratios = args.ratios or ([0.0, 0.5, 1.0, 1.5, 2.0] if formal else [0.0, 1.0])
    args.freqs = args.freqs or ([1.0e10, 1.0e14, 1.0e17] if formal else [1.0e14])
    args.output_dir = args.output_dir or ROOT / "paper" / "source_data" / "benchmarks" / "angular_sampling_full"
    figure_dir = args.figure_dir or ROOT / "paper" / "figures" / "benchmarks" / "angular_sampling_full"

    times = np.logspace(np.log10(args.time_min), np.log10(args.time_max), int(args.num_times))
    freqs = np.asarray(args.freqs, dtype=float)
    ratios = np.asarray(args.ratios, dtype=float)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    reference_flux = np.empty((ratios.size, freqs.size, times.size), dtype=float)
    candidate_flux = np.empty_like(reference_flux)
    candidate_patch_counts: list[int] = []
    candidate_theta_counts: list[int] = []
    candidate_phi_counts: list[int] = []

    reference_start = time.perf_counter()
    for i_ratio, ratio in enumerate(ratios):
        reference_flux[i_ratio] = _python_surface_axisym_reference(
            float(ratio * args.theta_c),
            times,
            freqs,
            gamma0=args.gamma0,
            theta_c=args.theta_c,
            theta_max=args.theta_max,
            patch_theta=ref_theta,
            patch_phi=ref_phi,
        )
    reference_elapsed = time.perf_counter() - reference_start

    candidate_start = time.perf_counter()
    for i_ratio, ratio in enumerate(ratios):
        model = _make_model(
            float(ratio * args.theta_c),
            gamma0=args.gamma0,
            theta_c=args.theta_c,
            theta_max=args.theta_max,
            backend="python_patch",
            sampling=args.candidate_sampling,
            patch_theta=cand_theta,
            patch_phi=cand_phi,
            beaming_factor=args.beaming_factor,
            beaming_resolution=args.beaming_resolution,
        )
        candidate, grid = patch_flux_grid(
            model,
            times,
            freqs,
            sampling=args.candidate_sampling,
            patch_theta=cand_theta,
            patch_phi=cand_phi,
        )
        candidate_flux[i_ratio] = np.asarray(candidate, dtype=float)
        candidate_patch_counts.append(int(grid.domega.size))
        candidate_theta_counts.append(int(grid.theta_centers.size))
        candidate_phi_counts.append(int(grid.phi_centers.size))
    candidate_elapsed = time.perf_counter() - candidate_start

    rows = []
    for i_ratio, ratio in enumerate(ratios):
        row = {"theta_obs_over_theta_c": float(ratio), **_metrics(candidate_flux[i_ratio], reference_flux[i_ratio])}
        rows.append(row)
        print(json.dumps(row, ensure_ascii=False))
    summary = {"theta_obs_over_theta_c": "all", **_metrics(candidate_flux, reference_flux)}
    print(json.dumps(summary, ensure_ascii=False))

    safe_sampling = args.candidate_sampling.replace("_", "-")
    safe_gamma0 = f"g{args.gamma0:g}".replace(".", "p")
    safe_theta_c = f"tc{args.theta_c:g}".replace(".", "p")
    stem = f"fullref{ref_theta}x{ref_phi}_vs_{safe_sampling}_{cand_theta}x{cand_phi}_{safe_gamma0}_{safe_theta_c}"
    metrics_csv = args.output_dir / f"{stem}_metrics.csv"
    lightcurve_csv = args.output_dir / f"{stem}_lightcurves.csv"
    npz_path = args.output_dir / f"{stem}.npz"
    lightcurve_path = figure_dir / f"{stem}_lightcurve_compare"
    relative_path = figure_dir / f"{stem}_relative_error"
    _write_metrics_csv(metrics_csv, rows)
    _write_lightcurve_csv(lightcurve_csv, ratios, times, freqs, reference_flux, candidate_flux)
    np.savez(
        npz_path,
        ratios=ratios,
        times=times,
        freqs=freqs,
        reference_flux=reference_flux,
        candidate_flux=candidate_flux,
        reference_elapsed_s=reference_elapsed,
        candidate_elapsed_s=candidate_elapsed,
        beaming_factor=args.beaming_factor,
        beaming_resolution=args.beaming_resolution,
        candidate_patch_counts=np.asarray(candidate_patch_counts, dtype=int),
        candidate_theta_counts=np.asarray(candidate_theta_counts, dtype=int),
        candidate_phi_counts=np.asarray(candidate_phi_counts, dtype=int),
    )
    candidate_label = (
        f"adaptive{candidate_theta_counts[0]}x{candidate_phi_counts[0]}"
        if len(set(candidate_theta_counts)) == 1 and len(set(candidate_phi_counts)) == 1
        else "adaptive-auto"
    )
    lightcurves = _plot_lightcurves(lightcurve_path, ratios, times, freqs, reference_flux, candidate_flux, candidate_label)
    relatives = _plot_relative(relative_path, ratios, times, freqs, reference_flux, candidate_flux)
    print(json.dumps(
        {
            "reference": f"python_surface_axisym uniform {ref_theta}x{ref_phi}",
            "candidate": f"python_patch {args.candidate_sampling} patch_theta={cand_theta} patch_phi_min={cand_phi}",
            "gamma0": args.gamma0,
            "theta_c": args.theta_c,
            "theta_max": args.theta_max,
            "beaming_factor": args.beaming_factor,
            "beaming_resolution": args.beaming_resolution,
            "reference_elapsed_s": reference_elapsed,
            "candidate_elapsed_s": candidate_elapsed,
            "patch_reduction_by_ratio": [1.0 - count / (ref_theta * ref_phi) for count in candidate_patch_counts],
            "candidate_patch_counts": candidate_patch_counts,
            "candidate_theta_counts": candidate_theta_counts,
            "candidate_phi_counts": candidate_phi_counts,
            "metrics_csv": str(metrics_csv),
            "lightcurve_csv": str(lightcurve_csv),
            "npz": str(npz_path),
            "lightcurve_figures": [str(path) for path in lightcurves],
            "relative_figures": [str(path) for path in relatives],
        },
        ensure_ascii=False,
        indent=2,
    ))


if __name__ == "__main__":
    main()
