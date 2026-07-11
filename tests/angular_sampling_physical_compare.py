from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from historical_angular_patch import make_gaussian_model, patch_flux_grid
from scripts.benchmarks.benchmark_common import PALETTE, plot_style, save_figure


def _make_model(theta_obs: float, *, sampling: str, patch_theta: int, patch_phi: int):
    del sampling, patch_theta, patch_phi
    return make_gaussian_model(theta_obs)


def _relative_metrics(adaptive: np.ndarray, uniform: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(adaptive) & np.isfinite(uniform) & (uniform > 0.0)
    if not np.any(mask):
        raise AssertionError("no positive finite uniform flux samples")
    rel = adaptive[mask] / uniform[mask] - 1.0
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
    uniform_flux: np.ndarray,
    adaptive_flux: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["theta_obs_over_theta_c", "time_s", "frequency_hz", "uniform", "adaptive", "relative"])
        for i_ratio, ratio in enumerate(ratios):
            for i_freq, freq in enumerate(freqs):
                for i_time, time_s in enumerate(times):
                    uniform = float(uniform_flux[i_ratio, i_freq, i_time])
                    adaptive = float(adaptive_flux[i_ratio, i_freq, i_time])
                    writer.writerow([float(ratio), float(time_s), float(freq), uniform, adaptive, adaptive / uniform - 1.0])


def _plot_lightcurves(
    path: Path,
    ratios: np.ndarray,
    times: np.ndarray,
    freqs: np.ndarray,
    uniform_flux: np.ndarray,
    adaptive_flux: np.ndarray,
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
            uniform = uniform_flux[i_ratio, i_freq]
            adaptive = adaptive_flux[i_ratio, i_freq]
            ax.loglog(times, uniform, color=PALETTE["blue"], marker="o", markevery=4, label="uniform")
            ax.loglog(times, adaptive, color=PALETTE["vermillion"], ls="--", marker="s", markevery=4, label="adaptive")
            ax.set_ylim(*_positive_ylim([uniform, adaptive]))
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


def _plot_relative_errors(
    path: Path,
    ratios: np.ndarray,
    times: np.ndarray,
    freqs: np.ndarray,
    uniform_flux: np.ndarray,
    adaptive_flux: np.ndarray,
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
            rel = adaptive_flux[i_ratio, i_freq] / uniform_flux[i_ratio, i_freq] - 1.0
            ax.semilogx(times, rel, color=PALETTE["green"], marker="o", markevery=4)
            ax.axhline(0.0, color="0.4", lw=0.8, ls=":")
            pad = max(float(np.max(np.abs(rel))) * 1.15, 1.0e-3)
            ax.set_ylim(-pad, pad)
            if i_ratio == 0:
                ax.set_title(f"{freq:.0e} Hz")
            if i_freq == 0:
                ax.set_ylabel(rf"$\theta_{{obs}}/\theta_c={ratio:g}$" + "\n" + r"$F_{ad}/F_{uni}-1$")
            if i_ratio == len(ratios) - 1:
                ax.set_xlabel("observer time [s]")
    outputs = save_figure(fig, path)
    plt.close(fig)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("quick", "formal"), default="formal")
    parser.add_argument("--patch-theta", type=int)
    parser.add_argument("--patch-phi", type=int)
    parser.add_argument("--theta-c", type=float, default=0.08)
    parser.add_argument("--ratios", type=float, nargs="+")
    parser.add_argument("--times", type=float, nargs="+", default=None)
    parser.add_argument("--time-min", type=float, default=1.0e3)
    parser.add_argument("--time-max", type=float, default=1.0e6)
    parser.add_argument("--num-times", type=int)
    parser.add_argument("--freqs", type=float, nargs="+")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--figure-dir", type=Path)
    args = parser.parse_args()

    formal = args.mode == "formal"
    args.patch_theta = args.patch_theta or (20 if formal else 8)
    args.patch_phi = args.patch_phi or (15 if formal else 6)
    args.num_times = args.num_times or (25 if formal else 5)
    args.ratios = args.ratios or ([0.0, 0.5, 1.0, 1.5, 2.0] if formal else [0.0, 1.0])
    args.freqs = args.freqs or ([1.0e10, 1.0e14, 1.0e17] if formal else [1.0e14])
    args.output_dir = args.output_dir or ROOT / "paper" / "source_data" / "benchmarks" / "angular_sampling_physical"
    figure_dir = args.figure_dir or ROOT / "paper" / "figures" / "benchmarks" / "angular_sampling_physical"

    times = (
        np.asarray(args.times, dtype=float)
        if args.times is not None
        else np.logspace(np.log10(args.time_min), np.log10(args.time_max), int(args.num_times))
    )
    freqs = np.asarray(args.freqs, dtype=float)
    ratios = np.asarray(args.ratios, dtype=float)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    uniform_flux = np.empty((ratios.size, freqs.size, times.size), dtype=float)
    adaptive_flux = np.empty_like(uniform_flux)
    for i_ratio, ratio in enumerate(ratios):
        theta_obs = float(ratio * args.theta_c)
        model = _make_model(theta_obs, sampling="uniform", patch_theta=args.patch_theta, patch_phi=args.patch_phi)
        uniform, _uniform_grid = patch_flux_grid(
            model, times, freqs, sampling="uniform", patch_theta=args.patch_theta, patch_phi=args.patch_phi
        )
        adaptive, _adaptive_grid = patch_flux_grid(
            model,
            times,
            freqs,
            sampling="dominant_region_ioka_v1",
            patch_theta=args.patch_theta,
            patch_phi=args.patch_phi,
        )
        uniform_flux[i_ratio] = np.asarray(uniform, dtype=float)
        adaptive_flux[i_ratio] = np.asarray(adaptive, dtype=float)
        metrics = _relative_metrics(adaptive_flux[i_ratio], uniform_flux[i_ratio])
        rows.append({"theta_obs_over_theta_c": float(ratio), **metrics})
        print(json.dumps(rows[-1], ensure_ascii=False))

    stem = f"gaussian_patch{args.patch_theta}x{args.patch_phi}"
    metrics_csv = args.output_dir / f"{stem}_metrics.csv"
    lightcurve_csv = args.output_dir / f"{stem}_lightcurves.csv"
    npz_path = args.output_dir / f"{stem}.npz"
    lightcurve_path = figure_dir / f"{stem}_lightcurve_compare"
    relative_path = figure_dir / f"{stem}_relative_error"
    _write_metrics_csv(metrics_csv, rows)
    _write_lightcurve_csv(lightcurve_csv, ratios, times, freqs, uniform_flux, adaptive_flux)
    np.savez(
        npz_path,
        ratios=ratios,
        times=times,
        freqs=freqs,
        uniform_flux=uniform_flux,
        adaptive_flux=adaptive_flux,
    )
    lightcurves = _plot_lightcurves(lightcurve_path, ratios, times, freqs, uniform_flux, adaptive_flux)
    relatives = _plot_relative_errors(relative_path, ratios, times, freqs, uniform_flux, adaptive_flux)
    print(json.dumps({"summary": rows}, ensure_ascii=False, indent=2))
    print(json.dumps(
        {
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
