from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _repo_path import ensure_repo_root_on_path


ensure_repo_root_on_path()

from fullhide_2d_smoke_bench import _thin_shell_projection_inputs
from src import Interpolation


def _top_hat_arrays(num_theta: int, num_phi: int) -> tuple[np.ndarray, np.ndarray]:
    boundary, radius, gamma, r_tobs, seed, obs, tobs, source = _thin_shell_projection_inputs()
    boundary[9] = boundary[8]
    projected = Interpolation.sed_interpolation(
        boundary,
        r_tobs,
        gamma,
        radius,
        source,
        seed,
        obs,
        tobs,
        num_theta,
        num_phi,
        1,
    )
    return tobs, projected[1]


def _structured_axis_arrays(num_theta: int, num_phi: int) -> tuple[np.ndarray, np.ndarray]:
    boundary, radius, gamma, r_tobs, seed, obs, tobs, source = _thin_shell_projection_inputs()
    boundary[9] = boundary[8]
    r_tobs_theta = np.repeat(r_tobs[:, None], num_theta, axis=1)
    gamma_theta = np.repeat(gamma[:, None], num_theta, axis=1)
    radius_theta = np.repeat(radius[:, None], num_theta, axis=1)
    source_theta = np.repeat(source[:, :, None], num_theta, axis=2)
    projected = Interpolation.sed_interpolation_structured(
        boundary,
        0.0,
        r_tobs_theta,
        gamma_theta,
        radius_theta,
        source_theta,
        seed,
        obs,
        tobs,
        num_phi,
        1,
    )
    return tobs, projected[1]


def _structured_phi_arrays(num_theta: int, num_phi: int) -> tuple[np.ndarray, np.ndarray]:
    boundary, radius, gamma, r_tobs, seed, obs, tobs, source = _thin_shell_projection_inputs()
    boundary[9] = boundary[8]
    r_tobs_phi = np.repeat(r_tobs[:, None, None], num_theta, axis=1)
    r_tobs_phi = np.repeat(r_tobs_phi, num_phi, axis=2)
    gamma_phi = np.repeat(gamma[:, None, None], num_theta, axis=1)
    gamma_phi = np.repeat(gamma_phi, num_phi, axis=2)
    radius_phi = np.repeat(radius[:, None, None], num_theta, axis=1)
    radius_phi = np.repeat(radius_phi, num_phi, axis=2)
    source_phi = np.repeat(source[:, :, None, None], num_theta, axis=2)
    source_phi = np.repeat(source_phi, num_phi, axis=3)
    projected = Interpolation.sed_interpolation_structured_phi(
        boundary,
        r_tobs_phi,
        gamma_phi,
        radius_phi,
        source_phi,
        seed,
        obs,
        tobs,
        1,
    )
    return tobs, projected[1]


def _compute() -> dict[str, np.ndarray]:
    tobs, top_hat_low = _top_hat_arrays(12, 12)
    tobs_ref, top_hat_reference = _top_hat_arrays(96, 72)
    tobs_axis, structured_axis_low = _structured_axis_arrays(12, 12)
    tobs_axis_ref, structured_axis_reference = _structured_axis_arrays(96, 72)
    tobs_phi, structured_phi_low = _structured_phi_arrays(12, 12)
    tobs_phi_ref, structured_phi_reference = _structured_phi_arrays(96, 72)
    np.testing.assert_allclose(tobs, tobs_ref, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(tobs, tobs_axis, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(tobs, tobs_axis_ref, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(tobs, tobs_phi, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(tobs, tobs_phi_ref, rtol=0.0, atol=0.0)
    return {
        "tobs": tobs,
        "top_hat_low": top_hat_low,
        "top_hat_reference": top_hat_reference,
        "structured_axis_low": structured_axis_low,
        "structured_axis_reference": structured_axis_reference,
        "structured_phi_low": structured_phi_low,
        "structured_phi_reference": structured_phi_reference,
    }


def _load_baseline(path: Path) -> dict[str, np.ndarray]:
    with np.load(path) as data:
        return {key: data[key] for key in data.files}


def _write_plot(baseline: dict[str, np.ndarray], current: dict[str, np.ndarray], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    time = current["tobs"]
    fig, (ax_flux, ax_err) = plt.subplots(
        2,
        1,
        figsize=(7.6, 6.4),
        sharex=True,
        height_ratios=(2.1, 1.0),
    )
    styles = (
        ("top_hat", "top-hat", "#0173b2"),
        ("structured_axis", "structured axis", "#029e73"),
        ("structured_phi", "structured phi", "#de8f05"),
    )
    for key, label, color in styles:
        ax_flux.loglog(
            time,
            current[f"{key}_reference"],
            color=color,
            lw=1.5,
            label=f"{label} after high grid",
        )
        ax_flux.loglog(
            time,
            current[f"{key}_low"],
            color=color,
            lw=1.1,
            ls="--",
            label=f"{label} after low grid",
        )
        old_rel = baseline[f"{key}_low"] / baseline[f"{key}_reference"] - 1.0
        new_rel = current[f"{key}_low"] / current[f"{key}_reference"] - 1.0
        ax_err.semilogx(time, old_rel, color=color, lw=1.0, ls=":", label=f"{label} before")
        ax_err.semilogx(time, new_rel, color=color, lw=1.4, label=f"{label} after")
    ax_flux.set_ylabel(r"$F_\nu$ (arb.)")
    ax_flux.legend(frameon=False, fontsize=7, ncol=2)
    ax_flux.grid(True, which="both", alpha=0.25)

    ax_err.axhline(0.0, color="black", lw=0.8)
    ax_err.set_xlabel(r"$t_{\rm obs}$ (s)")
    ax_err.set_ylabel("low/high - 1")
    ax_err.legend(frameon=False, fontsize=7, ncol=3)
    ax_err.grid(True, which="both", alpha=0.25)
    fig.suptitle("EATS angular projection: before/after result comparison")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--save-baseline", type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    current = _compute()
    if args.save_baseline is not None:
        args.save_baseline.parent.mkdir(parents=True, exist_ok=True)
        np.savez(args.save_baseline, **current)
    if args.baseline is not None:
        if args.output is None:
            raise SystemExit("--output is required with --baseline")
        _write_plot(_load_baseline(args.baseline), current, args.output)


if __name__ == "__main__":
    main()
