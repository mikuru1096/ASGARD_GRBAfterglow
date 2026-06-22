from __future__ import annotations

import argparse
import csv
from pathlib import Path
from time import perf_counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from _repo_path import ensure_repo_root_on_path


ensure_repo_root_on_path()

from asgard_core import Observer, top_hat_jet
from asgard_core.api_model import _direct_tophat_patch_config, _solve_patch_state
from asgard_core.asgard_state import project_flux_grid
from tests.public_api_builders import numerics, observer_grid, radiation, solver_options, top_hat_model


THETA_J_RAD = 0.1
THETA_MULTIPLES = np.array([2.0, 3.0, 4.0, 5.0], dtype=float)
THETA_VALUES = THETA_J_RAD * THETA_MULTIPLES
LC_TIMES_S = np.geomspace(1.0e2, 1.0e7, 88)
LC_FREQS_HZ = np.array([1.0e10, 1.0e14, 1.0e18], dtype=float)
SPEC_TIME_S = np.array([1.0e4], dtype=float)
SPEC_FREQS_HZ = np.geomspace(1.0e8, 1.0e20, 112)
REQUESTED_FREQS_HZ = np.unique(np.concatenate((LC_FREQS_HZ, SPEC_FREQS_HZ)))
SOLVE_TIMES_S = np.unique(np.concatenate((LC_TIMES_S, SPEC_TIME_S)))

EPSILON_B = 1.0e-3
EPSILON_B_FLOOR = 1.0e-5
MAGNETIC_DECAY_ALPHA_T = -0.4
MAGNETIC_DECAY_T0_S = 1.0

COLORS = ("#0072b2", "#d55e00", "#009e73", "#cc79a7")
BAND_LABELS = {
    1.0e10: r"$10^{10}$ Hz",
    1.0e14: r"$10^{14}$ Hz",
    1.0e18: r"$10^{18}$ Hz",
}


def benchmark_radiation(*, magnetic_decay: bool):
    params = dict(
        epsilon_e=0.1,
        epsilon_B=EPSILON_B,
        p=2.3,
        accelerated_electron_fraction=0.1,
        include_ssc=False,
        include_kn_correction=False,
    )
    if magnetic_decay:
        params |= dict(
            epsilon_b_floor=EPSILON_B_FLOOR,
            magnetic_decay_alpha_t=MAGNETIC_DECAY_ALPHA_T,
            magnetic_decay_t0_s=MAGNETIC_DECAY_T0_S,
        )
    return radiation(**params)


def build_model(solver: str, theta_v: float):
    is_2d = solver == "fullhide_2d"
    return top_hat_model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=300.0,
            opening_angle_rad=THETA_J_RAD,
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        observer=Observer(
            z=0.1,
            viewing_angle_rad=theta_v,
            viewing_azimuth_rad=0.0,
            luminosity_distance_cm=1.0e26,
        ),
        fwd_rad=benchmark_radiation(magnetic_decay=is_2d),
        rvs_rad=None,
        numerics=numerics(
            num_radius=88,
            eats_num_theta=36,
            eats_num_phi=24,
            num_observer_time=88,
            num_electron_gamma=61,
            num_photon_frequency=72,
            downstream_num_chi=12 if is_2d else None,
            num_threads=1,
            electron_adaptive_substeps=True,
            electron_substep_rtol=0.02,
            electron_substep_min=16,
            electron_substep_max=256,
            initial_radius_cm=1.0e14,
        ),
        observer_grid=observer_grid(time_min_s=float(LC_TIMES_S[0]), time_max_s=float(LC_TIMES_S[-1])),
        solver_options=solver_options(
            electron_solver=solver,
            geometry_projection="chi_eats_2d" if is_2d else "sed_legacy",
            ssc_cooling_mode="none",
        ),
    )


def run_1d_angle(theta_v: float) -> tuple[np.ndarray, np.ndarray, float]:
    model = build_model("fullhide_1d", theta_v)
    started = perf_counter()
    lc_flux = model.flux_density_grid(LC_TIMES_S, LC_FREQS_HZ, projection_kind="lightcurve").total
    spec_flux = model.flux_density_grid(SPEC_TIME_S, SPEC_FREQS_HZ, projection_kind="sed").total[:, 0]
    elapsed = perf_counter() - started
    return lc_flux, spec_flux, elapsed


def solve_2d_state():
    model = build_model("fullhide_2d", float(THETA_VALUES[0]))
    config = _direct_tophat_patch_config(model)
    timings: dict[str, float] = {}
    started = perf_counter()
    state = _solve_patch_state(
        model,
        config,
        SOLVE_TIMES_S,
        REQUESTED_FREQS_HZ,
        timings=timings,
        solve_reference_times_s=SOLVE_TIMES_S,
    )
    elapsed = perf_counter() - started
    return state, elapsed


def project_2d_angle(state, theta_v: float) -> tuple[np.ndarray, np.ndarray, float]:
    state.config.theta_v = float(theta_v)
    started = perf_counter()
    lc_state = project_flux_grid(state, LC_TIMES_S, LC_FREQS_HZ, projection_kind="lightcurve")
    spec_state = project_flux_grid(state, SPEC_TIME_S, SPEC_FREQS_HZ, projection_kind="sed")
    elapsed = perf_counter() - started
    return (
        np.asarray(lc_state.components["total"], dtype=float),
        np.asarray(spec_state.components["total"], dtype=float)[:, 0],
        elapsed,
    )


def compute() -> dict[str, np.ndarray]:
    flux_1d_lc = np.empty((THETA_VALUES.size, LC_FREQS_HZ.size, LC_TIMES_S.size), dtype=float)
    flux_2d_lc = np.empty_like(flux_1d_lc)
    flux_1d_spec = np.empty((THETA_VALUES.size, SPEC_FREQS_HZ.size), dtype=float)
    flux_2d_spec = np.empty_like(flux_1d_spec)
    one_d_elapsed = np.empty(THETA_VALUES.size, dtype=float)
    two_d_project_elapsed = np.empty(THETA_VALUES.size, dtype=float)

    for i, theta_v in enumerate(THETA_VALUES):
        print(f"running 1d theta/theta_j={THETA_MULTIPLES[i]:.0f}", flush=True)
        flux_1d_lc[i], flux_1d_spec[i], one_d_elapsed[i] = run_1d_angle(float(theta_v))

    print("solving 2d q-shell state with magnetic decay", flush=True)
    state, two_d_solve_elapsed = solve_2d_state()
    for i, theta_v in enumerate(THETA_VALUES):
        print(f"projecting 2d theta/theta_j={THETA_MULTIPLES[i]:.0f}", flush=True)
        flux_2d_lc[i], flux_2d_spec[i], two_d_project_elapsed[i] = project_2d_angle(state, float(theta_v))

    return {
        "theta_j_rad": np.array([THETA_J_RAD], dtype=float),
        "theta_multiples": THETA_MULTIPLES,
        "theta_values": THETA_VALUES,
        "lc_times_s": LC_TIMES_S,
        "lc_freqs_hz": LC_FREQS_HZ,
        "spec_time_s": SPEC_TIME_S,
        "spec_freqs_hz": SPEC_FREQS_HZ,
        "epsilon_B": np.array([EPSILON_B], dtype=float),
        "epsilon_b_floor": np.array([EPSILON_B_FLOOR], dtype=float),
        "magnetic_decay_alpha_t": np.array([MAGNETIC_DECAY_ALPHA_T], dtype=float),
        "magnetic_decay_t0_s": np.array([MAGNETIC_DECAY_T0_S], dtype=float),
        "flux_1d_lc": flux_1d_lc,
        "flux_2d_lc": flux_2d_lc,
        "flux_1d_spec": flux_1d_spec,
        "flux_2d_spec": flux_2d_spec,
        "one_d_elapsed_s": one_d_elapsed,
        "two_d_solve_once_s": np.array([two_d_solve_elapsed], dtype=float),
        "two_d_project_elapsed_s": two_d_project_elapsed,
    }


def ratio_stats(numerator: np.ndarray, denominator: np.ndarray) -> dict[str, float | int]:
    mask = np.isfinite(numerator) & np.isfinite(denominator) & (numerator > 0.0) & (denominator > 0.0)
    ratio = numerator[mask] / denominator[mask]
    log_ratio = np.log10(ratio)
    return {
        "num_points": int(ratio.size),
        "ratio_min": float(np.min(ratio)),
        "ratio_median": float(np.median(ratio)),
        "ratio_max": float(np.max(ratio)),
        "log10_ratio_rms": float(np.sqrt(np.mean(log_ratio * log_ratio))),
    }


def write_compare_summary(data: dict[str, np.ndarray], output: Path) -> None:
    rows: list[dict[str, object]] = []
    for i, theta_v in enumerate(THETA_VALUES):
        common = {
            "theta_v_rad": float(theta_v),
            "theta_v_over_theta_j": float(THETA_MULTIPLES[i]),
            "epsilon_B": EPSILON_B,
            "epsilon_b_floor": EPSILON_B_FLOOR,
            "magnetic_decay_alpha_t": MAGNETIC_DECAY_ALPHA_T,
            "magnetic_decay_t0_s": MAGNETIC_DECAY_T0_S,
            "one_d_elapsed_s": float(data["one_d_elapsed_s"][i]),
            "two_d_solve_once_s": float(data["two_d_solve_once_s"][0]),
            "two_d_project_elapsed_s": float(data["two_d_project_elapsed_s"][i]),
        }
        for j, nu in enumerate(LC_FREQS_HZ):
            rows.append(
                {
                    "kind": "lightcurve",
                    **common,
                    "label": f"nu={nu:.3e} Hz",
                    **ratio_stats(data["flux_2d_lc"][i, j], data["flux_1d_lc"][i, j]),
                }
            )
        rows.append(
            {
                "kind": "spectrum",
                **common,
                "label": f"t={SPEC_TIME_S[0]:.3e} s",
                **ratio_stats(data["flux_2d_spec"][i], data["flux_1d_spec"][i]),
            }
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_peak_summary(data: dict[str, np.ndarray], output: Path) -> None:
    rows: list[dict[str, object]] = []
    for i, theta_v in enumerate(THETA_VALUES):
        for j, nu in enumerate(LC_FREQS_HZ):
            f1 = data["flux_1d_lc"][i, j]
            f2 = data["flux_2d_lc"][i, j]
            idx1 = int(np.nanargmax(f1))
            idx2 = int(np.nanargmax(f2))
            rows.append(
                {
                    "theta_v_rad": float(theta_v),
                    "theta_v_over_theta_j": float(THETA_MULTIPLES[i]),
                    "nu_hz": float(nu),
                    "peak_time_1d_s": float(LC_TIMES_S[idx1]),
                    "peak_time_2d_s": float(LC_TIMES_S[idx2]),
                    "peak_flux_1d": float(f1[idx1]),
                    "peak_flux_2d": float(f2[idx2]),
                    "peak_flux_2d_over_1d": float(f2[idx2] / f1[idx1]),
                    "peak_time_2d_over_1d": float(LC_TIMES_S[idx2] / LC_TIMES_S[idx1]),
                }
            )
    with output.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def positive_ylim(values: list[np.ndarray], floor_fraction: float = 1.0e-8) -> tuple[float, float]:
    merged = np.concatenate([np.asarray(item, dtype=float).reshape(-1) for item in values])
    positive = merged[np.isfinite(merged) & (merged > 0.0)]
    ymax = float(np.max(positive))
    shown = positive[positive >= ymax * floor_fraction]
    ymin = float(np.min(shown))
    return (
        10.0 ** (np.floor(np.log10(ymin)) - 0.25),
        10.0 ** (np.ceil(np.log10(ymax)) + 0.25),
    )


def plot_lightcurves(data: dict[str, np.ndarray], output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.6), sharex=True)
    fig.suptitle(r"1D thin shell vs 2D q-shell with $\epsilon_B$ decay")
    for j, nu in enumerate(LC_FREQS_HZ):
        ax = axes[j]
        panel_values = []
        for i, theta_mult in enumerate(THETA_MULTIPLES):
            color = COLORS[i]
            y1 = data["flux_1d_lc"][i, j]
            y2 = data["flux_2d_lc"][i, j]
            panel_values.extend([y1, y2])
            ax.loglog(LC_TIMES_S, y1, color=color, ls="--", lw=1.5, alpha=0.85)
            ax.loglog(LC_TIMES_S, y2, color=color, ls="-", lw=1.9, marker="o", markevery=12, ms=3.0, mfc="none")
        ax.set_title(BAND_LABELS[float(nu)])
        ax.set_xlabel(r"$t_{\rm obs}$ (s)")
        ax.set_ylim(*positive_ylim(panel_values))
        ax.grid(True, which="both", alpha=0.18)
        ax.tick_params(which="both", direction="in", top=True, right=True)
    axes[0].set_ylabel(r"$F_\nu$ (erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$)")
    angle_handles = [
        Line2D([0], [0], color=COLORS[i], lw=2.0, label=rf"$\theta_v/\theta_j={theta_mult:.0f}$")
        for i, theta_mult in enumerate(THETA_MULTIPLES)
    ]
    solver_handles = [
        Line2D([0], [0], color="black", lw=1.8, ls="--", label="1D thin shell"),
        Line2D([0], [0], color="black", lw=2.0, ls="-", marker="o", mfc="none", ms=3.2, label="2D q-shell B decay"),
    ]
    axes[0].legend(handles=angle_handles, frameon=False, fontsize=8, loc="lower left")
    axes[1].legend(handles=solver_handles, frameon=False, fontsize=8, loc="lower left")
    axes[2].text(
        0.04,
        0.06,
        rf"$\alpha_t={MAGNETIC_DECAY_ALPHA_T:.1f}$, $\epsilon_{{B,\rm floor}}={EPSILON_B_FLOOR:.0e}$",
        transform=axes[2].transAxes,
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(output, dpi=180)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def plot_spectra(data: dict[str, np.ndarray], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    for i, theta_mult in enumerate(THETA_MULTIPLES):
        color = COLORS[i]
        y1 = SPEC_FREQS_HZ * data["flux_1d_spec"][i]
        y2 = SPEC_FREQS_HZ * data["flux_2d_spec"][i]
        ax.loglog(SPEC_FREQS_HZ, y1, color=color, ls="--", lw=1.5, alpha=0.85)
        ax.loglog(SPEC_FREQS_HZ, y2, color=color, ls="-", lw=1.9, marker="s", markevery=14, ms=3.0, mfc="none")
    ax.set_xlabel(r"$\nu$ (Hz)")
    ax.set_ylabel(r"$\nu F_\nu$ (erg cm$^{-2}$ s$^{-1}$)")
    ax.set_title(r"Spectrum at $t_{\rm obs}=10^4$ s")
    ax.set_ylim(*positive_ylim([SPEC_FREQS_HZ * data["flux_1d_spec"], SPEC_FREQS_HZ * data["flux_2d_spec"]]))
    ax.grid(True, which="both", alpha=0.18)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    angle_handles = [
        Line2D([0], [0], color=COLORS[i], lw=2.0, label=rf"$\theta_v/\theta_j={theta_mult:.0f}$")
        for i, theta_mult in enumerate(THETA_MULTIPLES)
    ]
    solver_handles = [
        Line2D([0], [0], color="black", lw=1.8, ls="--", label="1D thin shell"),
        Line2D([0], [0], color="black", lw=2.0, ls="-", marker="s", mfc="none", ms=3.2, label="2D q-shell B decay"),
    ]
    leg1 = ax.legend(handles=angle_handles, frameon=False, fontsize=8, loc="lower left")
    ax.add_artist(leg1)
    ax.legend(handles=solver_handles, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout(rect=(0, 0, 0.80, 1))
    fig.savefig(output, dpi=180)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def plot_multiband(data: dict[str, np.ndarray], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.2), sharex=True)
    for i, ax in enumerate(axes.flat):
        panel_values = []
        for j, nu in enumerate(LC_FREQS_HZ):
            color = ("#0072b2", "#009e73", "#d55e00")[j]
            y1 = data["flux_1d_lc"][i, j]
            y2 = data["flux_2d_lc"][i, j]
            panel_values.extend([y1, y2])
            ax.loglog(LC_TIMES_S, y1, color=color, ls="--", lw=1.35, alpha=0.85)
            ax.loglog(LC_TIMES_S, y2, color=color, ls="-", lw=1.75, marker="o", markevery=12, ms=2.8, mfc="none")
        ax.set_title(rf"$\theta_v/\theta_j={THETA_MULTIPLES[i]:.0f}$")
        ax.set_ylim(*positive_ylim(panel_values))
        ax.grid(True, which="both", alpha=0.18)
        ax.tick_params(which="both", direction="in", top=True, right=True)
    for ax in axes[1]:
        ax.set_xlabel(r"$t_{\rm obs}$ (s)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$F_\nu$ (erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$)")
    band_handles = [
        Line2D([0], [0], color=color, lw=2.0, label=label)
        for color, label in zip(("#0072b2", "#009e73", "#d55e00"), (BAND_LABELS[1.0e10], BAND_LABELS[1.0e14], BAND_LABELS[1.0e18]))
    ]
    solver_handles = [
        Line2D([0], [0], color="black", lw=1.8, ls="--", label="1D"),
        Line2D([0], [0], color="black", lw=2.0, ls="-", marker="o", mfc="none", ms=3.0, label="2D B decay"),
    ]
    axes[0, 0].legend(handles=band_handles, frameon=False, fontsize=8, loc="lower left")
    axes[0, 1].legend(handles=solver_handles, frameon=False, fontsize=8, loc="lower left")
    fig.suptitle(r"Band-by-band light curves")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output, dpi=180)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def write_outputs(data: dict[str, np.ndarray], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez(output_dir / "theta_j_multiples_bdecay_compare_data.npz", **data)
    write_compare_summary(data, output_dir / "theta_j_multiples_bdecay_compare_summary.csv")
    write_peak_summary(data, output_dir / "theta_j_multiples_bdecay_peak_summary.csv")
    plot_lightcurves(data, output_dir / "theta_j_multiples_bdecay_lightcurve_1d_vs_qshell.png")
    plot_spectra(data, output_dir / "theta_j_multiples_bdecay_spectrum_1d_vs_qshell.png")
    plot_multiband(data, output_dir / "theta_j_multiples_bdecay_multiband_lightcurve_1d_vs_qshell.png")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/benchmark_1d_vs_qshell_theta_j_multiples_bdecay_alpha04"),
    )
    args = parser.parse_args()
    data = compute()
    write_outputs(data, args.output_dir)
    print(args.output_dir)


if __name__ == "__main__":
    main()
