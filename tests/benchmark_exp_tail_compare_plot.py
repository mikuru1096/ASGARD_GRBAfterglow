from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _repo_path import ensure_repo_root_on_path


ensure_repo_root_on_path()

from asgard_core import Observer, UniformMedium, WindMedium, top_hat_jet
from tests.public_api_builders import numerics, observer_grid, radiation, solver_options, top_hat_model


SOLVERS = ("fullhide_1d", "fullhide_2d", "slc1_1d", "charint_1d", "charint_2d")
TIMES_S = np.array([1.0e3, 1.0e5, 1.0e7], dtype=float)
FREQUENCIES_HZ = np.logspace(7.0, 30.0, 180)
TIME_COLORS = {1.0e3: "#0072b2", 1.0e5: "#d55e00", 1.0e7: "#009e73"}
SOLVER_STYLES = {
    "fullhide_1d": "-",
    "fullhide_2d": ":",
    "slc1_1d": "--",
    "charint_1d": "-.",
    "charint_2d": (0, (5, 2)),
}
SOLVER_MARKERS = {
    "fullhide_1d": "o",
    "fullhide_2d": "s",
    "slc1_1d": "^",
    "charint_1d": "D",
    "charint_2d": "x",
}
MEDIA = (
    ("ISM Forward Shock", UniformMedium(number_density_cm3=1.0)),
    ("WIND Forward Shock", WindMedium(a_star=0.15, density_floor_cm3=1.0e-4, density_cap_cm3=30.0)),
)


def build_model(solver: str, medium):
    return top_hat_model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=300.0,
            opening_angle_rad=0.10,
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
            accelerated_electron_fraction=0.1,
            include_ssc=True,
            include_kn_correction=False,
        ),
        rvs_rad=None,
        numerics=numerics(
            num_radius=120,
            num_theta=120,
            num_phi=1,
            num_observer_time=120,
            num_electron_gamma=121,
            num_photon_frequency=121,
            num_chi=8 if solver.endswith("_2d") else None,
            num_threads=1,
            electron_adaptive_substeps=True,
            electron_substep_rtol=0.02,
            electron_substep_min=16,
            electron_substep_max=256,
            initial_radius_cm=1.0e14,
        ),
        observer_grid=observer_grid(time_min_s=float(TIMES_S[0]), time_max_s=float(TIMES_S[-1])),
        solver_options=solver_options(electron_solver=solver, ssc_cooling_mode="nakar_y_thomson"),
    )


def run_solver(solver: str, medium) -> dict[str, np.ndarray]:
    model = build_model(solver, medium)
    flux = model.flux_density_grid(TIMES_S, FREQUENCIES_HZ).total
    details = model.details()
    gamma = np.asarray(details.fwd.gamma_e, dtype=float)
    dnde_all = np.asarray(details.fwd.dN_dgamma_e, dtype=float)
    t_track = np.asarray(details.fwd.t_obs, dtype=float)
    indices = np.array([int(np.argmin(np.abs(t_track - t))) for t in TIMES_S], dtype=int)
    return {
        "sed": FREQUENCIES_HZ[:, None] * flux,
        "gamma": gamma,
        "dnde": dnde_all[:, indices],
        "track_time": t_track[indices],
    }


def compute() -> dict[str, dict[str, dict[str, np.ndarray]]]:
    results = {}
    for medium_name, medium in MEDIA:
        medium_results = {}
        for solver in SOLVERS:
            print(f"running {medium_name}: {solver}", flush=True)
            medium_results[solver] = run_solver(solver, medium)
        results[medium_name] = medium_results
    return results


def write_diagnostics(results: dict[str, dict[str, dict[str, np.ndarray]]], output: Path) -> None:
    rows = []
    for medium_name, _medium in MEDIA:
        ref = results[medium_name]["slc1_1d"]
        ref_gamma = ref["gamma"]
        ref_tail = ref_gamma * ref_gamma * ref["dnde"][:, -1]
        ref_mask = (ref_gamma >= 1.0e3) & np.isfinite(ref_tail) & (ref_tail > 0.0)
        ref_x = np.log10(ref_gamma[ref_mask])
        ref_y = np.log10(ref_tail[ref_mask])
        for solver in SOLVERS:
            bundle = results[medium_name][solver]
            gamma = bundle["gamma"]
            final_tail = gamma * gamma * bundle["dnde"][:, -1]
            ratio_mask = (gamma >= 1.0e3) & np.isfinite(final_tail) & (final_tail > 0.0)
            ratio_mask &= (gamma >= ref_gamma[ref_mask][0]) & (gamma <= ref_gamma[ref_mask][-1])
            if np.any(ratio_mask):
                ratio = final_tail[ratio_mask] / (10.0 ** np.interp(np.log10(gamma[ratio_mask]), ref_x, ref_y))
                max_ratio = float(np.nanmax(ratio))
            else:
                max_ratio = np.nan
            for time_index, time_s in enumerate(TIMES_S):
                dnde = bundle["dnde"][:, time_index]
                positive = np.isfinite(dnde) & (dnde > 0.0)
                peak_index = int(np.nanargmax(dnde))
                zero_after = np.flatnonzero((np.arange(dnde.size) > peak_index) & (~positive))
                last_positive = int(np.flatnonzero(positive)[-1]) if np.any(positive) else -1
                rows.append(
                    {
                        "medium": medium_name.split()[0],
                        "solver": solver,
                        "t_obs_s": f"{time_s:.8e}",
                        "nearest_track_time_s": f"{bundle['track_time'][time_index]:.8e}",
                        "peak_dnde": f"{np.nanmax(dnde):.8e}",
                        "last_positive_gamma": f"{gamma[last_positive]:.8e}" if last_positive >= 0 else "nan",
                        "first_zero_after_peak_gamma": (
                            f"{gamma[int(zero_after[0])]:.8e}" if zero_after.size else "nan"
                        ),
                        "max_final_tail_ratio_to_slc1": f"{max_ratio:.8e}",
                    }
                )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def tail_limits(values: list[np.ndarray], ymax: float) -> tuple[float, float] | None:
    if not values or ymax <= 0.0:
        return None
    merged = np.concatenate(values)
    mask = np.isfinite(merged) & (merged > 0.0) & (merged >= ymax * 1.0e-14)
    if not np.any(mask):
        return None
    ymin = 10.0 ** (np.floor(np.log10(float(np.nanmin(merged[mask])))) - 0.2)
    upper = 10.0 ** (np.ceil(np.log10(ymax)) + 0.2)
    return ymin, upper


def plot_ratio_to_slc1(ax, medium_results: dict[str, dict[str, np.ndarray]]) -> None:
    ref = medium_results["slc1_1d"]
    ref_gamma = ref["gamma"]
    ref_tail = ref_gamma * ref_gamma * ref["dnde"][:, -1]
    ref_mask = (ref_gamma >= 1.0e3) & np.isfinite(ref_tail) & (ref_tail > 0.0)
    ref_x = np.log10(ref_gamma[ref_mask])
    ref_y = np.log10(ref_tail[ref_mask])
    for solver in SOLVERS:
        bundle = medium_results[solver]
        gamma = bundle["gamma"]
        tail = gamma * gamma * bundle["dnde"][:, -1]
        mask = (gamma >= 1.0e3) & np.isfinite(tail) & (tail > 0.0)
        mask &= (gamma >= ref_gamma[ref_mask][0]) & (gamma <= ref_gamma[ref_mask][-1])
        if not np.any(mask):
            continue
        ratio = tail[mask] / (10.0 ** np.interp(np.log10(gamma[mask]), ref_x, ref_y))
        ax.loglog(
            gamma[mask],
            ratio,
            color="black",
            linestyle=SOLVER_STYLES[solver],
            lw=1.8,
            marker=SOLVER_MARKERS[solver],
            markevery=14,
            ms=3.2,
            mfc="none",
        )


def write_plot(results: dict[str, dict[str, dict[str, np.ndarray]]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 2, figsize=(15.2, 12.2), sharex="row")
    fig.suptitle("Forward-shock SED and high-energy electron-tail comparison", fontsize=16)
    for col, (medium_name, _medium) in enumerate(MEDIA):
        ax_sed = axes[0, col]
        ax_tail = axes[1, col]
        ax_ratio = axes[2, col]
        medium_results = results[medium_name]
        ax_sed.set_title(medium_name, fontsize=14)
        ax_sed.text(0.04, 0.92, f"({chr(97 + col)})", transform=ax_sed.transAxes, fontsize=13)
        ax_tail.text(0.04, 0.90, f"({chr(99 + col)})", transform=ax_tail.transAxes, fontsize=13)
        ax_ratio.text(0.04, 0.90, f"({chr(101 + col)})", transform=ax_ratio.transAxes, fontsize=13)
        ymax = 0.0
        positive_values = []
        for solver in SOLVERS:
            bundle = medium_results[solver]
            for time_index, time_s in enumerate(TIMES_S):
                style = SOLVER_STYLES[solver]
                marker = SOLVER_MARKERS[solver]
                color = TIME_COLORS[float(time_s)]
                ax_sed.loglog(FREQUENCIES_HZ, bundle["sed"][:, time_index], color=color, linestyle=style, lw=1.7)
                gamma = bundle["gamma"]
                tail = gamma * gamma * bundle["dnde"][:, time_index]
                ax_tail.loglog(
                    gamma,
                    tail,
                    color=color,
                    linestyle=style,
                    lw=1.55,
                    marker=marker,
                    markevery=18,
                    ms=3.2,
                    mec=color,
                    mfc="none",
                )
                mask = (gamma >= 1.0e2) & np.isfinite(tail) & (tail > 0.0)
                if np.any(mask):
                    ymax = max(ymax, float(np.nanmax(tail[mask])))
                    positive_values.append(tail[mask])
        plot_ratio_to_slc1(ax_ratio, medium_results)
        ax_sed.set_xlim(FREQUENCIES_HZ[0], FREQUENCIES_HZ[-1])
        ax_tail.set_xlim(1.0e2, 4.0e10)
        ax_ratio.set_xlim(1.0e3, 4.0e10)
        limits = tail_limits(positive_values, ymax)
        if limits is not None:
            ax_tail.set_ylim(*limits)
        ax_ratio.axhline(1.0, color="0.55", lw=1.0)
        ax_ratio.set_ylim(1.0e-2, 1.0e2)
        for ax in (ax_sed, ax_tail, ax_ratio):
            ax.grid(True, which="both", alpha=0.18)
            ax.tick_params(which="both", direction="in", top=True, right=True)
        ax_sed.set_xlabel("Observed frequency nu (Hz)")
        ax_tail.set_xlabel("Electron Lorentz factor gamma_e")
        ax_ratio.set_xlabel("Electron Lorentz factor gamma_e")

    axes[0, 0].set_ylabel("nu F_nu (erg cm^-2 s^-1)")
    axes[1, 0].set_ylabel("gamma_e^2 dN_e/dgamma_e")
    axes[2, 0].set_ylabel("tail ratio to slc1_1d")
    from matplotlib.lines import Line2D

    time_handles = [
        Line2D([0], [0], color=TIME_COLORS[float(time_s)], lw=2.0, label=f"t_obs={time_s:.0e} s")
        for time_s in TIMES_S
    ]
    solver_handles = [
        Line2D(
            [0],
            [0],
            color="black",
            lw=2.0,
            linestyle=SOLVER_STYLES[solver],
            marker=SOLVER_MARKERS[solver],
            ms=4.0,
            mfc="none",
            label=solver,
        )
        for solver in SOLVERS
    ]
    axes[0, 0].legend(handles=time_handles, title="Time", frameon=False, fontsize=9, title_fontsize=9, loc="lower left")
    axes[0, 1].legend(handles=solver_handles, title="Solver", frameon=False, fontsize=9, title_fontsize=9, loc="lower left")
    axes[2, 0].text(0.04, 0.12, "ratio panel uses t_obs=1e7 s", transform=axes[2, 0].transAxes, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output, dpi=180)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("output/benchmark_exp_tail"))
    args = parser.parse_args()
    results = compute()
    write_diagnostics(results, args.output_dir / "tail_diagnostics.csv")
    write_plot(results, args.output_dir / "spectrum_compare.png")
    print(args.output_dir / "spectrum_compare.png")


if __name__ == "__main__":
    main()
