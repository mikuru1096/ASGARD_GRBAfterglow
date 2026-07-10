from __future__ import annotations

from pathlib import Path
import argparse
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core import Model, Observer, UniformMedium, top_hat_jet
from asgard_core.asgard_paths import DOC_ROOT
from tests.public_api_builders import hadronic, numerics, observer_grid, radiation, reverse_shock, solver_options


OUTPUT_DIR = DOC_ROOT / "reverse_density_jump_tests"
OUTPUT_PATH = OUTPUT_DIR / "triple_density_jump_rs_fs_tophat.png"
JUMP_RADII_CM = (1.0e15, 1.0e16, 1.0e17)
JUMP_FACTOR = 1.0e2
JUMP_WIDTH_REL = 1.0e-1
REVERSE_SIGMA = 1.0e-1
BANDS_HZ = np.array([1.0e9, 1.0e14, 1.0e18], dtype=float)
BAND_LABELS = ("1 GHz radio", "1e14 Hz optical", "1e18 Hz X-ray")
MODE_GRIDS = {
    "quick": {"num_r": 60, "num_theta": 32, "num_tobs": 32, "num_gam_e": 41, "num_nu": 31, "times": 80},
    "formal": {"num_r": 120, "num_theta": 80, "num_tobs": 48, "num_gam_e": 81, "num_nu": 49, "times": 160},
}
PALETTE = {
    "total": "#0F4D92",
    "forward": "#42949E",
    "reverse": "#B64342",
    "baseline": "#767676",
    "density": "#272727",
    "jump": "#B64342",
}

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": True,
        "axes.spines.top": True,
        "axes.linewidth": 0.7,
        "legend.frameon": False,
    }
)


def _model(*, with_jumps: bool, grid: dict[str, int]) -> Model:
    model = Model(
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
        fwd_rad=radiation(
            epsilon_e=3.0e-1,
            epsilon_B=1.0e-5,
            p=2.3,
            accelerated_electron_fraction=1.0e-1,
            include_ssc=False,
        ),
        rvs_rad=radiation(
            epsilon_e=3.0e-1,
            epsilon_B=1.0e-5,
            p=2.4,
            accelerated_electron_fraction=1.0e-1,
            include_ssc=False,
        ),
        numerics=numerics(
            num_radius=grid["num_r"],
            eats_num_theta=grid["num_theta"],
            eats_num_phi=1,
            num_observer_time=grid["num_tobs"],
            num_electron_gamma=grid["num_gam_e"],
            num_photon_frequency=grid["num_nu"],
            num_threads=1,
            initial_radius_cm=1.0e13,
        ),
        observer_grid=observer_grid(time_min_s=1.0e0, time_max_s=1.0e8),
        solver_options=solver_options(
            electron_solver="fullhide_1d",
            ssc_cooling_mode="none",
        ),
        reverse_shock=reverse_shock(
            enabled=True,
            shell_duration_s=10.0,
            upstream_sigma=REVERSE_SIGMA,
            include_ssc=False,
            include_cross_zone_ic=False,
        ),
        hadronic=hadronic(),
    )
    if with_jumps:
        model.setups.jump_r_cm = tuple(float(value) for value in JUMP_RADII_CM)
        model.setups.jump_factor = (float(JUMP_FACTOR),) * len(JUMP_RADII_CM)
        model.setups.jump_width = (float(JUMP_WIDTH_REL),) * len(JUMP_RADII_CM)
    return model


def _density_enhancement(radius_cm: np.ndarray) -> np.ndarray:
    enhancement = np.ones_like(radius_cm, dtype=float)
    for radius_j in JUMP_RADII_CM:
        width_cm = JUMP_WIDTH_REL * radius_j
        enhancement = enhancement + (JUMP_FACTOR - 1.0) * np.exp(
            -((radius_cm - radius_j) ** 2) / (2.0 * width_cm**2)
        )
    return enhancement


def _positive_for_log(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float).copy()
    arr[(arr <= 0.0) | ~np.isfinite(arr)] = np.nan
    return arr


def _ylim(curves: list[np.ndarray]) -> tuple[float, float]:
    positive = [curve[np.isfinite(curve) & (curve > 0.0)] for curve in curves]
    positive = [curve for curve in positive if curve.size > 0]
    if not positive:
        raise RuntimeError("light-curve panel has no positive finite flux values")
    values = np.concatenate(positive)
    ymax = float(np.nanmax(values))
    shown = values[values >= ymax * 1.0e-8]
    ymin = float(np.nanmin(shown if shown.size else values))
    return 10.0 ** (np.log10(ymin) - 0.25), 10.0 ** (np.log10(ymax) + 0.25)


def _plot_light_curve_panel(
    ax,
    times_s: np.ndarray,
    *,
    band_index: int,
    no_jump_total: np.ndarray,
    jump_total: np.ndarray,
    jump_fwd: np.ndarray,
    jump_rev: np.ndarray,
) -> None:
    curves = [no_jump_total, jump_total, jump_fwd, jump_rev]
    ax.loglog(times_s, _positive_for_log(no_jump_total), color=PALETTE["baseline"], lw=1.0, ls="--", label="no-jump total")
    ax.loglog(times_s, _positive_for_log(jump_total), color=PALETTE["total"], lw=1.55, label="triple-jump total")
    ax.loglog(times_s, _positive_for_log(jump_fwd), color=PALETTE["forward"], lw=1.1, ls="-.", label="triple-jump FS")
    ax.loglog(times_s, _positive_for_log(jump_rev), color=PALETTE["reverse"], lw=1.1, ls=":", label="triple-jump RS")
    ax.set_ylim(*_ylim(curves))
    ax.set_title(BAND_LABELS[band_index])
    ax.grid(True, which="major", axis="y", alpha=0.16, linestyle="-", linewidth=0.5)
    ax.tick_params(axis="both", which="major", length=3.0, width=0.7)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.7)


def _add_panel_label(ax, label: str, *, x: float = -0.055, y: float = 1.05) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=8, fontweight="bold", ha="left", va="bottom")


def _save_figure(fig: plt.Figure, output: Path) -> list[Path]:
    output.parent.mkdir(parents=True, exist_ok=True)
    stem = output.with_suffix("")
    paths = [
        output,
        stem.with_suffix(".svg"),
        stem.with_suffix(".pdf"),
        stem.with_suffix(".tiff"),
    ]
    fig.savefig(paths[0], dpi=600, bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight")
    _strip_svg_trailing_whitespace(paths[1])
    fig.savefig(paths[2], bbox_inches="tight")
    fig.savefig(paths[3], dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    return paths


def _strip_svg_trailing_whitespace(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")


def _write_secondary_event_csv(details, output: Path) -> Path:
    if details.rev is None or details.rev.secondary_rs_event_active is None:
        raise RuntimeError("triple density-jump benchmark requires secondary RS event diagnostics")
    csv_path = output.with_suffix("").with_name(output.with_suffix("").name + "_secondary_rs_events.csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "jump_index",
            "jump_radius_cm",
            "jump_width_rel",
            "event_active",
            "start_radius_cm",
            "shock_end_radius_cm",
            "start_tobs_axis_s",
            "shock_end_tobs_axis_s",
            "cooling_tail_shells",
            "cooling_tail_max_B_g",
            "cooling_tail_max_luminosity",
            "cooling_tail_max_log10_B_jump",
        ])
        for i_jump, radius_j in enumerate(JUMP_RADII_CM):
            tail = _secondary_cooling_tail_metrics(details.rev, i_jump)
            writer.writerow([
                i_jump,
                radius_j,
                JUMP_WIDTH_REL,
                bool(details.rev.secondary_rs_event_active[i_jump]),
                float(details.rev.secondary_rs_start_radius[i_jump]),
                float(details.rev.secondary_rs_shock_end_radius[i_jump]),
                float(details.rev.secondary_rs_start_tobs_axis[i_jump]),
                float(details.rev.secondary_rs_shock_end_tobs_axis[i_jump]),
                tail["shells"],
                tail["max_B_g"],
                tail["max_luminosity"],
                tail["max_log10_B_jump"],
            ])
    return csv_path


def _write_secondary_energy_csv(details, model: Model, output: Path) -> Path:
    if (
        details.rev is None
        or details.rev.secondary_rs_dissipated_energy_erg is None
        or details.rev.secondary_rs_electron_injected_energy_erg is None
    ):
        raise RuntimeError("triple density-jump benchmark requires secondary RS energy diagnostics")
    dissipated = np.asarray(details.rev.secondary_rs_dissipated_energy_erg, dtype=float)
    electron = np.asarray(details.rev.secondary_rs_electron_injected_energy_erg, dtype=float)
    total_dissipated = float(np.sum(dissipated))
    total_electron = float(np.sum(electron))
    if total_dissipated <= 0.0:
        raise RuntimeError("secondary RS dissipated energy must be positive in the triple-jump benchmark")
    expected_electron = float(model.rvs_rad.eps_e) * total_dissipated
    csv_path = output.with_suffix("").with_name(output.with_suffix("").name + "_secondary_rs_energy.csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "secondary_rs_dissipated_energy_erg",
            "secondary_rs_electron_injected_energy_erg",
            "expected_electron_energy_erg",
            "electron_energy_fractional_error",
        ])
        writer.writerow([
            total_dissipated,
            total_electron,
            expected_electron,
            (total_electron - expected_electron) / expected_electron,
        ])
    return csv_path


def _write_adaptive_convergence_csv(model: Model, output: Path, times: np.ndarray, direct_total: np.ndarray) -> Path:
    adaptive = model.flux_density_grid_adaptive(times, BANDS_HZ)
    adaptive_total = np.asarray(adaptive.total, dtype=float)
    adaptive_times = np.asarray(adaptive.time_s, dtype=float)
    csv_path = output.with_suffix("").with_name(output.with_suffix("").name + "_adaptive_convergence.csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "band_hz",
            "user_time_count",
            "adaptive_time_count",
            "direct_peak_time_s",
            "adaptive_peak_time_s",
            "direct_peak_flux",
            "adaptive_peak_flux",
            "peak_time_ratio",
            "peak_flux_fractional_difference",
            "direct_integral_flux_dt",
            "adaptive_integral_flux_dt",
            "integral_fractional_difference",
        ])
        for i_band, nu_hz in enumerate(BANDS_HZ):
            direct_peak = int(np.nanargmax(direct_total[i_band]))
            adaptive_peak = int(np.nanargmax(adaptive_total[i_band]))
            direct_integral = float(np.trapezoid(direct_total[i_band], times))
            adaptive_integral = float(np.trapezoid(adaptive_total[i_band], adaptive_times))
            writer.writerow([
                float(nu_hz),
                int(times.size),
                int(adaptive_times.size),
                float(times[direct_peak]),
                float(adaptive_times[adaptive_peak]),
                float(direct_total[i_band, direct_peak]),
                float(adaptive_total[i_band, adaptive_peak]),
                float(adaptive_times[adaptive_peak] / times[direct_peak]),
                float(
                    (adaptive_total[i_band, adaptive_peak] - direct_total[i_band, direct_peak])
                    / adaptive_total[i_band, adaptive_peak]
                ),
                direct_integral,
                adaptive_integral,
                float((adaptive_integral - direct_integral) / adaptive_integral),
            ])
    return csv_path


def _secondary_cooling_tail_metrics(track, i_jump: int) -> dict[str, float | int]:
    radius = np.asarray(track.radius, dtype=float)
    branch_b = np.asarray(track.secondary_rs_branch_B, dtype=float)
    branch_lum = np.asarray(track.secondary_rs_branch_luminosity_syn, dtype=float)
    end_radius = float(track.secondary_rs_shock_end_radius[i_jump])
    mask = (radius > end_radius) & (branch_b[i_jump] > 0.0)
    if not np.any(mask):
        return {"shells": 0, "max_B_g": 0.0, "max_luminosity": 0.0, "max_log10_B_jump": 0.0}
    tail_b = branch_b[i_jump, mask]
    positive_b = tail_b[np.isfinite(tail_b) & (tail_b > 0.0)]
    if positive_b.size > 1:
        max_jump = float(np.max(np.abs(np.diff(np.log10(positive_b)))))
    else:
        max_jump = 0.0
    return {
        "shells": int(np.count_nonzero(mask)),
        "max_B_g": float(np.max(positive_b)),
        "max_luminosity": float(np.max(branch_lum[i_jump, :, mask])),
        "max_log10_B_jump": max_jump,
    }


def build_plot(*, mode: str, output: Path, times_count: int | None = None) -> Path:
    grid = dict(MODE_GRIDS[mode])
    if times_count is not None:
        grid["times"] = int(times_count)
    times = np.logspace(0.0, 8.0, grid["times"])
    no_jump_model = _model(with_jumps=False, grid=grid)
    triple_jump_model = _model(with_jumps=True, grid=grid)
    no_jump = no_jump_model.flux_density_grid(times, BANDS_HZ)
    triple_jump = triple_jump_model.flux_density_grid(times, BANDS_HZ)

    no_jump_total = np.asarray(no_jump.total, dtype=float)
    jump_total = np.asarray(triple_jump.total, dtype=float)
    jump_fwd = np.asarray(triple_jump.fwd.sync, dtype=float)
    jump_rev = np.asarray(triple_jump.rev.sync, dtype=float)

    fig = plt.figure(figsize=(7.2, 7.8), dpi=300)
    grid_spec = fig.add_gridspec(4, 1, height_ratios=(0.9, 1.25, 1.25, 1.25), hspace=0.50)
    density_ax = fig.add_subplot(grid_spec[0, 0])
    radius = np.logspace(14.5, 17.5, 800)
    enhancement = _density_enhancement(radius)
    density_ax.fill_between(radius, 1.0, enhancement, color="#D8E6F3", alpha=0.85, linewidth=0.0)
    density_ax.semilogx(radius, enhancement, color=PALETTE["density"], lw=1.1)
    for radius_j in JUMP_RADII_CM:
        density_ax.axvline(radius_j, color=PALETTE["jump"], lw=0.7, ls=":")
    density_ax.set_yscale("log")
    density_ax.set_ylim(0.8, 160.0)
    density_ax.set_yticks([1.0, 10.0, 100.0])
    density_ax.set_ylabel(r"$n/n_0$")
    density_ax.set_xlabel("Radius [cm]", labelpad=4.0)
    density_ax.set_title("Ambient density structure")
    density_ax.grid(True, which="major", axis="y", alpha=0.16, linestyle="-", linewidth=0.5)
    density_ax.tick_params(axis="both", which="major", length=3.0, width=0.7)
    for spine in density_ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.7)
    _add_panel_label(density_ax, "a", x=-0.06, y=1.08)

    density_ax.text(
        0.985,
        0.92,
        r"$\Gamma_0=100,\ \sigma=0.1$" "\n"
        r"$f_{\rm jump}=100,\ \sigma_R/R_{\rm jump}=0.1$",
        transform=density_ax.transAxes,
        fontsize=6.6,
        ha="right",
        va="top",
        linespacing=1.35,
    )

    axes = [fig.add_subplot(grid_spec[i + 1, 0]) for i in range(BANDS_HZ.size)]
    for i_band, ax in enumerate(axes):
        _plot_light_curve_panel(
            ax,
            times,
            band_index=i_band,
            no_jump_total=no_jump_total[i_band],
            jump_total=jump_total[i_band],
            jump_fwd=jump_fwd[i_band],
            jump_rev=jump_rev[i_band],
        )
        _add_panel_label(ax, chr(ord("b") + i_band))
        if i_band < BANDS_HZ.size - 1:
            ax.tick_params(axis="x", labelbottom=False)
        else:
            ax.set_xlabel("Observer time [s]")
    axes[0].set_ylabel(r"Flux density [erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$]")
    handles, labels = axes[-1].get_legend_handles_labels()
    legend = fig.legend(handles, labels, fontsize=6.4, loc="upper center", ncol=4, frameon=False, handlelength=2.5, bbox_to_anchor=(0.5, 0.998))
    for line in legend.get_lines():
        line.set_linewidth(1.6)
    fig.subplots_adjust(top=0.93, bottom=0.075, left=0.12, right=0.985)
    paths = _save_figure(fig, output)
    plt.close(fig)
    details = triple_jump_model.details(float(times[0]), float(times[-1]))
    event_csv = _write_secondary_event_csv(details, output)
    energy_csv = _write_secondary_energy_csv(details, triple_jump_model, output)
    adaptive_csv = _write_adaptive_convergence_csv(triple_jump_model, output, times, jump_total)

    for i_band, nu_hz in enumerate(BANDS_HZ):
        peak_index = int(np.nanargmax(jump_total[i_band]))
        print(
            f"band={nu_hz:.3e}Hz peak_time={times[peak_index]:.3e}s "
            f"peak_total={jump_total[i_band, peak_index]:.3e}"
        )
    for path in paths:
        print(f"wrote {path}")
    print(f"wrote {event_csv}")
    print(f"wrote {energy_csv}")
    print(f"wrote {adaptive_csv}")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate the ASGARD triple-density-jump RS+FS top-hat benchmark.")
    parser.add_argument("--mode", choices=sorted(MODE_GRIDS), default="formal")
    parser.add_argument("--times", type=int, default=None)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    build_plot(mode=args.mode, output=args.output, times_count=args.times)


if __name__ == "__main__":
    main()
