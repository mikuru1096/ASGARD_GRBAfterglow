from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from prompt.eats import EATSNumerics
from prompt.internal_shock import InternalShockNumerics, InternalShockShell, simulate_internal_shock
from prompt.radiation import InternalShockMicrophysics, RadiationNumerics, compute_prompt_observed_flux


OUTPUT_DIR = Path(__file__).resolve().parent / "results"
OKABE_ITO = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")


def main() -> None:
    slow = InternalShockShell(gamma=200.0, total_energy_iso_erg=4.0e52, duration_s=0.2, sigma=0.01)
    fast = InternalShockShell(gamma=600.0, total_energy_iso_erg=4.0e52, duration_s=0.2, sigma=0.03)
    microphysics = InternalShockMicrophysics(epsilon_e=0.1, epsilon_b=0.01, electron_index_p=2.3)
    solution = simulate_internal_shock(
        slow,
        fast,
        engine_gap_s=0.2,
        redshift=0.5,
        luminosity_distance_cm=1.0e28,
        opening_angle_rad=0.1,
        epsilon_e=microphysics.epsilon_e,
        epsilon_b=microphysics.epsilon_b,
        numerics=InternalShockNumerics(num_branch_steps=160),
    )
    observer_time = np.linspace(1.0e-4, 2.2, 384)
    observer_frequency = np.logspace(16.0, 25.0, 96)
    flux = compute_prompt_observed_flux(
        solution,
        observer_frequency,
        observer_time,
        microphysics=microphysics,
        radiation_numerics=RadiationNumerics(
            num_electron_gamma=201,
            num_photon_frequency=241,
            num_threads=8,
        ),
        eats_numerics=EATSNumerics(
            num_theta=48,
            num_phi=1,
            num_threads=8,
            adaptive_rtol=3.0e-3,
            adaptive_max_depth=8,
        ),
    )
    OUTPUT_DIR.mkdir(exist_ok=True)
    _write_flux_npz(observer_time, observer_frequency, flux)
    _write_diagnostics(solution, flux)
    _write_summary_table(observer_time, observer_frequency, flux.total)
    _plot_lightcurves(observer_time, observer_frequency, flux.total)
    _plot_spectra(observer_time, observer_frequency, flux.total)
    _plot_components(observer_time, observer_frequency, flux)


def _write_flux_npz(observer_time: np.ndarray, observer_frequency: np.ndarray, flux) -> None:
    np.savez(
        OUTPUT_DIR / "formal_flux.npz",
        observer_time_s=observer_time,
        observer_frequency_hz=observer_frequency,
        fs_sync=flux.fs_sync,
        fs_ssc=flux.fs_ssc,
        rs_sync=flux.rs_sync,
        rs_ssc=flux.rs_ssc,
        total=flux.total,
    )


def _write_diagnostics(solution, flux) -> None:
    diagnostics = {
        "radius_collision_cm": float(solution.radius_collision_cm),
        "gamma_contact": float(solution.gamma_contact),
        "fs_valid_shock": bool(solution.fs.valid_shock),
        "rs_valid_shock": bool(solution.rs.valid_shock),
        "fs_crossing_time_lab_s": float(solution.fs.jump.crossing_time_lab_s),
        "rs_crossing_time_lab_s": float(solution.rs.jump.crossing_time_lab_s),
        "sigma_slow": float(solution.slow_shell.sigma),
        "sigma_fast": float(solution.fast_shell.sigma),
        "slow_baryonic_mass_g": float(solution.slow_baryonic_mass_g),
        "fast_baryonic_mass_g": float(solution.fast_baryonic_mass_g),
        "fs_compression": float(solution.fs.jump.compression),
        "rs_compression": float(solution.rs.jump.compression),
        "fs_specific_internal": float(solution.fs.jump.specific_internal_energy),
        "rs_specific_internal": float(solution.rs.jump.specific_internal_energy),
        "fs_ordered_b_g_mean": float(np.mean(solution.fs.ordered_b_g)),
        "rs_ordered_b_g_mean": float(np.mean(solution.rs.ordered_b_g)),
        "fs_turbulent_b_g_mean": float(np.mean(solution.fs.turbulent_b_g)),
        "rs_turbulent_b_g_mean": float(np.mean(solution.rs.turbulent_b_g)),
        "fs_total_b_g_mean": float(np.mean(solution.fs.total_b_g)),
        "rs_total_b_g_mean": float(np.mean(solution.rs.total_b_g)),
        "fs_gamma_m_median": float(np.median(flux.fs_radiation.gamma_m[flux.fs_radiation.gamma_m > 0.0])),
        "rs_gamma_m_median": float(np.median(flux.rs_radiation.gamma_m[flux.rs_radiation.gamma_m > 0.0])),
        "fs_gamma_c_median": float(np.median(flux.fs_radiation.gamma_c[flux.fs_radiation.gamma_c > 0.0])),
        "rs_gamma_c_median": float(np.median(flux.rs_radiation.gamma_c[flux.rs_radiation.gamma_c > 0.0])),
        "fs_gamma_max_median": float(np.median(flux.fs_radiation.gamma_max[flux.fs_radiation.gamma_max > 0.0])),
        "rs_gamma_max_median": float(np.median(flux.rs_radiation.gamma_max[flux.rs_radiation.gamma_max > 0.0])),
        "fs_efficiency_median": float(np.median(flux.fs_radiation.efficiency[flux.fs_radiation.efficiency > 0.0])),
        "rs_efficiency_median": float(np.median(flux.rs_radiation.efficiency[flux.rs_radiation.efficiency > 0.0])),
        "fs_compactness_median": float(np.median(flux.fs_radiation.compactness[flux.fs_radiation.compactness > 0.0])),
        "rs_compactness_median": float(np.median(flux.rs_radiation.compactness[flux.rs_radiation.compactness > 0.0])),
        "fs_gamma_gamma_min_absorption": float(np.min(flux.fs_radiation.gamma_gamma_absorption)),
        "rs_gamma_gamma_min_absorption": float(np.min(flux.rs_radiation.gamma_gamma_absorption)),
    }
    (OUTPUT_DIR / "formal_diagnostics.json").write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")


def _write_summary_table(observer_time: np.ndarray, observer_frequency: np.ndarray, total: np.ndarray) -> None:
    bands = (("10 keV", 2.418e18), ("100 keV", 2.418e19), ("1 MeV", 2.418e20), ("100 MeV", 2.418e22))
    bolometric = np.trapezoid(total * observer_frequency[:, None], np.log(observer_frequency), axis=0)
    rows = [
        ("bolometric_grid", float(observer_time[np.argmax(bolometric)]), float(np.max(bolometric))),
    ]
    for label, frequency in bands:
        curve = frequency * _interp_log_flux(observer_frequency, total, frequency)
        rows.append((label, float(observer_time[np.argmax(curve)]), float(np.max(curve))))
    with (OUTPUT_DIR / "formal_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("band", "peak_time_s", "peak_nuFnu_erg_cm2_s"))
        writer.writerows(rows)


def _plot_lightcurves(observer_time: np.ndarray, observer_frequency: np.ndarray, total: np.ndarray) -> None:
    _set_plot_style()
    bands = (("10 keV", 2.418e18), ("100 keV", 2.418e19), ("1 MeV", 2.418e20), ("100 MeV", 2.418e22), ("1 GeV", 2.418e23))
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    for idx, (label, frequency) in enumerate(bands):
        curve = frequency * _interp_log_flux(observer_frequency, total, frequency)
        ax.plot(observer_time, curve, lw=2.0, color=OKABE_ITO[idx], label=label)
    ax.set_yscale("log")
    ax.set_xlabel("Observer time (s)")
    ax.set_ylabel(r"$\nu F_\nu$ (erg cm$^{-2}$ s$^{-1}$)")
    ax.set_title("Formal prompt internal-shock light curves")
    ax.legend(frameon=False, ncol=2)
    ax.grid(alpha=0.25, lw=0.6)
    _set_positive_ylim(ax, [line.get_ydata() for line in ax.lines])
    fig.savefig(OUTPUT_DIR / "formal_lightcurves.png", dpi=300)
    plt.close(fig)


def _plot_spectra(observer_time: np.ndarray, observer_frequency: np.ndarray, total: np.ndarray) -> None:
    _set_plot_style()
    bolometric = np.trapezoid(total * observer_frequency[:, None], np.log(observer_frequency), axis=0)
    active = np.flatnonzero(bolometric > np.max(bolometric) * 1.0e-4)
    indices = np.unique(
        np.array(
            [
                active[0],
                active[active.size // 4],
                int(np.argmax(bolometric)),
                active[(3 * active.size) // 4],
                active[-1],
            ],
            dtype=int,
        )
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    for idx, time_index in enumerate(indices):
        ax.plot(
            observer_frequency,
            observer_frequency * total[:, time_index],
            lw=2.0,
            color=OKABE_ITO[idx],
            label=f"t={observer_time[time_index]:.3g} s",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Observer frequency (Hz)")
    ax.set_ylabel(r"$\nu F_\nu$ (erg cm$^{-2}$ s$^{-1}$)")
    ax.set_title("Formal prompt internal-shock spectra")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25, lw=0.6, which="both")
    _set_positive_ylim(ax, [line.get_ydata() for line in ax.lines])
    fig.savefig(OUTPUT_DIR / "formal_spectra.png", dpi=300)
    plt.close(fig)


def _plot_components(observer_time: np.ndarray, observer_frequency: np.ndarray, flux) -> None:
    _set_plot_style()
    panels = (("100 keV", 2.418e19), ("1 MeV", 2.418e20))
    components = (
        ("FS sync", flux.fs_sync, "-", OKABE_ITO[0]),
        ("FS SSC", flux.fs_ssc, "--", OKABE_ITO[0]),
        ("RS sync", flux.rs_sync, "-", OKABE_ITO[1]),
        ("RS SSC", flux.rs_ssc, "--", OKABE_ITO[1]),
    )
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True, constrained_layout=True)
    for ax, (panel_label, frequency) in zip(axes, panels):
        plotted = []
        for label, matrix, linestyle, color in components:
            curve = frequency * _interp_log_flux(observer_frequency, matrix, frequency)
            plotted.append(curve)
            ax.plot(observer_time, curve, lw=1.8, ls=linestyle, color=color, label=label)
        ax.set_yscale("log")
        ax.set_ylabel(rf"{panel_label} $\nu F_\nu$")
        ax.grid(alpha=0.25, lw=0.6)
        ax.legend(frameon=False, ncol=2)
        _set_positive_ylim(ax, plotted)
    axes[-1].set_xlabel("Observer time (s)")
    fig.savefig(OUTPUT_DIR / "formal_components.png", dpi=300)
    plt.close(fig)


def _interp_log_flux(frequency_grid: np.ndarray, matrix: np.ndarray, target_frequency: float) -> np.ndarray:
    log_frequency = np.log(frequency_grid)
    values = np.zeros(matrix.shape[1], dtype=float)
    for time_index in range(matrix.shape[1]):
        column = matrix[:, time_index]
        positive = column > 0.0
        if np.count_nonzero(positive) >= 2:
            values[time_index] = np.exp(np.interp(np.log(target_frequency), log_frequency[positive], np.log(column[positive])))
    return values


def _set_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.linewidth": 0.9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.handlelength": 2.8,
        }
    )


def _set_positive_ylim(ax, curves: list[np.ndarray]) -> None:
    positive = np.concatenate([np.asarray(curve, dtype=float)[np.asarray(curve, dtype=float) > 0.0] for curve in curves])
    ymax = float(np.max(positive))
    informative = positive[positive > ymax * 1.0e-6]
    ymin = float(np.min(informative))
    ax.set_ylim(ymin / 2.0, ymax * 2.0)


if __name__ == "__main__":
    main()
