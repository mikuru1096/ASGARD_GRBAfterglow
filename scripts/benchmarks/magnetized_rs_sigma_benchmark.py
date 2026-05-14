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

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet
from asgard_core.asgard_config import FitConfig, ReverseShockConfig
from asgard_core.asgard_state import make_query_setup, solve_state_from_setup


OUTPUT_DIR = ROOT / "output" / "asgard_doc" / "magnetized_rs_sigma_benchmark"
SIGMAS = (0.0, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0)
COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#000000")
LINESTYLES = ("-", "--", "-.", ":", (0, (5, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)))
MARKERS = ("o", "s", "^", "D", "P", "X")


def _fit_config(sigma: float) -> FitConfig:
    return FitConfig(
        index_y=0,
        index_syn_integr=2,
        num_threads=1,
        num_gam_e=81,
        num_nu=61,
        num_r=120,
        num_theta=32,
        num_tobs=40,
        reverse=True,
        plot_lc=False,
        show_plots=False,
        reverse_shock=ReverseShockConfig(
            enabled=True,
            delta_t_s=10.0,
            epsilon_e=0.1,
            epsilon_b=1.0e-2,
            p=2.4,
            f_e=1.0,
            sigma=sigma,
        ),
    )


def _model(sigma: float) -> Model:
    return Model(
        jet=TophatJet(E_iso=1.0e53, Gamma0=1.0e2, theta_j=1.0e-1, duration=10.0),
        medium=ISM(n_ism=1.0e-1),
        observer=Observer(z=0.4, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.5, xi_N=1.0e-1, ssc=False),
        rvs_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-2, p=2.4, xi_N=1.0, ssc=False),
        setups=Setups(
            rvs_shock=True,
            reverse_delta_t_s=10.0,
            reverse_sigma=sigma,
            fwd_ssc=False,
            rvs_ssc=False,
            num_threads=1,
            num_gam_e=81,
            num_nu=81,
            num_r=140,
            num_theta=32,
            num_tobs=48,
            observer_time_min_s=1.0e0,
            observer_time_max_s=1.0e7,
            electron_solver="fullhide_1d",
        ),
    )


def _style(index: int) -> dict[str, object]:
    return {
        "color": COLORS[index % len(COLORS)],
        "linestyle": LINESTYLES[index % len(LINESTYLES)],
        "marker": MARKERS[index % len(MARKERS)],
        "markevery": 14,
        "markersize": 4.0,
        "linewidth": 1.8,
    }


def _positive(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr) & (arr > 0.0)]


def _display_ylim(curves: list[np.ndarray], *, rel_floor: float = 1.0e-10, pad_dex: float = 0.35) -> tuple[float, float]:
    positives = [_positive(curve) for curve in curves]
    positives = [curve for curve in positives if curve.size]
    if not positives:
        return 1.0e-40, 1.0
    panel = np.concatenate(positives)
    ymax = float(np.nanmax(panel))
    floor = ymax * rel_floor
    shown = panel[panel >= floor]
    if shown.size == 0:
        shown = panel
    ymin = float(np.nanmin(shown))
    return 10.0 ** (np.log10(ymin) - pad_dex), 10.0 ** (np.log10(ymax) + pad_dex)


def _mask_for_display(values: np.ndarray, *, rel_floor: float = 1.0e-10) -> np.ndarray:
    arr = np.asarray(values, dtype=float).copy()
    positive = _positive(arr)
    if positive.size == 0:
        return arr
    floor = float(np.nanmax(positive)) * rel_floor
    arr[(arr > 0.0) & (arr < floor)] = np.nan
    return arr


def _logjump(values: np.ndarray) -> float:
    positive = _positive(values)
    if positive.size < 2:
        return float("nan")
    return float(np.nanmax(np.abs(np.diff(np.log10(positive)))))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _solve_dynamics_scan(sigmas: tuple[float, ...]) -> tuple[list[dict[str, object]], list[tuple[float, np.ndarray, object, np.ndarray, np.ndarray]]]:
    rows: list[dict[str, object]] = []
    series = []
    for sigma in sigmas:
        cfg = _fit_config(sigma)
        setup = make_query_setup(cfg, np.logspace(2.0, 6.0, 10), np.array([1.0e9, 1.0e14]))
        state = solve_state_from_setup(cfg, setup)
        rs = state.dynamics.reverse_shock
        if rs is None:
            raise RuntimeError("reverse shock state missing")
        t_obs = np.asarray(state.dynamics.r_tobs, dtype=float)
        e3 = np.asarray(rs.internal_energy_erg, dtype=float) / np.asarray(rs.comoving_volume_cm3, dtype=float)
        u3_per_m3 = np.asarray(rs.internal_energy_erg, dtype=float) / np.asarray(rs.swept_mass_g, dtype=float)
        post = np.asarray(state.dynamics.radius, dtype=float) >= float(rs.r_cross)
        post_u = u3_per_m3[post]
        post_mono = bool(np.all(np.diff(post_u) <= 0.0)) if post_u.size > 2 else True
        rows.append(
            {
                "sigma": sigma,
                "t_cross_s": float(rs.t_cross),
                "r_cross_cm": float(rs.r_cross),
                "ordered_B_cross_G": float(rs.ordered_magnetic_cross_g),
                "max_B3_G": float(np.nanmax(rs.magnetic_field_g)),
                "max_gamma34": float(np.nanmax(rs.gamma34)),
                "max_e3_erg_cm3": float(np.nanmax(e3)),
                "logjump_B3": _logjump(rs.magnetic_field_g),
                "logjump_gamma34": _logjump(rs.gamma34),
                "logjump_e3": _logjump(e3),
                "logjump_U3_per_M3": _logjump(u3_per_m3),
                "post_U3_per_M3_monotone_decreasing": post_mono,
            }
        )
        series.append((sigma, t_obs, rs, e3, u3_per_m3))
    return rows, series


def _solve_observable_scan(sigmas: tuple[float, ...]):
    lc_times = np.logspace(0.0, 7.0, 96)
    lc_bands = np.array([1.0e9, 1.0e14, 1.0e18], dtype=float)
    sed_times = np.array([1.0e2, 1.0e3, 1.0e4, 1.0e5], dtype=float)
    sed_freqs = np.logspace(7.0, 22.0, 180)
    lc_data = {}
    sed_data = {}
    lc_rows: list[dict[str, object]] = []
    sed_rows: list[dict[str, object]] = []
    for sigma in sigmas:
        model = _model(sigma)
        lc = model.flux_density_grid(lc_times, lc_bands)
        sed = model.flux_density_grid(sed_times, sed_freqs)
        lc_total = np.asarray(lc.total, dtype=float)
        lc_rs = np.asarray(lc.rev.sync, dtype=float)
        sed_total = np.asarray(sed.total, dtype=float)
        sed_rs = np.asarray(sed.rev.sync, dtype=float)
        lc_data[sigma] = (lc_total, lc_rs)
        sed_data[sigma] = (sed_total, sed_rs)
        for i_band, nu_hz in enumerate(lc_bands):
            total = lc_total[i_band]
            rs = lc_rs[i_band]
            peak_idx = int(np.nanargmax(total))
            rs_peak_idx = int(np.nanargmax(rs)) if np.nanmax(rs) > 0.0 else 0
            lc_rows.append(
                {
                    "sigma": sigma,
                    "nu_hz": float(nu_hz),
                    "total_peak_time_s": float(lc_times[peak_idx]),
                    "total_peak_fnu_cgs": float(total[peak_idx]),
                    "rs_peak_time_s": float(lc_times[rs_peak_idx]),
                    "rs_peak_fnu_cgs": float(rs[rs_peak_idx]),
                    "rs_to_total_at_total_peak": float(rs[peak_idx] / total[peak_idx]) if total[peak_idx] > 0.0 else float("nan"),
                }
            )
        for i_time, time_s in enumerate(sed_times):
            nu_fnu = sed_freqs * sed_total[:, i_time]
            rs_nu_fnu = sed_freqs * sed_rs[:, i_time]
            peak_idx = int(np.nanargmax(nu_fnu))
            rs_peak_idx = int(np.nanargmax(rs_nu_fnu)) if np.nanmax(rs_nu_fnu) > 0.0 else 0
            sed_rows.append(
                {
                    "sigma": sigma,
                    "time_s": float(time_s),
                    "total_peak_nu_hz": float(sed_freqs[peak_idx]),
                    "total_peak_nu_fnu_cgs": float(nu_fnu[peak_idx]),
                    "rs_peak_nu_hz": float(sed_freqs[rs_peak_idx]),
                    "rs_peak_nu_fnu_cgs": float(rs_nu_fnu[rs_peak_idx]),
                }
            )
    return lc_times, lc_bands, lc_data, lc_rows, sed_times, sed_freqs, sed_data, sed_rows


def _plot_dynamics(path: Path, series) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    x_min = 1.0
    x_max = max(float(item[1][-1]) for item in series)
    panels = [
        ("B3 total [G]", lambda item: item[2].magnetic_field_g, 1.0e-12),
        ("gamma34", lambda item: item[2].gamma34, 1.0e-12),
        ("U3/V3 [erg cm^-3]", lambda item: item[3], 1.0e-10),
        ("U3/M3 [erg g^-1]", lambda item: item[4], 1.0e-10),
    ]
    for ax, (ylabel, getter, rel_floor) in zip(axes.flat, panels):
        curves = [np.asarray(getter(item), dtype=float) for item in series]
        for i, item in enumerate(series):
            sigma, t_obs = item[0], item[1]
            y_plot = _mask_for_display(getter(item), rel_floor=rel_floor)
            y_plot = np.where(t_obs >= x_min, y_plot, np.nan)
            ax.loglog(t_obs, y_plot, label=f"sigma={sigma:g}", **_style(i))
        ax.set_ylabel(ylabel)
        ax.set_ylim(*_display_ylim(curves, rel_floor=rel_floor))
        ax.grid(True, which="both", alpha=0.25)
    axes[0, 0].legend(fontsize=8, ncol=2)
    for ax in axes[-1, :]:
        ax.set_xlabel("observer time [s]")
    axes[0, 0].set_xlim(x_min, x_max)
    fig.suptitle("ASGARD magnetized RS sigma scan: dynamics", y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_lightcurves(path: Path, lc_times: np.ndarray, lc_bands: np.ndarray, lc_data: dict[float, tuple[np.ndarray, np.ndarray]]) -> None:
    labels = ("1 GHz", "1e14 Hz", "1e18 Hz")
    fig, axes = plt.subplots(3, 2, figsize=(13, 12), sharex=True)
    for i_band, label in enumerate(labels):
        total_curves = [lc_data[sigma][0][i_band] for sigma in SIGMAS]
        rs_curves = [lc_data[sigma][1][i_band] for sigma in SIGMAS]
        for i_sigma, sigma in enumerate(SIGMAS):
            total, rs = lc_data[sigma]
            axes[i_band, 0].loglog(lc_times, _mask_for_display(total[i_band], rel_floor=1.0e-8), label=f"sigma={sigma:g}" if i_band == 0 else None, **_style(i_sigma))
            axes[i_band, 1].loglog(lc_times, _mask_for_display(rs[i_band], rel_floor=1.0e-8), **_style(i_sigma))
        axes[i_band, 0].set_ylabel(f"total Fnu {label}")
        axes[i_band, 1].set_ylabel(f"RS sync Fnu {label}")
        axes[i_band, 0].set_ylim(*_display_ylim(total_curves, rel_floor=1.0e-8))
        axes[i_band, 1].set_ylim(*_display_ylim(rs_curves, rel_floor=1.0e-8))
        for ax in axes[i_band, :]:
            ax.grid(True, which="both", alpha=0.25)
    axes[0, 0].legend(fontsize=8, ncol=2)
    for ax in axes[-1, :]:
        ax.set_xlabel("observer time [s]")
    axes[0, 0].set_xlim(lc_times[0], lc_times[-1])
    fig.suptitle("ASGARD magnetized RS sigma scan: light curves", y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_sed(path_total: Path, path_rs: Path, sed_times: np.ndarray, sed_freqs: np.ndarray, sed_data: dict[float, tuple[np.ndarray, np.ndarray]]) -> None:
    for component, path, title, rel_floor in (
        ("total", path_total, "ASGARD magnetized RS sigma scan: total SED", 1.0e-10),
        ("rs", path_rs, "ASGARD magnetized RS sigma scan: RS synch SED", 1.0e-8),
    ):
        fig, axes = plt.subplots(2, 2, figsize=(13, 10), sharex=True)
        for i_time, (ax, time_s) in enumerate(zip(axes.flat, sed_times)):
            curves = []
            for sigma in SIGMAS:
                total, rs = sed_data[sigma]
                matrix = total if component == "total" else rs
                curves.append(sed_freqs * matrix[:, i_time])
            for i_sigma, sigma in enumerate(SIGMAS):
                total, rs = sed_data[sigma]
                matrix = total if component == "total" else rs
                y = sed_freqs * matrix[:, i_time]
                ax.loglog(sed_freqs, _mask_for_display(y, rel_floor=rel_floor), label=f"sigma={sigma:g}" if i_time == 0 else None, **_style(i_sigma))
            ax.set_title(f"t = {time_s:g} s")
            ax.set_ylim(*_display_ylim(curves, rel_floor=rel_floor))
            ax.grid(True, which="both", alpha=0.25)
        for ax in axes[:, 0]:
            ax.set_ylabel("nu Fnu [cgs]")
        for ax in axes[-1, :]:
            ax.set_xlabel("frequency [Hz]")
        axes[0, 0].legend(fontsize=8, ncol=2)
        axes[0, 0].set_xlim(sed_freqs[0], sed_freqs[-1])
        fig.suptitle(title, y=0.995)
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    dyn_rows, dyn_series = _solve_dynamics_scan(SIGMAS)
    _write_csv(outdir / "sigma_scan_summary.csv", dyn_rows)
    _plot_dynamics(outdir / "magnetized_rs_sigma_dynamics_readable.png", dyn_series)

    lc_times, lc_bands, lc_data, lc_rows, sed_times, sed_freqs, sed_data, sed_rows = _solve_observable_scan(SIGMAS)
    _write_csv(outdir / "sigma_scan_lightcurve_summary.csv", lc_rows)
    _write_csv(outdir / "sigma_scan_sed_summary.csv", sed_rows)
    _plot_lightcurves(outdir / "magnetized_rs_sigma_lightcurves_readable.png", lc_times, lc_bands, lc_data)
    _plot_sed(
        outdir / "magnetized_rs_sigma_sed_total_readable.png",
        outdir / "magnetized_rs_sigma_sed_rs_readable.png",
        sed_times,
        sed_freqs,
        sed_data,
    )

    print(f"wrote {outdir}")


if __name__ == "__main__":
    main()
