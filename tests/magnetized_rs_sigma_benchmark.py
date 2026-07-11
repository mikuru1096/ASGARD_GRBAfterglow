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

from asgard_core import Model, Observer, UniformMedium, WindMedium, top_hat_jet
from asgard_core.asgard_config import ReverseShockConfig, RuntimeConfig
from asgard_core.asgard_setup import build_setup
from asgard_core.asgard_state import query_cfg, solve_setup
from tests.public_api_builders import hadronic, numerics, observer_grid, radiation, reverse_shock, solver_options
from scripts.benchmarks.benchmark_common import DATA_ROOT, FIGURE_ROOT, environment, plot_style, save_figure, write_json


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIRS = {medium: FIGURE_ROOT / "magnetized_rs_sigma" / medium for medium in ("ism", "wind")}
SIGMAS = (0.0, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0)
PLOT_SIGMAS = (0.0, 1.0e-3, 1.0e-1, 1.0, 10.0, 1000.0)
COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#000000")
MARKERS = ("o", "s", "^", "D", "P", "X")
SCENARIO = {
    "e_iso": 1.0e54,
    "gamma0": 150.0,
    "duration": 500.0,
    "theta_j": 1.0e-1,
    "z": 0.4,
    "fwd_eps_e": 1.0e-1,
    "fwd_eps_b": 1.0e-3,
    "fwd_p": 2.5,
    "fwd_xi": 1.0e-1,
    "rs_eps_e": 0.3,
    "rs_eps_b": 1.0e-1,
    "rs_p": 2.4,
    "rs_xi": 1.0,
}
MODE_GRIDS = {
    "quick": {
        "dyn_r": 80,
        "dyn_theta": 20,
        "model_gam": 41,
        "model_nu": 51,
        "model_r": 80,
        "model_theta": 20,
        "model_tobs": 32,
        "lc_times": 160,
        "sed_freqs": 96,
    },
    "formal": {
        "dyn_r": 120,
        "dyn_theta": 32,
        "model_gam": 81,
        "model_nu": 81,
        "model_r": 140,
        "model_theta": 32,
        "model_tobs": 48,
        "lc_times": 160,
        "sed_freqs": 180,
    },
}

plt.rcParams.update(plot_style())


def _cooling_index(cooling_mode: str) -> int:
    if cooling_mode == "none":
        return 0
    if cooling_mode == "numeric_ic_kn":
        return 1
    if cooling_mode == "nakar_y_thomson":
        return 2
    raise ValueError("cooling mode must be none, numeric_ic_kn, or nakar_y_thomson")


def _medium_config(medium: str) -> dict[str, float]:
    if medium == "ism":
        return {"d_ne": 10.0, "a_star": -1.0, "r0": 0.0}
    if medium == "wind":
        return {"d_ne": 1.0e-15, "a_star": 1.0e-1, "r0": 0.0}
    raise ValueError(f"unknown medium {medium!r}")


def _medium_model(medium: str):
    if medium == "ism":
        return UniformMedium(number_density_cm3=10.0)
    if medium == "wind":
        return WindMedium(a_star=1.0e-1, density_floor_cm3=1.0e-15, density_cap_cm3=None)
    raise ValueError(f"unknown medium {medium!r}")


def _runtime_config(sigma: float, medium: str, grid: dict[str, int], cooling_mode: str) -> RuntimeConfig:
    medium_kwargs = _medium_config(medium)
    return RuntimeConfig(
        index_y=_cooling_index(cooling_mode),
        index_syn_integr=2,
        electron_solver="fullhide_1d",
        num_threads=1,
        num_gam_e=grid["model_gam"],
        num_nu=61,
        num_r=grid["dyn_r"],
        eats_num_theta=grid["dyn_theta"],
        num_tobs=40,
        t_obs_min_log10=-2.0,
        t_obs_max_log10=7.0,
        reverse=True,
        include_forward_ssc=False,
        e_iso=SCENARIO["e_iso"],
        eta_0=SCENARIO["gamma0"],
        opening_angle_jet=SCENARIO["theta_j"],
        z=SCENARIO["z"],
        epsilon_e=SCENARIO["fwd_eps_e"],
        epsilon_b=SCENARIO["fwd_eps_b"],
        p=SCENARIO["fwd_p"],
        f_e=SCENARIO["fwd_xi"],
        **medium_kwargs,
        reverse_shock=ReverseShockConfig(
            enabled=True,
            delta_t_s=SCENARIO["duration"],
            epsilon_e=SCENARIO["rs_eps_e"],
            epsilon_b=SCENARIO["rs_eps_b"],
            p=SCENARIO["rs_p"],
            f_e=SCENARIO["rs_xi"],
            sigma=sigma,
        ),
    )


def _model(sigma: float, medium: str, grid: dict[str, int], cooling_mode: str) -> Model:
    return Model(
        jet=top_hat_jet(
            energy_iso_erg=SCENARIO["e_iso"],
            initial_lorentz_factor=SCENARIO["gamma0"],
            opening_angle_rad=SCENARIO["theta_j"],
            shell_duration_s=SCENARIO["duration"],
            magnetar=None,
            spreading=False,
        ),
        medium=_medium_model(medium),
        observer=Observer(
            z=SCENARIO["z"],
            viewing_angle_rad=0.0,
            viewing_azimuth_rad=0.0,
            luminosity_distance_cm=None,
        ),
        fwd_rad=radiation(
            epsilon_e=SCENARIO["fwd_eps_e"],
            epsilon_B=SCENARIO["fwd_eps_b"],
            p=SCENARIO["fwd_p"],
            accelerated_electron_fraction=SCENARIO["fwd_xi"],
            include_ssc=False,
        ),
        rvs_rad=radiation(
            epsilon_e=SCENARIO["rs_eps_e"],
            epsilon_B=SCENARIO["rs_eps_b"],
            p=SCENARIO["rs_p"],
            accelerated_electron_fraction=SCENARIO["rs_xi"],
            include_ssc=False,
        ),
        numerics=numerics(
            num_radius=grid["model_r"],
            eats_num_theta=grid["model_theta"],
            eats_num_phi=1,
            num_observer_time=grid["model_tobs"],
            num_electron_gamma=grid["model_gam"],
            num_photon_frequency=grid["model_nu"],
            num_threads=1,
            initial_radius_cm=1.0e14,
        ),
        observer_grid=observer_grid(time_min_s=1.0e-2, time_max_s=1.0e7),
        solver_options=solver_options(
            electron_solver="fullhide_1d",
            ssc_cooling_mode=cooling_mode,
        ),
        reverse_shock=reverse_shock(
            enabled=True,
            shell_duration_s=SCENARIO["duration"],
            upstream_sigma=sigma,
            include_ssc=False,
            include_cross_zone_ic=False,
        ),
        hadronic=hadronic(),
    )


def _style(index: int) -> dict[str, object]:
    return {
        "color": COLORS[index % len(COLORS)],
        "linestyle": "-" if index < len(COLORS) else "--",
        "marker": MARKERS[index % len(MARKERS)],
        "markevery": 14,
        "markersize": 4.0,
        "linewidth": 1.8,
    }


def _positive(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr) & (arr > 0.0)]


def _positive_for_log(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float).copy()
    arr[(arr <= 0.0) | ~np.isfinite(arr)] = np.nan
    return arr


def _logjump(values: np.ndarray) -> float:
    positive = _positive(values)
    if positive.size < 2:
        return float("nan")
    return float(np.nanmax(np.abs(np.diff(np.log10(positive)))))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _solve_dynamics_scan(
    sigmas: tuple[float, ...],
    medium: str,
    grid: dict[str, int],
    cooling_mode: str,
) -> tuple[list[dict[str, object]], list[tuple[float, np.ndarray, object, np.ndarray, np.ndarray]]]:
    rows: list[dict[str, object]] = []
    series = []
    for sigma in sigmas:
        config = _runtime_config(sigma, medium, grid, cooling_mode)
        observer_time = np.logspace(-2.0, 7.0, 12)
        query_config = query_cfg(config, observer_time)
        setup = build_setup(query_config, np.array([1.0e9, 1.0e14]), observer_time_s=observer_time)
        state = solve_setup(query_config, setup)
        rs = state.dynamics.reverse_shock
        if rs is None:
            raise RuntimeError("reverse shock state missing")
        causality = rs.causality
        if causality is None:
            raise RuntimeError("reverse shock causality diagnostics missing")
        t_obs = np.asarray(state.dynamics.r_tobs, dtype=float)
        internal_energy = np.asarray(rs.internal_energy_erg, dtype=float)
        swept_mass = np.asarray(rs.swept_mass_g, dtype=float)
        e3 = internal_energy / np.asarray(rs.comoving_volume_cm3, dtype=float)
        u3_per_m3 = np.divide(internal_energy, swept_mass, out=np.zeros_like(internal_energy), where=swept_mass > 0.0)
        post = np.asarray(state.dynamics.radius, dtype=float) >= float(rs.r_cross)
        post_u = u3_per_m3[post]
        post_mono = bool(np.all(np.diff(post_u) <= 0.0)) if post_u.size > 2 else True
        rows.append(
            {
                "medium": medium,
                "cooling_mode": cooling_mode,
                "sigma": sigma,
                "t_cross_s": float(rs.t_cross),
                "r_cross_cm": float(rs.r_cross),
                "global_rs_allowed": bool(causality.global_reverse_shock_allowed),
                "pressure_balance_condition_seen": bool(causality.pressure_balance_condition_seen),
                "local_fast_condition_seen": bool(causality.local_fast_condition_seen),
                "actual_rs_started": bool(causality.reverse_shock_started),
                "criteria_agree": bool(causality.criteria_agree),
                "contact_radius_cm": float(causality.contact_radius_cm),
                "reference_crossing_radius_cm": float(causality.reference_crossing_radius_cm),
                "pressure_balance_start_radius_cm": float(causality.pressure_balance_start_radius_cm),
                "pressure_balance_start_tobs_s": float(causality.pressure_balance_start_tobs_s),
                "pressure_balance_start_ratio": float(causality.pressure_balance_start_ratio),
                "pressure_balance_start_contact_fraction": float(causality.pressure_balance_start_contact_fraction),
                "fast_wave_crossing_radius_cm": float(causality.fast_wave_crossing_radius_cm),
                "fast_wave_crossing_tobs_s": float(causality.fast_wave_crossing_tobs_s),
                "local_start_radius_cm": float(causality.local_start_radius_cm),
                "local_start_tobs_s": float(causality.local_start_tobs_s),
                "actual_start_radius_cm": float(causality.actual_start_radius_cm),
                "actual_start_tobs_s": float(causality.actual_start_tobs_s),
                "actual_start_pressure_ratio": float(causality.actual_start_pressure_ratio),
                "actual_start_contact_fraction": float(causality.actual_start_contact_fraction),
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


def _solve_observable_scan(sigmas: tuple[float, ...], medium: str, grid: dict[str, int], cooling_mode: str):
    lc_times = np.logspace(-2.0, 7.0, grid["lc_times"])
    lc_bands = np.array([1.0e9, 1.0e14, 1.0e18], dtype=float)
    sed_times = np.array([1.0e2, 1.0e3, 1.0e4, 1.0e5], dtype=float)
    sed_freqs = np.logspace(7.0, 22.0, grid["sed_freqs"])
    lc_data = {}
    sed_data = {}
    lc_rows: list[dict[str, object]] = []
    sed_rows: list[dict[str, object]] = []
    for sigma in sigmas:
        model = _model(sigma, medium, grid, cooling_mode)
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
                    "medium": medium,
                    "cooling_mode": cooling_mode,
                    "sigma": sigma,
                    "nu_hz": float(nu_hz),
                    "total_peak_time_s": float(lc_times[peak_idx]),
                    "total_peak_fnu_cgs": float(total[peak_idx]),
                    "rs_peak_time_s": float(lc_times[rs_peak_idx]),
                    "rs_peak_fnu_cgs": float(rs[rs_peak_idx]),
                    "rs_to_total_at_total_peak": (
                        float(rs[peak_idx] / total[peak_idx]) if total[peak_idx] > 0.0 else float("nan")
                    ),
                }
            )
        for i_time, time_s in enumerate(sed_times):
            nu_fnu = sed_freqs * sed_total[:, i_time]
            rs_nu_fnu = sed_freqs * sed_rs[:, i_time]
            peak_idx = int(np.nanargmax(nu_fnu))
            rs_peak_idx = int(np.nanargmax(rs_nu_fnu)) if np.nanmax(rs_nu_fnu) > 0.0 else 0
            sed_rows.append(
                {
                    "medium": medium,
                    "cooling_mode": cooling_mode,
                    "sigma": sigma,
                    "time_s": float(time_s),
                    "total_peak_nu_hz": float(sed_freqs[peak_idx]),
                    "total_peak_nu_fnu_cgs": float(nu_fnu[peak_idx]),
                    "rs_peak_nu_hz": float(sed_freqs[rs_peak_idx]),
                    "rs_peak_nu_fnu_cgs": float(rs_nu_fnu[rs_peak_idx]),
                }
            )
    return lc_times, lc_bands, lc_data, lc_rows, sed_times, sed_freqs, sed_data, sed_rows


def _plot_dynamics(path: Path, series, medium: str, cooling_mode: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    x_min = 1.0e-2
    x_max = max(float(item[1][-1]) for item in series)
    panels = [
        ("B3 total [G]", lambda item: item[2].magnetic_field_g, 1.0e-12),
        ("gamma34", lambda item: item[2].gamma34, 1.0e-12),
        ("U3/V3 [erg cm^-3]", lambda item: item[3], 1.0e-10),
        ("U3/M3 [erg g^-1]", lambda item: item[4], 1.0e-10),
    ]
    for ax, (ylabel, getter, rel_floor) in zip(axes.flat, panels):
        for i, item in enumerate(series):
            sigma, t_obs = item[0], item[1]
            y_plot = _positive_for_log(getter(item))
            y_plot = np.where(t_obs >= x_min, y_plot, np.nan)
            ax.loglog(t_obs, y_plot, label=f"sigma={sigma:g}", **_style(i))
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.25)
    axes[0, 0].legend(fontsize=8, ncol=2)
    for ax in axes[-1, :]:
        ax.set_xlabel("observer time [s]")
    axes[0, 0].set_xlim(x_min, x_max)
    fig.suptitle(f"ASGARD magnetized RS sigma scan ({medium}, {cooling_mode}): dynamics", y=0.995)
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)


def _plot_lightcurves(
    path: Path,
    sigmas: tuple[float, ...],
    lc_times: np.ndarray,
    lc_bands: np.ndarray,
    lc_data: dict[float, tuple[np.ndarray, np.ndarray]],
    medium: str,
    cooling_mode: str,
) -> None:
    labels = ("1 GHz", "1e14 Hz", "1e18 Hz")
    fig, axes = plt.subplots(3, 2, figsize=(13, 12), sharex=True)
    for i_band, label in enumerate(labels):
        for i_sigma, sigma in enumerate(sigmas):
            total, rs = lc_data[sigma]
            axes[i_band, 0].loglog(
                lc_times,
                _positive_for_log(total[i_band]),
                label=f"sigma={sigma:g}" if i_band == 0 else None,
                **_style(i_sigma),
            )
            axes[i_band, 1].loglog(
                lc_times,
                _positive_for_log(rs[i_band]),
                **_style(i_sigma),
            )
        axes[i_band, 0].set_ylabel(f"total Fnu {label}")
        axes[i_band, 1].set_ylabel(f"RS sync Fnu {label}")
        for ax in axes[i_band, :]:
            ax.grid(True, which="both", alpha=0.25)
    axes[0, 0].legend(fontsize=8, ncol=2)
    for ax in axes[-1, :]:
        ax.set_xlabel("observer time [s]")
    axes[0, 0].set_xlim(lc_times[0], lc_times[-1])
    fig.suptitle(f"ASGARD magnetized RS sigma scan ({medium}, {cooling_mode}): light curves", y=0.995)
    fig.tight_layout()
    save_figure(fig, path)
    plt.close(fig)


def _plot_sed(
    path_total: Path,
    path_rs: Path,
    sigmas: tuple[float, ...],
    sed_times: np.ndarray,
    sed_freqs: np.ndarray,
    sed_data: dict[float, tuple[np.ndarray, np.ndarray]],
    medium: str,
    cooling_mode: str,
) -> None:
    panels = (
        ("total", path_total, f"ASGARD magnetized RS sigma scan ({medium}, {cooling_mode}): total SED", 1.0e-10),
        ("rs", path_rs, f"ASGARD magnetized RS sigma scan ({medium}, {cooling_mode}): RS synch SED", 1.0e-8),
    )
    for component, path, title, rel_floor in panels:
        fig, axes = plt.subplots(2, 2, figsize=(13, 10), sharex=True)
        for i_time, (ax, time_s) in enumerate(zip(axes.flat, sed_times)):
            for i_sigma, sigma in enumerate(sigmas):
                total, rs = sed_data[sigma]
                matrix = total if component == "total" else rs
                y = sed_freqs * matrix[:, i_time]
                ax.loglog(
                    sed_freqs,
                    _positive_for_log(y),
                    label=f"sigma={sigma:g}" if i_time == 0 else None,
                    **_style(i_sigma),
                )
            ax.set_title(f"t = {time_s:g} s")
            ax.grid(True, which="both", alpha=0.25)
        for ax in axes[:, 0]:
            ax.set_ylabel("nu Fnu [cgs]")
        for ax in axes[-1, :]:
            ax.set_xlabel("frequency [Hz]")
        axes[0, 0].legend(fontsize=8, ncol=2)
        axes[0, 0].set_xlim(sed_freqs[0], sed_freqs[-1])
        fig.suptitle(title, y=0.995)
        fig.tight_layout()
        save_figure(fig, path)
        plt.close(fig)


def run_medium(medium: str, mode: str, cooling_mode: str, output_dir: Path | None, data_root: Path) -> Path:
    grid = MODE_GRIDS[mode]
    outdir = OUTPUT_DIRS[medium] if output_dir is None else output_dir
    outdir.mkdir(parents=True, exist_ok=True)
    dyn_rows, dyn_series = _solve_dynamics_scan(SIGMAS, medium, grid, cooling_mode)
    data_dir = data_root / medium
    data_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(data_dir / "sigma_scan_summary.csv", dyn_rows)
    raw_dyn = []
    for sigma, time, reverse, energy_density, energy_mass in dyn_series:
        for values in zip(time, reverse.magnetic_field_g, reverse.gamma34, energy_density, energy_mass):
            raw_dyn.append(dict(zip(("time_s", "magnetic_field_g", "gamma34", "energy_density_erg_cm3", "specific_energy_erg_g"), values)) | {"sigma": sigma})
    _write_csv(data_dir / "sigma_scan_dynamics_raw.csv", raw_dyn)
    plot_series = [item for item in dyn_series if item[0] in PLOT_SIGMAS]
    _plot_dynamics(outdir / "magnetized_rs_sigma_dynamics", plot_series, medium, cooling_mode)

    lc_times, lc_bands, lc_data, lc_rows, sed_times, sed_freqs, sed_data, sed_rows = _solve_observable_scan(
        SIGMAS,
        medium,
        grid,
        cooling_mode,
    )
    _write_csv(data_dir / "sigma_scan_lightcurve_summary.csv", lc_rows)
    _write_csv(data_dir / "sigma_scan_sed_summary.csv", sed_rows)
    raw_lc = []
    raw_sed = []
    for sigma in SIGMAS:
        total, reverse = lc_data[sigma]
        for band, time, total_value, reverse_value in zip(np.repeat(lc_bands, lc_times.size), np.tile(lc_times, lc_bands.size), total.ravel(), reverse.ravel()):
            raw_lc.append({"sigma": sigma, "band_hz": band, "time_s": time, "total_fnu_cgs": total_value, "reverse_sync_fnu_cgs": reverse_value})
        total, reverse = sed_data[sigma]
        for time, frequency, total_value, reverse_value in zip(np.repeat(sed_times, sed_freqs.size), np.tile(sed_freqs, sed_times.size), total.T.ravel(), reverse.T.ravel()):
            raw_sed.append({"sigma": sigma, "time_s": time, "frequency_hz": frequency, "total_fnu_cgs": total_value, "reverse_sync_fnu_cgs": reverse_value})
    _write_csv(data_dir / "sigma_scan_lightcurve_raw.csv", raw_lc)
    _write_csv(data_dir / "sigma_scan_sed_raw.csv", raw_sed)
    write_json(data_dir / "metadata.json", environment(mode, threads=1, grid=grid, repeats=1) | {"medium": medium, "cooling_mode": cooling_mode, "sigmas": SIGMAS, "plotted_sigmas": PLOT_SIGMAS})
    _plot_lightcurves(
        outdir / "magnetized_rs_sigma_lightcurves",
        PLOT_SIGMAS,
        lc_times,
        lc_bands,
        lc_data,
        medium,
        cooling_mode,
    )
    _plot_sed(
        outdir / "magnetized_rs_sigma_sed_total",
        outdir / "magnetized_rs_sigma_sed_rs",
        PLOT_SIGMAS,
        sed_times,
        sed_freqs,
        sed_data,
        medium,
        cooling_mode,
    )
    return outdir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--medium", choices=("ism", "wind", "both"), default="both")
    parser.add_argument("--mode", choices=tuple(MODE_GRIDS), default="formal")
    parser.add_argument(
        "--cooling-mode",
        choices=("none", "numeric_ic_kn", "nakar_y_thomson"),
        default="nakar_y_thomson",
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--data-dir", type=Path, default=DATA_ROOT / "magnetized_rs_sigma")
    args = parser.parse_args()
    if args.output_dir is not None and args.medium == "both":
        raise ValueError("--output-dir is only valid when --medium is ism or wind")
    media = ("ism", "wind") if args.medium == "both" else (args.medium,)
    for medium in media:
        outdir = run_medium(medium, args.mode, args.cooling_mode, args.output_dir, args.data_dir)
        print(f"wrote {outdir}")


if __name__ == "__main__":
    main()
