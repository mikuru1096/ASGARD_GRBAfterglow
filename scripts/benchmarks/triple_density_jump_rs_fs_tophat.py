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

from asgard_core import Model, Observer, TabulatedMedium, UniformMedium, WindMedium, top_hat_jet
from tests.public_api_builders import hadronic, numerics, observer_grid, radiation, reverse_shock, solver_options
from scripts.benchmarks.benchmark_common import DATA_ROOT, FIGURE_ROOT, environment, plot_style, save_figure, write_json
from asgard_core.asgard_physics_utils import ambient_density


OUTPUT_PATH = FIGURE_ROOT / "triple_density_jump" / "triple_density_jump_rs_fs_tophat"
JUMP_RADII_CM = (1.0e15, 1.0e16, 1.0e17)
JUMP_FACTOR = 1.0e2
JUMP_WIDTH_REL = 1.0e-1
REVERSE_SIGMA = 1.0e-1
BANDS_HZ = np.array([1.0e9, 1.0e14, 1.0e18], dtype=float)
BAND_LABELS = ("1 GHz radio", "1e14 Hz optical", "1e18 Hz X-ray")
PROFILE_DOMAIN_CM = (1.0e13, 3.0e17)
PROFILE_POINTS = 96
PROFILE_BASE_DENSITY_CM3 = 1.0
CSM_DIR = ROOT / "paper" / "source_data" / "csm"
PION_RAW_PROFILE = CSM_DIR / "pion_eta_car_sph1d_n2048_external_raw.csv"
PION_INTERFACE_PROFILE = CSM_DIR / "pion_eta_car_sph1d_n2048_external_96knots.csv"
PION_INJECTION_RADIUS_CM = 1.32e16
PION_OBSERVER_TIME_MIN_S = 1.0e2
PION_OBSERVER_TIME_MAX_S = 5.0e5
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

plt.rcParams.update(plot_style())


def _model(
    *,
    with_jumps: bool,
    grid: dict[str, int],
    medium=None,
    initial_radius_cm: float = 1.0e13,
    observer_time_min_s: float = 1.0e0,
    observer_time_max_s: float = 1.0e8,
    reverse_enabled: bool = True,
) -> Model:
    model = Model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=100.0,
            opening_angle_rad=0.1,
            shell_duration_s=10.0,
            magnetar=None,
            spreading=False,
        ),
        medium=UniformMedium(number_density_cm3=1.0) if medium is None else medium,
        observer=Observer(z=0.1, viewing_angle_rad=0.0, viewing_azimuth_rad=0.0, luminosity_distance_cm=1.0e26),
        fwd_rad=radiation(
            epsilon_e=3.0e-1,
            epsilon_B=1.0e-5,
            p=2.3,
            accelerated_electron_fraction=1.0e-1,
            include_ssc=False,
        ),
        rvs_rad=(radiation(
            epsilon_e=3.0e-1,
            epsilon_B=1.0e-5,
            p=2.4,
            accelerated_electron_fraction=1.0e-1,
            include_ssc=False,
        ) if reverse_enabled else None),
        numerics=numerics(
            num_radius=grid["num_r"],
            eats_num_theta=grid["num_theta"],
            eats_num_phi=1,
            num_observer_time=grid["num_tobs"],
            num_electron_gamma=grid["num_gam_e"],
            num_photon_frequency=grid["num_nu"],
            num_threads=1,
            initial_radius_cm=initial_radius_cm,
        ),
        observer_grid=observer_grid(time_min_s=observer_time_min_s, time_max_s=observer_time_max_s),
        solver_options=solver_options(
            electron_solver="fullhide_1d",
            ssc_cooling_mode="none",
        ),
        reverse_shock=reverse_shock(
            enabled=reverse_enabled,
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


def _write_secondary_event_csv(details, output: Path) -> Path:
    if details.rev is None or details.rev.secondary_rs_event_active is None:
        raise RuntimeError("triple density-jump benchmark requires secondary RS event diagnostics")
    csv_path = output.with_suffix("").with_name(output.with_suffix("").name + "_secondary_rs_events.csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
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


def _event_energy_rows(track, epsilon_e: float) -> list[tuple[int, float, float, float, float]]:
    dissipated = np.asarray(track.secondary_rs_dissipated_energy_erg, dtype=float)
    electron = np.asarray(track.secondary_rs_electron_injected_energy_erg, dtype=float)
    branch_gm = np.asarray(track.secondary_rs_branch_gamma_m, dtype=float)
    if branch_gm.shape[1] != dissipated.size:
        raise RuntimeError("secondary RS branch and total energy grids do not match")
    owner = branch_gm > 1.0
    if np.any(np.count_nonzero(owner, axis=0) > 1):
        raise RuntimeError("event-resolved energy requires non-overlapping dissipative branches")
    if np.any((dissipated > 0.0) & ~np.any(owner, axis=0)):
        raise RuntimeError("positive secondary RS dissipation has no owning branch")
    rows = []
    for index in range(branch_gm.shape[0]):
        diss_j = float(np.sum(dissipated[owner[index]]))
        electron_j = float(np.sum(electron[owner[index]]))
        expected_j = epsilon_e * diss_j
        error_j = 0.0 if expected_j == 0.0 else (electron_j - expected_j) / expected_j
        rows.append((index, diss_j, electron_j, expected_j, error_j))
    return rows


def _write_secondary_energy_csv(details, model: Model, output: Path) -> Path:
    if (details.rev is None or details.rev.secondary_rs_dissipated_energy_erg is None
            or details.rev.secondary_rs_electron_injected_energy_erg is None
            or details.rev.secondary_rs_branch_gamma_m is None):
        raise RuntimeError("density-structure benchmark requires secondary RS energy diagnostics")
    rows = _event_energy_rows(details.rev, float(model.rvs_rad.eps_e))
    total_dissipated = float(sum(row[1] for row in rows))
    total_electron = float(sum(row[2] for row in rows))
    if total_dissipated <= 0.0:
        raise RuntimeError("secondary RS dissipated energy must be positive")
    expected_electron = float(model.rvs_rad.eps_e) * total_dissipated
    csv_path = output.with_suffix("").with_name(output.with_suffix("").name + "_secondary_rs_energy.csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("event_index", "secondary_rs_dissipated_energy_erg",
                         "secondary_rs_electron_injected_energy_erg", "expected_electron_energy_erg",
                         "electron_energy_fractional_error"))
        writer.writerows(rows)
        writer.writerow(("total", total_dissipated, total_electron, expected_electron,
                         (total_electron - expected_electron) / expected_electron))
    return csv_path


def _write_adaptive_convergence_csv(model: Model, output: Path, times: np.ndarray, direct_total: np.ndarray) -> Path:
    adaptive = model.flux_density_grid_adaptive(times, BANDS_HZ)
    adaptive_total = np.asarray(adaptive.total, dtype=float)
    adaptive_times = np.asarray(adaptive.time_s, dtype=float)
    csv_path = output.with_suffix("").with_name(output.with_suffix("").name + "_adaptive_convergence.csv")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
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



def _pion_lbv_profile() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with PION_RAW_PROFILE.open(encoding="utf-8", newline="") as stream:
        raw = list(csv.DictReader(stream))
    with PION_INTERFACE_PROFILE.open(encoding="utf-8", newline="") as stream:
        interface = list(csv.DictReader(stream))
    raw_index = np.asarray([int(row["raw_index"]) for row in raw], dtype=int)
    raw_r = np.asarray([float(row["r_cm"]) for row in raw], dtype=float)
    raw_n = np.asarray([float(row["proton_equivalent_n_cm3"]) for row in raw], dtype=float)
    interface_index = np.asarray([int(row["raw_index"]) for row in interface], dtype=int)
    interface_r = np.asarray([float(row["r_cm"]) for row in interface], dtype=float)
    interface_n = np.asarray([float(row["proton_equivalent_n_cm3"]) for row in interface], dtype=float)
    if raw_r.size != 1968 or interface_r.size != PROFILE_POINTS:
        raise RuntimeError("PION CSM source rows do not match the tracked 2048-cell external domain")
    if raw_r[0] < PION_INJECTION_RADIUS_CM or interface_r[0] != raw_r[0]:
        raise RuntimeError("PION CSM interface must begin at the first cell center outside the wind injection sphere")
    if interface_r[-1] != raw_r[-1] or np.any(np.diff(raw_r) <= 0.0) or np.any(np.diff(interface_r) <= 0.0):
        raise RuntimeError("PION CSM radius grids must share ordered endpoints")
    if not np.all(np.isfinite(raw_n)) or not np.all(raw_n > 0.0) or not np.all(np.isfinite(interface_n)) or not np.all(interface_n > 0.0):
        raise RuntimeError("PION proton-equivalent density must be finite and positive")
    positions = np.searchsorted(raw_index, interface_index)
    if np.any(positions >= raw_index.size) or not np.array_equal(raw_index[positions], interface_index):
        raise RuntimeError("PION 96-point interface contains an index absent from the raw profile")
    if not np.array_equal(raw_r[positions], interface_r) or not np.array_equal(raw_n[positions], interface_n):
        raise RuntimeError("PION interface points must be unchanged raw simulation cells")
    return raw_r, raw_n, interface_r, interface_n


def _tabulated_profile(scenario: str) -> tuple[np.ndarray, np.ndarray, tuple[float, ...]]:
    edge_points = []
    for center in JUMP_RADII_CM:
        half_width = 4.0 * JUMP_WIDTH_REL * center
        edge_points.extend((center - half_width, center, center + half_width))
    base_count = PROFILE_POINTS - len(edge_points)
    radius = np.unique(np.concatenate((
        np.logspace(np.log10(PROFILE_DOMAIN_CM[0]), np.log10(PROFILE_DOMAIN_CM[1]), base_count),
        np.asarray(edge_points, dtype=float),
    )))
    if radius.size != PROFILE_POINTS:
        raise RuntimeError("tabulated density profile grid does not contain 96 unique points")
    density = np.full_like(radius, PROFILE_BASE_DENSITY_CM3)

    def add_clump(values: np.ndarray, center: float) -> None:
        half_width = 4.0 * JUMP_WIDTH_REL * center
        x = (radius - center) / half_width
        inside = np.abs(x) < 1.0
        bump = np.zeros_like(radius)
        bump[inside] = np.exp(1.0 - 1.0 / (1.0 - x[inside] ** 2))
        values += PROFILE_BASE_DENSITY_CM3 * (JUMP_FACTOR - 1.0) * bump

    if scenario == "fig10":
        termination, shell_center = JUMP_RADII_CM[:2]
        density = np.where(radius < termination,
                           PROFILE_BASE_DENSITY_CM3 * (radius / termination) ** -2.0,
                           PROFILE_BASE_DENSITY_CM3)
        add_clump(density, shell_center)
        return radius, density, (shell_center,)
    if scenario == "fig12":
        for center in JUMP_RADII_CM:
            add_clump(density, center)
        return radius, density, JUMP_RADII_CM
    raise ValueError(f"unknown density-structure scenario: {scenario}")


def _write_profile_csv(path: Path, radius: np.ndarray, density: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("radius_cm", "density_cm3"))
        writer.writerows(zip(radius, density))


def _write_profile_dynamics_csv(track, path: Path) -> None:
    mass = np.asarray(track.secondary_rs_branch_swept_mass_g, dtype=float)
    energy = np.asarray(track.secondary_rs_branch_internal_energy_erg, dtype=float)
    volume = np.asarray(track.secondary_rs_branch_comoving_volume_cm3, dtype=float)
    bfield = np.asarray(track.secondary_rs_branch_B, dtype=float)
    gamma_m = np.asarray(track.secondary_rs_branch_gamma_m, dtype=float)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("event_index", "radius_index", "observer_time_s", "radius_cm", "gamma_bulk",
                         "forward_B_g", "branch_swept_mass_g", "branch_internal_energy_erg",
                         "branch_comoving_volume_cm3", "branch_B_g", "branch_gamma_m"))
        for event in range(mass.shape[0]):
            for index in range(mass.shape[1]):
                writer.writerow((event, index, track.t_obs[index], track.radius[index], track.Gamma[index],
                                 track.B_comv[index], mass[event, index], energy[event, index],
                                 volume[event, index], bfield[event, index], gamma_m[event, index]))


def _write_profile_event_csv(track, centers: tuple[float, ...], path: Path) -> None:
    active = np.asarray(track.secondary_rs_event_active, dtype=bool)
    if active.size != len(centers):
        raise RuntimeError("tabulated profile event count does not match the designed compression centers")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("event_index", "input_center_cm", "event_active", "start_radius_cm",
                         "shock_end_radius_cm", "start_tobs_axis_s", "shock_end_tobs_axis_s"))
        for index, center in enumerate(centers):
            writer.writerow((index, center, bool(active[index]), track.secondary_rs_start_radius[index],
                             track.secondary_rs_shock_end_radius[index],
                             track.secondary_rs_start_tobs_axis[index],
                             track.secondary_rs_shock_end_tobs_axis[index]))


def _write_profile_flux_csv(times: np.ndarray, flux, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("band_hz", "time_s", "total_flux_cgs", "forward_sync_flux_cgs",
                         "reverse_sync_flux_cgs"))
        for band_index, band in enumerate(BANDS_HZ):
            for time_index, time in enumerate(times):
                writer.writerow((band, time, flux.total[band_index, time_index],
                                 flux.fwd.sync[band_index, time_index],
                                 flux.rev.sync[band_index, time_index]))


def build_density_structure(*, scenario: str, mode: str, data_dir: Path,
                            times_count: int | None = None) -> Path:
    grid = dict(MODE_GRIDS[mode])
    if times_count is not None:
        grid["times"] = int(times_count)
    radius, density, centers = _tabulated_profile(scenario)
    model = _model(with_jumps=False, grid=grid, medium=TabulatedMedium(radius, density, scenario))
    times = np.logspace(0.0, 8.0, grid["times"])
    flux = model.flux_density_grid(times, BANDS_HZ)
    details = model.details(float(times[0]), float(times[-1]))
    if details.rev is None or details.rev.secondary_rs_event_active is None:
        raise RuntimeError("tabulated density benchmark requires secondary RS diagnostics")
    name = {"fig10": "fig10_wind_termination_shell", "fig12": "fig12_smooth_clumpy_medium"}[scenario]
    stem = data_dir / name
    profile_path = stem.with_name(name + "_profile.csv")
    dynamics_path = stem.with_name(name + "_dynamics.csv")
    event_path = stem.with_name(name + "_events.csv")
    flux_path = stem.with_name(name + "_flux.csv")
    _write_profile_csv(profile_path, radius, density)
    _write_profile_dynamics_csv(details.rev, dynamics_path)
    _write_profile_event_csv(details.rev, centers, event_path)
    energy_path = _write_secondary_energy_csv(details, model, stem)
    _write_profile_flux_csv(times, flux, flux_path)
    meta = environment(mode, threads=1, grid=grid, repeats=1)
    meta.update({
        "scenario": scenario,
        "physical_assumptions": ("wind R^-2 to n0 plateau, finite compact-support shell, then n0 plateau"
                                 if scenario == "fig10"
                                 else "three continuous compact-support clumps on a constant n0 baseline"),
        "decision_value": ("isolate the secondary-shock response of a wind termination plus finite shell"
                           if scenario == "fig10"
                           else "test continuity and event attribution for three smooth tabulated clumps"),
        "profile_domain_cm": PROFILE_DOMAIN_CM,
        "profile_points": PROFILE_POINTS,
        "baseline_density_cm3": PROFILE_BASE_DENSITY_CM3,
        "common_event_centers_cm": JUMP_RADII_CM,
        "active_compression_centers_cm": centers,
        "density_contrast": JUMP_FACTOR,
        "relative_width": JUMP_WIDTH_REL,
        "event_rule": "contiguous dn/dR>0 intervals are candidates; reverse_contact source must be positive",
        "source_data": [str(profile_path), str(dynamics_path), str(event_path),
                        str(energy_path), str(flux_path)],
    })
    metadata_path = stem.with_name(name + "_metadata.json")
    write_json(metadata_path, meta)
    for output in (profile_path, dynamics_path, event_path, energy_path, flux_path, metadata_path):
        print(f"wrote {output}")
    return stem


def build_external_media(*, mode: str, data_dir: Path, times_count: int | None = None) -> Path:
    grid = dict(MODE_GRIDS[mode])
    if times_count is not None:
        grid["times"] = int(times_count)
    raw_r, raw_n, profile_r, profile_n = _pion_lbv_profile()
    times = np.logspace(
        np.log10(PION_OBSERVER_TIME_MIN_S),
        np.log10(PION_OBSERVER_TIME_MAX_S),
        grid["times"],
    )
    common = {
        "grid": grid,
        "initial_radius_cm": float(profile_r[0]),
        "observer_time_min_s": PION_OBSERVER_TIME_MIN_S,
        "observer_time_max_s": PION_OBSERVER_TIME_MAX_S,
        "reverse_enabled": False,
    }
    models = {
        "ism": _model(with_jumps=False, **common),
        "wind_r2": _model(
            with_jumps=False,
            medium=WindMedium(a_star=1.0e-2, density_floor_cm3=1.0, density_cap_cm3=1.0e4),
            **common,
        ),
        "explicit_jump": _model(with_jumps=True, **common),
        "pion_lbv": _model(
            with_jumps=False,
            medium=TabulatedMedium(profile_r, profile_n, "PION LBV eruption CSM"),
            **common,
        ),
    }
    jump_model = models["explicit_jump"]
    jump_model.setups.jump_r_cm = (JUMP_RADII_CM[2],)
    jump_model.setups.jump_factor = (JUMP_FACTOR,)
    jump_model.setups.jump_width = (JUMP_WIDTH_REL,)

    csv_path = data_dir / "fig2_external_media_response.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("medium", "kind", "radius_index", "radius_cm", "density_cm3",
                         "observer_time_s", "gamma_bulk", "bfield_g", "frequency_hz", "flux_cgs"))
        for label, model in models.items():
            flux = model.flux_density_grid(times, np.array([1.0e14]))
            track = model.details(float(times[0]), float(times[-1])).fwd
            density = np.asarray(ambient_density(track.radius, model.setups), dtype=float)
            state_value = np.column_stack(
                (track.radius, density, track.t_obs, track.Gamma, track.B_comv)
            )
            flux_value = np.asarray(flux.total, dtype=float)
            if not np.all(np.isfinite(state_value)) or not np.all(state_value > 0.0):
                raise RuntimeError(f"non-finite or non-positive {label} forward-shock state")
            if not np.all(np.isfinite(flux_value)) or not np.all(flux_value > 0.0):
                raise RuntimeError(f"non-finite or non-positive {label} forward-shock flux")
            if label == "pion_lbv":
                if np.min(track.radius) < raw_r[0] or np.max(track.radius) > raw_r[-1]:
                    raise RuntimeError(
                        "ASGARD PION-CSM trajectory leaves the tracked raw simulation domain: "
                        f"[{track.radius[0]:.12e}, {track.radius[-1]:.12e}] versus "
                        f"[{raw_r[0]:.12e}, {raw_r[-1]:.12e}]"
                    )
                for index, radius in enumerate(raw_r):
                    writer.writerow((label, "profile_raw", index, radius, raw_n[index],
                                     "", "", "", "", ""))
                for index, radius in enumerate(profile_r):
                    writer.writerow((label, "profile_interface", index, radius, profile_n[index],
                                     "", "", "", "", ""))
            for index, radius in enumerate(track.radius):
                writer.writerow((label, "state", index, radius, density[index], track.t_obs[index],
                                 track.Gamma[index], track.B_comv[index], "", ""))
            for index, time in enumerate(times):
                writer.writerow((label, "flux", "", "", "", time, "", "", 1.0e14, flux.total[0, index]))

    meta = environment(mode, threads=1, grid=grid, repeats=1)
    meta.update({
        "physical_assumptions": (
            "FS-only comparison with the same top-hat ejecta and microphysics for a uniform ISM, "
            "R^-2 wind, one explicit Gaussian density jump, and the PION Eta Car spherical-reduction "
            "late-eruption CSM; primary reverse shock disabled in all four columns"
        ),
        "decision_value": "trace n(R) into Gamma(R), B'(R), and the 1e14 Hz forward observer light curve",
        "frequency_hz": 1.0e14,
        "explicit_jump_radius_cm": JUMP_RADII_CM[2],
        "common_initial_radius_cm": float(profile_r[0]),
        "observer_time_min_s": PION_OBSERVER_TIME_MIN_S,
        "observer_time_domain_s": [float(times[0]), float(times[-1])],
        "pion_raw_radius_domain_cm": [float(raw_r[0]), float(raw_r[-1])],
        "pion_injection_radius_cm": PION_INJECTION_RADIUS_CM,
        "pion_density_column": "proton_equivalent_n_cm3 = rho_g_cm3 / m_p",
        "pion_hydrogen_diagnostic": "nH_X033_cm3 = 0.33 rho_g_cm3 / m_p; not passed to ASGARD dynamics",
        "pion_resolution": "1D spherical, 2048 cells; external-domain raw rows=1968; interface knots=96",
        "pion_convergence": (
            "structure positions and shell mass converge from 1024 to 2048, but cellwise "
            "log-density L1=0.0632 dex exceeds the 0.05-dex threshold; not fully converged"
        ),
        "pion_formal_angle_patch_dependency": False,
        "source_data": [
            str(csv_path),
            str(PION_RAW_PROFILE.relative_to(ROOT)),
            str(PION_INTERFACE_PROFILE.relative_to(ROOT)),
        ],
    })
    metadata_path = data_dir / "fig2_external_media_response.json"
    write_json(metadata_path, meta)
    print(f"wrote {csv_path}")
    print(f"wrote {metadata_path}")
    return csv_path

def build_plot(*, mode: str, output: Path, data_dir: Path, times_count: int | None = None) -> Path:
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
        fontsize=7.0,
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
    legend = fig.legend(handles, labels, fontsize=7.0, loc="upper center", ncol=4, frameon=False, handlelength=2.5, bbox_to_anchor=(0.5, 0.998))
    for line in legend.get_lines():
        line.set_linewidth(1.6)
    fig.subplots_adjust(top=0.93, bottom=0.075, left=0.12, right=0.985)
    paths = save_figure(fig, output)
    plt.close(fig)
    details = triple_jump_model.details(float(times[0]), float(times[-1]))
    data_stem = data_dir / "triple_density_jump"
    event_csv = _write_secondary_event_csv(details, data_stem)
    energy_csv = _write_secondary_energy_csv(details, triple_jump_model, data_stem)
    adaptive_csv = _write_adaptive_convergence_csv(triple_jump_model, data_stem, times, jump_total)
    raw_path = data_stem.with_name("triple_density_jump_flux.csv")
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    with raw_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("band_hz", "time_s", "no_jump_total", "jump_total", "jump_forward_sync", "jump_reverse_sync"))
        for band, time, base, total, fwd, rev in zip(np.repeat(BANDS_HZ, times.size), np.tile(times, BANDS_HZ.size), no_jump_total.ravel(), jump_total.ravel(), jump_fwd.ravel(), jump_rev.ravel()):
            writer.writerow((band, time, base, total, fwd, rev))
    meta = environment(mode, threads=1, grid=grid, repeats=1)
    meta.update({"jump_radii_cm": JUMP_RADII_CM, "jump_factor": JUMP_FACTOR, "jump_width_relative": JUMP_WIDTH_REL, "reverse_sigma": REVERSE_SIGMA, "plot_display_floor_peak_fraction": 1.0e-8})
    write_json(data_dir / "metadata.json", meta)

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
    parser = argparse.ArgumentParser(description="Generate density-structure RS+FS top-hat benchmarks.")
    parser.add_argument("--mode", choices=sorted(MODE_GRIDS), default="formal")
    parser.add_argument("--times", type=int, default=None)
    parser.add_argument("--scenario", choices=("triple", "fig2", "fig10", "fig12", "density_structure"),
                        default="triple")
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--data-dir", type=Path, default=DATA_ROOT / "triple_density_jump")
    args = parser.parse_args()
    if args.scenario == "triple":
        build_plot(mode=args.mode, output=args.output, data_dir=args.data_dir, times_count=args.times)
        return
    if args.scenario == "fig2":
        build_external_media(mode=args.mode, data_dir=args.data_dir, times_count=args.times)
        return
    scenarios = ("fig10", "fig12") if args.scenario == "density_structure" else (args.scenario,)
    for scenario in scenarios:
        build_density_structure(scenario=scenario, mode=args.mode, data_dir=args.data_dir,
                                times_count=args.times)


if __name__ == "__main__":
    main()
