from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

from _repo_path import ensure_repo_root_on_path


ensure_repo_root_on_path()

from asgard_core.api_model import _direct_tophat_patch_config, _solve_patch_state
from src import constants
from benchmark_theta_j_multiples_magnetic_decay import (
    EPSILON_B,
    EPSILON_B_FLOOR,
    MAGNETIC_DECAY_ALPHA_T,
    MAGNETIC_DECAY_T0_S,
    THETA_J_RAD,
    build_model,
)


NU_OBS_HZ = 1.0e10
THETA_MULTIPLES = np.array([2.0, 3.0, 4.0, 5.0], dtype=float)
SKYMAP_THETA_MULTIPLE = 3.0
CENTROID_TIMES_S = np.geomspace(1.0e4, 1.0e7, 64)
SKYMAP_TIMES_S = np.array([1.0e6, 3.0e6, 8.0e6], dtype=float)
SOLVE_TIMES_S = np.unique(np.concatenate((CENTROID_TIMES_S, SKYMAP_TIMES_S)))
MAP_NUM_THETA = 48
MAP_NUM_PHI_FULL = 96
SKYMAP_NPIX = 150
RAD_TO_MAS = 180.0 / np.pi * 3600.0 * 1.0e3
COLORS = ("#0072b2", "#d55e00", "#009e73", "#cc79a7")


@dataclass
class MomentSeries:
    time_s: np.ndarray
    flux: np.ndarray
    x_rad: np.ndarray
    y_rad: np.ndarray
    offset_rad: np.ndarray


@dataclass
class MapProduct:
    image_flux: np.ndarray
    x_mas: np.ndarray
    y_mas: np.ndarray
    flux: np.ndarray
    x_centroid_rad: np.ndarray
    y_centroid_rad: np.ndarray


def solve_state_for_solver(solver: str):
    model = build_model(solver, THETA_J_RAD * SKYMAP_THETA_MULTIPLE)
    config = _direct_tophat_patch_config(model)
    return model, _solve_patch_state(
        model,
        config,
        SOLVE_TIMES_S,
        np.array([NU_OBS_HZ], dtype=float),
        solve_reference_times_s=SOLVE_TIMES_S,
    )


def angular_grid(theta_j: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta_edges = np.linspace(0.0, theta_j, MAP_NUM_THETA + 1)
    theta = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    phi_edges = np.linspace(0.0, 2.0 * np.pi, MAP_NUM_PHI_FULL + 1)
    phi = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    domega_theta = np.cos(theta_edges[:-1]) - np.cos(theta_edges[1:])
    domega = domega_theta[:, None] * (phi_edges[1] - phi_edges[0])
    return theta, phi, domega


def sky_basis(theta_v: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sightline = direction_vector(theta_v, 0.0)
    trial = np.array([0.0, 0.0, 1.0], dtype=float)
    sky_x = np.cross(trial, sightline)
    if np.linalg.norm(sky_x) < 1.0e-12:
        sky_x = np.array([1.0, 0.0, 0.0], dtype=float)
    else:
        sky_x = sky_x / np.linalg.norm(sky_x)
    sky_y = np.cross(sightline, sky_x)
    return sightline, sky_x, sky_y / np.linalg.norm(sky_y)


def direction_vector(theta: float | np.ndarray, phi: float | np.ndarray) -> np.ndarray:
    return np.stack(
        (
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta) + np.zeros_like(np.asarray(phi, dtype=float)),
        ),
        axis=-1,
    )


def sky_coordinates_mas(radius_cm: np.ndarray, direction: np.ndarray, sightline, sky_x, sky_y, d_ang_cm: float):
    position = np.asarray(radius_cm, dtype=float)[:, None] * direction[None, :]
    los = np.sum(position * sightline[None, :], axis=1)
    transverse = position - los[:, None] * sightline[None, :]
    x_rad = transverse @ sky_x / d_ang_cm
    y_rad = transverse @ sky_y / d_ang_cm
    return x_rad * RAD_TO_MAS, y_rad * RAD_TO_MAS, x_rad, y_rad


def interp_sed_values(seed_log: np.ndarray, source: np.ndarray, shell_i: np.ndarray, shell_ratio: np.ndarray, target_log: np.ndarray) -> np.ndarray:
    values = np.zeros(target_log.size, dtype=float)
    inside = (target_log > seed_log[0]) & (target_log <= seed_log[-1])
    idx = np.flatnonzero(inside)
    i_hi = np.searchsorted(seed_log, target_log[idx], side="left")
    i_lo = i_hi - 1
    freq_ratio = (target_log[idx] - seed_log[i_lo]) / (seed_log[i_hi] - seed_log[i_lo])
    k = shell_i[idx]
    r = shell_ratio[idx]
    y_lo = (1.0 - r) * source[i_lo, k] + r * source[i_lo, k + 1]
    y_hi = (1.0 - r) * source[i_hi, k] + r * source[i_hi, k + 1]
    positive = (y_lo > 0.0) & (y_hi > 0.0)
    local = np.empty(idx.size, dtype=float)
    local[positive] = np.exp(np.log(y_lo[positive]) + freq_ratio[positive] * (np.log(y_hi[positive]) - np.log(y_lo[positive])))
    local[~positive] = (1.0 - freq_ratio[~positive]) * y_lo[~positive] + freq_ratio[~positive] * y_hi[~positive]
    values[idx] = local
    return values


def chi_escape(tau_front: np.ndarray, tau_cell: np.ndarray) -> np.ndarray:
    small = tau_cell <= 1.0e-6
    out = np.empty_like(tau_cell)
    out[small] = np.exp(-tau_front[small]) * (1.0 - 0.5 * tau_cell[small] + tau_cell[small] * tau_cell[small] / 6.0)
    out[~small] = np.exp(-tau_front[~small]) * (1.0 - np.exp(-tau_cell[~small])) / tau_cell[~small]
    return out


def accumulate_cell(
    *,
    times_s: np.ndarray,
    seed_log: np.ndarray,
    nu_obs_log: float,
    z: float,
    direction: np.ndarray,
    sightline,
    sky_x,
    sky_y,
    d_ang_cm: float,
    domega: float,
    arrival: np.ndarray,
    radius: np.ndarray,
    gamma: np.ndarray,
    source: np.ndarray,
    flux: np.ndarray,
    x_moment: np.ndarray,
    y_moment: np.ndarray,
    samples: list[tuple[int, float, float, float]] | None,
) -> None:
    if not np.all(arrival[1:] > arrival[:-1]):
        raise ValueError("skymap diagnostic requires monotonic EATS arrival time for each emitting cell.")
    shell_i = np.searchsorted(arrival, times_s, side="right") - 1
    valid = (shell_i >= 0) & (shell_i < arrival.size - 1)
    time_i = np.flatnonzero(valid)
    shell_i = shell_i[valid]
    ratio = (times_s[valid] - arrival[shell_i]) / (arrival[shell_i + 1] - arrival[shell_i])
    gamma_seg = np.exp(np.log(gamma[shell_i]) + ratio * (np.log(gamma[shell_i + 1]) - np.log(gamma[shell_i])))
    radius_seg = (1.0 - ratio) * radius[shell_i] + ratio * radius[shell_i + 1]
    beta = np.sqrt(1.0 - gamma_seg**-2)
    mu = float(direction @ sightline)
    doppler_inv = gamma_seg * (1.0 - beta * mu)
    weights = domega / (4.0 * np.pi) * doppler_inv**-3
    target_log = nu_obs_log + np.log(doppler_inv) + np.log1p(z)
    contribution = weights * interp_sed_values(seed_log, source, shell_i, ratio, target_log)
    x_mas, y_mas, x_rad, y_rad = sky_coordinates_mas(radius_seg, direction, sightline, sky_x, sky_y, d_ang_cm)
    np.add.at(flux, time_i, contribution)
    np.add.at(x_moment, time_i, contribution * x_rad)
    np.add.at(y_moment, time_i, contribution * y_rad)
    if samples is not None:
        positive = contribution > 0.0
        samples.extend(
            (int(time_i[j]), float(x_mas[j]), float(y_mas[j]), float(contribution[j]))
            for j in np.flatnonzero(positive)
        )


def state_moments_and_samples(state, theta_v: float, times_s: np.ndarray, *, collect_samples: bool = False):
    theta, phi, domega = angular_grid(THETA_J_RAD)
    z = float(state.config.z)
    d_ang_cm = float(state.setup.luminosity_distance_cm) / (1.0 + z) ** 2
    sightline, sky_x, sky_y = sky_basis(theta_v)
    seed_log = np.log(np.asarray(state.setup.seed_frequency_hz, dtype=float))
    nu_obs_log = float(np.log(NU_OBS_HZ))
    flux = np.zeros(times_s.size, dtype=float)
    x_moment = np.zeros_like(flux)
    y_moment = np.zeros_like(flux)
    samples: list[tuple[int, float, float, float]] | None = [] if collect_samples else None

    is_2d = state.electron.l_syn_spec_chi is not None and state.electron.tau_syn_chi is not None
    for i_theta, theta_center in enumerate(theta):
        for phi_center in phi:
            direction = direction_vector(float(theta_center), float(phi_center))
            omega = float(domega[i_theta, 0])
            if is_2d:
                accumulate_2d_direction(state, times_s, seed_log, nu_obs_log, z, direction, sightline, sky_x, sky_y, d_ang_cm, omega, flux, x_moment, y_moment, samples)
            else:
                accumulate_1d_direction(state, times_s, seed_log, nu_obs_log, z, direction, sightline, sky_x, sky_y, d_ang_cm, omega, flux, x_moment, y_moment, samples)

    x_rad = np.zeros_like(flux)
    y_rad = np.zeros_like(flux)
    valid = flux > 0.0
    x_rad[valid] = x_moment[valid] / flux[valid]
    y_rad[valid] = y_moment[valid] / flux[valid]
    moments = MomentSeries(times_s, flux, x_rad, y_rad, np.sqrt(x_rad * x_rad + y_rad * y_rad))
    return moments, samples or []


def accumulate_1d_direction(state, times_s, seed_log, nu_obs_log, z, direction, sightline, sky_x, sky_y, d_ang_cm, omega, flux, x_moment, y_moment, samples) -> None:
    radius = np.asarray(state.components.fwd.radius_cm, dtype=float)
    gamma = np.asarray(state.components.fwd.gamma, dtype=float)
    arrival = np.asarray(state.components.fwd.characteristic_time_s, dtype=float) + radius * (1.0 - float(direction @ sightline)) * (1.0 + z) / constants.para_c
    source = np.asarray(state.components.fwd_sync, dtype=float)
    accumulate_cell(
        times_s=times_s,
        seed_log=seed_log,
        nu_obs_log=nu_obs_log,
        z=z,
        direction=direction,
        sightline=sightline,
        sky_x=sky_x,
        sky_y=sky_y,
        d_ang_cm=d_ang_cm,
        domega=omega,
        arrival=arrival,
        radius=radius,
        gamma=gamma,
        source=source,
        flux=flux,
        x_moment=x_moment,
        y_moment=y_moment,
        samples=samples,
    )


def accumulate_2d_direction(state, times_s, seed_log, nu_obs_log, z, direction, sightline, sky_x, sky_y, d_ang_cm, omega, flux, x_moment, y_moment, samples) -> None:
    radius_front = np.asarray(state.components.fwd.radius_cm, dtype=float)
    r_tobs = np.asarray(state.components.fwd.characteristic_time_s, dtype=float)
    radius_chi = np.asarray(state.electron.chi_radius_cm, dtype=float)
    gamma_chi = np.asarray(state.electron.chi_gamma_bulk, dtype=float)
    source_chi = np.asarray(state.electron.l_syn_spec_chi, dtype=float) * np.asarray(state.observer.prefactor, dtype=float)[:, None, :]
    tau_chi = np.asarray(state.electron.tau_syn_chi, dtype=float)
    chi_weight = np.asarray(state.electron.chi_dvolume_weight, dtype=float)
    tau_prefix = np.concatenate((np.zeros((tau_chi.shape[0], 1, tau_chi.shape[2]), dtype=float), np.cumsum(tau_chi, axis=1)), axis=1)
    mu = float(direction @ sightline)
    for i_chi in range(radius_chi.shape[0]):
        escape = chi_escape(tau_prefix[:, i_chi, :], tau_chi[:, i_chi, :])
        source = source_chi[:, i_chi, :] * chi_weight[i_chi, :][None, :] * escape
        arrival = r_tobs + (1.0 + z) * (radius_front - radius_chi[i_chi, :] * mu) / constants.para_c
        accumulate_cell(
            times_s=times_s,
            seed_log=seed_log,
            nu_obs_log=nu_obs_log,
            z=z,
            direction=direction,
            sightline=sightline,
            sky_x=sky_x,
            sky_y=sky_y,
            d_ang_cm=d_ang_cm,
            domega=omega,
            arrival=arrival,
            radius=radius_chi[i_chi, :],
            gamma=gamma_chi[i_chi, :],
            source=source,
            flux=flux,
            x_moment=x_moment,
            y_moment=y_moment,
            samples=samples,
        )


def render_samples(samples: list[tuple[int, float, float, float]], times_s: np.ndarray, half_extent_mas: float) -> MapProduct:
    edges = np.linspace(-half_extent_mas, half_extent_mas, SKYMAP_NPIX + 1)
    pixel_area = (edges[1] - edges[0]) ** 2
    image = np.zeros((times_s.size, SKYMAP_NPIX, SKYMAP_NPIX), dtype=float)
    flux = np.zeros(times_s.size, dtype=float)
    x_moment = np.zeros(times_s.size, dtype=float)
    y_moment = np.zeros(times_s.size, dtype=float)
    sample_array = np.asarray(samples, dtype=float)
    time_i = sample_array[:, 0].astype(int)
    x_mas = sample_array[:, 1]
    y_mas = sample_array[:, 2]
    value = sample_array[:, 3]
    x_i = np.searchsorted(edges, x_mas, side="right") - 1
    y_i = np.searchsorted(edges, y_mas, side="right") - 1
    inside = (x_i >= 0) & (x_i < SKYMAP_NPIX) & (y_i >= 0) & (y_i < SKYMAP_NPIX)
    np.add.at(image, (time_i[inside], x_i[inside], y_i[inside]), value[inside] / pixel_area)
    np.add.at(flux, time_i, value)
    np.add.at(x_moment, time_i, value * x_mas / RAD_TO_MAS)
    np.add.at(y_moment, time_i, value * y_mas / RAD_TO_MAS)
    x_centroid = np.zeros(times_s.size, dtype=float)
    y_centroid = np.zeros(times_s.size, dtype=float)
    valid = flux > 0.0
    x_centroid[valid] = x_moment[valid] / flux[valid]
    y_centroid[valid] = y_moment[valid] / flux[valid]
    centers = 0.5 * (edges[:-1] + edges[1:])
    return MapProduct(image, centers, centers, flux, x_centroid, y_centroid)


def apparent_beta(time_s: np.ndarray, x_rad: np.ndarray, y_rad: np.ndarray, d_ang_cm: float, z: float) -> tuple[np.ndarray, np.ndarray]:
    dt = np.diff(time_s)
    dx = np.diff(x_rad) * d_ang_cm
    dy = np.diff(y_rad) * d_ang_cm
    beta = (1.0 + z) * np.sqrt(dx * dx + dy * dy) / dt / constants.para_c
    return np.sqrt(time_s[:-1] * time_s[1:]), beta


def compute_products() -> tuple[dict[str, dict[float, MomentSeries]], dict[str, MapProduct], dict[str, np.ndarray]]:
    model_1d, state_1d = solve_state_for_solver("fullhide_1d")
    model_2d, state_2d = solve_state_for_solver("fullhide_2d")
    states = {"1d": state_1d, "2d": state_2d}
    d_ang_cm = float(model_2d.observer.lumi_dist_cm) / (1.0 + float(model_2d.observer.z)) ** 2
    z = float(model_2d.observer.z)

    moments: dict[str, dict[float, MomentSeries]] = {"1d": {}, "2d": {}}
    for label, state in states.items():
        for theta_mult in THETA_MULTIPLES:
            theta_v = float(theta_mult * THETA_J_RAD)
            print(f"centroid {label} theta/theta_j={theta_mult:.0f}", flush=True)
            moments[label][float(theta_mult)], _ = state_moments_and_samples(state, theta_v, CENTROID_TIMES_S)

    sample_outputs = {}
    all_samples = {}
    theta_v_map = SKYMAP_THETA_MULTIPLE * THETA_J_RAD
    for label, state in states.items():
        print(f"skymap {label} theta/theta_j={SKYMAP_THETA_MULTIPLE:.0f}", flush=True)
        _map_moments, samples = state_moments_and_samples(state, theta_v_map, SKYMAP_TIMES_S, collect_samples=True)
        all_samples[label] = samples

    coord_max = 0.0
    for samples in all_samples.values():
        for _time_i, x_mas, y_mas, value in samples:
            if value > 0.0:
                coord_max = max(coord_max, abs(x_mas), abs(y_mas))
    half_extent = 1.08 * coord_max
    for label, samples in all_samples.items():
        sample_outputs[label] = render_samples(samples, SKYMAP_TIMES_S, half_extent)

    meta = {
        "d_ang_cm": np.array([d_ang_cm], dtype=float),
        "z": np.array([z], dtype=float),
        "theta_j_rad": np.array([THETA_J_RAD], dtype=float),
        "nu_obs_hz": np.array([NU_OBS_HZ], dtype=float),
        "epsilon_B": np.array([EPSILON_B], dtype=float),
        "epsilon_b_floor": np.array([EPSILON_B_FLOOR], dtype=float),
        "magnetic_decay_alpha_t": np.array([MAGNETIC_DECAY_ALPHA_T], dtype=float),
        "magnetic_decay_t0_s": np.array([MAGNETIC_DECAY_T0_S], dtype=float),
    }
    return moments, sample_outputs, meta


def write_tables(moments: dict[str, dict[float, MomentSeries]], meta: dict[str, np.ndarray], output_dir: Path) -> None:
    rows = []
    d_ang_cm = float(meta["d_ang_cm"][0])
    z = float(meta["z"][0])
    for label, by_angle in moments.items():
        for theta_mult, series in by_angle.items():
            beta_t, beta = apparent_beta(series.time_s, series.x_rad, series.y_rad, d_ang_cm, z)
            beta_lookup = {float(t): float(b) for t, b in zip(beta_t, beta)}
            for i_time, time_s in enumerate(series.time_s):
                rows.append(
                    {
                        "model": label,
                        "theta_v_over_theta_j": theta_mult,
                        "time_s": float(time_s),
                        "flux_fnu": float(series.flux[i_time]),
                        "x_centroid_mas": float(series.x_rad[i_time] * RAD_TO_MAS),
                        "y_centroid_mas": float(series.y_rad[i_time] * RAD_TO_MAS),
                        "centroid_offset_mas": float(series.offset_rad[i_time] * RAD_TO_MAS),
                        "beta_app_at_midpoint": "",
                    }
                )
            for time_mid, beta_value in beta_lookup.items():
                rows.append(
                    {
                        "model": label,
                        "theta_v_over_theta_j": theta_mult,
                        "time_s": time_mid,
                        "flux_fnu": "",
                        "x_centroid_mas": "",
                        "y_centroid_mas": "",
                        "centroid_offset_mas": "",
                        "beta_app_at_midpoint": beta_value,
                    }
                )
    with (output_dir / "skymap_centroid_motion_summary.csv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_skymaps(maps: dict[str, MapProduct], output: Path) -> None:
    fig, axes = plt.subplots(2, SKYMAP_TIMES_S.size, figsize=(15.5, 7.2), sharex=True, sharey=True, layout="constrained")
    fig.suptitle(rf"Sky maps at $\nu={NU_OBS_HZ:.0e}$ Hz, $\theta_v/\theta_j={SKYMAP_THETA_MULTIPLE:.0f}$")
    extent = [
        maps["1d"].x_mas[0],
        maps["1d"].x_mas[-1],
        maps["1d"].y_mas[0],
        maps["1d"].y_mas[-1],
    ]
    for i_time, time_s in enumerate(SKYMAP_TIMES_S):
        column_values = np.concatenate((maps["1d"].image_flux[i_time].ravel(), maps["2d"].image_flux[i_time].ravel()))
        positive = column_values[column_values > 0.0]
        norm = LogNorm(vmin=float(np.min(positive)), vmax=float(np.max(positive)))
        for i_row, label in enumerate(("1d", "2d")):
            ax = axes[i_row, i_time]
            im = ax.imshow(
                maps[label].image_flux[i_time].T,
                origin="lower",
                extent=extent,
                norm=norm,
                cmap="magma",
                interpolation="nearest",
                aspect="equal",
            )
            ax.scatter(0.0, 0.0, marker="+", s=36, color="white", lw=0.9)
            ax.scatter(
                maps[label].x_centroid_rad[i_time] * RAD_TO_MAS,
                maps[label].y_centroid_rad[i_time] * RAD_TO_MAS,
                marker="x",
                s=34,
                color="#00ffff",
                lw=1.0,
            )
            if i_row == 0:
                ax.set_title(rf"$t={time_s / 86400.0:.1f}$ d")
            ax.tick_params(which="both", direction="in", top=True, right=True)
        cbar = fig.colorbar(im, ax=axes[:, i_time], fraction=0.045, pad=0.015, shrink=0.82)
        cbar.set_label(r"$F_\nu$ per mas$^2$")
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$x$ (mas)")
    axes[0, 0].set_ylabel("1D thin shell\n" + r"$y$ (mas)")
    axes[1, 0].set_ylabel("2D q-shell B decay\n" + r"$y$ (mas)")
    fig.savefig(output, dpi=180)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def plot_centroid_offset(moments: dict[str, dict[float, MomentSeries]], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    for i, theta_mult in enumerate(THETA_MULTIPLES):
        color = COLORS[i]
        for label, linestyle in (("1d", "--"), ("2d", "-")):
            series = moments[label][float(theta_mult)]
            ax.loglog(series.time_s / 86400.0, series.offset_rad * RAD_TO_MAS, color=color, ls=linestyle, lw=1.8)
    ax.set_xlabel(r"$t_{\rm obs}$ (d)")
    ax.set_ylabel("centroid offset (mas)")
    ax.set_title(rf"Flux-centroid offset at $\nu={NU_OBS_HZ:.0e}$ Hz")
    ax.grid(True, which="both", alpha=0.18)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    angle_handles = [
        Line2D([0], [0], color=COLORS[i], lw=2.0, label=rf"$\theta_v/\theta_j={theta_mult:.0f}$")
        for i, theta_mult in enumerate(THETA_MULTIPLES)
    ]
    model_handles = [
        Line2D([0], [0], color="black", lw=1.8, ls="--", label="1D"),
        Line2D([0], [0], color="black", lw=2.0, ls="-", label="2D B decay"),
    ]
    leg1 = ax.legend(handles=angle_handles, frameon=False, fontsize=8, loc="lower right")
    ax.add_artist(leg1)
    ax.legend(handles=model_handles, frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def plot_apparent_speed(moments: dict[str, dict[float, MomentSeries]], meta: dict[str, np.ndarray], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    d_ang_cm = float(meta["d_ang_cm"][0])
    z = float(meta["z"][0])
    for i, theta_mult in enumerate(THETA_MULTIPLES):
        color = COLORS[i]
        for label, linestyle in (("1d", "--"), ("2d", "-")):
            series = moments[label][float(theta_mult)]
            time_mid, beta = apparent_beta(series.time_s, series.x_rad, series.y_rad, d_ang_cm, z)
            ax.semilogx(time_mid / 86400.0, beta, color=color, ls=linestyle, lw=1.8)
    ax.axhline(1.0, color="0.5", lw=1.0, ls=":")
    ax.set_xlabel(r"$t_{\rm obs}$ (d)")
    ax.set_ylabel(r"$\beta_{\rm app}$")
    ax.set_title(r"Apparent centroid speed")
    ax.grid(True, which="both", alpha=0.18)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    angle_handles = [
        Line2D([0], [0], color=COLORS[i], lw=2.0, label=rf"$\theta_v/\theta_j={theta_mult:.0f}$")
        for i, theta_mult in enumerate(THETA_MULTIPLES)
    ]
    model_handles = [
        Line2D([0], [0], color="black", lw=1.8, ls="--", label="1D"),
        Line2D([0], [0], color="black", lw=2.0, ls="-", label="2D B decay"),
    ]
    leg1 = ax.legend(handles=angle_handles, frameon=False, fontsize=8, loc="upper right")
    ax.add_artist(leg1)
    ax.legend(handles=model_handles, frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def write_npz(moments: dict[str, dict[float, MomentSeries]], maps: dict[str, MapProduct], meta: dict[str, np.ndarray], output_dir: Path) -> None:
    payload = dict(meta)
    payload["centroid_times_s"] = CENTROID_TIMES_S
    payload["skymap_times_s"] = SKYMAP_TIMES_S
    payload["theta_multiples"] = THETA_MULTIPLES
    for label, by_angle in moments.items():
        for theta_mult, series in by_angle.items():
            key = f"{label}_theta{theta_mult:.0f}"
            payload[f"{key}_flux"] = series.flux
            payload[f"{key}_x_rad"] = series.x_rad
            payload[f"{key}_y_rad"] = series.y_rad
            payload[f"{key}_offset_rad"] = series.offset_rad
    for label, product in maps.items():
        payload[f"{label}_skymap_image_flux"] = product.image_flux
        payload[f"{label}_skymap_x_mas"] = product.x_mas
        payload[f"{label}_skymap_y_mas"] = product.y_mas
        payload[f"{label}_skymap_flux"] = product.flux
        payload[f"{label}_skymap_x_centroid_rad"] = product.x_centroid_rad
        payload[f"{label}_skymap_y_centroid_rad"] = product.y_centroid_rad
    np.savez(output_dir / "skymap_centroid_motion_data.npz", **payload)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("output/benchmark_1d_vs_qshell_skymap_motion_bdecay_alpha04"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    moments, maps, meta = compute_products()
    write_npz(moments, maps, meta, args.output_dir)
    write_tables(moments, meta, args.output_dir)
    plot_skymaps(maps, args.output_dir / "skymap_1d_vs_2d_bdecay.png")
    plot_centroid_offset(moments, args.output_dir / "centroid_offset_1d_vs_2d_bdecay.png")
    plot_apparent_speed(moments, meta, args.output_dir / "apparent_speed_1d_vs_2d_bdecay.png")
    print(args.output_dir)


if __name__ == "__main__":
    main()
