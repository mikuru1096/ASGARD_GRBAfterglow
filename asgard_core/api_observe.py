from __future__ import annotations

from typing import Any

import numpy as np

from .api_fit import Param
from .api_model import (
    Hadronic,
    Model,
    Numerics,
    Observer,
    ObserverGrid,
    PolarizationResult,
    Radiation,
    ReverseShock,
    SkyImage,
    SolverOptions,
)
from src import constants


def _polarization(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    *,
    magnetic_geometry: str,
    local_emissivity: str,
) -> PolarizationResult:
    times_s = np.asarray(times_s, dtype=float)
    nu_hz = np.asarray(nu_hz, dtype=float)
    if times_s.ndim != 1 or nu_hz.ndim != 1:
        raise ValueError("polarization() requires 1D times_s and nu_hz grids.")
    if times_s.size == 0 or nu_hz.size == 0:
        raise ValueError("polarization() grids must be non-empty.")
    if np.any(times_s <= 0.0) or np.any(nu_hz <= 0.0):
        raise ValueError("polarization() times and frequencies must be positive.")
    if magnetic_geometry not in {"shock_random", "toroidal"}:
        raise ValueError("magnetic_geometry must be 'shock_random' or 'toroidal'.")
    if local_emissivity not in {"analytic", "analytic_then_kernel"}:
        raise ValueError("local_emissivity must be 'analytic' or 'analytic_then_kernel'.")
    if magnetic_geometry == "toroidal" and model.jet.kind != "tophat":
        raise NotImplementedError("toroidal polarization currently requires an axisymmetric top-hat jet.")
    state = _observe_state(model, times_s, nu_hz)
    theta, phi, omega = _anglegrid(state)
    sight, sx, sy = _skybasis(float(state.config.theta_v), float(model.observer.phi_obs))
    dirs = _directions(theta, phi)
    comp: dict[str, dict[str, np.ndarray]] = {}
    total_i = np.zeros((nu_hz.size, times_s.size), dtype=float)
    total_q = np.zeros_like(total_i)
    total_u = np.zeros_like(total_i)

    for name, branch, source, pval in _polsources(model, state):
        stokes = _stokesgrid(
            state,
            times_s,
            nu_hz,
            branch,
            source,
            dirs,
            omega,
            sight,
            sx,
            sy,
            magnetic_geometry,
            float(pval),
        )
        comp[name] = {"I": stokes[0], "Q": stokes[1], "U": stokes[2]}
        total_i = total_i + stokes[0]
        total_q = total_q + stokes[1]
        total_u = total_u + stokes[2]

    pol = np.zeros_like(total_i)
    ang = np.zeros_like(total_i)
    active = total_i > 0.0
    pol[active] = np.sqrt(total_q[active] * total_q[active] + total_u[active] * total_u[active]) / total_i[active]
    ang[active] = 0.5 * np.arctan2(total_u[active], total_q[active])
    return PolarizationResult(
        I_sync=total_i,
        Q=total_q,
        U=total_u,
        linear_polarization=pol,
        polarization_angle_rad=ang,
        components=comp,
    )


def _skyimage(model: Model, times_s: np.ndarray, nu_obs: float, fov: float, npixel: int):
    if model.observer.lumi_dist_cm is None or model.observer.lumi_dist_cm <= 0.0:
        raise ValueError("Observer.luminosity_distance_cm must be set for sky_image().")
    state = _observe_state(model, times_s, np.array([nu_obs], dtype=float))
    theta, phi, omega = _anglegrid(state)
    sight, sx, sy = _skybasis(float(state.config.theta_v), float(model.observer.phi_obs))
    dirs = _directions(theta, phi)
    edges = np.linspace(-0.5 * fov, 0.5 * fov, npixel + 1)
    pix = edges[1] - edges[0]
    image = np.zeros((times_s.size, npixel, npixel), dtype=float)
    flux = np.zeros(times_s.size, dtype=float)
    xmom = np.zeros_like(flux)
    ymom = np.zeros_like(flux)
    _accumimage(
        state,
        times_s,
        float(nu_obs),
        state.components.fwd,
        np.asarray(state.components.total, dtype=float),
        dirs,
        omega,
        sight,
        sx,
        sy,
        edges,
        image,
        flux,
        xmom,
        ymom,
    )
    centroid_x = np.zeros_like(flux)
    centroid_y = np.zeros_like(flux)
    active = flux > 0.0
    centroid_x[active] = xmom[active] / flux[active]
    centroid_y[active] = ymom[active] / flux[active]
    direct = _directflux(state, times_s, float(nu_obs))
    rendered = image.sum(axis=(1, 2)) * pix * pix
    return SkyImage(
        image=image,
        extent=np.array([edges[0], edges[-1], edges[0], edges[-1]], dtype=float),
        pixel_solid_angle=float(pix * pix),
        pixel_size=float(pix),
        direct_flux=direct,
        rendered_flux=rendered,
        normalization_scale=np.ones_like(rendered),
        x_centroid=centroid_x,
        y_centroid=centroid_y,
    )


def _observe_state(model: Model, times_s: np.ndarray, nu_hz: np.ndarray):
    if model.jet.kind != "tophat" or not model._supports_direct_kernel():
        raise NotImplementedError("sky-image and polarization products currently support the direct top-hat backend.")
    from .api_model import _currentconfig, _direct_tophat_patch_config, _solve_patch_state

    ref = np.unique(np.concatenate((np.asarray(times_s, dtype=float), model.default_detail_times())))
    config = _direct_tophat_patch_config(model, baseconfig=_currentconfig(model))
    return _solve_patch_state(
        model,
        config,
        ref,
        np.asarray(nu_hz, dtype=float),
        solve_reference_times_s=ref,
    )


def _directflux(state, times_s: np.ndarray, nu_obs: float) -> np.ndarray:
    from .api_adaptive import _observe_parts

    observed = _observe_parts(state, times_s, np.array([nu_obs], dtype=float), projection_kind="lightcurve")
    return np.asarray(observed.total[0, :], dtype=float)


def _anglegrid(state) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta_n = int(state.config.eats_num_theta)
    phi_n = int(state.config.eats_num_phi)
    theta_e = np.linspace(0.0, float(state.config.opening_angle_jet), theta_n + 1)
    phi_e = np.linspace(0.0, 2.0 * np.pi, 2 * phi_n + 1)
    theta = 0.5 * (theta_e[:-1] + theta_e[1:])
    phi = 0.5 * (phi_e[:-1] + phi_e[1:])
    omega = (np.cos(theta_e[:-1]) - np.cos(theta_e[1:]))[:, None] * (phi_e[1] - phi_e[0])
    return theta, phi, omega


def _directions(theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    tt, pp = np.meshgrid(theta, phi, indexing="ij")
    return np.stack((np.sin(tt) * np.cos(pp), np.sin(tt) * np.sin(pp), np.cos(tt)), axis=-1)


def _skybasis(theta_v: float, phi_v: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sight = np.array(
        [np.sin(theta_v) * np.cos(phi_v), np.sin(theta_v) * np.sin(phi_v), np.cos(theta_v)],
        dtype=float,
    )
    trial = np.array([0.0, 0.0, 1.0], dtype=float)
    sx = np.cross(trial, sight)
    if np.linalg.norm(sx) == 0.0:
        sx = np.array([1.0, 0.0, 0.0], dtype=float)
    else:
        sx = sx / np.linalg.norm(sx)
    sy = np.cross(sight, sx)
    return sight, sx, sy / np.linalg.norm(sy)


def _interp(seed_log: np.ndarray, source: np.ndarray, shell: np.ndarray, ratio: np.ndarray, target_log: np.ndarray) -> np.ndarray:
    values = np.zeros(target_log.size, dtype=float)
    inside = (target_log > seed_log[0]) & (target_log <= seed_log[-1])
    idx = np.flatnonzero(inside)
    if idx.size == 0:
        return values
    hi = np.searchsorted(seed_log, target_log[idx], side="left")
    lo = hi - 1
    freq = (target_log[idx] - seed_log[lo]) / (seed_log[hi] - seed_log[lo])
    col = shell[idx]
    rr = ratio[idx]
    ylo = (1.0 - rr) * source[lo, col] + rr * source[lo, col + 1]
    yhi = (1.0 - rr) * source[hi, col] + rr * source[hi, col + 1]
    pos = (ylo > 0.0) & (yhi > 0.0)
    local = np.empty(idx.size, dtype=float)
    local[pos] = np.exp(np.log(ylo[pos]) + freq[pos] * (np.log(yhi[pos]) - np.log(ylo[pos])))
    local[~pos] = (1.0 - freq[~pos]) * ylo[~pos] + freq[~pos] * yhi[~pos]
    values[idx] = local
    return values


def _cellflux(
    state,
    times_s: np.ndarray,
    nu_hz: float,
    branch,
    source: np.ndarray,
    direction: np.ndarray,
    omega: float,
    sight: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    radius = np.asarray(branch.radius_cm, dtype=float)
    gamma = np.asarray(branch.gamma, dtype=float)
    arrival = np.asarray(branch.characteristic_time_s, dtype=float) + radius * (1.0 - float(direction @ sight)) * (1.0 + float(state.config.z)) / constants.para_c
    shell = np.searchsorted(arrival, times_s, side="right") - 1
    valid = (shell >= 0) & (shell < arrival.size - 1)
    tidx = np.flatnonzero(valid)
    shell = shell[valid]
    ratio = (times_s[valid] - arrival[shell]) / (arrival[shell + 1] - arrival[shell])
    gam = np.exp(np.log(gamma[shell]) + ratio * (np.log(gamma[shell + 1]) - np.log(gamma[shell])))
    rad = (1.0 - ratio) * radius[shell] + ratio * radius[shell + 1]
    beta = np.sqrt(1.0 - gam**-2)
    mu = float(direction @ sight)
    dinv = gam * (1.0 - beta * mu)
    weight = omega / (4.0 * np.pi) * dinv**-3
    target = np.log(nu_hz) + np.log(dinv) + np.log1p(float(state.config.z))
    vals = weight * _interp(np.log(np.asarray(state.setup.seed_frequency_hz, dtype=float)), source, shell, ratio, target)
    return tidx, vals, rad

def _skypos(radius: np.ndarray, direction: np.ndarray, sight: np.ndarray, sx: np.ndarray, sy: np.ndarray, d_ang: float) -> tuple[np.ndarray, np.ndarray]:
    pos = radius[:, None] * direction[None, :]
    los = np.sum(pos * sight[None, :], axis=1)
    trans = pos - los[:, None] * sight[None, :]
    return trans @ sx / d_ang, trans @ sy / d_ang


def _accumimage(state, times_s, nu_obs, branch, source, dirs, omega, sight, sx, sy, edges, image, flux, xmom, ymom) -> None:
    d_ang = float(state.setup.luminosity_distance_cm) / (1.0 + float(state.config.z)) ** 2
    for i in range(dirs.shape[0]):
        for j in range(dirs.shape[1]):
            direction = dirs[i, j]
            tidx, vals, rad = _cellflux(state, times_s, nu_obs, branch, source, direction, float(omega[i, 0]), sight)
            x, y = _skypos(rad, direction, sight, sx, sy, d_ang)
            xi = np.searchsorted(edges, x, side="right") - 1
            yi = np.searchsorted(edges, y, side="right") - 1
            inside = (xi >= 0) & (xi < image.shape[1]) & (yi >= 0) & (yi < image.shape[2])
            pixarea = (edges[1] - edges[0]) ** 2
            np.add.at(image, (tidx[inside], xi[inside], yi[inside]), vals[inside] / pixarea)
            np.add.at(flux, tidx, vals)
            np.add.at(xmom, tidx, vals * x)
            np.add.at(ymom, tidx, vals * y)


def _polsources(model: Model, state) -> list[tuple[str, Any, np.ndarray, float]]:
    out = [("fwd_sync", state.components.fwd, np.asarray(state.components.fwd_sync, dtype=float), float(model.fwd_rad.p))]
    if state.components.rev is not None and state.components.rev_sync is not None:
        prad = model.rvs_rad if model.rvs_rad is not None else model.fwd_rad
        out.append(("rev_sync", state.components.rev, np.asarray(state.components.rev_sync, dtype=float), float(prad.p)))
    return out


def _evpa(direction: np.ndarray, sight: np.ndarray, sx: np.ndarray, sy: np.ndarray, geom: str) -> float:
    if geom == "toroidal":
        axis = np.array([0.0, 0.0, 1.0], dtype=float)
        bvec = np.cross(axis, direction)
    else:
        bvec = direction - float(direction @ sight) * sight
    proj = bvec - float(bvec @ sight) * sight
    bx = float(proj @ sx)
    by = float(proj @ sy)
    return np.arctan2(by, bx) + 0.5 * np.pi


def _stokesgrid(state, times_s, nu_hz, branch, source, dirs, omega, sight, sx, sy, geom, pval):
    pi0 = (pval + 1.0) / (pval + 7.0 / 3.0)
    igrid = np.zeros((nu_hz.size, times_s.size), dtype=float)
    qgrid = np.zeros_like(igrid)
    ugrid = np.zeros_like(igrid)
    for inu, nu in enumerate(nu_hz):
        for i in range(dirs.shape[0]):
            for j in range(dirs.shape[1]):
                direction = dirs[i, j]
                tidx, vals, _rad = _cellflux(state, times_s, float(nu), branch, source, direction, float(omega[i, 0]), sight)
                chi = _evpa(direction, sight, sx, sy, geom)
                np.add.at(igrid[inu], tidx, vals)
                np.add.at(qgrid[inu], tidx, pi0 * vals * np.cos(2.0 * chi))
                np.add.at(ugrid[inu], tidx, pi0 * vals * np.sin(2.0 * chi))
    return igrid, qgrid, ugrid


def _parampath(model: Model, param: Param) -> str:
    if param.path is not None:
        return param.path
    name = param.name
    path_map = {
        "energy_iso_erg": "jet.energy_iso_erg",
        "log10_energy_iso_erg": "jet.energy_iso_erg",
        "initial_lorentz_factor": "jet.initial_lorentz_factor",
        "log10_initial_lorentz_factor": "jet.initial_lorentz_factor",
        "core_angle_rad": "jet.opening_angle_rad" if model.jet.kind == "tophat" else "jet.core_angle_rad",
        "opening_angle_rad": "jet.opening_angle_rad",
        "outer_angle_rad": "jet.outer_angle_rad",
        "energy_index": "jet.energy_index",
        "lorentz_index": "jet.lorentz_index",
        "shell_duration_s": "jet.shell_duration_s",
        "viewing_angle_rad": "observer.viewing_angle_rad",
        "viewing_azimuth_rad": "observer.viewing_azimuth_rad",
        "z": "observer.z",
        "luminosity_distance_cm": "observer.luminosity_distance_cm",
        "epsilon_e": "fwd_rad.epsilon_e",
        "epsilon_B": "fwd_rad.epsilon_B",
        "p": "fwd_rad.p",
        "proton_energy_fraction": "fwd_rad.proton_energy_fraction",
        "epsilon_b_floor": "fwd_rad.epsilon_b_floor",
        "magnetic_decay_alpha_t": "fwd_rad.magnetic_decay_alpha_t",
        "magnetic_decay_t0_s": "fwd_rad.magnetic_decay_t0_s",
        "accelerated_electron_fraction": "fwd_rad.accelerated_electron_fraction",
        "acceleration_efficiency": "fwd_rad.acceleration_efficiency",
        "reverse_proton_energy_fraction": "fwd_rad.reverse_proton_energy_fraction",
        "number_density_cm3": "medium.number_density_cm3",
        "a_star": "medium.A_star",
        "density_floor_cm3": "medium.density_floor_cm3",
        "density_cap_cm3": "medium.density_cap_cm3",
        "reverse_shock_enabled": "reverse_shock.enabled",
        "reverse_shock_shell_duration_s": "reverse_shock.shell_duration_s",
        "reverse_shock_upstream_sigma": "reverse_shock.upstream_sigma",
        "hadronic_num_proton_gamma": "hadronic.num_proton_gamma",
        "hadronic_num_neutrino_frequency": "hadronic.num_neutrino_frequency",
        "hadronic_pair_cascade_iterations": "hadronic.pair_cascade_iterations",
    }
    if name not in path_map:
        raise KeyError(f"Cannot infer canonical parameter path for {param.name}; pass Param(name, path, lower, upper, scale).")
    return path_map[name]


def _asmodel(cfg: Any) -> Model:
    if isinstance(cfg, Model):
        return cfg
    if cfg is None:
        raise ValueError("Either a Model or cfg must be provided.")
    if not isinstance(cfg, dict):
        raise TypeError("cfg must be a Model or a dictionary containing explicit public objects.")

    def build(key: str, cls):
        value = cfg[key]
        return value if isinstance(value, cls) else cls(**value)

    return Model(
        medium=cfg["medium"],
        jet=cfg["jet"],
        observer=build("observer", Observer),
        fwd_rad=build("fwd_rad", Radiation),
        rvs_rad=None if cfg.get("rvs_rad") is None else build("rvs_rad", Radiation),
        numerics=build("numerics", Numerics),
        observer_grid=build("observer_grid", ObserverGrid),
        solver_options=build("solver_options", SolverOptions),
        reverse_shock=build("reverse_shock", ReverseShock),
        hadronic=build("hadronic", Hadronic),
    )
