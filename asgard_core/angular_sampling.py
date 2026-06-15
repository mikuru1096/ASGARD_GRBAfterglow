from __future__ import annotations

from dataclasses import dataclass

import numpy as np


SUPPORTED_PATCH_SAMPLING = ("uniform", "dominant_region_ioka_v1", "dominant_region_ioka_time_v1")


@dataclass(frozen=True)
class PatchGrid:
    theta_centers: np.ndarray
    theta_edges: np.ndarray
    phi_centers: np.ndarray
    phi_edges: np.ndarray
    domega: np.ndarray
    patch_half_angle: np.ndarray


def build_patch_grid(model, observer_time_s: np.ndarray | None = None) -> PatchGrid:
    sampling = str(getattr(model.setups, "patch_sampling", "uniform")).lower()
    if sampling == "uniform":
        return _build_uniform_grid(model)
    if sampling == "dominant_region_ioka_v1":
        return _build_dominant_region_ioka_grid(model)
    if sampling == "dominant_region_ioka_time_v1":
        if observer_time_s is None:
            raise ValueError("dominant_region_ioka_time_v1 requires observer_time_s for true-Gamma weighting.")
        return _build_dominant_region_ioka_grid(model, observer_time_s=np.asarray(observer_time_s, dtype=float))
    raise ValueError(
        f"Unknown patch_sampling={sampling!r}; expected one of {SUPPORTED_PATCH_SAMPLING}."
    )


def _patch_counts(model) -> tuple[int, int]:
    n_theta = int(model.setups.patch_theta)
    n_phi = int(model.setups.patch_phi)
    if n_theta <= 0 or n_phi <= 0:
        raise ValueError("patch_theta and patch_phi must be positive for patch sampling.")
    return n_theta, n_phi


def _build_uniform_grid(model) -> PatchGrid:
    n_theta, n_phi = _patch_counts(model)
    theta_edges = np.linspace(0.0, float(model.jet.theta_max), n_theta + 1)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, n_phi + 1)
    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    domega = _solid_angle_cells(theta_edges, phi_edges)
    return PatchGrid(
        theta_centers=theta_centers,
        theta_edges=theta_edges,
        phi_centers=phi_centers,
        phi_edges=phi_edges,
        domega=domega,
        patch_half_angle=np.sqrt(np.maximum(domega, 1.0e-12) / np.pi),
    )


def _build_dominant_region_ioka_grid(
    model,
    observer_time_s: np.ndarray | None = None,
) -> PatchGrid:
    n_theta, n_phi = _patch_counts(model)
    theta_max = float(model.jet.theta_max)
    theta_obs = float(model.observer.theta_obs)
    phi_obs = float(model.observer.phi_obs)
    theta_scan = np.linspace(0.0, theta_max, max(4 * n_theta, 256))
    phi_scan = np.linspace(0.0, 2.0 * np.pi, max(4 * n_phi, 256), endpoint=False)
    if observer_time_s is None:
        scan_weight = _dominant_weight(model, theta_scan, phi_scan, theta_obs, phi_obs)
        gamma_time = None
    else:
        sampling_times = _sampling_times(model, observer_time_s)
        gamma_time = _pilot_gamma_theta_time(model, theta_scan, sampling_times)
        scan_weight = _time_dependent_dominant_weight(
            model,
            theta_scan,
            phi_scan,
            theta_obs,
            phi_obs,
            gamma_time,
        )
    theta_weight = np.sin(theta_scan) * np.sqrt(np.mean(scan_weight, axis=1))
    theta_count = n_theta
    if gamma_time is not None:
        gamma_envelope = np.max(gamma_time, axis=1)
        theta_weight = np.maximum(
            _normalized_density(theta_scan, theta_weight),
            _normalized_density(theta_scan, gamma_envelope),
        )
        theta_count = _beaming_resolved_theta_count(model, theta_scan, gamma_envelope, n_theta)
    theta_centers, theta_edges = _weighted_centers_edges(theta_scan, theta_weight, theta_count)

    if _axisymmetric_jet(model):
        if theta_obs == 0.0:
            phi_count = n_phi
        elif gamma_time is not None:
            phi_count = _beaming_resolved_phi_count(model, theta_scan, gamma_time, n_phi)
        else:
            phi_count = n_phi
        phi_centers, phi_edges, phi_weights = _axisymmetric_phi_quadrature(phi_count, phi_obs)
        domega = _solid_angle_cells_from_phi_weights(theta_edges, phi_weights)
    else:
        phi_edges = np.linspace(0.0, 2.0 * np.pi, n_phi + 1)
        phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
        domega = _solid_angle_cells(theta_edges, phi_edges)
    return PatchGrid(
        theta_centers=theta_centers,
        theta_edges=theta_edges,
        phi_centers=phi_centers,
        phi_edges=phi_edges,
        domega=domega,
        patch_half_angle=np.sqrt(np.maximum(domega, 1.0e-12) / np.pi),
    )


def _dominant_weight(
    model,
    theta: np.ndarray,
    phi: np.ndarray,
    theta_obs: float,
    phi_obs: float,
) -> np.ndarray:
    theta_mesh = theta[:, None]
    phi_mesh = phi[None, :]
    energy = _evaluate_jet_function(model.jet.energy_iso, phi_mesh, theta_mesh)
    gamma = _evaluate_jet_function(model.jet.gamma0, phi_mesh, theta_mesh)
    separation = _angular_separation(theta_mesh, phi_mesh, theta_obs, phi_obs)
    weight = energy * _doppler_factor(gamma, separation) ** 3
    return _with_structure_floor(weight)


def _axisymmetric_jet(model) -> bool:
    return str(getattr(model.jet, "kind", "")).lower() in {
        "gaussian",
        "powerlaw",
        "twocomponent",
        "steppowerlaw",
        "tophat",
    }


def _axisymmetric_phi_quadrature(n_phi: int, phi_obs: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if n_phi < 2:
        raise ValueError("axisymmetric half-plane phi quadrature requires patch_phi >= 2.")
    half_centers = np.linspace(0.0, np.pi, n_phi)
    phi_centers = np.mod(float(phi_obs) + half_centers, 2.0 * np.pi)
    phi_edges = np.linspace(float(phi_obs), float(phi_obs) + np.pi, n_phi + 1)
    step = np.pi / float(n_phi - 1)
    phi_weights = np.full(n_phi, 2.0 * step, dtype=float)
    phi_weights[0] = step
    phi_weights[-1] = step
    return phi_centers, phi_edges, phi_weights


def _time_dependent_dominant_weight(
    model,
    theta: np.ndarray,
    phi: np.ndarray,
    theta_obs: float,
    phi_obs: float,
    gamma_time: np.ndarray,
) -> np.ndarray:
    theta_mesh = theta[:, None]
    phi_mesh = phi[None, :]
    energy = _evaluate_jet_function(model.jet.energy_iso, phi_mesh, theta_mesh)
    separation = _angular_separation(theta_mesh, phi_mesh, theta_obs, phi_obs)
    weighted = np.zeros_like(energy, dtype=float)
    for i_time in range(gamma_time.shape[1]):
        weighted += _doppler_factor(gamma_time[:, i_time][:, None], separation) ** 3
    weighted = energy * weighted / float(gamma_time.shape[1])
    return _with_structure_floor(weighted)


def _beaming_resolved_phi_count(
    model,
    theta: np.ndarray,
    gamma_time: np.ndarray,
    minimum_count: int,
) -> int:
    factor = float(getattr(model.setups, "patch_sampling_beaming_factor", 3.0))
    resolution = float(getattr(model.setups, "patch_sampling_beaming_resolution", 8.0))
    if factor <= 0.0 or resolution <= 0.0:
        raise ValueError("patch_sampling_beaming_factor and patch_sampling_beaming_resolution must be positive.")
    angular_frequency = float(np.max(np.asarray(gamma_time, dtype=float) * np.sin(theta)[:, None]))
    required = int(np.ceil(np.pi * resolution * angular_frequency / factor)) + 1
    return max(int(minimum_count), required, 2)


def _beaming_resolved_theta_count(
    model,
    theta: np.ndarray,
    gamma_envelope: np.ndarray,
    minimum_count: int,
) -> int:
    factor = float(getattr(model.setups, "patch_sampling_beaming_factor", 3.0))
    resolution = float(getattr(model.setups, "patch_sampling_beaming_resolution", 8.0))
    if factor <= 0.0 or resolution <= 0.0:
        raise ValueError("patch_sampling_beaming_factor and patch_sampling_beaming_resolution must be positive.")
    angular_frequency = float(np.max(np.asarray(gamma_envelope, dtype=float)))
    required = int(np.ceil(float(model.jet.theta_max) * resolution * angular_frequency / factor)) + 1
    return max(int(minimum_count), required)


def _normalized_density(x: np.ndarray, density: np.ndarray) -> np.ndarray:
    positive = np.asarray(density, dtype=float)
    total = float(np.trapezoid(positive, np.asarray(x, dtype=float)))
    if total <= 0.0:
        raise ValueError("adaptive theta sampling requires positive density.")
    return positive / total


def _sampling_times(model, observer_time_s: np.ndarray) -> np.ndarray:
    observer_time_s = np.unique(np.asarray(observer_time_s, dtype=float))
    count = int(getattr(model.setups, "patch_sampling_num_times", 12))
    if observer_time_s.size <= count:
        return observer_time_s
    log_times = np.linspace(np.log(observer_time_s[0]), np.log(observer_time_s[-1]), count)
    return np.exp(log_times)


def _pilot_gamma_theta_time(model, theta: np.ndarray, observer_time_s: np.ndarray) -> np.ndarray:
    from concurrent.futures import ThreadPoolExecutor

    from ASGARD.api_model import _build_fit_config_for_patch
    from asgard_core.asgard_setup import build_simulation_setup
    from asgard_core.asgard_runtime import solve_dynamics
    from asgard_core.asgard_state import make_query_cfg

    observer_time_s = np.asarray(observer_time_s, dtype=float)
    sample_count = int(getattr(model.setups, "patch_sampling_pilot_theta", 0))
    if sample_count <= 0:
        sample_count = max(2 * int(model.setups.patch_theta), 48)
    pilot_theta = np.linspace(0.0, float(model.jet.theta_max), sample_count)
    def solve_pilot_theta(theta_center: float) -> np.ndarray:
        e_iso = model.jet.energy_iso(0.0, float(theta_center))
        gamma0 = model.jet.gamma0(0.0, float(theta_center))
        if e_iso <= 0.0 or gamma0 <= 1.0:
            return np.ones(observer_time_s.shape, dtype=float)
        config = _build_fit_config_for_patch(
            model,
            phi_center=0.0,
            theta_v=0.0,
            opening_angle_jet=float(model.jet.theta_max / sample_count),
            e_iso=float(e_iso),
            gamma0=float(gamma0),
            theta_center=float(theta_center),
        )
        query_config = make_query_cfg(config, observer_time_s)
        setup = build_simulation_setup(query_config)
        dynamics = solve_dynamics(setup.boundary, query_config)
        return np.interp(
            np.log(observer_time_s),
            np.log(np.asarray(dynamics.r_tobs, dtype=float)),
            np.asarray(dynamics.r_gamma, dtype=float),
        )

    worker_count = int(getattr(model.setups, "structured_outer_threads", 0) or getattr(model.setups, "num_threads", 1))
    if worker_count > 1:
        with ThreadPoolExecutor(max_workers=min(worker_count, pilot_theta.size)) as executor:
            pilot_rows = list(executor.map(solve_pilot_theta, (float(value) for value in pilot_theta)))
        pilot_gamma = np.asarray(pilot_rows, dtype=float)
    else:
        pilot_gamma = np.asarray([solve_pilot_theta(float(value)) for value in pilot_theta], dtype=float)
    gamma_time = np.empty((theta.size, observer_time_s.size), dtype=float)
    for i_time in range(observer_time_s.size):
        gamma_time[:, i_time] = np.interp(theta, pilot_theta, pilot_gamma[:, i_time])
    return gamma_time


def _evaluate_jet_function(function, phi: np.ndarray, theta: np.ndarray) -> np.ndarray:
    phi_b, theta_b = np.broadcast_arrays(phi, theta)
    try:
        values = np.asarray(function(phi_b, theta_b), dtype=float)
        if values.shape == phi_b.shape:
            return values
    except (TypeError, ValueError):
        pass
    values = np.empty(phi_b.shape, dtype=float)
    for index in np.ndindex(phi_b.shape):
        values[index] = function(float(phi_b[index]), float(theta_b[index]))
    return values


def _doppler_factor(gamma: np.ndarray, angle: np.ndarray) -> np.ndarray:
    gamma = np.asarray(gamma, dtype=float)
    beta = np.sqrt(np.maximum(gamma * gamma - 1.0, 0.0)) / gamma
    return 1.0 / (gamma * (1.0 - beta * np.cos(angle)))


def _angular_separation(theta: float, phi: np.ndarray, theta_obs: float, phi_obs: float) -> np.ndarray:
    cos_sep = (
        np.cos(theta) * np.cos(theta_obs)
        + np.sin(theta) * np.sin(theta_obs) * np.cos(phi - phi_obs)
    )
    return np.arccos(np.clip(cos_sep, -1.0, 1.0))


def _weighted_centers_edges(x: np.ndarray, weight: np.ndarray, count: int) -> tuple[np.ndarray, np.ndarray]:
    pdf = np.asarray(weight, dtype=float)
    increments = 0.5 * (pdf[:-1] + pdf[1:]) * np.diff(x)
    total = float(np.sum(increments))
    if total <= 0.0:
        raise ValueError("Weighted patch sampling requires positive total angular weight.")
    cdf = np.zeros_like(x, dtype=float)
    cdf[1:] = np.cumsum(increments) / total
    cdf[-1] = 1.0
    edges = np.interp(np.linspace(0.0, 1.0, count + 1), cdf, x)
    centers = np.interp((np.arange(count, dtype=float) + 0.5) / count, cdf, x)
    return centers, edges


def _with_structure_floor(weight: np.ndarray) -> np.ndarray:
    peak = float(np.max(weight))
    if peak <= 0.0:
        raise ValueError("dominant_region_ioka_v1 found no positive angular emission weight.")
    return weight + 0.01 * peak


def _weighted_periodic_centers_edges(
    phi: np.ndarray, weight: np.ndarray, count: int
) -> tuple[np.ndarray, np.ndarray]:
    phi_periodic = np.concatenate((phi, [2.0 * np.pi]))
    weight_periodic = np.concatenate((weight, [weight[0]]))
    centers, _ = _weighted_centers_edges(phi_periodic, weight_periodic, count)
    edges = np.empty(count + 1, dtype=float)
    edges[0] = 0.0
    edges[-1] = 2.0 * np.pi
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    return centers, edges


def _solid_angle_cells(theta_edges: np.ndarray, phi_edges: np.ndarray) -> np.ndarray:
    dcos = np.cos(theta_edges[:-1]) - np.cos(theta_edges[1:])
    dphi = np.diff(phi_edges)
    return dcos[:, None] * dphi[None, :]


def _solid_angle_cells_from_phi_weights(theta_edges: np.ndarray, phi_weights: np.ndarray) -> np.ndarray:
    dcos = np.cos(theta_edges[:-1]) - np.cos(theta_edges[1:])
    return dcos[:, None] * np.asarray(phi_weights, dtype=float)[None, :]
