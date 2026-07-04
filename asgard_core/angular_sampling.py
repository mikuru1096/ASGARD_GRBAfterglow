from __future__ import annotations

from dataclasses import dataclass

import numpy as np


SUPPORTED_PATCH_SAMPLING = ("uniform", "dominant_region_ioka_v1", "dominant_region_ioka_time_v1")
AXISYMMETRIC_JET_KINDS = frozenset(("tophat", "gaussian", "powerlaw", "twocomponent", "steppowerlaw", "tabulated"))


def is_axisymmetric_jet(jet) -> bool:
    return str(getattr(jet, "kind", "")).lower() in AXISYMMETRIC_JET_KINDS


@dataclass(frozen=True)
class PatchGrid:
    theta_centers: np.ndarray
    theta_edges: np.ndarray
    phi_centers: np.ndarray
    phi_edges: np.ndarray
    domega: np.ndarray
    patch_half_angle: np.ndarray


def build_patch_grid(
    model,
    observer_time_s: np.ndarray | None = None,
    theta_count: int | None = None,
    phi_count: int | None = None,
) -> PatchGrid:
    sampling = str(model.setups.patch_sampling).lower()
    if sampling == "uniform":
        return _build_uniform_grid(model, theta_count=theta_count, phi_count=phi_count)
    if sampling == "dominant_region_ioka_v1":
        return _build_dominant_grid(model, theta_count=theta_count, phi_count=phi_count)
    if sampling == "dominant_region_ioka_time_v1":
        if observer_time_s is None:
            raise ValueError("dominant_region_ioka_time_v1 requires observer_time_s.")
        return _build_dominant_grid(
            model,
            observer_time_s=np.asarray(observer_time_s, dtype=float),
            theta_count=theta_count,
            phi_count=phi_count,
        )
    raise ValueError(f"unsupported patch_sampling: {sampling}")


def _patch_counts(model, theta_count: int | None = None, phi_count: int | None = None) -> tuple[int, int]:
    n_theta = int(model.setups.structured_num_theta if theta_count is None else theta_count)
    n_phi = int(model.setups.structured_num_phi if phi_count is None else phi_count)
    if n_theta <= 0 or n_phi <= 0:
        raise ValueError("structured_num_theta and structured_num_phi must be positive for patch sampling.")
    return n_theta, n_phi


def _build_uniform_grid(model, theta_count: int | None = None, phi_count: int | None = None) -> PatchGrid:
    n_theta, n_phi = _patch_counts(model, theta_count=theta_count, phi_count=phi_count)
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
        patch_half_angle=np.sqrt(domega / np.pi),
    )


def _build_dominant_grid(
    model,
    observer_time_s: np.ndarray | None = None,
    theta_count: int | None = None,
    phi_count: int | None = None,
) -> PatchGrid:
    n_theta, n_phi = _patch_counts(model, theta_count=theta_count, phi_count=phi_count)
    theta_scan = np.linspace(0.0, float(model.jet.theta_max), max(4 * n_theta, 256))
    phi_scan = np.linspace(0.0, 2.0 * np.pi, max(4 * n_phi, 256), endpoint=False)
    if observer_time_s is None:
        scan_weight = _dominant_weight(model, theta_scan, phi_scan)
        gamma_time = None
    else:
        sample_time = _sampling_times(model, observer_time_s)
        gamma_time = _pilot_gamma_theta_time(model, theta_scan, sample_time, n_theta)
        scan_weight = _time_weight(model, theta_scan, phi_scan, gamma_time)
    theta_weight = np.sin(theta_scan) * np.sqrt(np.mean(scan_weight, axis=1))
    if gamma_time is not None:
        gamma_envelope = np.max(gamma_time, axis=1)
        theta_weight = np.maximum(
            _normalized_density(theta_scan, theta_weight),
            _normalized_density(theta_scan, gamma_envelope),
        )
        n_theta = _beaming_theta_count(model, theta_scan, gamma_envelope, n_theta)
    theta_centers, theta_edges = _weighted_centers_edges(theta_scan, theta_weight, n_theta)
    if gamma_time is not None:
        n_phi = _beaming_phi_count(model, theta_scan, gamma_time, n_phi)
    phi_centers, phi_edges, phi_weights = _axisymmetric_phi_quadrature(n_phi, model.observer.phi_obs)
    domega = _solid_angle_cells_from_phi_weights(theta_edges, phi_weights)
    return PatchGrid(
        theta_centers=theta_centers,
        theta_edges=theta_edges,
        phi_centers=phi_centers,
        phi_edges=phi_edges,
        domega=domega,
        patch_half_angle=np.sqrt(domega / np.pi),
    )


def _dominant_weight(model, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    theta_mesh, phi_mesh = theta[:, None], phi[None, :]
    energy = model.jet.energy_iso(phi_mesh, theta_mesh)
    gamma = model.jet.gamma0(phi_mesh, theta_mesh)
    sep = angular_separation(theta_mesh, phi_mesh, model.observer.theta_obs, model.observer.phi_obs)
    return _with_structure_floor(energy * _doppler_factor(gamma, sep) ** 3)


def _time_weight(model, theta: np.ndarray, phi: np.ndarray, gamma_time: np.ndarray) -> np.ndarray:
    theta_mesh, phi_mesh = theta[:, None], phi[None, :]
    energy = model.jet.energy_iso(phi_mesh, theta_mesh)
    sep = angular_separation(theta_mesh, phi_mesh, model.observer.theta_obs, model.observer.phi_obs)
    weighted = np.zeros_like(sep, dtype=float)
    for i_time in range(gamma_time.shape[1]):
        weighted += _doppler_factor(gamma_time[:, i_time][:, None], sep) ** 3
    return _with_structure_floor(energy * weighted / float(gamma_time.shape[1]))


def _sampling_times(model, observer_time_s: np.ndarray) -> np.ndarray:
    observer_time_s = np.unique(np.asarray(observer_time_s, dtype=float))
    count = int(model.setups.patch_sampling_num_times)
    if observer_time_s.size <= count:
        return observer_time_s
    return np.exp(np.linspace(np.log(observer_time_s[0]), np.log(observer_time_s[-1]), count))


def _pilot_gamma_theta_time(
    model,
    theta: np.ndarray,
    observer_time_s: np.ndarray,
    theta_count: int,
) -> np.ndarray:
    from asgard_core.api_model import _build_fit_config_for_patch
    from asgard_core.asgard_runtime import solve_dynamics
    from asgard_core.asgard_setup import build_setup
    from asgard_core.asgard_state import make_query_cfg

    sample_count = int(model.setups.patch_sampling_pilot_theta)
    if sample_count <= 0:
        sample_count = max(2 * theta_count, 48)
    pilot_theta = np.linspace(0.0, float(model.jet.theta_max), sample_count)
    pilot_gamma = []
    for theta_center in pilot_theta:
        config = _build_fit_config_for_patch(
            model,
            theta_v=0.0,
            opening_angle_jet=float(model.jet.theta_max / sample_count),
            e_iso=float(model.jet.energy_iso(0.0, float(theta_center))),
            gamma0=float(model.jet.gamma0(0.0, float(theta_center))),
            theta_center=float(theta_center),
        )
        query = make_query_cfg(config, observer_time_s)
        dynamics = solve_dynamics(build_setup(query).boundary, query)
        pilot_gamma.append(np.interp(np.log(observer_time_s), np.log(dynamics.r_tobs), dynamics.r_gamma))
    pilot_gamma = np.asarray(pilot_gamma, dtype=float)
    gamma_time = np.empty((theta.size, observer_time_s.size), dtype=float)
    for i_time in range(observer_time_s.size):
        gamma_time[:, i_time] = np.interp(theta, pilot_theta, pilot_gamma[:, i_time])
    return gamma_time


def _axisymmetric_phi_quadrature(n_phi: int, phi_obs: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if int(n_phi) < 2:
        raise ValueError("dominant-region patch sampling requires structured_num_phi >= 2.")
    centers = np.linspace(0.0, np.pi, int(n_phi))
    phi_centers = np.mod(float(phi_obs) + centers, 2.0 * np.pi)
    phi_edges = np.linspace(float(phi_obs), float(phi_obs) + np.pi, int(n_phi) + 1)
    step = np.pi / float(int(n_phi) - 1)
    weights = np.full(int(n_phi), 2.0 * step, dtype=float)
    weights[0] = step
    weights[-1] = step
    return phi_centers, phi_edges, weights


def _beaming_theta_count(model, theta: np.ndarray, gamma: np.ndarray, minimum: int) -> int:
    factor = float(model.setups.patch_sampling_beaming_factor)
    resolution = float(model.setups.patch_sampling_beaming_resolution)
    required = int(np.ceil(float(model.jet.theta_max) * resolution * float(np.max(gamma)) / factor)) + 1
    return max(int(minimum), required)


def _beaming_phi_count(model, theta: np.ndarray, gamma_time: np.ndarray, minimum: int) -> int:
    factor = float(model.setups.patch_sampling_beaming_factor)
    resolution = float(model.setups.patch_sampling_beaming_resolution)
    angular_frequency = float(np.max(gamma_time * np.sin(theta)[:, None]))
    required = int(np.ceil(np.pi * resolution * angular_frequency / factor)) + 1
    return max(int(minimum), required, 2)


def _normalized_density(x: np.ndarray, density: np.ndarray) -> np.ndarray:
    total = float(np.trapezoid(np.asarray(density, dtype=float), np.asarray(x, dtype=float)))
    if total <= 0.0:
        raise ValueError("adaptive theta sampling requires positive density.")
    return np.asarray(density, dtype=float) / total


def _weighted_centers_edges(x: np.ndarray, weight: np.ndarray, count: int) -> tuple[np.ndarray, np.ndarray]:
    pdf = np.asarray(weight, dtype=float)
    increments = 0.5 * (pdf[:-1] + pdf[1:]) * np.diff(x)
    total = float(np.sum(increments))
    if total <= 0.0:
        raise ValueError("weighted patch sampling requires positive total angular weight.")
    cdf = np.zeros_like(x, dtype=float)
    cdf[1:] = np.cumsum(increments) / total
    cdf[-1] = 1.0
    return (
        np.interp((np.arange(int(count), dtype=float) + 0.5) / float(count), cdf, x),
        np.interp(np.linspace(0.0, 1.0, int(count) + 1), cdf, x),
    )


def _with_structure_floor(weight: np.ndarray) -> np.ndarray:
    peak = float(np.max(weight))
    if peak <= 0.0:
        raise ValueError("dominant-region sampling found no positive angular weight.")
    return weight + 0.01 * peak


def _doppler_factor(gamma: np.ndarray, angle: np.ndarray) -> np.ndarray:
    gamma = np.asarray(gamma, dtype=float)
    beta = np.sqrt(gamma * gamma - 1.0) / gamma
    return 1.0 / (gamma * (1.0 - beta * np.cos(angle)))


def angular_separation(
    theta: float | np.ndarray,
    phi: float | np.ndarray,
    theta_obs: float,
    phi_obs: float,
) -> np.ndarray:
    cos_sep = (
        np.cos(theta) * np.cos(theta_obs)
        + np.sin(theta) * np.sin(theta_obs) * np.cos(phi - phi_obs)
    )
    return np.arccos(np.clip(cos_sep, -1.0, 1.0))


def _solid_angle_cells(theta_edges: np.ndarray, phi_edges: np.ndarray) -> np.ndarray:
    return (np.cos(theta_edges[:-1]) - np.cos(theta_edges[1:]))[:, None] * np.diff(phi_edges)[None, :]


def _solid_angle_cells_from_phi_weights(theta_edges: np.ndarray, phi_weights: np.ndarray) -> np.ndarray:
    return (
        (np.cos(theta_edges[:-1]) - np.cos(theta_edges[1:]))[:, None]
        * np.asarray(phi_weights, dtype=float)[None, :]
    )
