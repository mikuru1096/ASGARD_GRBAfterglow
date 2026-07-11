from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from asgard_core import (
    Model,
    Numerics,
    Observer,
    ObserverGrid,
    Radiation,
    ReverseShock,
    SolverOptions,
    UniformMedium,
    gaussian_jet,
)
from asgard_core.api_model import _build_fit_config_for_patch, _solve_patch_state
from asgard_core.asgard_state import project_flux


@dataclass(frozen=True)
class PatchGrid:
    theta_centers: np.ndarray
    theta_edges: np.ndarray
    phi_centers: np.ndarray
    phi_edges: np.ndarray
    domega: np.ndarray


def make_gaussian_model(
    theta_obs: float,
    *,
    gamma0: float = 120.0,
    theta_c: float = 0.08,
    theta_max: float = 0.24,
    beaming_factor: float = 3.0,
    beaming_resolution: float = 8.0,
    pilot_theta: int = 0,
    sampling_num_times: int = 12,
) -> Model:
    return Model(
        jet=gaussian_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=gamma0,
            core_angle_rad=theta_c,
            outer_angle_rad=theta_max,
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        medium=UniformMedium(number_density_cm3=1.0),
        observer=Observer(
            z=0.1,
            viewing_angle_rad=theta_obs,
            viewing_azimuth_rad=0.0,
            luminosity_distance_cm=1.0e26,
        ),
        fwd_rad=Radiation(
            epsilon_e=0.1,
            epsilon_B=1.0e-3,
            p=2.3,
            proton_energy_fraction=0.0,
            epsilon_b_floor=None,
            magnetic_decay_alpha_t=0.0,
            magnetic_decay_t0_s=1.0,
            accelerated_electron_fraction=0.1,
            thermal_electrons=False,
            include_ssc=False,
            proton_synch=True,
            include_pgamma=False,
            bethe_heitler=False,
            hadronic_inverse_compton=False,
            pp=False,
            neutrino=False,
            acceleration_efficiency=1.0,
            reverse_proton_energy_fraction=0.0,
            pgamma_scheme="disabled",
            pair_production=False,
        ),
        rvs_rad=None,
        numerics=Numerics(
            structured_num_theta=12,
            structured_num_phi=24,
            num_radius=16,
            eats_num_theta=12,
            eats_num_phi=12,
            downstream_num_chi=None,
            num_observer_time=12,
            num_electron_gamma=21,
            num_photon_frequency=17,
            num_threads=1,
            electron_adaptive_substeps=False,
            electron_substep_rtol=0.02,
            electron_substep_min=100,
            electron_substep_max=1000,
            initial_radius_cm=1.0e14,
        ),
        observer_grid=ObserverGrid(time_min_s=1.0e3, time_max_s=1.0e6),
        solver_options=SolverOptions(
            electron_solver="fullhide_1d",
            dynamics_solver="forward_legacy",
            geometry_projection="sed_legacy",
            electron_photon_coupling="separated",
            ssc_cooling_mode="none",
            synchrotron_integration="fixed_grid",
            cooling_kernel="legacy",
            structured_backend="fortran_1d",
            patch_sampling="uniform",
            patch_sampling_pilot_theta=pilot_theta,
            patch_sampling_num_times=sampling_num_times,
            patch_sampling_beaming_factor=beaming_factor,
            patch_sampling_beaming_resolution=beaming_resolution,
            structured_parallel_mode="outer",
            structured_outer_threads=None,
            structured_inner_threads=None,
            fullhide2d_transport_model="legacy",
            fullhide2d_stochastic_accel_norm=0.0,
            fullhide2d_escape_mode="closed",
        ),
        reverse_shock=ReverseShock(
            enabled=False,
            shell_duration_s=10.0,
            upstream_sigma=0.0,
            include_cross_zone_ic=False,
            include_ssc=False,
        ),
        hadronic=_disabled_hadronic(),
    )


def _disabled_hadronic():
    from asgard_core import Hadronic

    return Hadronic(
        enabled=False,
        solver="legacy_1d",
        num_proton_gamma=161,
        num_neutrino_frequency=121,
        pgamma_scheme="disabled",
        pair_cascade_iterations=1,
    )


def patch_flux_grid(
    model: Model,
    times: np.ndarray,
    freqs: np.ndarray,
    *,
    sampling: str,
    patch_theta: int,
    patch_phi: int,
    observer_time_s: np.ndarray | None = None,
) -> tuple[np.ndarray, PatchGrid]:
    grid = build_historical_grid(
        model,
        sampling=sampling,
        patch_theta=patch_theta,
        patch_phi=patch_phi,
        observer_time_s=times if observer_time_s is None else observer_time_s,
    )
    total = np.zeros((len(freqs), len(times)), dtype=float)
    for i_theta, theta_center in enumerate(grid.theta_centers):
        e_iso = float(model.jet.energy_iso(0.0, float(theta_center)))
        gamma0 = float(model.jet.gamma0(0.0, float(theta_center)))
        if e_iso <= 0.0 or gamma0 <= 1.0:
            continue
        representative = float(np.mean(grid.domega[i_theta]))
        state = _solve_patch_state(
            model,
            _build_fit_config_for_patch(
                model,
                theta_v=0.0,
                opening_angle_jet=float(np.sqrt(representative / np.pi)),
                e_iso=e_iso,
                gamma0=gamma0,
                theta_center=float(theta_center),
            ),
            np.asarray(times, dtype=float),
            np.asarray(freqs, dtype=float),
        )
        for i_phi, phi_center in enumerate(grid.phi_centers):
            state.config.theta_v = float(_angular_separation(theta_center, phi_center, model.observer.theta_obs, model.observer.phi_obs))
            state.config.opening_angle_jet = float(np.sqrt(float(grid.domega[i_theta, i_phi]) / np.pi))
            observed = project_flux(state, times, freqs, projection_kind="lightcurve")
            total += np.asarray(observed.components["total"], dtype=float)
    return total, grid


def build_historical_grid(
    model: Model,
    *,
    sampling: str,
    patch_theta: int,
    patch_phi: int,
    observer_time_s: np.ndarray | None,
) -> PatchGrid:
    sampling = str(sampling).lower()
    if sampling == "uniform":
        return _uniform_grid(model, patch_theta, patch_phi)
    if sampling == "dominant_region_ioka_v1":
        return _dominant_grid(model, patch_theta, patch_phi, observer_time_s=None)
    if sampling == "dominant_region_ioka_time_v1":
        if observer_time_s is None:
            raise ValueError("dominant_region_ioka_time_v1 requires observer_time_s.")
        return _dominant_grid(model, patch_theta, patch_phi, observer_time_s=np.asarray(observer_time_s, dtype=float))
    raise ValueError(f"unsupported historical sampling: {sampling}")


def _uniform_grid(model: Model, patch_theta: int, patch_phi: int) -> PatchGrid:
    theta_edges = np.linspace(0.0, float(model.jet.theta_max), int(patch_theta) + 1)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, int(patch_phi) + 1)
    return PatchGrid(
        theta_centers=0.5 * (theta_edges[:-1] + theta_edges[1:]),
        theta_edges=theta_edges,
        phi_centers=0.5 * (phi_edges[:-1] + phi_edges[1:]),
        phi_edges=phi_edges,
        domega=_solid_angle_cells(theta_edges, phi_edges),
    )


def _dominant_grid(
    model: Model,
    patch_theta: int,
    patch_phi: int,
    observer_time_s: np.ndarray | None,
) -> PatchGrid:
    theta_max = float(model.jet.theta_max)
    theta_scan = np.linspace(0.0, theta_max, max(4 * int(patch_theta), 256))
    phi_scan = np.linspace(0.0, 2.0 * np.pi, max(4 * int(patch_phi), 256), endpoint=False)
    if observer_time_s is None:
        scan_weight = _dominant_weight(model, theta_scan, phi_scan)
        gamma_time = None
    else:
        sampling_times = _sampling_times(model, observer_time_s)
        gamma_time = _pilot_gamma_theta_time(model, theta_scan, sampling_times, int(patch_theta))
        scan_weight = _time_weight(model, theta_scan, phi_scan, gamma_time)
    theta_weight = np.sin(theta_scan) * np.sqrt(np.mean(scan_weight, axis=1))
    theta_count = int(patch_theta)
    if gamma_time is not None:
        gamma_envelope = np.max(gamma_time, axis=1)
        theta_weight = np.maximum(_normalized_density(theta_scan, theta_weight), _normalized_density(theta_scan, gamma_envelope))
        theta_count = _beaming_theta_count(model, theta_scan, gamma_envelope, int(patch_theta))
    theta_centers, theta_edges = _weighted_centers_edges(theta_scan, theta_weight, theta_count)
    phi_count = int(patch_phi)
    if gamma_time is not None:
        phi_count = _beaming_phi_count(model, theta_scan, gamma_time, int(patch_phi))
    phi_centers, phi_edges, phi_weights = _axisymmetric_phi_quadrature(phi_count, model.observer.phi_obs)
    return PatchGrid(
        theta_centers=theta_centers,
        theta_edges=theta_edges,
        phi_centers=phi_centers,
        phi_edges=phi_edges,
        domega=_solid_angle_cells_from_phi_weights(theta_edges, phi_weights),
    )


def _dominant_weight(model: Model, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    theta_mesh, phi_mesh = theta[:, None], phi[None, :]
    energy = model.jet.energy_iso(phi_mesh, theta_mesh)
    gamma = model.jet.gamma0(phi_mesh, theta_mesh)
    sep = _angular_separation(theta_mesh, phi_mesh, model.observer.theta_obs, model.observer.phi_obs)
    return _with_structure_floor(energy * _doppler_factor(gamma, sep) ** 3)


def _time_weight(model: Model, theta: np.ndarray, phi: np.ndarray, gamma_time: np.ndarray) -> np.ndarray:
    theta_mesh, phi_mesh = theta[:, None], phi[None, :]
    energy = model.jet.energy_iso(phi_mesh, theta_mesh)
    sep = _angular_separation(theta_mesh, phi_mesh, model.observer.theta_obs, model.observer.phi_obs)
    weighted = np.zeros_like(sep, dtype=float)
    for i_time in range(gamma_time.shape[1]):
        weighted += _doppler_factor(gamma_time[:, i_time][:, None], sep) ** 3
    return _with_structure_floor(energy * weighted / float(gamma_time.shape[1]))


def _sampling_times(model: Model, observer_time_s: np.ndarray) -> np.ndarray:
    observer_time_s = np.unique(np.asarray(observer_time_s, dtype=float))
    count = int(model.setups.patch_sampling_num_times)
    if observer_time_s.size <= count:
        return observer_time_s
    return np.exp(np.linspace(np.log(observer_time_s[0]), np.log(observer_time_s[-1]), count))


def _pilot_gamma_theta_time(model: Model, theta: np.ndarray, observer_time_s: np.ndarray, patch_theta: int) -> np.ndarray:
    from asgard_core.asgard_runtime import solve_dynamics
    from asgard_core.asgard_setup import build_setup
    from asgard_core.asgard_state import query_cfg

    sample_count = int(model.setups.patch_sampling_pilot_theta)
    if sample_count <= 0:
        sample_count = max(2 * patch_theta, 48)
    pilot_theta = np.linspace(0.0, float(model.jet.theta_max), sample_count)
    pilot_gamma = []
    for theta_center in pilot_theta:
        e_iso = float(model.jet.energy_iso(0.0, float(theta_center)))
        gamma0 = float(model.jet.gamma0(0.0, float(theta_center)))
        config = _build_fit_config_for_patch(
            model,
            theta_v=0.0,
            opening_angle_jet=float(model.jet.theta_max / sample_count),
            e_iso=e_iso,
            gamma0=gamma0,
            theta_center=float(theta_center),
        )
        query = query_cfg(config, observer_time_s)
        dynamics = solve_dynamics(build_setup(query).boundary, query)
        pilot_gamma.append(np.interp(np.log(observer_time_s), np.log(dynamics.r_tobs), dynamics.r_gamma))
    pilot_gamma = np.asarray(pilot_gamma, dtype=float)
    gamma_time = np.empty((theta.size, observer_time_s.size), dtype=float)
    for i_time in range(observer_time_s.size):
        gamma_time[:, i_time] = np.interp(theta, pilot_theta, pilot_gamma[:, i_time])
    return gamma_time


def _axisymmetric_phi_quadrature(n_phi: int, phi_obs: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    half_centers = np.linspace(0.0, np.pi, int(n_phi))
    phi_centers = np.mod(float(phi_obs) + half_centers, 2.0 * np.pi)
    phi_edges = np.linspace(float(phi_obs), float(phi_obs) + np.pi, int(n_phi) + 1)
    step = np.pi / float(int(n_phi) - 1)
    weights = np.full(int(n_phi), 2.0 * step, dtype=float)
    weights[0] = step
    weights[-1] = step
    return phi_centers, phi_edges, weights


def _beaming_theta_count(model: Model, theta: np.ndarray, gamma: np.ndarray, minimum: int) -> int:
    factor = float(model.setups.patch_sampling_beaming_factor)
    resolution = float(model.setups.patch_sampling_beaming_resolution)
    required = int(np.ceil(float(model.jet.theta_max) * resolution * float(np.max(gamma)) / factor)) + 1
    return max(int(minimum), required)


def _beaming_phi_count(model: Model, theta: np.ndarray, gamma_time: np.ndarray, minimum: int) -> int:
    factor = float(model.setups.patch_sampling_beaming_factor)
    resolution = float(model.setups.patch_sampling_beaming_resolution)
    angular_frequency = float(np.max(gamma_time * np.sin(theta)[:, None]))
    required = int(np.ceil(np.pi * resolution * angular_frequency / factor)) + 1
    return max(int(minimum), required, 2)


def _normalized_density(x: np.ndarray, density: np.ndarray) -> np.ndarray:
    total = float(np.trapezoid(np.asarray(density, dtype=float), np.asarray(x, dtype=float)))
    if total <= 0.0:
        raise ValueError("historical adaptive theta sampling requires positive density.")
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
    beta = np.sqrt(np.maximum(gamma * gamma - 1.0, 0.0)) / gamma
    return 1.0 / (gamma * (1.0 - beta * np.cos(angle)))


def _angular_separation(theta: float | np.ndarray, phi: float | np.ndarray, theta_obs: float, phi_obs: float) -> np.ndarray:
    cos_sep = np.cos(theta) * np.cos(theta_obs) + np.sin(theta) * np.sin(theta_obs) * np.cos(phi - phi_obs)
    return np.arccos(np.clip(cos_sep, -1.0, 1.0))


def _solid_angle_cells(theta_edges: np.ndarray, phi_edges: np.ndarray) -> np.ndarray:
    return (np.cos(theta_edges[:-1]) - np.cos(theta_edges[1:]))[:, None] * np.diff(phi_edges)[None, :]


def _solid_angle_cells_from_phi_weights(theta_edges: np.ndarray, phi_weights: np.ndarray) -> np.ndarray:
    return (np.cos(theta_edges[:-1]) - np.cos(theta_edges[1:]))[:, None] * np.asarray(phi_weights, dtype=float)[None, :]
