from __future__ import annotations

from typing import Any

import numpy as np
from scipy.special import erf

from asgard_core.angular_sampling import is_axisymmetric_jet
from asgard_core.asgard_state import (
    SolveState,
    _build_observer_setup_from_state,
    _forward_synchrotron_absorption_transfer,
)
from src.Electron.electron_radiation import electron_radiation_kernel as electron_radiation_module
from asgard_core.asgard_config import RuntimeConfig, SpectrumOutputConfig
from asgard_core.asgard_postprocess import (
    OUTPUT_BANDS,
    build_spectrum_dataset_names,
    build_spectrum_frequency_grid,
    build_multiband_observer_frequencies,
    combine_multiband_flux,
    compute_light_curve_redchi,
    compute_spectrum_redchi,
    select_spectrum_time_index,
)
from src import constants
from .api_adaptive import _observe_parts, _observe_total
from .api_fit import FitResult, Param
from .api_model import (
    Model,
    Numerics,
    Observer,
    ObserverGrid,
    PolarizationResult,
    Radiation,
    ReverseShock,
    SolverOptions,
    TabulatedMedium,
    Hadronic,
    SkyImage,
    UniformMedium,
    WindMedium,
    top_hat_jet,
    _direct_tophat_patch_config,
    _iter_patch_elements,
    _iter_solved_patch_elements,
    _project_surface_element,
    _solve_patch_state,
)


def _total_matrix(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    timings: dict[str, float] | None = None,
    projection_kind: str = "lightcurve",
) -> np.ndarray:
    times_s = np.asarray(times_s, dtype=float)
    nu_hz = np.asarray(nu_hz, dtype=float)
    if model.jet.kind == "tophat" and model._supports_direct_kernel():
        config = _direct_tophat_patch_config(model)
        state = _solve_patch_state(model, config, times_s, nu_hz, timings=timings)
        return _observe_total(state, times_s, nu_hz, timings=timings, projection_kind=projection_kind)
    total = np.zeros((nu_hz.shape[0], times_s.shape[0]), dtype=float)
    for _patch, state in _iter_solved_patch_elements(
        model,
        times_s,
        nu_hz,
        _iter_patch_elements(model),
        timings=timings,
    ):
        total += _observe_total(state, times_s, nu_hz, timings=timings, projection_kind=projection_kind)
    return total


def _compute_polarization(
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
        raise ValueError("polarization() requires one-dimensional times_s and nu_hz grids.")
    if times_s.size == 0 or nu_hz.size == 0:
        raise ValueError("polarization() grids must be non-empty.")
    if np.any(times_s <= 0.0) or np.any(nu_hz <= 0.0):
        raise ValueError("polarization() times and frequencies must be positive.")
    if magnetic_geometry not in {"shock_random", "toroidal"}:
        raise ValueError("magnetic_geometry must be 'shock_random' or 'toroidal'.")
    if local_emissivity not in {"analytic", "analytic_then_kernel"}:
        raise ValueError("local_emissivity must be 'analytic' or 'analytic_then_kernel'.")
    if magnetic_geometry == "toroidal" and not is_axisymmetric_jet(model.jet):
        raise NotImplementedError("toroidal polarization currently requires an axisymmetric jet.")

    total_i = np.zeros((nu_hz.shape[0], times_s.shape[0]), dtype=float)
    total_q, total_u = np.zeros((2, *total_i.shape), dtype=float)
    components = {
        "fwd_sync": _empty_stokes(total_i),
        "rev_sync": _empty_stokes(total_i),
        "hadronic_sync": _empty_stokes(total_i),
        "rev_hadronic_sync": _empty_stokes(total_i),
    }
    sightline, sky_x_axis, sky_y_axis = _sky_basis(model.observer)
    active_patch_found = False

    for patch, state in _iter_solved_patch_elements(model, times_s, nu_hz, _iter_patch_elements(model)):
        active_patch_found = True
        patch_direction = _direction_vector(patch.theta_center, patch.phi_center)
        cos2pa, sin2pa = _patch_polarization_angle_factors(
            magnetic_geometry,
            patch_direction,
            sightline,
            sky_x_axis,
            sky_y_axis,
        )
        _accumulate_patch_polarization(
            components["fwd_sync"],
            state,
            state.components.fwd_sync,
            times_s,
            nu_hz,
            _component_polarization_fraction(
                state,
                model.fwd_rad.p,
                local_emissivity,
                reverse=False,
            ),
            cos2pa,
            sin2pa,
            _shock_random_anisotropy(magnetic_geometry, patch.theta_v, state.components.fwd.gamma),
            patch.domega,
        )
        if state.reverse_emission is not None and model.rvs_rad is not None:
            rev_sync = np.asarray(state.reverse_emission.l_syn_spec, dtype=float) * np.asarray(state.observer.prefactor, dtype=float)
            _accumulate_patch_polarization(
                components["rev_sync"],
                state,
                rev_sync,
                times_s,
                nu_hz,
                _component_polarization_fraction(
                    state,
                    model.rvs_rad.p,
                    local_emissivity,
                    reverse=True,
                ),
                cos2pa,
                sin2pa,
                _shock_random_anisotropy(magnetic_geometry, patch.theta_v, _reverse_shock_gamma_array(state)),
                patch.domega,
            )
            rev_hadronic_sync = _patch_reverse_hadronic_synchrotron_component(state)
            if rev_hadronic_sync is not None:
                rev_hadronic_source, rev_hadronic_pi = rev_hadronic_sync
                _accumulate_patch_polarization(
                    components["rev_hadronic_sync"],
                    state,
                    rev_hadronic_source,
                    times_s,
                    nu_hz,
                    rev_hadronic_pi,
                    cos2pa,
                    sin2pa,
                    _shock_random_anisotropy(magnetic_geometry, patch.theta_v, _reverse_shock_gamma_array(state)),
                    patch.domega,
                )
        hadronic_sync = _patch_hadronic_synchrotron_component(state)
        if hadronic_sync is not None:
            hadronic_source, hadronic_pi = hadronic_sync
            _accumulate_patch_polarization(
                components["hadronic_sync"],
                state,
                hadronic_source,
                times_s,
                nu_hz,
                hadronic_pi,
                cos2pa,
                sin2pa,
                _shock_random_anisotropy(magnetic_geometry, patch.theta_v, state.components.fwd.gamma),
                patch.domega,
            )

    if not active_patch_found:
        raise ValueError("No active jet patches were found for polarization().")
    for component in components.values():
        total_i += component["I"]
        total_q += component["Q"]
        total_u += component["U"]
    linear = np.zeros_like(total_i)
    positive = total_i > 0.0
    linear[positive] = np.sqrt(total_q[positive] * total_q[positive] + total_u[positive] * total_u[positive]) / total_i[positive]
    pa = 0.5 * np.arctan2(total_u, total_q)
    return PolarizationResult(
        I_sync=total_i,
        Q=total_q,
        U=total_u,
        linear_polarization=linear,
        polarization_angle_rad=pa,
        components=components,
    )


def _empty_stokes(template: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "I": np.zeros_like(template),
        "Q": np.zeros_like(template),
        "U": np.zeros_like(template),
    }


def _synchrotron_intrinsic_polarization(p_index: float) -> float:
    p = float(p_index)
    return (p + 1.0) / (p + 7.0 / 3.0)


def _component_polarization_fraction(
    state: SolveState,
    p_index: float,
    local_emissivity: str,
    *,
    reverse: bool,
) -> float | np.ndarray:
    if local_emissivity == "analytic":
        return _synchrotron_intrinsic_polarization(p_index)
    if reverse:
        if state.dynamics.reverse_shock is None or state.reverse_emission is None:
            return _synchrotron_intrinsic_polarization(p_index)
        gamma_grid = np.asarray(state.dynamics.reverse_shock.gam_e, dtype=float)
        distribution = np.asarray(state.dynamics.reverse_shock.d_n_gam_e, dtype=float)
        magnetic_field = np.asarray(state.dynamics.reverse_shock.magnetic_field_g, dtype=float)
        radius_cm = np.full(magnetic_field.size, float(state.dynamics.reverse_shock.r_cross))
    else:
        gamma_grid = np.asarray(state.electron.gam_e, dtype=float)
        distribution = np.asarray(state.electron.d_n_gam_e, dtype=float)
        magnetic_field = np.asarray(state.components.fwd.magnetic_field_g, dtype=float)
        radius_cm = np.asarray(state.components.fwd.radius_cm, dtype=float)
    frequency = np.asarray(state.setup.seed_frequency_hz, dtype=float)
    pi_grid = np.zeros((frequency.size, magnetic_field.size), dtype=float)
    for i_shell in range(magnetic_field.size):
        _, _, pi_nu = electron_radiation_module.get_syn_polarization_selected(
            int(state.config.index_syn_integr),
            float(radius_cm[i_shell]),
            float(magnetic_field[i_shell]),
            int(state.config.num_threads),
            gamma_grid,
            distribution[:, i_shell],
            frequency,
            float(p_index),
        )
        pi_grid[:, i_shell] = np.asarray(pi_nu, dtype=float)
    return pi_grid


def _shock_random_anisotropy(magnetic_geometry: str, theta_v: float, gamma: np.ndarray) -> np.ndarray:
    if magnetic_geometry != "shock_random":
        return np.ones_like(np.asarray(gamma, dtype=float))
    gamma_arr = np.asarray(gamma, dtype=float)
    beta = np.sqrt(1.0 - gamma_arr ** (-2))
    mu = np.cos(float(theta_v))
    mu_prime = (mu - beta) / (1.0 - beta * mu)
    sin2_prime = 1.0 - mu_prime * mu_prime
    return sin2_prime / (1.0 + mu_prime * mu_prime)


def _reverse_shock_gamma_array(state: SolveState) -> np.ndarray:
    return np.full_like(
        np.asarray(state.dynamics.radius, dtype=float),
        float(state.dynamics.reverse_shock.gam20),
    )


def _patch_polarization_angle_factors(
    magnetic_geometry: str,
    patch_direction: np.ndarray,
    sightline: np.ndarray,
    sky_x_axis: np.ndarray,
    sky_y_axis: np.ndarray,
) -> tuple[float, float]:
    if magnetic_geometry == "shock_random":
        e_vector = patch_direction - np.dot(patch_direction, sightline) * sightline
    else:
        jet_axis = np.array([0.0, 0.0, 1.0], dtype=float)
        b_vector = np.cross(jet_axis, patch_direction)
        e_vector = np.cross(sightline, b_vector)
    e_x = float(np.dot(e_vector, sky_x_axis))
    e_y = float(np.dot(e_vector, sky_y_axis))
    norm = np.hypot(e_x, e_y)
    if norm <= 0.0:
        return 0.0, 0.0
    e_x = e_x / norm
    e_y = e_y / norm
    return e_x * e_x - e_y * e_y, 2.0 * e_x * e_y


def _accumulate_patch_polarization(
    accumulator: dict[str, np.ndarray],
    state: SolveState,
    source_component: np.ndarray,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    intrinsic_pi: float | np.ndarray,
    cos2pa: float,
    sin2pa: float,
    shell_anisotropy: np.ndarray,
    patch_solid_angle: float,
) -> None:
    source = np.asarray(source_component, dtype=float)
    if not np.any(source):
        return
    pi_local = np.asarray(intrinsic_pi, dtype=float)
    if pi_local.ndim == 0:
        polarized_source = source * float(pi_local)
    else:
        polarized_source = source * pi_local
    polarized_source = polarized_source * np.asarray(shell_anisotropy, dtype=float)[None, :]
    observer_setup = _build_observer_setup_from_state(state, times_s)
    i_obs = _project_surface_element(
        observer_setup,
        state.components.fwd.characteristic_time_s,
        state.components.fwd.gamma,
        state.components.fwd.radius_cm,
        state.setup.seed_frequency_hz,
        source,
        nu_hz,
        patch_solid_angle,
    )
    p_obs = _project_surface_element(
        observer_setup,
        state.components.fwd.characteristic_time_s,
        state.components.fwd.gamma,
        state.components.fwd.radius_cm,
        state.setup.seed_frequency_hz,
        polarized_source,
        nu_hz,
        patch_solid_angle,
    )
    accumulator["I"] += i_obs
    accumulator["Q"] += p_obs * cos2pa
    accumulator["U"] += p_obs * sin2pa


def _patch_hadronic_synchrotron_component(state: SolveState) -> tuple[np.ndarray, np.ndarray] | None:
    if state.hadronic is None:
        return None
    import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module

    luminosity = np.asarray(state.hadronic.l_had_syn_spec, dtype=float)
    weighted_pi_luminosity = luminosity * _hadronic_synchrotron_pi_grid(
        state,
        hadronic_fortran_module,
        state.hadronic.gam_p * constants.para_m_p_gev,
        state.hadronic.d_n_gam_p / constants.para_m_p_gev,
        constants.para_m_p_gev,
        luminosity,
    )
    if state.hadronic.l_had_pion_synch is not None:
        pion_luminosity = np.asarray(state.hadronic.l_had_pion_synch, dtype=float)
        pion_density = (
            np.asarray(state.hadronic.d_n_gam_pi_plus, dtype=float)
            + np.asarray(state.hadronic.d_n_gam_pi_minus, dtype=float)
        ) / constants.para_m_pi_charged_gev
        luminosity = luminosity + pion_luminosity
        weighted_pi_luminosity = weighted_pi_luminosity + pion_luminosity * _hadronic_synchrotron_pi_grid(
            state,
            hadronic_fortran_module,
            state.hadronic.gam_secondary * constants.para_m_pi_charged_gev,
            pion_density,
            constants.para_m_pi_charged_gev,
            pion_luminosity,
        )
    if state.hadronic.l_had_muon_synch is not None:
        muon_luminosity = np.asarray(state.hadronic.l_had_muon_synch, dtype=float)
        muon_density = (
            np.asarray(state.hadronic.d_n_gam_mu_minus_left, dtype=float)
            + np.asarray(state.hadronic.d_n_gam_mu_minus_right, dtype=float)
            + np.asarray(state.hadronic.d_n_gam_mu_plus_left, dtype=float)
            + np.asarray(state.hadronic.d_n_gam_mu_plus_right, dtype=float)
        ) / constants.para_m_mu_gev
        luminosity = luminosity + muon_luminosity
        weighted_pi_luminosity = weighted_pi_luminosity + muon_luminosity * _hadronic_synchrotron_pi_grid(
            state,
            hadronic_fortran_module,
            state.hadronic.gam_secondary * constants.para_m_mu_gev,
            muon_density,
            constants.para_m_mu_gev,
            muon_luminosity,
        )
    transfer = _forward_synchrotron_absorption_transfer(
        electron=state.electron,
        radius_cm=state.dynamics.radius,
        magnetic_field_g=state.components.fwd.magnetic_field_g,
        seed_frequency_hz=state.setup.seed_frequency_hz,
        config=state.config,
    )
    source = luminosity * transfer * np.asarray(state.observer.prefactor, dtype=float)
    pi_grid = np.zeros_like(luminosity)
    positive = luminosity > 0.0
    pi_grid[positive] = weighted_pi_luminosity[positive] / luminosity[positive]
    return source, pi_grid


def _patch_reverse_hadronic_synchrotron_component(state: SolveState) -> tuple[np.ndarray, np.ndarray] | None:
    if state.reverse_emission is None or state.reverse_emission.rs_hadronic is None:
        return None
    import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module

    rs_hadronic = state.reverse_emission.rs_hadronic
    luminosity = np.asarray(rs_hadronic.l_had_syn_spec, dtype=float)
    source = luminosity * np.asarray(state.observer.prefactor, dtype=float)
    pi_grid = _generic_hadronic_synchrotron_pi_grid(
        hadronic_fortran_module,
        rs_hadronic.gam_p * constants.para_m_p_gev,
        np.asarray(rs_hadronic.d_n_gam_p, dtype=float) / constants.para_m_p_gev,
        np.asarray(state.setup.seed_frequency_hz, dtype=float),
        np.asarray(state.dynamics.reverse_shock.magnetic_field_g, dtype=float),
        constants.para_m_p_gev,
        float(state.config.hadronic.p_p),
        luminosity,
    )
    return source, pi_grid


def _hadronic_synchrotron_pi_grid(
    state: SolveState,
    hadronic_fortran_module,
    hadron_energy_gev: np.ndarray,
    density_per_gev: np.ndarray,
    particle_mass_gev: float,
    luminosity_grid: np.ndarray,
) -> np.ndarray:
    energy = np.asarray(hadron_energy_gev, dtype=float)
    density = np.asarray(density_per_gev, dtype=float)
    frequency = np.asarray(state.setup.seed_frequency_hz, dtype=float)
    magnetic_field = np.asarray(state.components.fwd.magnetic_field_g, dtype=float)
    return _generic_hadronic_synchrotron_pi_grid(
        hadronic_fortran_module,
        energy,
        density,
        frequency,
        magnetic_field,
        particle_mass_gev,
        float(state.config.hadronic.p_p),
        luminosity_grid,
    )


def _generic_hadronic_synchrotron_pi_grid(
    hadronic_fortran_module,
    hadron_energy_gev: np.ndarray,
    density_per_gev: np.ndarray,
    frequency_hz: np.ndarray,
    magnetic_field_g: np.ndarray,
    particle_mass_gev: float,
    p_index: float,
    luminosity_grid: np.ndarray,
) -> np.ndarray:
    energy = np.asarray(hadron_energy_gev, dtype=float)
    density = np.asarray(density_per_gev, dtype=float)
    frequency = np.asarray(frequency_hz, dtype=float)
    magnetic_field = np.asarray(magnetic_field_g, dtype=float)
    luminosity = np.asarray(luminosity_grid, dtype=float)
    analytic_pi = _synchrotron_intrinsic_polarization(p_index)
    pi_grid = np.zeros((frequency.size, magnetic_field.size), dtype=float)
    for i_shell in range(magnetic_field.size):
        if not np.any(luminosity[:, i_shell] > 0.0):
            pi_grid[:, i_shell] = analytic_pi
            continue
        if magnetic_field[i_shell] <= 0.0:
            raise RuntimeError("positive hadronic synchrotron luminosity requires magnetic_field_g > 0.")
        pi_grid[:, i_shell] = np.asarray(
            hadronic_fortran_module.fs_hadronic_syn_polarization_shell(
                energy,
                density[:, i_shell],
                frequency,
                float(particle_mass_gev),
                float(magnetic_field[i_shell]),
                float(p_index),
            ),
            dtype=float,
        )
    return pi_grid


def _render_sky_image(model: Model, times_s: np.ndarray, nu_obs: float, fov: float, npixel: int) -> SkyImage:
    base_pixel_size = float(fov) / float(npixel)
    patch_cache: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]] = []
    if model.observer.lumi_dist_cm is None or model.observer.lumi_dist_cm <= 0.0:
        raise ValueError("Observer.luminosity_distance_cm must be set for sky_image().")
    angular_diameter_distance_cm = model.observer.lumi_dist_cm / (1.0 + model.observer.z) ** 2
    sightline, sky_x_axis, sky_y_axis = _sky_basis(model.observer)
    frequencies_hz = np.array([nu_obs], dtype=float)
    required_half_fov = 0.5 * float(fov)
    total_solid_angle = 2.0 * np.pi * (1.0 - np.cos(float(model.jet.theta_max)))
    collapse_phi = model.jet.kind == "tophat" and abs(float(model.observer.theta_obs)) <= 1.0e-12
    img_patches = _iter_img_patches(model, npixel, collapse_phi=collapse_phi)
    for patch, state in _iter_solved_patch_elements(
        model,
        times_s,
        frequencies_hz,
        img_patches,
        opening_angle_jet=model.jet.theta_max,
    ):
        observed = _observe_parts(state, times_s, frequencies_hz)
        patch_flux = np.asarray(observed.total[0, :], dtype=float)
        radius_cm = np.exp(
            np.interp(
                np.log(times_s),
                np.log(np.asarray(state.components.fwd.characteristic_time_s, dtype=float)),
                np.log(np.asarray(state.components.fwd.radius_cm, dtype=float)),
                left=np.log(float(state.components.fwd.radius_cm[0])),
                right=np.log(float(state.components.fwd.radius_cm[-1])),
            )
        )
        if not np.any(np.isfinite(patch_flux) & (patch_flux > 0.0)):
            continue

        x_center, y_center = _project_patch_to_sky(
            radius_cm,
            patch.theta_center,
            patch.phi_center,
            sightline,
            sky_x_axis,
            sky_y_axis,
            angular_diameter_distance_cm,
        )
        sigma = np.maximum(
            radius_cm * np.sin(max(patch.half_angle, 1.0e-12)) / angular_diameter_distance_cm / 2.0,
            0.5 * base_pixel_size,
        )
        patch_weight = patch.domega / total_solid_angle if total_solid_angle > 0.0 else 0.0
        patch_cache.append((patch_flux, x_center, y_center, sigma, patch_weight))
        span = np.max(np.maximum(np.abs(x_center), np.abs(y_center)) + 8.0 * sigma)
        required_half_fov = max(required_half_fov, float(span))

    fov_eff = max(float(fov), 2.0 * required_half_fov)
    pixel_size = fov_eff / float(npixel)
    extent = np.array([-0.5 * fov_eff, 0.5 * fov_eff, -0.5 * fov_eff, 0.5 * fov_eff], dtype=float)
    pixel_axis = np.linspace(extent[0] + 0.5 * pixel_size, extent[1] - 0.5 * pixel_size, npixel)
    x_grid, y_grid = np.meshgrid(pixel_axis, pixel_axis, indexing="ij")
    image = np.zeros((times_s.shape[0], npixel, npixel), dtype=float)
    x_axis = x_grid[:, 0]
    y_axis = y_grid[0, :]
    for patch_flux, x_center, y_center, sigma, patch_weight in patch_cache:
        image += (
            patch_weight
            * patch_flux[:, None, None]
            * _gaussian_splat_stack(x_axis, y_axis, pixel_size, x_center, y_center, sigma)
        )

    direct_total = np.asarray(model.flux_density_grid(times_s, frequencies_hz).total[0, :], dtype=float)
    raw_total = image.sum(axis=(1, 2)) * pixel_size * pixel_size
    scale = np.where(raw_total > 0.0, direct_total / raw_total, 0.0)
    image *= scale[:, None, None]
    rendered_total = image.sum(axis=(1, 2)) * pixel_size * pixel_size
    x_centroid, y_centroid = np.zeros((2, times_s.shape[0]), dtype=float)
    x_weights = image.sum(axis=2)
    y_weights = image.sum(axis=1)
    total_brightness = image.sum(axis=(1, 2))
    valid = total_brightness > 0.0
    if np.any(valid):
        x_centroid[valid] = (
            np.sum(x_weights[valid, :] * pixel_axis[None, :], axis=1) / total_brightness[valid]
        )
        y_centroid[valid] = (
            np.sum(y_weights[valid, :] * pixel_axis[None, :], axis=1) / total_brightness[valid]
        )

    return SkyImage(
        image=image,
        extent=extent,
        pixel_solid_angle=pixel_size * pixel_size,
        pixel_size=pixel_size,
        direct_flux=direct_total,
        rendered_flux=rendered_total,
        normalization_scale=scale,
        x_centroid=x_centroid,
        y_centroid=y_centroid,
    )


def _iter_img_patches(model: Model, npixel: int, *, collapse_phi: bool = False):
    theta_bins = min(model.setups.patch_theta, max(2, int(np.ceil(np.sqrt(npixel) / 6.0))))
    phi_bins = min(model.setups.patch_phi, max(12, 6 * theta_bins))
    theta_edges = np.linspace(0.0, model.jet.theta_max, theta_bins + 1)
    for i_theta in range(theta_bins):
        theta1 = theta_edges[i_theta]
        theta2 = theta_edges[i_theta + 1]
        theta_center = 0.5 * (theta1 + theta2)
        if collapse_phi:
            domega = (np.cos(theta1) - np.cos(theta2)) * (2.0 * np.pi)
            patch_half_angle = np.sqrt(max(domega, 1.0e-12) / np.pi)
            yield 0.0, theta_center, patch_half_angle, domega
            continue
        phi_edges = np.linspace(0.0, 2.0 * np.pi, phi_bins + 1)
        for i_phi in range(phi_bins):
            phi1 = phi_edges[i_phi]
            phi2 = phi_edges[i_phi + 1]
            phi_center = 0.5 * (phi1 + phi2)
            domega = (np.cos(theta1) - np.cos(theta2)) * (phi2 - phi1)
            patch_half_angle = np.sqrt(max(domega, 1.0e-12) / np.pi)
            yield phi_center, theta_center, patch_half_angle, domega


def _sky_basis(observer: Observer) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sightline = _direction_vector(observer.theta_obs, observer.phi_obs)
    trial = np.array([0.0, 0.0, 1.0], dtype=float)
    sky_x = np.cross(trial, sightline)
    if np.linalg.norm(sky_x) < 1.0e-12:
        sky_x = np.array([1.0, 0.0, 0.0], dtype=float)
    else:
        sky_x = sky_x / np.linalg.norm(sky_x)
    sky_y = np.cross(sightline, sky_x)
    sky_y = sky_y / np.linalg.norm(sky_y)
    return sightline, sky_x, sky_y


def _direction_vector(theta: float, phi: float) -> np.ndarray:
    return np.array(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ],
        dtype=float,
    )


def _project_patch_to_sky(
    radius_cm: np.ndarray,
    theta_center: float,
    phi_center: float,
    sightline: np.ndarray,
    sky_x_axis: np.ndarray,
    sky_y_axis: np.ndarray,
    angular_diameter_distance_cm: float,
) -> tuple[np.ndarray, np.ndarray]:
    direction = _direction_vector(theta_center, phi_center)
    position = radius_cm[:, None] * direction[None, :]
    line_of_sight = np.sum(position * sightline[None, :], axis=1)
    transverse = position - line_of_sight[:, None] * sightline[None, :]
    x_center = transverse @ sky_x_axis / angular_diameter_distance_cm
    y_center = transverse @ sky_y_axis / angular_diameter_distance_cm
    return x_center, y_center


def _gaussian_splat_stack(
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    pixel_size: float,
    x_center: np.ndarray,
    y_center: np.ndarray,
    sigma: np.ndarray,
) -> np.ndarray:
    sigma = np.maximum(np.asarray(sigma, dtype=float), 1.0e-30)
    x_axis = np.asarray(x_axis, dtype=float)
    y_axis = np.asarray(y_axis, dtype=float)
    x_edges = np.concatenate((x_axis - 0.5 * pixel_size, np.array([x_axis[-1] + 0.5 * pixel_size], dtype=float)))
    y_edges = np.concatenate((y_axis - 0.5 * pixel_size, np.array([y_axis[-1] + 0.5 * pixel_size], dtype=float)))
    sigma_view = np.sqrt(2.0) * sigma[:, None]
    x_frac = 0.5 * (
        erf((x_edges[1:][None, :] - np.asarray(x_center, dtype=float)[:, None]) / sigma_view)
        - erf((x_edges[:-1][None, :] - np.asarray(x_center, dtype=float)[:, None]) / sigma_view)
    )
    y_frac = 0.5 * (
        erf((y_edges[1:][None, :] - np.asarray(y_center, dtype=float)[:, None]) / sigma_view)
        - erf((y_edges[:-1][None, :] - np.asarray(y_center, dtype=float)[:, None]) / sigma_view)
    )
    return x_frac[:, :, None] * y_frac[:, None, :] / (pixel_size * pixel_size)


def _param_path(model: Model, param: Param) -> str:
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


def _as_model(cfg: Any) -> Model:
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


def observe(
    model: Model,
    config: RuntimeConfig | None = None,
    spectrum_output: SpectrumOutputConfig | None = None,
):
    from asgard_core.asgard_config import FitResult as _PhysicalFitResult

    setups = model.setups
    t_obs_s = np.logspace(
        np.log10(setups.observer_time_min_s),
        np.log10(setups.observer_time_max_s),
        setups.num_tobs,
    )

    n_xrt, all_freqs = build_multiband_observer_frequencies()
    grid_result = model.flux_density_grid(t_obs_s, all_freqs, projection_kind="lightcurve")
    bands_flux = combine_multiband_flux(grid_result.total, all_freqs, n_xrt)

    details = model.details()
    char_time = details.fwd.t_obs

    redchi = 0.0
    if config is not None:
        redchi = compute_light_curve_redchi(bands_flux, t_obs_s, config)

    spectrum_freq_hz = None
    spectrum_fnu = None
    spectrum_redchi = None
    spectrum_time_s = None
    if spectrum_output is not None and spectrum_output.enabled:
        spec_freqs = build_spectrum_frequency_grid(
            type('_Cfg', (), {'spectrum_output': spectrum_output})()
        )
        spec_grid = model.flux_density_grid(t_obs_s, spec_freqs, projection_kind="sed")
        spectrum_freq_hz = spec_freqs
        spectrum_fnu = spec_grid.total
        if spectrum_output.time_s is not None or spectrum_output.dataset_names:
            spec_idx = select_spectrum_time_index(t_obs_s, spectrum_output.time_s)
            spectrum_time_s = float(t_obs_s[spec_idx])
            spectrum_redchi = compute_spectrum_redchi(
                spectrum_fnu[:, spec_idx],
                spectrum_freq_hz,
                spectrum_names=spectrum_output.dataset_names or build_spectrum_dataset_names(),
            )

    return _PhysicalFitResult(
        t_obs_s=t_obs_s,
        characteristic_time_s=char_time,
        bands=OUTPUT_BANDS,
        bands_flux=bands_flux,
        redchi=redchi,
        spectrum_redchi=spectrum_redchi,
        spectrum_time_s=spectrum_time_s,
        spectrum_freq_hz=spectrum_freq_hz,
        spectrum_fnu=spectrum_fnu,
    )


def _build_model_from_fit_config(config: RuntimeConfig) -> Model:
    ssc_enabled = config.index_y != 0
    kn_enabled = config.index_y == 1
    reverse = config.reverse_shock
    hadronic = config.hadronic
    reverse_enabled = bool(reverse.enabled)

    if config.density_profile_radius_cm or config.density_profile_n_cm3:
        medium = TabulatedMedium(config.density_profile_radius_cm, config.density_profile_n_cm3, label="runtime_density_profile")
    elif config.a_star > 0.0:
        medium = WindMedium(a_star=config.a_star, density_floor_cm3=config.d_ne, density_cap_cm3=None)
    else:
        medium = UniformMedium(number_density_cm3=config.d_ne)

    reverse_delta_t_s = 10.0 if reverse.delta_t_s is None else float(reverse.delta_t_s)

    rvs_rad = None
    if reverse_enabled:
        rvs_rad = Radiation(
            epsilon_e=reverse.epsilon_e if reverse.epsilon_e is not None else config.epsilon_e,
            epsilon_B=reverse.epsilon_b if reverse.epsilon_b is not None else config.epsilon_b,
            p=reverse.p if reverse.p is not None else config.p,
            proton_energy_fraction=0.0,
            epsilon_b_floor=None,
            magnetic_decay_alpha_t=config.magnetic_decay_alpha_t,
            magnetic_decay_t0_s=config.magnetic_decay_t0_s,
            accelerated_electron_fraction=reverse.f_e if reverse.f_e is not None else config.f_e,
            thermal_electrons=False,
            include_ssc=bool(reverse.include_ssc),
            include_kn_correction=kn_enabled,
            proton_synch=True,
            include_pgamma=False,
            bethe_heitler=False,
            hadronic_inverse_compton=False,
            pp=False,
            neutrino=False,
            acceleration_efficiency=hadronic.eta_acc,
            reverse_proton_energy_fraction=hadronic.reverse_epsilon_p,
            pgamma_scheme=hadronic.pgamma_scheme,
            pair_production=False,
        )

    return Model(
        jet=top_hat_jet(
            energy_iso_erg=config.e_iso,
            initial_lorentz_factor=config.eta_0,
            opening_angle_rad=config.opening_angle_jet,
            shell_duration_s=reverse_delta_t_s if reverse_enabled else None,
            magnetar=None,
            spreading=False,
        ),
        medium=medium,
        observer=Observer(
            z=config.z,
            viewing_angle_rad=config.theta_v,
            viewing_azimuth_rad=0.0,
            luminosity_distance_cm=config.luminosity_distance_cm_override,
        ),
        fwd_rad=Radiation(
            epsilon_e=config.epsilon_e,
            epsilon_B=config.epsilon_b,
            epsilon_b_floor=config.epsilon_b_floor,
            magnetic_decay_alpha_t=config.magnetic_decay_alpha_t,
            magnetic_decay_t0_s=config.magnetic_decay_t0_s,
            p=config.p,
            proton_energy_fraction=hadronic.epsilon_p,
            accelerated_electron_fraction=config.f_e,
            thermal_electrons=config.thermal_electrons,
            include_ssc=ssc_enabled,
            include_kn_correction=kn_enabled,
            proton_synch=hadronic.include_proton_synch,
            include_pgamma=hadronic.include_pg,
            bethe_heitler=hadronic.include_bethe_heitler,
            hadronic_inverse_compton=hadronic.include_hadronic_inverse_compton,
            pp=hadronic.include_pp,
            neutrino=hadronic.include_neutrino,
            acceleration_efficiency=hadronic.eta_acc,
            reverse_proton_energy_fraction=hadronic.reverse_epsilon_p,
            pgamma_scheme=hadronic.pgamma_scheme,
            pair_production=hadronic.include_pair_production,
        ),
        rvs_rad=rvs_rad,
        numerics=Numerics(
            num_radius=config.num_r,
            num_theta=config.num_theta,
            num_phi=config.num_phi,
            num_observer_time=config.num_tobs,
            num_electron_gamma=config.num_gam_e,
            num_photon_frequency=config.num_nu,
            num_chi=config.num_chi,
            num_threads=config.num_threads,
            electron_adaptive_substeps=config.electron_adaptive_substeps,
            electron_substep_rtol=config.electron_substep_rtol,
            electron_substep_min=config.electron_substep_min,
            electron_substep_max=config.electron_substep_max,
            initial_radius_cm=config.initial_radius_cm,
        ),
        observer_grid=ObserverGrid(
            time_min_s=10 ** config.t_obs_min_log10,
            time_max_s=10 ** config.t_obs_max_log10,
        ),
        solver_options=SolverOptions(
            electron_solver=config.electron_solver,
            dynamics_solver=config.dynamics_kernel,
            geometry_projection=config.geometry_kernel,
            electron_photon_coupling=config.electron_photon_coupling,
            ssc_cooling_mode="none" if config.index_y == 0 else "numeric_ic_kn" if config.index_y == 1 else "nakar_y_thomson",
            synchrotron_integration="fixed_grid",
            cooling_kernel=config.cooling_kernel,
            radiation_kernel=config.radiation_kernel,
            structured_backend=config.structured_backend,
            patch_sampling=config.patch_sampling,
            patch_projection=config.patch_projection,
            patch_sampling_pilot_theta=config.patch_sampling_pilot_theta,
            patch_sampling_num_times=config.patch_sampling_num_times,
            patch_sampling_beaming_factor=config.patch_sampling_beaming_factor,
            patch_sampling_beaming_resolution=config.patch_sampling_beaming_resolution,
            structured_parallel_mode=config.structured_parallel_mode,
            structured_outer_threads=config.structured_outer_threads,
            structured_inner_threads=config.structured_inner_threads,
            projection_adaptive_rtol=config.projection_adaptive_rtol,
            projection_adaptive_max_depth=config.projection_adaptive_max_depth,
            fullhide2d_transport_model=config.fullhide2d_transport_model,
            fullhide2d_stochastic_accel_norm=config.fullhide2d_stochastic_accel_norm,
            fullhide2d_escape_mode=config.fullhide2d_escape_mode,
            nu_callback=config.nu_callback,
        ),
        reverse_shock=ReverseShock(
            enabled=reverse_enabled,
            shell_duration_s=reverse_delta_t_s,
            upstream_sigma=float(reverse.sigma) if reverse_enabled else 0.0,
            include_cross_zone_ic=bool(reverse.include_cross_zone_ic) if reverse_enabled else False,
            include_ssc=bool(reverse.include_ssc) if reverse_enabled else False,
        ),
        hadronic=Hadronic(
            enabled=hadronic.enabled,
            solver=hadronic.solver,
            num_proton_gamma=hadronic.num_gam_p,
            num_neutrino_frequency=hadronic.num_nu_nu,
            pgamma_scheme=hadronic.pgamma_scheme,
            pair_cascade_iterations=hadronic.pair_cascade_iterations,
        ),
    )


def run_fit(config: RuntimeConfig) -> FitResult:
    return observe(_build_model_from_fit_config(config), config=config, spectrum_output=config.spectrum_output)
