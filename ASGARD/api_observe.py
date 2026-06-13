from __future__ import annotations

from copy import deepcopy
from typing import Any, Optional

import numpy as np

from asgard_core.asgard_state import (
    SolveState,
    _build_observer_setup_from_state,
    _forward_synchrotron_absorption_transfer,
    _project_component,
    project_flux_grid,
)
from src.Electron.electron_radiation import electron_radiation_kernel as electron_radiation_module
from asgard_core.asgard_config import FitConfig, SpectrumOutputConfig
from asgard_core.asgard_observables import OUTPUT_BANDS, build_multiband_observer_frequencies, combine_multiband_flux
from asgard_core.asgard_postprocess import (
    build_spectrum_dataset_names,
    build_spectrum_frequency_grid,
    compute_light_curve_redchi,
    compute_spectrum_redchi,
    select_spectrum_time_index,
)
from src import Interpolation, constants

from .api_adaptive import _observe_parts, _observe_total
from .api_fit import FitResult
from .api_model import (
    GaussianJet,
    ISM,
    Jet,
    Magnetar,
    Model,
    Observer,
    PolarizationResult,
    PowerLawJet,
    Radiation,
    Setups,
    SkyImage,
    StepPowerLawJet,
    TophatJet,
    TwoComponentJet,
    Wind,
    _angular_separation,
    _build_fit_config_for_patch,
    _extract_pair_flux,
    _iter_patch_elements,
    _iter_patches,
    _solve_patch_state,
)


def _direct_total(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    timings: Optional[dict[str, float]] = None,
    projection_kind: str = "lightcurve",
) -> np.ndarray:
    config = _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )
    state = _solve_patch_state(model, config, times_s, nu_hz, timings=timings)
    return _observe_total(state, times_s, nu_hz, timings=timings, projection_kind=projection_kind)


def _patch_total(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    timings: Optional[dict[str, float]] = None,
    projection_kind: str = "lightcurve",
) -> np.ndarray:
    total = np.zeros((nu_hz.shape[0], times_s.shape[0]), dtype=float)
    for phi_center, theta_center, patch_half_angle in _iter_patches(model):
        e_iso = model.jet.energy_iso(phi_center, theta_center)
        gamma0 = model.jet.gamma0(phi_center, theta_center)
        if e_iso <= 0.0 or gamma0 <= 1.0:
            continue
        theta_v = _angular_separation(theta_center, phi_center, model.observer.theta_obs, model.observer.phi_obs)
        config = _build_fit_config_for_patch(
            model,
            phi_center=phi_center,
            theta_v=theta_v,
            opening_angle_jet=patch_half_angle,
            e_iso=e_iso,
            gamma0=gamma0,
            theta_center=theta_center,
        )
        state = _solve_patch_state(model, config, times_s, nu_hz, timings=timings)
        total += _observe_total(state, times_s, nu_hz, timings=timings, projection_kind=projection_kind)
    return total


def _total_matrix(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    timings: Optional[dict[str, float]] = None,
    projection_kind: str = "lightcurve",
) -> np.ndarray:
    times_s = np.asarray(times_s, dtype=float)
    nu_hz = np.asarray(nu_hz, dtype=float)
    if model.jet.kind == "tophat" and model._supports_direct_kernel():
        return _direct_total(model, times_s, nu_hz, timings=timings, projection_kind=projection_kind)
    return _patch_total(model, times_s, nu_hz, timings=timings, projection_kind=projection_kind)


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
    if magnetic_geometry == "toroidal" and not _jet_is_axisymmetric(model.jet):
        raise NotImplementedError("toroidal polarization currently requires an axisymmetric jet.")

    total_i = np.zeros((nu_hz.shape[0], times_s.shape[0]), dtype=float)
    total_q = np.zeros_like(total_i)
    total_u = np.zeros_like(total_i)
    components = {
        "fwd_sync": _empty_stokes(total_i),
        "rev_sync": _empty_stokes(total_i),
        "hadronic_sync": _empty_stokes(total_i),
        "rev_hadronic_sync": _empty_stokes(total_i),
    }
    sightline, sky_x_axis, sky_y_axis = _sky_basis(model.observer)
    active_patch_found = False

    for phi_center, theta_center, patch_half_angle, patch_solid_angle in _iter_patch_elements(model):
        e_iso = model.jet.energy_iso(phi_center, theta_center)
        gamma0 = model.jet.gamma0(phi_center, theta_center)
        if e_iso <= 0.0 or gamma0 <= 1.0:
            continue
        active_patch_found = True
        patch_direction = _direction_vector(theta_center, phi_center)
        theta_v = _angular_separation(theta_center, phi_center, model.observer.theta_obs, model.observer.phi_obs)
        config = _build_fit_config_for_patch(
            model,
            phi_center=phi_center,
            theta_v=theta_v,
            opening_angle_jet=patch_half_angle,
            e_iso=e_iso,
            gamma0=gamma0,
            theta_center=theta_center,
        )
        state = _solve_patch_state(model, config, times_s, nu_hz)
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
            _shock_random_anisotropy(magnetic_geometry, theta_v, state.components.fwd.gamma),
            patch_solid_angle,
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
                _shock_random_anisotropy(magnetic_geometry, theta_v, _reverse_shock_gamma_array(state)),
                patch_solid_angle,
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
                    _shock_random_anisotropy(magnetic_geometry, theta_v, _reverse_shock_gamma_array(state)),
                    patch_solid_angle,
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
                _shock_random_anisotropy(magnetic_geometry, theta_v, state.components.fwd.gamma),
                patch_solid_angle,
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
        np.asarray(state.reverse_emission.magnetic_field_g, dtype=float),
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


def _project_surface_element(
    setup,
    characteristic_time_s: np.ndarray,
    gamma: np.ndarray,
    radius_cm: np.ndarray,
    seed_frequency_hz: np.ndarray,
    absorbed_spectral_flux: np.ndarray,
    frequencies_hz: np.ndarray,
    patch_solid_angle: float,
) -> np.ndarray:
    """按真实球面面元dOmega投影一个patch中心代表的同步辐射。"""
    if not np.any(absorbed_spectral_flux):
        return np.zeros((frequencies_hz.shape[0], setup.observer_time_s.shape[0]), dtype=float)
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    order = np.argsort(frequencies_hz)
    sorted_frequencies = frequencies_hz[order]
    flux_sorted = Interpolation.sed_interpolation_surface_element(
        setup.boundary,
        characteristic_time_s,
        gamma,
        radius_cm,
        absorbed_spectral_flux,
        seed_frequency_hz,
        sorted_frequencies,
        setup.observer_time_s,
        float(patch_solid_angle),
    )
    if np.array_equal(order, np.arange(order.shape[0])):
        return flux_sorted
    flux_matrix = np.empty_like(flux_sorted)
    flux_matrix[order] = flux_sorted
    return flux_matrix


def _patch_hadronic_synchrotron_component(state: SolveState) -> tuple[np.ndarray, np.ndarray] | None:
    if state.hadronic is None:
        return None
    from src import constants
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
    from src import constants
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
    angular_diameter_distance_cm = _angular_diameter_distance_cm(model.observer)
    sightline, sky_x_axis, sky_y_axis = _sky_basis(model.observer)
    frequencies_hz = np.array([nu_obs], dtype=float)
    required_half_fov = 0.5 * float(fov)
    total_solid_angle = 2.0 * np.pi * (1.0 - np.cos(float(model.jet.theta_max)))
    collapse_phi = _can_collapse_sky_image_phi(model)
    for phi_center, theta_center, patch_half_angle, domega in _iter_img_patches(
        model,
        npixel,
        collapse_phi=collapse_phi,
    ):
        e_iso = model.jet.energy_iso(phi_center, theta_center)
        gamma0 = model.jet.gamma0(phi_center, theta_center)
        if e_iso <= 0.0 or gamma0 <= 1.0:
            continue
        theta_v = _angular_separation(theta_center, phi_center, model.observer.theta_obs, model.observer.phi_obs)
        config = _build_fit_config_for_patch(
            model,
            phi_center=phi_center,
            theta_v=theta_v,
            opening_angle_jet=model.jet.theta_max,
            e_iso=e_iso,
            gamma0=gamma0,
            theta_center=theta_center,
        )
        state = _solve_patch_state(model, config, times_s, frequencies_hz)
        observed = _observe_parts(state, times_s, frequencies_hz)
        patch_flux = np.asarray(observed.total[0, :], dtype=float)
        radius_cm = _interpolate_positive_series(
            state.components.fwd.characteristic_time_s,
            state.components.fwd.radius_cm,
            times_s,
        )
        if not np.any(np.isfinite(patch_flux) & (patch_flux > 0.0)):
            continue

        x_center, y_center = _project_patch_to_sky(
            radius_cm,
            theta_center,
            phi_center,
            sightline,
            sky_x_axis,
            sky_y_axis,
            angular_diameter_distance_cm,
        )
        sigma = np.maximum(
            radius_cm * np.sin(max(patch_half_angle, 1.0e-12)) / angular_diameter_distance_cm / 2.0,
            0.5 * base_pixel_size,
        )
        patch_weight = domega / total_solid_angle if total_solid_angle > 0.0 else 0.0
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
    x_centroid = np.zeros(times_s.shape[0], dtype=float)
    y_centroid = np.zeros(times_s.shape[0], dtype=float)
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


def _can_collapse_sky_image_phi(model: Model) -> bool:
    return model.jet.kind == "tophat" and abs(float(model.observer.theta_obs)) <= 1.0e-12


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


def _jet_is_axisymmetric(jet: Jet) -> bool:
    return jet.kind in ("tophat", "gaussian", "powerlaw", "twocomponent", "steppowerlaw")


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


def _angular_diameter_distance_cm(observer: Observer) -> float:
    if observer.lumi_dist_cm is None or observer.lumi_dist_cm <= 0.0:
        raise ValueError("Observer.lumi_dist_cm must be set for sky_image().")
    return observer.lumi_dist_cm / (1.0 + observer.z) ** 2


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
    try:
        from scipy.special import erf as _erf  # type: ignore
    except Exception:
        from math import erf as _scalar_erf

        _erf = np.vectorize(_scalar_erf)

    sigma = np.maximum(np.asarray(sigma, dtype=float), 1.0e-30)
    x_axis = np.asarray(x_axis, dtype=float)
    y_axis = np.asarray(y_axis, dtype=float)
    x_edges = np.concatenate((x_axis - 0.5 * pixel_size, np.array([x_axis[-1] + 0.5 * pixel_size], dtype=float)))
    y_edges = np.concatenate((y_axis - 0.5 * pixel_size, np.array([y_axis[-1] + 0.5 * pixel_size], dtype=float)))
    sigma_view = np.sqrt(2.0) * sigma[:, None]
    x_frac = 0.5 * (
        _erf((x_edges[1:][None, :] - np.asarray(x_center, dtype=float)[:, None]) / sigma_view)
        - _erf((x_edges[:-1][None, :] - np.asarray(x_center, dtype=float)[:, None]) / sigma_view)
    )
    y_frac = 0.5 * (
        _erf((y_edges[1:][None, :] - np.asarray(y_center, dtype=float)[:, None]) / sigma_view)
        - _erf((y_edges[:-1][None, :] - np.asarray(y_center, dtype=float)[:, None]) / sigma_view)
    )
    return x_frac[:, :, None] * y_frac[:, None, :] / (pixel_size * pixel_size)


def _interpolate_positive_series(source_t: np.ndarray, source_y: np.ndarray, target_t: np.ndarray) -> np.ndarray:
    source_t = np.asarray(source_t, dtype=float)
    source_y = np.asarray(source_y, dtype=float)
    target_t = np.asarray(target_t, dtype=float)
    if np.all(source_t > 0.0) and np.all(source_y > 0.0) and np.all(target_t > 0.0):
        return np.exp(
            np.interp(
                np.log(target_t),
                np.log(source_t),
                np.log(source_y),
                left=np.log(source_y[0]),
                right=np.log(source_y[-1]),
            )
        )
    return np.interp(target_t, source_t, source_y, left=source_y[0], right=source_y[-1])


def _set_dotted_attr(obj: Any, path: str, value: Any) -> None:
    target = obj
    parts = path.split(".")
    for name in parts[:-1]:
        target = getattr(target, name)
    setattr(target, parts[-1], value)


def _evaluate_flux_observations(model: Model, times_s: np.ndarray, frequencies_hz: np.ndarray) -> np.ndarray:
    times_s = np.asarray(times_s, dtype=float)
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    if frequencies_hz.ndim != 1 or times_s.ndim != 1:
        raise ValueError("Flux-density observations must be one-dimensional arrays.")
    if frequencies_hz.shape == times_s.shape:
        unique_freqs, inverse = np.unique(frequencies_hz, return_inverse=True)
        grid = model.flux_density_grid(times_s, unique_freqs).total
        return grid[inverse, np.arange(times_s.shape[0])]
    return model.flux_density_grid(times_s, frequencies_hz).total


def _param_path(model: Model, param: ParamDef) -> str:
    if param.path is not None:
        return param.path
    name = param.name.lower()
    alias_map = {
        "e_iso": "jet.E_iso",
        "log10_eiso": "jet.E_iso",
        "log10_e_iso": "jet.E_iso",
        "gamma0": "jet.lf",
        "log10_gamma0": "jet.lf",
        "theta_c": "jet.theta_j" if model.jet.kind == "tophat" else "jet.theta_c",
        "theta_j": "jet.theta_j",
        "theta_obs": "observer.theta_obs",
        "z": "observer.z",
        "lumi_dist": "observer.lumi_dist_cm",
        "lumi_dist_cm": "observer.lumi_dist_cm",
        "eps_e": "fwd_rad.eps_e",
        "epsilon_e": "fwd_rad.eps_e",
        "eps_b": "fwd_rad.eps_B",
        "epsilon_b": "fwd_rad.eps_B",
        "p": "fwd_rad.p",
        "xi_n": "fwd_rad.xi_N",
        "f_e": "fwd_rad.xi_N",
        "n0": "medium.n_ism",
        "n_ism": "medium.n_ism",
        "d_ne": "medium.n_ism",
        "a_star": "medium.A_star",
        "astar": "medium.A_star",
        "e_iso_c": "jet.E_iso_c",
        "e_iso_n": "jet.E_iso_n",
        "e_iso_outer": "jet.E_iso_w",
        "e_iso_w": "jet.E_iso_w",
        "gamma0_c": "jet.lf_c",
        "gamma0_n": "jet.lf_n",
        "gamma0_outer": "jet.lf_w",
        "gamma0_w": "jet.lf_w",
        "theta_n": "jet.theta_n",
        "theta_o": "jet.theta_w",
        "theta_w": "jet.theta_w",
        "k": "jet.k_e",
        "k_e": "jet.k_e",
        "k_g": "jet.k_g",
        "duration": "jet.duration",
        "l0": "jet.magnetar.L0",
        "t0": "jet.magnetar.t0",
        "q": "jet.magnetar.q",
        "eps_e_r": "rvs_rad.eps_e",
        "epsilon_e_r": "rvs_rad.eps_e",
        "eps_b_r": "rvs_rad.eps_B",
        "epsilon_b_r": "rvs_rad.eps_B",
        "p_r": "rvs_rad.p",
        "xi_n_r": "rvs_rad.xi_N",
        "f_e_r": "rvs_rad.xi_N",
        "reverse_sigma": "setups.reverse_sigma",
    }
    if name not in alias_map:
        raise KeyError(f"Cannot infer parameter path for {param.name}.")
    return alias_map[name]


def _as_model(cfg: Any) -> Model:
    if isinstance(cfg, Model):
        return cfg
    if cfg is None:
        raise ValueError("Either a Model or cfg must be provided.")
    if isinstance(cfg, Setups):
        cfg = cfg.__dict__.copy()
    if not isinstance(cfg, dict):
        raise TypeError("cfg must be a Model or a dictionary of model options.")
    setups_source = cfg.get("setups", Setups())
    setups = deepcopy(setups_source if isinstance(setups_source, Setups) else Setups(**setups_source))
    if "reverse_sigma" in cfg:
        setups.reverse_sigma = float(cfg["reverse_sigma"])
    observer = cfg.get(
        "observer",
        Observer(
            z=cfg.get("z", setups.z),
            theta_obs=cfg.get("theta_obs", setups.theta_obs),
            phi_obs=cfg.get("phi_obs", setups.phi_obs),
            lumi_dist=cfg.get("lumi_dist", setups.lumi_dist),
            lumi_dist_cm=cfg.get("lumi_dist_cm"),
        ),
    )
    medium = cfg.get("medium")
    medium_kind = str(cfg.get("medium_type", cfg.get("medium_name", setups.medium or "ism"))).lower()
    if isinstance(medium, str):
        medium_kind = medium.lower()
        medium = None
    if medium is None:
        if medium_kind == "wind":
            medium = Wind(A_star=cfg.get("A_star", cfg.get("Astar", 1.0)), n0=cfg.get("n0", cfg.get("n_ism", 0.1)))
        else:
            medium = ISM(n0=cfg.get("n0", cfg.get("n_ism", 0.1)))

    jet = cfg.get("jet")
    jet_kind = str(cfg.get("jet_type", cfg.get("jet_name", setups.jet or "tophat"))).lower()
    if isinstance(jet, str):
        jet_kind = jet.lower()
        jet = None
    if jet is None:
        magnetar = None
        if "magnetar" in cfg and cfg["magnetar"] is not None:
            source = cfg["magnetar"]
            if isinstance(source, Magnetar):
                magnetar = source
            else:
                magnetar = Magnetar(L0=source["L0"], t0=source["t0"], q=source["q"])
        elif {"L0", "t0", "q"} <= set(cfg.keys()):
            magnetar = Magnetar(L0=cfg["L0"], t0=cfg["t0"], q=cfg["q"])
        if jet_kind == "gaussian":
            jet = GaussianJet(
                E_iso=cfg["E_iso"],
                Gamma0=cfg["Gamma0"],
                theta_c=cfg["theta_c"],
                theta_max=cfg.get("theta_max", 0.6),
                duration=cfg.get("duration"),
                magnetar=magnetar,
                spreading=cfg.get("spreading", False),
            )
        elif jet_kind == "powerlaw":
            jet = PowerLawJet(
                E_iso=cfg["E_iso"],
                Gamma0=cfg["Gamma0"],
                theta_c=cfg["theta_c"],
                k=cfg.get("k"),
                k_e=cfg.get("k_e"),
                k_g=cfg.get("k_g"),
                theta_max=cfg.get("theta_max", np.pi / 2.0),
                duration=cfg.get("duration"),
                magnetar=magnetar,
                spreading=cfg.get("spreading", False),
            )
        elif jet_kind == "twocomponent":
            jet = TwoComponentJet(
                E_iso_c=cfg["E_iso_c"],
                Gamma0_c=cfg["Gamma0_c"],
                theta_c=cfg["theta_c"],
                E_iso_outer=cfg["E_iso_outer"],
                Gamma0_outer=cfg["Gamma0_outer"],
                theta_o=cfg["theta_o"],
                duration=cfg.get("duration"),
                magnetar=magnetar,
                spreading=cfg.get("spreading", False),
            )
        elif jet_kind == "steppowerlaw":
            jet = StepPowerLawJet(
                E_iso_c=cfg["E_iso_c"],
                Gamma0_c=cfg["Gamma0_c"],
                theta_c=cfg["theta_c"],
                E_iso_w=cfg["E_iso_w"],
                Gamma0_w=cfg["Gamma0_w"],
                theta_w=cfg["theta_w"],
                k=cfg.get("k", 2.0),
                k_e=cfg.get("k_e"),
                k_g=cfg.get("k_g"),
                duration=cfg.get("duration"),
                magnetar=magnetar,
                spreading=cfg.get("spreading", False),
            )
        else:
            jet = TophatJet(
                E_iso=cfg["E_iso"],
                Gamma0=cfg["Gamma0"],
                theta_c=cfg["theta_c"],
                duration=cfg.get("duration"),
                magnetar=magnetar,
                spreading=cfg.get("spreading", False),
            )

    fwd_rad = cfg.get(
        "fwd_rad",
        Radiation(
            eps_e=cfg.get("eps_e", cfg.get("epsilon_e", 1.0e-1)),
            eps_B=cfg.get("eps_B", cfg.get("epsilon_B", 1.0e-3)),
            epsilon_b_floor=cfg.get("eps_B_floor", cfg.get("epsilon_B_floor")),
            magnetic_decay_alpha_t=cfg.get("magnetic_decay_alpha_t", 0.0),
            magnetic_decay_t0_s=cfg.get("magnetic_decay_t0_s", 1.0),
            p=cfg.get("p", 2.5),
            xi_N=cfg.get("xi_N", cfg.get("f_e", 1.0e-1)),
            thermal_electrons=cfg.get("thermal_electrons", False),
            ssc=cfg.get("ssc", setups.fwd_ssc),
            kn=cfg.get("kn", setups.kn),
        ),
    )
    rvs_rad = cfg.get("rvs_rad")
    if rvs_rad is None and cfg.get("reverse", setups.rvs_shock):
        rvs_rad = Radiation(
            eps_e=cfg.get("eps_e_r", cfg.get("eps_e", 1.0e-1)),
            eps_B=cfg.get("eps_B_r", cfg.get("eps_B", 1.0e-2)),
            p=cfg.get("p_r", cfg.get("p", 2.4)),
            xi_N=cfg.get("xi_N_r", cfg.get("f_e_r", 1.0)),
            thermal_electrons=cfg.get("thermal_electrons_r", False),
            ssc=cfg.get("rvs_ssc", setups.rvs_ssc),
            kn=cfg.get("kn", setups.kn),
        )

    resolutions = cfg.get("resolutions")
    return Model(medium=medium, jet=jet, observer=observer, fwd_rad=fwd_rad, rvs_rad=rvs_rad, setups=setups, resolutions=resolutions)


def observe(
    model: Model,
    config: Optional[FitConfig] = None,
    spectrum_output: Optional[SpectrumOutputConfig] = None,
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

    rs_nu_m = None
    rs_nu_c = None
    rs_nu_a = None
    if details.rev is not None:
        rs_nu_m = details.rev.nu_m
        rs_nu_c = details.rev.nu_c
        rs_nu_a = details.rev.nu_a

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
        if spectrum_output.time_s is not None or len(spectrum_output.dataset_names) > 0:
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
        nu_m=details.fwd.nu_m,
        nu_c=details.fwd.nu_c,
        nu_a=details.fwd.nu_a,
        rs_nu_m=rs_nu_m,
        rs_nu_c=rs_nu_c,
        rs_nu_a=rs_nu_a,
        spectrum_time_s=spectrum_time_s,
        spectrum_freq_hz=spectrum_freq_hz,
        spectrum_fnu=spectrum_fnu,
    )


def _build_model_from_fit_config(config: FitConfig) -> Model:
    ssc_enabled = config.index_y != 0
    kn_enabled = config.index_y == 1
    reverse = getattr(config, "reverse_shock", None)
    reverse_enabled = bool(reverse and reverse.enabled)

    if config.a_star > 0.0:
        medium = Wind(A_star=config.a_star, n0=config.d_ne)
    else:
        medium = ISM(n_ism=config.d_ne)

    reverse_delta_t_s = 10.0
    if reverse and reverse.delta_t_s is not None:
        reverse_delta_t_s = float(reverse.delta_t_s)

    rvs_rad = None
    if reverse_enabled:
        rvs_rad = Radiation(
            eps_e=reverse.epsilon_e if reverse.epsilon_e is not None else config.epsilon_e,
            eps_B=reverse.epsilon_b if reverse.epsilon_b is not None else config.epsilon_b,
            p=reverse.p if reverse.p is not None else config.p,
            xi_N=reverse.f_e if reverse.f_e is not None else config.f_e,
            ssc=bool(reverse.include_ssc),
            kn=kn_enabled,
        )

    return Model(
        jet=TophatJet(E_iso=config.e_iso, Gamma0=config.eta_0, theta_j=config.opening_angle_jet, duration=reverse_delta_t_s if reverse_enabled else None),
        medium=medium,
        observer=Observer(z=config.z, theta_obs=config.theta_v, lumi_dist_cm=config.luminosity_distance_cm_override),
        fwd_rad=Radiation(
            eps_e=config.epsilon_e,
            eps_B=config.epsilon_b,
            epsilon_b_floor=config.epsilon_b_floor,
            magnetic_decay_alpha_t=config.magnetic_decay_alpha_t,
            magnetic_decay_t0_s=config.magnetic_decay_t0_s,
            p=config.p,
            xi_N=config.f_e,
            thermal_electrons=config.thermal_electrons,
            ssc=ssc_enabled,
            kn=kn_enabled,
        ),
        rvs_rad=rvs_rad,
        setups=Setups(
            z=config.z,
            theta_obs=config.theta_v,
            rvs_shock=reverse_enabled,
            fwd_ssc=ssc_enabled,
            rvs_ssc=bool(reverse.include_ssc) if reverse_enabled else False,
            ssc_cooling=ssc_enabled,
            kn=kn_enabled,
            num_threads=config.num_threads,
            num_gam_e=config.num_gam_e,
            num_nu=config.num_nu,
            num_r=config.num_r,
            num_theta=config.num_theta,
            num_phi=config.num_phi,
            num_tobs=config.num_tobs,
            observer_time_min_s=10 ** config.t_obs_min_log10,
            observer_time_max_s=10 ** config.t_obs_max_log10,
            initial_radius_cm=config.initial_radius_cm,
            reverse_delta_t_s=reverse_delta_t_s,
            reverse_sigma=float(reverse.sigma) if reverse_enabled else 0.0,
            include_cross_zone_ic=bool(reverse.include_cross_zone_ic) if reverse_enabled else False,
            weno5=config.weno5,
            electron_solver=config.electron_solver,
            cooling_kernel=config.cooling_kernel,
            radiation_kernel=config.radiation_kernel,
            dynamics_kernel=config.dynamics_kernel,
            geometry_kernel=config.geometry_kernel,
            electron_photon_coupling=config.electron_photon_coupling,
            structured_backend=config.structured_backend,
            structured_parallel_mode=config.structured_parallel_mode,
            structured_outer_threads=config.structured_outer_threads,
            structured_inner_threads=config.structured_inner_threads,
            num_chi=config.num_chi,
            fullhide2d_transport_model=config.fullhide2d_transport_model,
            fullhide2d_stochastic_accel_norm=config.fullhide2d_stochastic_accel_norm,
            fullhide2d_escape_mode=config.fullhide2d_escape_mode,
            electron_adaptive_substeps=config.electron_adaptive_substeps,
            electron_substep_rtol=config.electron_substep_rtol,
            electron_substep_min=config.electron_substep_min,
            electron_substep_max=config.electron_substep_max,
            index_dyn=config.index_dyn,
            index_syn_integr=config.index_syn_integr,
        ),
    )


def run_fit(config: Optional[FitConfig] = None) -> FitResult:
    cfg = FitConfig() if config is None else config
    model = _build_model_from_fit_config(cfg)
    return observe(model, config=cfg, spectrum_output=cfg.spectrum_output)
