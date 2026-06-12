from __future__ import annotations

from copy import deepcopy
from typing import Any, Optional

import numpy as np

from asgard_core.asgard_state import (
    FluxComponents,
    SolveState,
    _build_observer_setup_from_state,
    make_query_setup,
    make_tgrid,
    _forward_synchrotron_absorption_transfer,
    _project_component,
    project_flux_grid,
    solve_state_from_setup,
)
from src.Electron.electron_radiation import electron_radiation_kernel as electron_radiation_module
from asgard_core.asgard_config import FitConfig, ReverseShockConfig, SpectrumOutputConfig
from asgard_core.asgard_config import HadronicConfig
from asgard_core.asgard_observables import OUTPUT_BANDS, build_multiband_observer_frequencies, combine_multiband_flux
from asgard_core.asgard_postprocess import (
    build_spectrum_dataset_names,
    build_spectrum_frequency_grid,
    compute_light_curve_redchi,
    compute_spectrum_redchi,
    select_spectrum_time_index,
)
from asgard_core.asgard_presets import build_baseline_config
from src import Interpolation, constants

from .api_adaptive import _observe_parts, _observe_total
from .api_fit import FitResult
from .api_model import (
    CharTrack,
    FluxPair,
    FluxResult,
    GaussianJet,
    ISM,
    Jet,
    Magnetar,
    Medium,
    Model,
    Observer,
    PolarizationResult,
    PowerLawJet,
    Radiation,
    Setups,
    SkyImage,
    StepPowerLawJet,
    TophatJet,
    TrackBundle,
    TwoComponentJet,
    Wind,
)

def _solve_patch_state(
    model: Model,
    config: FitConfig,
    times_s: np.ndarray,
    requested_frequencies_hz: np.ndarray | None,
    timings: Optional[dict[str, float]] = None,
    solve_reference_times_s: np.ndarray | None = None,
) -> SolveState:
    solve_reference = times_s if solve_reference_times_s is None else np.asarray(solve_reference_times_s, dtype=float)
    base_count = max(int(model.setups.num_tobs), int(np.unique(solve_reference).size))
    solve_times_s = make_tgrid(solve_reference, base_count)
    if solve_reference.size > 1:
        solve_t_min = min(float(model.setups.observer_time_min_s), float(np.min(solve_reference)))
        solve_t_max_requested = float(np.max(solve_reference))
        solve_count = max(base_count, model._detail_time_count(solve_t_min, solve_t_max_requested))
        if solve_t_max_requested <= solve_t_min:
            solve_t_max = max(float(model.setups.observer_time_max_s), solve_t_min * 10.0)
        else:
            log_t_min = np.log10(solve_t_min)
            log_t_max_requested = np.log10(solve_t_max_requested)
            if solve_count <= 2:
                log_t_max = log_t_max_requested
            else:
                log_step = (log_t_max_requested - log_t_min) / float(solve_count - 2)
                log_t_max = log_t_max_requested + log_step
            solve_t_max = 10.0**log_t_max
        solve_times_s = np.logspace(np.log10(solve_t_min), np.log10(solve_t_max), solve_count)
    setup = make_query_setup(config, solve_times_s, requested_frequencies_hz)
    state = solve_state_from_setup(
        config,
        setup,
        timings=timings,
        requested_frequencies_hz=requested_frequencies_hz,
    )
    return state


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


def _solve_direct_model(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    solve_reference_times_s: np.ndarray | None = None,
    projection_kind: str = "lightcurve",
) -> tuple[FluxResult, TrackBundle]:
    config = _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )
    state = _solve_patch_state(
        model,
        config,
        times_s,
        nu_hz,
        solve_reference_times_s=solve_reference_times_s,
    )
    observed = _observe_parts(state, times_s, nu_hz, projection_kind=projection_kind)
    details = _make_details(state.components, patches=[{"phi": 0.0, "theta": 0.0, "weight": 1.0}], state=state)
    return observed, details


def _evaluate_direct_details(model: Model, times_s: np.ndarray) -> TrackBundle:
    config = _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )
    state = _solve_patch_state(model, config, times_s, None)
    return _make_details(state.components, patches=[{"phi": 0.0, "theta": 0.0, "weight": 1.0}], state=state)


def _solve_patch_model(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    solve_reference_times_s: np.ndarray | None = None,
    projection_kind: str = "lightcurve",
) -> tuple[FluxResult, TrackBundle]:
    if str(getattr(model.setups, "structured_backend", "fortran_1d")).lower() != "python_patch":
        if solve_reference_times_s is not None:
            raise NotImplementedError("structured_backend='fortran_1d' does not yet support external solve_reference_times_s.")
        from asgard_core.structured_jet_kernel import solve_structured_jet_fortran

        return solve_structured_jet_fortran(model, times_s, nu_hz, _build_fit_config_for_patch)

    return _solve_patch_model_python(
        model,
        times_s,
        nu_hz,
        solve_reference_times_s=solve_reference_times_s,
        projection_kind=projection_kind,
    )


def _solve_patch_model_python(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    solve_reference_times_s: np.ndarray | None = None,
    projection_kind: str = "lightcurve",
) -> tuple[FluxResult, TrackBundle]:
    total = np.zeros((nu_hz.shape[0], times_s.shape[0]), dtype=float)
    fwd_sync_total = np.zeros_like(total)
    fwd_ssc_total = np.zeros_like(total)
    rev_sync_total = np.zeros_like(total)
    rev_ssc_total = np.zeros_like(total)
    cross_ic_total = np.zeros_like(total)
    patches_meta: list[dict[str, float]] = []
    details_fwd: Optional[CharTrack] = None
    details_rev: Optional[CharTrack] = None

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
        state = _solve_patch_state(
            model,
            config,
            times_s,
            nu_hz,
            solve_reference_times_s=solve_reference_times_s,
        )
        observed = _observe_parts(state, times_s, nu_hz, projection_kind=projection_kind)
        total += observed.total
        fwd_sync_total += observed.fwd.sync
        fwd_ssc_total += observed.fwd.ssc
        rev_sync_total += observed.rev.sync
        rev_ssc_total += observed.rev.ssc
        if observed.cross_ic is not None:
            cross_ic_total += observed.cross_ic
        patches_meta.append(
            {
                "phi": float(phi_center),
                "theta": float(theta_center),
                "theta_v": float(theta_v),
                "half_angle": float(patch_half_angle),
                "E_iso": float(e_iso),
                "Gamma0": float(gamma0),
            }
        )
        if details_fwd is None:
            details = _make_details(state.components, patches_meta, state=state)
            details_fwd = details.fwd
            details_rev = details.rev

    if details_fwd is None:
        raise ValueError("No active jet patches were found for the requested structured jet.")
    return (
        FluxResult(
            total=total,
            fwd=FluxPair(sync=fwd_sync_total, ssc=fwd_ssc_total),
            rev=FluxPair(sync=rev_sync_total, ssc=rev_ssc_total),
            cross_ic=cross_ic_total,
        ),
        TrackBundle(fwd=details_fwd, rev=details_rev, patches=patches_meta),
    )


def _patch_details(model: Model, times_s: np.ndarray) -> TrackBundle:
    patches_meta: list[dict[str, float]] = []
    first_component: Optional[FluxComponents] = None
    first_details: Optional[TrackBundle] = None

    for phi_center, theta_center, patch_half_angle in _iter_patches(model):
        e_iso = model.jet.energy_iso(phi_center, theta_center)
        gamma0 = model.jet.gamma0(phi_center, theta_center)
        if e_iso <= 0.0 or gamma0 <= 1.0:
            continue
        theta_v = _angular_separation(theta_center, phi_center, model.observer.theta_obs, model.observer.phi_obs)
        patches_meta.append(
            {
                "phi": float(phi_center),
                "theta": float(theta_center),
                "theta_v": float(theta_v),
                "half_angle": float(patch_half_angle),
                "E_iso": float(e_iso),
                "Gamma0": float(gamma0),
            }
        )
        if first_component is None:
            config = _build_fit_config_for_patch(
                model,
                phi_center=phi_center,
                theta_v=theta_v,
                opening_angle_jet=patch_half_angle,
                e_iso=e_iso,
                gamma0=gamma0,
                theta_center=theta_center,
            )
            first_state = _solve_patch_state(model, config, times_s, None)
            first_component = first_state.components
            first_details = _make_details(first_component, patches_meta, state=first_state)

    if first_component is None or first_details is None:
        raise ValueError("No active jet patches were found for the requested structured jet.")
    return first_details


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


def _build_fit_config_for_patch(
    model: Model,
    *,
    phi_center: float,
    theta_v: float,
    opening_angle_jet: float,
    e_iso: float,
    gamma0: float,
    theta_center: Optional[float] = None,
) -> FitConfig:
    if getattr(model.jet, "spreading", False):
        raise NotImplementedError("Jet spreading is not implemented in the current ASGARD backend.")
    reverse_delta_t = model.setups.reverse_delta_t_s
    if getattr(model.jet, "duration", None) is not None:
        reverse_delta_t = float(model.jet.duration)
    index_y = 0
    # `ssc` controls whether the SSC radiation component is emitted.
    # `ssc_cooling` controls whether IC cooling is included in the electron solver.
    if model.setups.ssc_cooling:
        index_y = 1 if model.fwd_rad.kn else 2
    config = build_baseline_config(
        num_threads=model.setups.num_threads,
        num_gam_e=model.setups.num_gam_e,
        num_nu=model.setups.num_nu,
        num_r=model.setups.num_r,
        num_theta=model.setups.num_theta,
        num_phi=model.setups.num_phi,
        num_tobs=model.setups.num_tobs,
        t_obs_min_log10=float(np.log10(model.setups.observer_time_min_s)),
        t_obs_max_log10=float(np.log10(model.setups.observer_time_max_s)),
        index_dyn=model.setups.index_dyn,
        index_syn_integr=model.setups.index_syn_integr,
        electron_solver=model.setups.electron_solver,
        cooling_kernel=model.setups.cooling_kernel,
        radiation_kernel=model.setups.radiation_kernel,
        dynamics_kernel=model.setups.dynamics_kernel,
        geometry_kernel=model.setups.geometry_kernel,
        electron_photon_coupling=model.setups.electron_photon_coupling,
        structured_backend=model.setups.structured_backend,
        structured_parallel_mode=model.setups.structured_parallel_mode,
        structured_outer_threads=model.setups.structured_outer_threads,
        structured_inner_threads=model.setups.structured_inner_threads,
        num_chi=model.setups.num_chi,
        fullhide2d_transport_model=model.setups.fullhide2d_transport_model,
        fullhide2d_stochastic_accel_norm=model.setups.fullhide2d_stochastic_accel_norm,
        fullhide2d_escape_mode=model.setups.fullhide2d_escape_mode,
        electron_adaptive_substeps=model.setups.electron_adaptive_substeps,
        electron_substep_rtol=model.setups.electron_substep_rtol,
        electron_substep_min=model.setups.electron_substep_min,
        electron_substep_max=model.setups.electron_substep_max,
        weno5=model.setups.weno5,
        z=model.observer.z,
        theta_v=theta_v,
        opening_angle_jet=opening_angle_jet,
        e_iso=e_iso,
        eta_0=gamma0,
        epsilon_e=model.fwd_rad.eps_e,
        epsilon_b=model.fwd_rad.eps_B,
        epsilon_b_floor=model.fwd_rad.epsilon_b_floor,
        magnetic_decay_alpha_t=model.fwd_rad.magnetic_decay_alpha_t,
        magnetic_decay_t0_s=model.fwd_rad.magnetic_decay_t0_s,
        p=model.fwd_rad.p,
        f_e=model.fwd_rad.xi_N,
        thermal_electrons=bool(model.fwd_rad.thermal_electrons),
        index_y=index_y,
        include_forward_ssc=model.fwd_rad.ssc,
        luminosity_distance_cm_override=model.observer.lumi_dist_cm,
        initial_radius_cm=model.setups.initial_radius_cm,
        spectrum_output=SpectrumOutputConfig(enabled=False),
    )
    config.hadronic = HadronicConfig(
        enabled=bool(model.setups.hadronic_enabled and model.fwd_rad.epsilon_p > 0.0),
        solver=str(model.setups.hadronic_solver),
        epsilon_p=float(model.fwd_rad.epsilon_p),
        p_p=float(model.fwd_rad.p),
        eta_acc=float(model.fwd_rad.eta_acc),
        num_gam_p=int(model.setups.num_gam_p),
        include_proton_synch=bool(model.fwd_rad.proton_synch),
        include_pg=bool(model.fwd_rad.pg),
        include_bethe_heitler=bool(model.fwd_rad.bethe_heitler),
        include_hadronic_inverse_compton=bool(model.fwd_rad.hadronic_inverse_compton),
        include_pp=bool(model.fwd_rad.pp),
        include_neutrino=bool(model.fwd_rad.neutrino),
        include_pair_production=bool(model.fwd_rad.pair_production),
        pgamma_scheme=str(model.setups.pgamma_scheme if model.setups.pgamma_scheme != "disabled" else model.fwd_rad.pgamma_scheme),
        num_nu_nu=int(model.setups.num_nu_nu),
        reverse_enabled=bool(model.setups.rvs_shock and model.fwd_rad.reverse_epsilon_p > 0.0),
        reverse_epsilon_p=float(model.fwd_rad.reverse_epsilon_p),
        pair_cascade_iterations=int(model.setups.pair_cascade_iterations),
        pp_model=int(getattr(model.fwd_rad, 'pp_model', -1)),
    )
    magnetar = getattr(model.jet, "magnetar", None)
    if magnetar is not None and _jet_magnetar_active(model.jet, 0.0 if theta_center is None else theta_center):
        config.l_inj_0 = float(magnetar.L0)
        config.e_inj_t2 = float(magnetar.t0)
        config.q_inj = float(magnetar.q)
    config.reverse = model.rvs_rad is not None
    if model.rvs_rad is not None:
        config.reverse_shock = ReverseShockConfig(
            enabled=True,
            delta_t_s=reverse_delta_t,
            epsilon_e=model.rvs_rad.eps_e,
            epsilon_b=model.rvs_rad.eps_B,
            p=model.rvs_rad.p,
            f_e=model.rvs_rad.xi_N,
            sigma=model.setups.reverse_sigma,
            include_ssc=model.rvs_rad.ssc,
            include_cross_zone_ic=model.setups.include_cross_zone_ic,
        )
    kernel_medium = _project_medium_to_kernel(
        model.medium,
        phi_center=phi_center,
        theta_center=0.0 if theta_center is None else theta_center,
    )
    for key, value in kernel_medium.items():
        setattr(config, key, value)
    if model.medium.kind == "ism" and float(model.setups.f_jump) != 1.0:
        config.r_tr = float(model.setups.r_tr)
        config.f_jump = float(model.setups.f_jump)
        config.f_wide = float(model.setups.f_wide)
    return config


def _make_details(
    components: FluxComponents,
    patches: list[dict[str, float]],
    state: Optional[SolveState] = None,
) -> TrackBundle:
    fwd_gamma_e = None if state is None else np.asarray(state.electron.gam_e, dtype=float)
    fwd_dnde = None if state is None else np.asarray(state.electron.d_n_gam_e, dtype=float)
    fwd_dnde_bh = None if state is None or state.electron.d_n_gam_e_bh is None else np.asarray(state.electron.d_n_gam_e_bh, dtype=float)
    if state is not None and state.hadronic is not None and state.hadronic.d_n_gam_e_bh is not None:
        fwd_dnde_bh = np.asarray(state.hadronic.d_n_gam_e_bh, dtype=float)
    fwd_dnde_chi = None if state is None or state.electron.d_n_gam_e_chi is None else np.asarray(state.electron.d_n_gam_e_chi, dtype=float)
    fwd_chi_grid = None if state is None or state.electron.chi_grid is None else np.asarray(state.electron.chi_grid, dtype=float)
    fwd_lsyn_chi = None if state is None or state.electron.l_syn_spec_chi is None else np.asarray(state.electron.l_syn_spec_chi, dtype=float)
    fwd_seed_chi = None if state is None or state.electron.seed_syn_chi is None else np.asarray(state.electron.seed_syn_chi, dtype=float)
    fwd_tau_chi = None if state is None or state.electron.tau_syn_chi is None else np.asarray(state.electron.tau_syn_chi, dtype=float)
    fwd_chi_radius = None if state is None or state.electron.chi_radius_cm is None else np.asarray(state.electron.chi_radius_cm, dtype=float)
    fwd_chi_gamma = None if state is None or state.electron.chi_gamma_bulk is None else np.asarray(state.electron.chi_gamma_bulk, dtype=float)
    fwd_chi_weight = None if state is None or state.electron.chi_dvolume_weight is None else np.asarray(state.electron.chi_dvolume_weight, dtype=float)
    fwd_gamma_p = None if state is None or state.hadronic is None else np.asarray(state.hadronic.gam_p, dtype=float)
    fwd_dndp = None if state is None or state.hadronic is None else np.asarray(state.hadronic.d_n_gam_p, dtype=float)
    fwd_gamma_secondary = None if state is None or state.hadronic is None or state.hadronic.gam_secondary is None else np.asarray(state.hadronic.gam_secondary, dtype=float)
    fwd_dndn = None if state is None or state.hadronic is None or state.hadronic.d_n_gam_n is None else np.asarray(state.hadronic.d_n_gam_n, dtype=float)
    fwd_dndpi_plus = None if state is None or state.hadronic is None or state.hadronic.d_n_gam_pi_plus is None else np.asarray(state.hadronic.d_n_gam_pi_plus, dtype=float)
    fwd_dndpi_minus = None if state is None or state.hadronic is None or state.hadronic.d_n_gam_pi_minus is None else np.asarray(state.hadronic.d_n_gam_pi_minus, dtype=float)
    fwd_dndmu_ml = None if state is None or state.hadronic is None or state.hadronic.d_n_gam_mu_minus_left is None else np.asarray(state.hadronic.d_n_gam_mu_minus_left, dtype=float)
    fwd_dndmu_mr = None if state is None or state.hadronic is None or state.hadronic.d_n_gam_mu_minus_right is None else np.asarray(state.hadronic.d_n_gam_mu_minus_right, dtype=float)
    fwd_dndmu_pl = None if state is None or state.hadronic is None or state.hadronic.d_n_gam_mu_plus_left is None else np.asarray(state.hadronic.d_n_gam_mu_plus_left, dtype=float)
    fwd_dndmu_pr = None if state is None or state.hadronic is None or state.hadronic.d_n_gam_mu_plus_right is None else np.asarray(state.hadronic.d_n_gam_mu_plus_right, dtype=float)
    fwd_had_gamma = None
    if state is not None and state.hadronic is not None:
        fwd_had_gamma = np.asarray(state.hadronic.l_had_syn_spec + state.hadronic.l_had_pg_gamma, dtype=float)
        if state.hadronic.l_had_pion_synch is not None:
            fwd_had_gamma = fwd_had_gamma + np.asarray(state.hadronic.l_had_pion_synch, dtype=float)
        if state.hadronic.l_had_muon_synch is not None:
            fwd_had_gamma = fwd_had_gamma + np.asarray(state.hadronic.l_had_muon_synch, dtype=float)
        if state.hadronic.l_had_pion_inverse_compton is not None:
            fwd_had_gamma = fwd_had_gamma + np.asarray(state.hadronic.l_had_pion_inverse_compton, dtype=float)
        if state.hadronic.l_had_muon_inverse_compton is not None:
            fwd_had_gamma = fwd_had_gamma + np.asarray(state.hadronic.l_had_muon_inverse_compton, dtype=float)
    fwd_had_pi_syn = None if state is None or state.hadronic is None or state.hadronic.l_had_pion_synch is None else np.asarray(state.hadronic.l_had_pion_synch, dtype=float)
    fwd_had_mu_syn = None if state is None or state.hadronic is None or state.hadronic.l_had_muon_synch is None else np.asarray(state.hadronic.l_had_muon_synch, dtype=float)
    fwd_had_pi_ic = None if state is None or state.hadronic is None or state.hadronic.l_had_pion_inverse_compton is None else np.asarray(state.hadronic.l_had_pion_inverse_compton, dtype=float)
    fwd_had_mu_ic = None if state is None or state.hadronic is None or state.hadronic.l_had_muon_inverse_compton is None else np.asarray(state.hadronic.l_had_muon_inverse_compton, dtype=float)
    fwd_nu_freq = None
    fwd_nu_lum = None
    if state is not None and state.hadronic is not None and state.config.hadronic.include_neutrino:
        fwd_nu_freq = np.asarray(state.hadronic.neutrino_frequency_hz, dtype=float)
        fwd_nu_lum = np.asarray(state.hadronic.neutrino_luminosity, dtype=float)
    fwd_had_syn = None
    fwd_had_pg_gamma = None
    fwd_had_bh = None
    fwd_had_hic = None
    fwd_am3_power = None
    fwd_tau_pg = None
    fwd_tau_bh = None
    fwd_pg_survival = None
    fwd_timings = None
    fwd_seed_freq = None
    if state is not None and state.hadronic is not None:
        fwd_seed_freq = np.asarray(state.photon_field.seed_frequency_hz, dtype=float)
        fwd_had_syn = np.asarray(state.hadronic.l_had_syn_spec, dtype=float)
        fwd_had_pg_gamma = np.asarray(state.hadronic.l_had_pg_gamma, dtype=float)
        if state.hadronic.l_had_bethe_heitler is not None:
            fwd_had_bh = np.asarray(state.hadronic.l_had_bethe_heitler, dtype=float)
        if state.hadronic.l_had_hadronic_inverse_compton is not None:
            fwd_had_hic = np.asarray(state.hadronic.l_had_hadronic_inverse_compton, dtype=float)
        if state.hadronic.am3_process_power is not None:
            fwd_am3_power = np.asarray(state.hadronic.am3_process_power, dtype=float)
        if state.hadronic.tau_pg is not None:
            fwd_tau_pg = np.asarray(state.hadronic.tau_pg, dtype=float)
        if state.hadronic.tau_bh is not None:
            fwd_tau_bh = np.asarray(state.hadronic.tau_bh, dtype=float)
        if state.hadronic.pg_photon_survival is not None:
            fwd_pg_survival = np.asarray(state.hadronic.pg_photon_survival, dtype=float)
        fwd_timings = dict(state.hadronic.timings) if state.hadronic.timings else {}
    rev_gamma_e = None
    rev_dnde = None
    if state is not None and state.dynamics.reverse_shock is not None:
        rev_gamma_e = np.asarray(state.dynamics.reverse_shock.gam_e, dtype=float)
        rev_dnde = np.asarray(state.dynamics.reverse_shock.d_n_gam_e, dtype=float)
    return TrackBundle(
        fwd=CharTrack(
            t_obs=components.fwd.characteristic_time_s,
            radius=components.fwd.radius_cm,
            Gamma=components.fwd.gamma,
            N_p=np.asarray(components.fwd.swept_mass_g, dtype=float) / constants.para_m_p,
            Doppler=components.fwd.doppler,
            B_comv=components.fwd.magnetic_field_g,
            nu_m=components.fwd.nu_m,
            nu_c=components.fwd.nu_c,
            nu_a=components.fwd.nu_a,
            nu_M=components.fwd.nu_M,
            cooling_timescale_s=components.fwd.cooling_timescale_s,
            dynamical_timescale_s=components.fwd.dynamical_timescale_s,
            gamma_e=fwd_gamma_e,
            dN_dgamma_e=fwd_dnde,
            dN_dgamma_e_bh=fwd_dnde_bh,
            dN_dgamma_e_chi=fwd_dnde_chi,
            chi_grid=fwd_chi_grid,
            l_syn_spec_chi=fwd_lsyn_chi,
            seed_syn_chi=fwd_seed_chi,
            tau_syn_chi=fwd_tau_chi,
            chi_radius_cm=fwd_chi_radius,
            chi_gamma_bulk=fwd_chi_gamma,
            chi_dvolume_weight=fwd_chi_weight,
            gamma_p=fwd_gamma_p,
            dN_dgamma_p=fwd_dndp,
            gamma_secondary=fwd_gamma_secondary,
            dN_dgamma_n=fwd_dndn,
            dN_dgamma_pi_plus=fwd_dndpi_plus,
            dN_dgamma_pi_minus=fwd_dndpi_minus,
            dN_dgamma_mu_minus_left=fwd_dndmu_ml,
            dN_dgamma_mu_minus_right=fwd_dndmu_mr,
            dN_dgamma_mu_plus_left=fwd_dndmu_pl,
            dN_dgamma_mu_plus_right=fwd_dndmu_pr,
            hadronic_gamma=fwd_had_gamma,
            hadronic_pion_synch=fwd_had_pi_syn,
            hadronic_muon_synch=fwd_had_mu_syn,
            hadronic_pion_inverse_compton=fwd_had_pi_ic,
            hadronic_muon_inverse_compton=fwd_had_mu_ic,
            neutrino_frequency_hz=fwd_nu_freq,
            neutrino_luminosity=fwd_nu_lum,
            seed_frequency_hz=fwd_seed_freq,
            l_had_syn_spec=fwd_had_syn,
            l_had_pg_gamma_spec=fwd_had_pg_gamma,
            l_had_bethe_heitler_spec=fwd_had_bh,
            l_had_hadronic_ic_spec=fwd_had_hic,
            am3_process_power=fwd_am3_power,
            tau_pg=fwd_tau_pg,
            tau_bh=fwd_tau_bh,
            pg_photon_survival=fwd_pg_survival,
            timings=fwd_timings,
        ),
        rev=None
        if components.rev is None
        else CharTrack(
            t_obs=components.rev.characteristic_time_s,
            radius=components.rev.radius_cm,
            Gamma=components.rev.gamma,
            N_p=np.asarray(components.rev.swept_mass_g, dtype=float) / constants.para_m_p,
            Doppler=components.rev.doppler,
            B_comv=components.rev.magnetic_field_g,
            nu_m=components.rev.nu_m,
            nu_c=components.rev.nu_c,
            nu_a=components.rev.nu_a,
            nu_M=components.rev.nu_M,
            cooling_timescale_s=components.rev.cooling_timescale_s,
            dynamical_timescale_s=components.rev.dynamical_timescale_s,
            gamma_e=rev_gamma_e,
            dN_dgamma_e=rev_dnde,
        ),
        patches=patches,
    )


def _iter_patches(model: Model):
    for phi_center, theta_center, patch_half_angle, _domega in _iter_patch_elements(model):
        yield phi_center, theta_center, patch_half_angle


def _iter_patch_elements(model: Model):
    """生成真实球面经纬面元：中心方向、等面积圆半角和dOmega。"""
    theta_edges = np.linspace(0.0, model.jet.theta_max, model.setups.patch_theta + 1)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, model.setups.patch_phi + 1)
    for i_theta in range(model.setups.patch_theta):
        theta1 = theta_edges[i_theta]
        theta2 = theta_edges[i_theta + 1]
        theta_center = 0.5 * (theta1 + theta2)
        for i_phi in range(model.setups.patch_phi):
            phi1 = phi_edges[i_phi]
            phi2 = phi_edges[i_phi + 1]
            phi_center = 0.5 * (phi1 + phi2)
            domega = (np.cos(theta1) - np.cos(theta2)) * (phi2 - phi1)
            patch_half_angle = np.sqrt(max(domega, 1.0e-12) / np.pi)
            yield phi_center, theta_center, patch_half_angle, domega


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


def _angular_separation(theta1: float, phi1: float, theta2: float, phi2: float) -> float:
    cos_alpha = (
        np.cos(theta1) * np.cos(theta2)
        + np.sin(theta1) * np.sin(theta2) * np.cos(phi1 - phi2)
    )
    return float(np.arccos(np.clip(cos_alpha, -1.0, 1.0)))

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


def _extract_pair_flux(grid: np.ndarray, times_s: np.ndarray, frequencies_hz: np.ndarray) -> np.ndarray:
    if grid.ndim != 2:
        return np.asarray(grid, dtype=float)
    if frequencies_hz.shape == times_s.shape:
        unique_freqs, inverse = np.unique(frequencies_hz, return_inverse=True)
        if unique_freqs.shape[0] == grid.shape[0]:
            return grid[inverse, np.arange(times_s.shape[0])]
    if grid.shape == (frequencies_hz.shape[0], times_s.shape[0]):
        return np.diag(grid)
    raise ValueError("Cannot extract pairwise flux from the provided grid.")


def _jet_magnetar_active(jet: JetProfile, theta_center: float) -> bool:
    if jet.kind in ("tophat", "gaussian", "powerlaw", "steppowerlaw"):
        return theta_center <= getattr(jet, "theta_c", getattr(jet, "theta_j", jet.theta_max))
    if jet.kind == "twocomponent":
        return theta_center <= jet.theta_n
    return True


def _project_medium_to_kernel(medium: Medium, *, phi_center: float, theta_center: float) -> dict[str, float]:
    if medium.kind in ("ism", "wind"):
        return medium.to_kernel_params()

    radius = np.logspace(9.0, 19.0, 256)
    density = np.asarray(medium.density(phi_center, theta_center, radius), dtype=float)
    if density.shape != radius.shape:
        raise ValueError("Custom medium callable must return a density array with the same shape as the radius grid.")
    if not np.all(np.isfinite(density)) or np.any(density <= 0.0):
        raise ValueError("Custom medium callable must return positive finite densities.")

    candidates = [
        _fit_ism_medium(radius, density),
        _fit_wind_medium(radius, density),
        _fit_jump_medium(radius, density),
    ]
    errors = [_medium_fit_error(radius, density, candidate) for candidate in candidates]
    return candidates[int(np.argmin(errors))]


def _fit_ism_medium(radius: np.ndarray, density: np.ndarray) -> dict[str, float]:
    return {
        "d_ne": float(np.exp(np.mean(np.log(density)))),
        "a_star": -1.0,
        "r0": 1.0e9,
        "r_tr": 1.0e18,
        "f_jump": 1.0,
        "f_wide": 0.1,
    }


def _fit_wind_medium(radius: np.ndarray, density: np.ndarray) -> dict[str, float]:
    logr = np.log(radius)
    logd = np.log(density)
    slopes = np.gradient(logd, logr)
    wind_mask = np.abs(slopes + 2.0) < 0.35
    if np.any(wind_mask):
        a_star = float(np.median(density[wind_mask] * radius[wind_mask] ** 2 / 3.0e35))
    else:
        a_star = 0.0

    tail_slope = float(np.median(slopes[-16:]))
    floor = 0.0 if abs(tail_slope + 2.0) < 0.35 else float(np.min(density))
    r0 = 0.0
    if a_star > 0.0:
        wind = a_star * 3.0e35 / radius**2
        head_slope = float(np.median(slopes[:16]))
        if abs(head_slope) < 0.35 and density[0] < 0.95 * wind[0]:
            plateau_mask = np.abs(slopes) < 0.35
            if np.any(plateau_mask):
                n0 = float(np.median(density[plateau_mask]))
                r0 = float(np.sqrt(a_star * 3.0e35 / n0))
    return {
        "d_ne": floor,
        "a_star": a_star,
        "r0": r0,
        "r_tr": 1.0e18,
        "f_jump": 1.0,
        "f_wide": 0.1,
    }


def _fit_jump_medium(radius: np.ndarray, density: np.ndarray) -> dict[str, float]:
    floor = float(np.min(density))
    peak_idx = int(np.argmax(density))
    peak = float(density[peak_idx])
    logr = np.log10(radius)
    excess = np.maximum(density - floor, 0.0)
    if peak <= floor or np.sum(excess) <= 0.0:
        return _fit_ism_medium(radius, density)
    center = logr[peak_idx]
    variance = np.sum(excess * (logr - center) ** 2) / np.sum(excess)
    width = float(np.sqrt(max(variance, 0.03**2)))
    return {
        "d_ne": floor,
        "a_star": -1.0,
        "r0": 1.0e9,
        "r_tr": float(radius[peak_idx]),
        "f_jump": float(max(peak / floor, 1.0)),
        "f_wide": width,
    }


def _medium_fit_error(radius: np.ndarray, density: np.ndarray, params: dict[str, float]) -> float:
    model = _evaluate_kernel_density(radius, params)
    return float(np.mean((np.log10(model) - np.log10(density)) ** 2))


def _evaluate_kernel_density(radius: np.ndarray, params: dict[str, float]) -> np.ndarray:
    radius = np.asarray(radius, dtype=float)
    if params["a_star"] > 0.0:
        wind = params["a_star"] * 3.0e35 / radius**2
        density = np.where(wind <= params["d_ne"] / 4.0, params["d_ne"], wind)
    else:
        density = params["d_ne"] * (
            1.0
            + (params["f_jump"] - 1.0)
            * np.exp(-(np.log10(radius) - np.log10(params["r_tr"])) ** 2 / (2.0 * params["f_wide"] ** 2))
        )
    if np.any(radius < params["r0"]) and params["a_star"] > 0.0:
        density = density.copy()
        density[radius < params["r0"]] = params["a_star"] * 3.0e35 / params["r0"] ** 2
    return density


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
