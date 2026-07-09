from __future__ import annotations

from functools import cache

import numpy as np
from astropy import units
from astropy.cosmology import FlatLambdaCDM

from asgard_core.asgard_config import RuntimeConfig, MAX_DENSITY_JUMPS, MAX_DENSITY_PROFILE_POINTS, SimulationSetup
from asgard_core.asgard_physics_utils import densityjumps, densityprofile
from src import constants


R0_INDEX = 26


@cache
def _lumdist(redshift: float) -> float:
    cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)
    return float(cosmo.luminosity_distance(redshift).to(units.cm).value)


def build_setup(
    config: RuntimeConfig,
    requested_frequencies_hz: np.ndarray | None = None,
) -> SimulationSetup:
    geometry_kernel = str(config.geometry_kernel).lower()
    if geometry_kernel not in {"sed_legacy", "sed_adaptive_theta", "chi_eats_2d"}:
        raise ValueError("geometry_kernel must be 'sed_legacy', 'sed_adaptive_theta', or 'chi_eats_2d'.")
    if geometry_kernel == "chi_eats_2d" and not str(config.electron_solver).lower().endswith("_2d"):
        raise ValueError("geometry_kernel='chi_eats_2d' requires a 2d electron solver.")
    if geometry_kernel == "sed_adaptive_theta":
        if float(config.projection_adaptive_rtol) <= 0.0:
            raise ValueError("projection_adaptive_rtol must be positive for sed_adaptive_theta.")
        if int(config.projection_adaptive_max_depth) < 0:
            raise ValueError("projection_adaptive_max_depth must be non-negative for sed_adaptive_theta.")
    if config.luminosity_distance_cm_override is None:
        luminosity_distance_cm = _lumdist(float(config.z))
    else:
        luminosity_distance_cm = config.luminosity_distance_cm_override
    return SimulationSetup(
        luminosity_distance_cm=luminosity_distance_cm,
        boundary=build_boundary(config, luminosity_distance_cm),
        seed_frequency_hz=seedgrid(config, requested_frequencies_hz),
        observer_time_s=np.logspace(config.t_obs_min_log10, config.t_obs_max_log10, config.num_tobs),
    )


def seedgrid(config: RuntimeConfig, requested_frequencies_hz: np.ndarray | None = None) -> np.ndarray:
    seed_min_hz = 1.0e-8 * constants.para_ev2hz
    seed_max_hz = 1.0e4 * constants.para_tev2hz
    if requested_frequencies_hz is not None:
        requested = np.asarray(requested_frequencies_hz, dtype=float)
        requested = requested[np.isfinite(requested) & (requested > 0.0)]
        if requested.size > 0:
            gamma_max = max(float(config.eta_0), 1.0)
            beta_max = np.sqrt(1.0 - gamma_max**-2) if gamma_max > 1.0 else 0.0
            min_doppler_den = gamma_max * (1.0 - beta_max)
            max_doppler_den = gamma_max * (1.0 + beta_max)
            redshift_factor = 1.0 + float(config.z)
            seed_min_hz = min(
                seed_min_hz,
                0.5 * float(np.min(requested)) * min_doppler_den * redshift_factor,
            )
            seed_max_hz = max(
                seed_max_hz,
                2.0 * float(np.max(requested)) * max_doppler_den * redshift_factor,
            )
    v_seed_min = np.log10(seed_min_hz)
    v_seed_max = np.log10(seed_max_hz)
    return np.logspace(v_seed_min, v_seed_max, config.num_nu)


def build_boundary(config: RuntimeConfig, luminosity_distance_cm: float) -> np.ndarray:
    t_end = config.t_obs_max_log10
    u_0 = 1.0e13
    r_0 = config.initial_radius_cm
    epsilon_b_floor = config.epsilon_b if config.epsilon_b_floor is None else config.epsilon_b_floor
    transport_model = str(config.fullhide2d_transport_model).lower()
    escape_mode = str(config.fullhide2d_escape_mode).lower()
    if transport_model not in {"legacy", "pwn_cr_v1"}:
        raise ValueError("fullhide2d_transport_model must be 'legacy' or 'pwn_cr_v1'.")
    if escape_mode not in {"closed", "free_outer"}:
        raise ValueError("fullhide2d_escape_mode must be 'closed' or 'free_outer'.")
    jump_r, jump_factor, jump_width = densityjumps(config)
    if jump_r.size > MAX_DENSITY_JUMPS:
        raise ValueError(f"At most {MAX_DENSITY_JUMPS} density jumps are supported.")
    profile_r, profile_n = densityprofile(config)
    if jump_r.size > 0 and profile_r.size > 0:
        raise ValueError("density profile and density jumps are mutually exclusive.")
    jump_r_pad = np.zeros(MAX_DENSITY_JUMPS, dtype=float)
    jump_factor_pad, jump_width_pad = np.ones((2, MAX_DENSITY_JUMPS), dtype=float)
    profile_r_pad = np.zeros(MAX_DENSITY_PROFILE_POINTS, dtype=float)
    profile_n_pad = np.ones(MAX_DENSITY_PROFILE_POINTS, dtype=float)
    if jump_r.size > 0:
        jump_r_pad[: jump_r.size] = jump_r
        jump_factor_pad[: jump_factor.size] = jump_factor
        jump_width_pad[: jump_width.size] = jump_width
    if profile_r.size > 0:
        profile_r_pad[: profile_r.size] = profile_r
        profile_n_pad[: profile_n.size] = profile_n

    boundary = [
        config.eta_0,
        0.0,
        u_0,
        r_0,
        config.epsilon_e,
        config.epsilon_b,
        config.p,
        config.z,
        config.opening_angle_jet,
        config.theta_v,
        config.d_ne,
        config.a_star,
        luminosity_distance_cm,
        config.e_iso,
        t_end,
        config.f_e,
        config.e_inj_t1,
        config.e_inj_t2,
        config.l_inj_0,
        config.q_inj,
        config.r_tr,
        config.f_jump,
        config.f_wide,
        epsilon_b_floor,
        config.magnetic_decay_alpha_t,
        config.magnetic_decay_t0_s,
        config.r0,
        float(jump_r.size),
        *jump_r_pad,
        *jump_factor_pad,
        *jump_width_pad,
        float(profile_r.size),
        *profile_r_pad,
        *profile_n_pad,
    ]
    transport_selector = 1.0 if transport_model == "pwn_cr_v1" else 0.0
    escape_selector = 1.0 if escape_mode == "free_outer" else 0.0
    stochastic_accel_norm = float(config.fullhide2d_stochastic_accel_norm)
    boundary.extend([transport_selector, stochastic_accel_norm, escape_selector])
    return np.array(boundary, dtype=float)
