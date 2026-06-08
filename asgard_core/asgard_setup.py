from __future__ import annotations

from functools import lru_cache

import numpy as np
from astropy import units
from astropy.cosmology import FlatLambdaCDM

from asgard_core.asgard_config import FitConfig, SimulationSetup
from src import constants


@lru_cache(maxsize=None)
def _luminosity_distance_cm(redshift: float) -> float:
    cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)
    return float(cosmo.luminosity_distance(redshift).to(units.cm).value)


def build_simulation_setup(
    config: FitConfig,
    requested_frequencies_hz: np.ndarray | None = None,
) -> SimulationSetup:
    if config.luminosity_distance_cm_override is None:
        luminosity_distance_cm = _luminosity_distance_cm(float(config.z))
    else:
        luminosity_distance_cm = config.luminosity_distance_cm_override
    return SimulationSetup(
        luminosity_distance_cm=luminosity_distance_cm,
        boundary=build_boundary(config, luminosity_distance_cm),
        seed_frequency_hz=build_seed_frequency_grid(config, requested_frequencies_hz),
        observer_time_s=build_observer_time_grid(config),
    )


def build_observer_time_grid(config: FitConfig) -> np.ndarray:
    return np.logspace(config.t_obs_min_log10, config.t_obs_max_log10, config.num_tobs)


def build_seed_frequency_grid(config: FitConfig, requested_frequencies_hz: np.ndarray | None = None) -> np.ndarray:
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


def build_boundary(config: FitConfig, luminosity_distance_cm: float) -> np.ndarray:
    t_end = config.t_obs_max_log10
    m_0 = 1.0e12
    u_0 = 1.0e13
    r_0 = config.initial_radius_cm
    epsilon_b_floor = config.epsilon_b if config.epsilon_b_floor is None else config.epsilon_b_floor
    transport_model = str(config.fullhide2d_transport_model).lower()
    escape_mode = str(config.fullhide2d_escape_mode).lower()
    if transport_model not in {"legacy", "pwn_cr_v1"}:
        raise ValueError("fullhide2d_transport_model must be 'legacy' or 'pwn_cr_v1'.")
    if escape_mode not in {"closed", "free_outer"}:
        raise ValueError("fullhide2d_escape_mode must be 'closed' or 'free_outer'.")
    boundary = [
        config.eta_0,
        m_0,
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
    ]
    transport_selector = 1.0 if transport_model == "pwn_cr_v1" else 0.0
    escape_selector = 1.0 if escape_mode == "free_outer" else 0.0
    stochastic_accel_norm = float(config.fullhide2d_stochastic_accel_norm)
    if transport_selector != 0.0 or stochastic_accel_norm != 0.0 or escape_selector != 0.0:
        boundary.extend([transport_selector, stochastic_accel_norm, escape_selector])
    return np.array(boundary, dtype=float)
