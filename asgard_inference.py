from __future__ import annotations

import numpy as np

from asgard_component_backend import solve_model_state_from_setup, extract_physical_solution_from_state
from asgard_models import FitConfig
from asgard_postprocess import compute_band_fluxes, compute_light_curve_redchi
from asgard_presets import build_baseline_config
from asgard_setup import build_simulation_setup


def build_inference_config(
    *,
    e_iso: float,
    p: float,
    eta_0: float,
    epsilon_e: float,
    epsilon_b: float,
    f_e: float,
    d_ne: float,
    opening_angle_jet: float,
    **overrides,
) -> FitConfig:
    return build_baseline_config(
        f_e=f_e,
        e_iso=e_iso,
        d_ne=d_ne,
        eta_0=eta_0,
        p=p,
        opening_angle_jet=opening_angle_jet,
        epsilon_e=epsilon_e,
        epsilon_b=epsilon_b,
        z=1.0,
        num_gam_e=101,
        **overrides,
    )


def build_log_inference_config(
    *,
    e_iso_log10: float,
    p: float,
    eta_0_log10: float,
    epsilon_e_log10: float,
    epsilon_b_log10: float,
    f_e_log10: float,
    d_ne_log10: float,
    opening_angle_jet_log10: float,
    **overrides,
) -> FitConfig:
    return build_inference_config(
        f_e=10.0**f_e_log10,
        e_iso=10.0**e_iso_log10,
        d_ne=10.0**d_ne_log10,
        eta_0=10.0**eta_0_log10,
        p=p,
        opening_angle_jet=10.0**opening_angle_jet_log10,
        epsilon_e=10.0**epsilon_e_log10,
        epsilon_b=10.0**epsilon_b_log10,
        **overrides,
    )


def evaluate_fit_loglike(config: FitConfig) -> float:
    setup = build_simulation_setup(config)
    state = solve_model_state_from_setup(config, setup)
    physical = extract_physical_solution_from_state(state)
    bands_flux = compute_band_fluxes(setup, physical, config)
    redchi = compute_light_curve_redchi(bands_flux, setup.observer_time_s, config)
    if np.isnan(redchi):
        redchi = np.inf
    return -0.5 * redchi
