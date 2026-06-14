"""
Shared physics utility functions for ASGARD.

This module consolidates physics calculations that were previously
duplicated across multiple modules (asgard_state, asgard_runtime,
asgard_ssc, asgard_coupling).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from src import constants

if TYPE_CHECKING:
    from asgard_core.asgard_config import FitConfig


def density_jump_arrays(config: FitConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    jump_r = np.asarray(config.jump_r_cm, dtype=float)
    jump_factor = np.asarray(config.jump_factor, dtype=float)
    jump_width = np.asarray(config.jump_width_log10, dtype=float)
    if jump_r.size == 0 and jump_factor.size == 0 and jump_width.size == 0:
        if float(config.f_jump) == 1.0:
            return np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float)
        return (
            np.array([float(config.r_tr)], dtype=float),
            np.array([float(config.f_jump)], dtype=float),
            np.array([float(config.f_wide)], dtype=float),
        )
    if jump_r.shape != jump_factor.shape or jump_r.shape != jump_width.shape:
        raise ValueError("jump_r_cm, jump_factor, and jump_width_log10 must have the same length.")
    if jump_r.size == 0:
        return jump_r, jump_factor, jump_width
    if not np.all(np.isfinite(jump_r)) or not np.all(np.isfinite(jump_factor)) or not np.all(np.isfinite(jump_width)):
        raise ValueError("density jump arrays must contain finite values.")
    if np.any(jump_r <= 0.0) or np.any(jump_factor <= 0.0) or np.any(jump_width <= 0.0):
        raise ValueError("density jump radii, factors, and widths must be positive.")
    return jump_r, jump_factor, jump_width


def ambient_density(radius_cm: np.ndarray | float, config: FitConfig) -> np.ndarray | float:
    """
    Calculate ambient medium density at given radius.

    Handles both ISM and wind profiles, plus density jumps.
    Works with both scalar and array inputs.
    """
    radius = np.asarray(radius_cm, dtype=float)
    scalar_input = radius.ndim == 0

    if config.a_star > 0.0:
        d_ne_wind = config.a_star * 3.0e35 / radius**2
        density = np.where(d_ne_wind <= config.d_ne / 4.0, config.d_ne, d_ne_wind)
    else:
        density = np.full_like(radius, float(config.d_ne), dtype=float)

    jump_r, jump_factor, jump_width = density_jump_arrays(config)
    if jump_r.size > 0:
        enhancement = np.ones_like(radius, dtype=float)
        for radius_j, factor_j, width_j in zip(jump_r, jump_factor, jump_width):
            width_cm = width_j * radius_j
            enhancement = enhancement + (factor_j - 1.0) * np.exp(
                -((radius - radius_j) ** 2) / (2.0 * width_cm**2)
            )
        density = density * enhancement

    # Apply inner boundary cutoff for wind
    if config.a_star > 0.0 and config.r0 > 0.0:
        if scalar_input:
            if radius < config.r0:
                density = config.a_star * 3.0e35 / config.r0**2
        else:
            mask = radius < config.r0
            if np.any(mask):
                density = density.copy()
                density[mask] = config.a_star * 3.0e35 / config.r0**2

    return float(density) if scalar_input else density


def compute_doppler(gamma: np.ndarray | float, redshift: float, theta_obs: float = 0.0) -> np.ndarray | float:
    """
    Calculate Doppler factor for given bulk Lorentz factor.

    For on-axis viewing (theta_obs=0), simplifies to 1/(gamma*(1-beta)*(1+z)).
    """
    gamma_arr = np.asarray(gamma, dtype=float)
    scalar_input = gamma_arr.ndim == 0

    beta = np.sqrt(1.0 - gamma_arr**-2)

    if theta_obs == 0.0:
        # On-axis simplification
        doppler = 1.0 / (gamma_arr * (1.0 - beta) * (1.0 + redshift))
    else:
        # General case
        cos_theta = np.cos(theta_obs)
        doppler = 1.0 / (gamma_arr * (1.0 - beta * cos_theta) * (1.0 + redshift))

    return float(doppler) if scalar_input else doppler


def doppler_denominator(gamma_bulk: float, redshift: float) -> float:
    """
    Calculate Doppler denominator: gamma*(1-beta)*(1+z).

    This is the inverse of the Doppler factor for on-axis viewing.
    """
    beta_bulk = np.sqrt(1.0 - gamma_bulk**-2)
    return gamma_bulk * (1.0 - beta_bulk) * (1.0 + redshift)


def compute_magnetic_field(
    gamma: np.ndarray | float,
    radius_cm: np.ndarray | float,
    config: FitConfig,
    swept_mass_g: np.ndarray | float | None = None,
) -> np.ndarray | float:
    """
    Calculate magnetic field strength.

    If swept_mass_g is provided, uses it directly.
    Otherwise, computes from ambient density.
    """
    if swept_mass_g is not None:
        # Direct calculation from swept mass
        gamma_arr = np.asarray(gamma, dtype=float)
        radius_arr = np.asarray(radius_cm, dtype=float)
        swept_arr = np.asarray(swept_mass_g, dtype=float)

        density = swept_arr / (4.0 / 3.0 * np.pi * radius_arr**3 * constants.para_m_p)
        b_field = 0.39 * np.sqrt(config.epsilon_b * density * gamma_arr * np.maximum(gamma_arr - 1.0, 0.0))
    else:
        # Calculate from ambient density
        density = ambient_density(radius_cm, config)
        gamma_arr = np.asarray(gamma, dtype=float)
        b_field = 0.39 * np.sqrt(config.epsilon_b * density * gamma_arr * np.maximum(gamma_arr - 1.0, 0.0))

    scalar_input = np.asarray(gamma).ndim == 0
    return float(b_field) if scalar_input else b_field


def compute_maximum_synchrotron_frequency(
    gamma: np.ndarray,
    radius_cm: np.ndarray,
    config: FitConfig,
) -> np.ndarray:
    """
    Calculate maximum synchrotron frequency nu_M.

    This is the frequency corresponding to the maximum electron Lorentz factor.
    """
    magnetic_field = compute_magnetic_field(gamma, radius_cm, config)
    doppler = compute_doppler(gamma, config.z, config.theta_v)
    gam_e_max = 3.0 * constants.para_m_energy / np.sqrt(8.0 * magnetic_field * constants.para_e**3)
    return 4.2e6 * magnetic_field * gam_e_max**2 * doppler
