"""Shared physical formulas used by ASGARD runtime glue."""

from __future__ import annotations

from decimal import Decimal, localcontext
from typing import TYPE_CHECKING

import numpy as np

from asgard_core.asgard_config import MAX_DENSITY_PROFILE_POINTS as PROFILE_MAX
from src import constants

if TYPE_CHECKING:
    from asgard_core.asgard_config import RuntimeConfig


def densityjumps(config: RuntimeConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    jump_r = np.asarray(config.jump_r_cm, dtype=float)
    jump_factor = np.asarray(config.jump_factor, dtype=float)
    jump_width = np.asarray(config.jump_width, dtype=float)
    if jump_r.size == 0 and jump_factor.size == 0 and jump_width.size == 0:
        if float(config.f_jump) == 1.0:
            return np.array([], dtype=float), np.array([], dtype=float), np.array([], dtype=float)
        return (
            np.array([float(config.r_tr)], dtype=float),
            np.array([float(config.f_jump)], dtype=float),
            np.array([float(config.f_wide)], dtype=float),
        )
    if jump_r.shape != jump_factor.shape or jump_r.shape != jump_width.shape:
        raise ValueError("jump_r_cm, jump_factor, and jump_width must have the same length.")
    if jump_r.size == 0:
        return jump_r, jump_factor, jump_width
    if not np.all(np.isfinite(jump_r)) or not np.all(np.isfinite(jump_factor)) or not np.all(np.isfinite(jump_width)):
        raise ValueError("density jump arrays must contain finite values.")
    if np.any(jump_r <= 0.0) or np.any(jump_factor <= 0.0) or np.any(jump_width <= 0.0):
        raise ValueError("density jump radii, factors, and widths must be positive.")
    return jump_r, jump_factor, jump_width


def densityprofile(config: RuntimeConfig) -> tuple[np.ndarray, np.ndarray]:
    profile_r = np.asarray(config.density_profile_radius_cm, dtype=float)
    profile_n = np.asarray(config.density_profile_n_cm3, dtype=float)
    if profile_r.size == 0 and profile_n.size == 0:
        return profile_r, profile_n
    if profile_r.shape != profile_n.shape:
        raise ValueError("density_profile_radius_cm and density_profile_n_cm3 must have the same length.")
    if profile_r.size < 2:
        raise ValueError("density_profile requires at least 2 radius-density points.")
    if profile_r.size > PROFILE_MAX:
        raise ValueError(f"At most {PROFILE_MAX} density profile points are supported.")
    if not np.all(np.isfinite(profile_r)) or not np.all(np.isfinite(profile_n)):
        raise ValueError("density profile arrays must contain finite values.")
    if np.any(profile_r <= 0.0) or np.any(profile_n <= 0.0):
        raise ValueError("density profile radii and densities must be positive.")
    if np.any(np.diff(profile_r) <= 0.0):
        raise ValueError("density profile radii must be strictly increasing.")
    return profile_r, profile_n


def loglog_interp(radius_cm, log_r, log_n) -> np.ndarray | float:
    """Interpolate logged profile nodes and extrapolate with their edge slopes."""
    radius = np.asarray(radius_cm, dtype=float)
    logx = np.log(radius)
    index = np.searchsorted(log_r[1:-1], logx, side="right")
    slope = (log_n[index + 1] - log_n[index]) / (log_r[index + 1] - log_r[index])
    density = np.exp(log_n[index] + slope * (logx - log_r[index]))
    return float(density) if radius.ndim == 0 else density


def profile_power(profile_r, profile_n) -> np.ndarray:
    """Return q=d ln(r^3 n)/d ln(r) without cancellation near q=0."""
    values = []
    with localcontext() as ctx:
        ctx.prec = 50
        for rleft, rright, nleft, nright in zip(profile_r[:-1], profile_r[1:], profile_n[:-1], profile_n[1:]):
            rscale = Decimal.from_float(float(rright)) / Decimal.from_float(float(rleft))
            nscale = Decimal.from_float(float(nright)) / Decimal.from_float(float(nleft))
            product = nscale * rscale**3
            values.append(float(product.ln() / rscale.ln()))
    return np.asarray(values, dtype=float)


def profile_moment(radius_cm, r0, profile_r, log_r, log_n, power) -> float:
    """Return the exact radial moment integral_0^R r^2 n(r) dr of a log-log profile."""
    radius = float(radius_cm)
    moment = 0.0
    left = 0.0
    if r0 > 0.0:
        core = min(radius, r0)
        if core > 0.0:
            index = np.searchsorted(profile_r[1:-1], r0, side="right")
            logcore = log_n[index] + (power[index] - 3.0) * (np.log(r0) - log_r[index])
            moment = np.exp(logcore + 3.0 * np.log(core) - np.log(3.0))
        if radius <= r0:
            return moment
        left = r0

    for index, qval in enumerate(power):
        lo = left if index == 0 else max(left, profile_r[index])
        hi = radius if index == power.size - 1 else min(radius, profile_r[index + 1])
        if hi <= lo:
            continue
        logw = log_n[index] + 3.0 * log_r[index]
        if lo == 0.0:
            moment += np.exp(logw + qval * (np.log(hi) - log_r[index]) - np.log(qval))
        elif qval == 0.0:
            moment += np.exp(logw + np.log(np.log(hi / lo)))
        else:
            span = qval * np.log(hi / lo)
            if span > 0.0:
                logedge = logw + qval * (np.log(hi) - log_r[index])
                factor = -np.expm1(-span)
            else:
                logedge = logw + qval * (np.log(lo) - log_r[index])
                factor = -np.expm1(span)
            moment += np.exp(logedge + np.log(factor) - np.log(abs(qval)))
    return float(moment)


def profile_crossing(target, r0, profile_r, profile_n) -> float:
    """Solve R integral_0^R r^2 n(r) dr = target for a tabulated profile."""
    from scipy.optimize import brentq

    log_r = np.log(profile_r)
    log_n = np.log(profile_n)
    power = profile_power(profile_r, profile_n)
    moment = lambda value: profile_moment(value, r0, profile_r, log_r, log_n, power)
    base = max(float(r0), float(profile_r[0]))
    edge = target / moment(base)
    if edge == base:
        return base
    lower, upper = sorted((base, edge))
    logtarget = np.log(target)
    root = brentq(
        lambda value: value + np.log(moment(np.exp(value))) - logtarget,
        np.log(lower),
        np.log(upper),
    )
    return float(np.exp(root))


def reverse_mass(config: RuntimeConfig) -> float:
    return config.e_iso / ((1.0 + config.reverse_shock.sigma) * config.eta_0 * constants.para_c**2)


def ambient_density(radius_cm: np.ndarray | float, config: RuntimeConfig) -> np.ndarray | float:
    """
    Calculate ambient medium density at given radius.

    Handles both ISM and wind profiles, plus density jumps.
    Works with both scalar and array inputs.
    """
    radius = np.asarray(radius_cm, dtype=float)
    scalar_input = radius.ndim == 0

    profile_r, profile_n = densityprofile(config)
    if profile_r.size > 0:
        eval_radius = np.maximum(radius, config.r0) if config.r0 > 0.0 else radius
        return loglog_interp(eval_radius, np.log(profile_r), np.log(profile_n))

    if config.a_star > 0.0:
        winddensity = config.a_star * 3.0e35 / radius**2
        density = np.where(winddensity <= config.d_ne / 4.0, config.d_ne, winddensity)
    else:
        density = np.full_like(radius, float(config.d_ne), dtype=float)

    jump_r, jump_factor, jump_width = densityjumps(config)
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


def doppler_factor(gamma: np.ndarray | float, redshift: float, theta_obs: float = 0.0) -> np.ndarray | float:
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


def magfield(
    gamma: np.ndarray | float,
    radius_cm: np.ndarray | float,
    config: RuntimeConfig,
    swept_mass: np.ndarray | float | None = None,
) -> np.ndarray | float:
    """
    Calculate magnetic field strength.

    If swept_mass is provided, uses it directly.
    Otherwise, computes from ambient density.
    """
    gamma_arr = np.asarray(gamma, dtype=float)
    if swept_mass is None:
        density = ambient_density(radius_cm, config)
    else:
        radius_arr = np.asarray(radius_cm, dtype=float)
        swept_arr = np.asarray(swept_mass, dtype=float)
        density = swept_arr / (4.0 / 3.0 * np.pi * radius_arr**3 * constants.para_m_p)

    b_field = 0.39 * np.sqrt(config.epsilon_b * density * gamma_arr * (gamma_arr - 1.0))
    scalar_input = np.asarray(b_field).ndim == 0
    return float(b_field) if scalar_input else b_field
