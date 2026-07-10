from __future__ import annotations

from dataclasses import dataclass
from math import pi
from typing import cast

import numpy as np
from scipy.optimize import brentq

from src import constants
from src.Dynamics import rs_prompt_jump


@dataclass(frozen=True)
class InternalShockShell:
    gamma: float
    total_energy_iso_erg: float
    duration_s: float
    sigma: float


@dataclass(frozen=True)
class InternalShockNumerics:
    num_branch_steps: int = 64


@dataclass(frozen=True)
class BranchJump:
    name: str
    valid_shock: bool
    gamma_rel: float
    beta_upstream: float
    beta_shock_lab: float | None
    axis_rate: float | None
    crossing_time_lab_s: float | None
    compression: float
    specific_internal_energy: float
    magnetic_fraction: float
    pressure_dyn_cm2: float


@dataclass(frozen=True)
class BranchHistory:
    name: str
    valid_shock: bool
    characteristic_time_s: np.ndarray
    gamma: np.ndarray
    radius_cm: np.ndarray
    shell_density_cm3: np.ndarray
    thermal_energy_density_erg_cm3: np.ndarray
    thermal_luminosity_comoving_erg_s: np.ndarray
    electron_luminosity_comoving_erg_s: np.ndarray
    ordered_b_g: np.ndarray
    turbulent_b_g: np.ndarray
    total_b_g: np.ndarray
    swept_mass_g: np.ndarray
    internal_energy_erg: np.ndarray
    comoving_volume_cm3: np.ndarray
    upstream_gamma: float
    upstream_baryonic_mass_g: float
    upstream_lab_width_cm: float
    jump: BranchJump


@dataclass(frozen=True)
class InternalShockSolution:
    slow_shell: InternalShockShell
    fast_shell: InternalShockShell
    radius_collision_cm: float
    gamma_contact: float
    beta_contact: float
    engine_gap_s: float
    redshift: float
    luminosity_distance_cm: float
    opening_angle_rad: float
    viewing_angle_rad: float
    slow_baryonic_mass_g: float
    fast_baryonic_mass_g: float
    fs: BranchHistory
    rs: BranchHistory


def simulate_internal_shock(
    slow: InternalShockShell,
    fast: InternalShockShell,
    *,
    engine_gap_s: float,
    redshift: float,
    luminosity_distance_cm: float,
    opening_angle_rad: float,
    viewing_angle_rad: float = 0.0,
    epsilon_e: float = 0.1,
    epsilon_b: float = 0.01,
    numerics: InternalShockNumerics | None = None,
) -> InternalShockSolution:
    _validate_shell("slow", slow)
    _validate_shell("fast", fast)
    if fast.gamma <= slow.gamma:
        raise ValueError("fast shell gamma must exceed slow shell gamma.")
    if engine_gap_s <= 0.0 or redshift < 0.0 or luminosity_distance_cm <= 0.0:
        raise ValueError("engine_gap_s, redshift, and luminosity_distance_cm define the observer boundary.")
    if opening_angle_rad <= 0.0 or viewing_angle_rad < 0.0:
        raise ValueError("opening/viewing angles must define a physical EATS geometry.")
    if epsilon_e < 0.0 or epsilon_b < 0.0:
        raise ValueError("epsilon_e and epsilon_b must be non-negative.")

    numerics = InternalShockNumerics() if numerics is None else numerics
    if numerics.num_branch_steps < 4:
        raise ValueError("num_branch_steps must be at least 4 for EATS interpolation.")

    beta_slow = beta_from_gamma(slow.gamma)
    beta_fast = beta_from_gamma(fast.gamma)
    speed_gap = beta_gap(slow.gamma, fast.gamma, beta_slow, beta_fast)
    radius_collision = constants.para_c * engine_gap_s * beta_slow * beta_fast / speed_gap
    slow_mass = baryonic_mass_g(slow)
    fast_mass = baryonic_mass_g(fast)
    slow_density = shell_proper_number_density(slow, radius_collision)
    fast_density = shell_proper_number_density(fast, radius_collision)
    gamma_contact = _solve_contact_gamma(
        slow.gamma,
        fast.gamma,
        slow_density,
        fast_density,
        slow.sigma,
        fast.sigma,
    )
    beta_contact = beta_from_gamma(gamma_contact)

    fs_jump = _build_jump(
        "fs",
        gamma_contact,
        beta_contact,
        upstream_shell=slow,
        upstream_density_cm3=slow_density,
        shock_direction=1.0,
    )
    rs_jump = _build_jump(
        "rs",
        gamma_contact,
        beta_contact,
        upstream_shell=fast,
        upstream_density_cm3=fast_density,
        shock_direction=-1.0,
    )
    fs = _build_branch_history(
        "fs",
        slow,
        fs_jump,
        gamma_contact,
        radius_collision,
        redshift,
        epsilon_e,
        epsilon_b,
        numerics.num_branch_steps,
    )
    rs = _build_branch_history(
        "rs",
        fast,
        rs_jump,
        gamma_contact,
        radius_collision,
        redshift,
        epsilon_e,
        epsilon_b,
        numerics.num_branch_steps,
    )
    return InternalShockSolution(
        slow_shell=slow,
        fast_shell=fast,
        radius_collision_cm=radius_collision,
        gamma_contact=gamma_contact,
        beta_contact=beta_contact,
        engine_gap_s=engine_gap_s,
        redshift=redshift,
        luminosity_distance_cm=luminosity_distance_cm,
        opening_angle_rad=opening_angle_rad,
        viewing_angle_rad=viewing_angle_rad,
        slow_baryonic_mass_g=slow_mass,
        fast_baryonic_mass_g=fast_mass,
        fs=fs,
        rs=rs,
    )


def fast_shock_allowed(gamma_rel: float, sigma: float) -> bool:
    return bool(rs_prompt_jump(gamma_rel, sigma)[3])


def baryonic_mass_g(shell: InternalShockShell) -> float:
    matter_fraction = 1.0 / (1.0 + shell.sigma)
    return shell.total_energy_iso_erg * matter_fraction / (shell.gamma * constants.para_c**2)


def beta_from_gamma(gamma: float) -> float:
    return float(np.sqrt((gamma - 1.0) * (gamma + 1.0)) / gamma)


def beta_gap(gamma_a: float, gamma_b: float, beta_a: float, beta_b: float) -> float:
    numerator = (gamma_b - gamma_a) * (gamma_b + gamma_a)
    return numerator / (gamma_a**2 * gamma_b**2 * (beta_a + beta_b))


def shell_lab_width_cm(shell: InternalShockShell) -> float:
    return constants.para_c * shell.duration_s


def shell_proper_number_density(shell: InternalShockShell, radius_cm: float) -> float:
    return baryonic_mass_g(shell) / (
        4.0 * pi * radius_cm**2 * shell.gamma * shell_lab_width_cm(shell) * constants.para_m_p
    )


def relative_gamma(gamma_a: float, gamma_b: float) -> float:
    beta_a = np.sqrt((gamma_a - 1.0) * (gamma_a + 1.0)) / gamma_a
    beta_b = np.sqrt((gamma_b - 1.0) * (gamma_b + 1.0)) / gamma_b
    u_rel = (gamma_b - gamma_a) * (gamma_b + gamma_a) / (gamma_a * gamma_b * (beta_a + beta_b))
    return float(np.hypot(1.0, u_rel))


def _validate_shell(name: str, shell: InternalShockShell) -> None:
    if shell.gamma <= 1.0:
        raise ValueError(f"{name} shell gamma must exceed 1.")
    if shell.total_energy_iso_erg <= 0.0 or shell.duration_s <= 0.0:
        raise ValueError(f"{name} shell energy and duration must be positive.")
    if shell.sigma < 0.0:
        raise ValueError(f"{name} shell sigma must be non-negative.")


def _solve_contact_gamma(
    gamma_slow: float,
    gamma_fast: float,
    n_slow: float,
    n_fast: float,
    sigma_slow: float,
    sigma_fast: float,
) -> float:
    lo = gamma_slow
    hi = gamma_fast

    def balance(gamma_contact: float) -> float:
        g_fs = relative_gamma(gamma_contact, gamma_slow)
        g_rs = relative_gamma(gamma_contact, gamma_fast)
        p_fs = _postshock_pressure(g_fs, n_slow, sigma_slow)
        p_rs = _postshock_pressure(g_rs, n_fast, sigma_fast)
        return p_fs - p_rs

    return float(brentq(balance, lo, hi))


def _postshock_pressure(gamma_rel: float, upstream_density_cm3: float, sigma: float) -> float:
    _, compression, specific_internal, _ = rs_prompt_jump(gamma_rel, sigma)
    return _jump_pressure(gamma_rel, upstream_density_cm3, sigma, compression, specific_internal)


def _jump_pressure(
    gamma_rel: float,
    upstream_n: float,
    sigma: float,
    compression: float,
    specific_internal: float,
) -> float:
    gammahat = 4.0 / 3.0 + 1.0 / (3.0 * gamma_rel)
    rest_energy_density = upstream_n * constants.para_m_p * constants.para_c**2
    thermal_pressure = (gammahat - 1.0) * compression * specific_internal * rest_energy_density
    ordered_magnetic_pressure = 0.5 * compression * compression * sigma * rest_energy_density
    return thermal_pressure + ordered_magnetic_pressure


def _build_jump(
    name: str,
    gamma_contact: float,
    beta_contact: float,
    *,
    upstream_shell: InternalShockShell,
    upstream_density_cm3: float,
    shock_direction: float,
) -> BranchJump:
    gamma_rel = relative_gamma(gamma_contact, upstream_shell.gamma)
    u_down, compression, specific_internal, shock_allowed = rs_prompt_jump(gamma_rel, upstream_shell.sigma)
    pressure = _jump_pressure(
        gamma_rel, upstream_density_cm3, upstream_shell.sigma, compression, specific_internal
    )
    upstream_beta = beta_from_gamma(upstream_shell.gamma)
    magnetic_fraction = upstream_shell.sigma / (1.0 + upstream_shell.sigma)
    if not shock_allowed:
        return BranchJump(
            name=name,
            valid_shock=False,
            gamma_rel=gamma_rel,
            beta_upstream=upstream_beta,
            beta_shock_lab=None,
            axis_rate=None,
            crossing_time_lab_s=None,
            compression=compression,
            specific_internal_energy=specific_internal,
            magnetic_fraction=magnetic_fraction,
            pressure_dyn_cm2=pressure,
        )
    shock_speed_cd = u_down / np.sqrt(1.0 + u_down**2)
    if shock_direction > 0.0:
        contact_gap = beta_gap(upstream_shell.gamma, gamma_contact, upstream_beta, beta_contact)
    else:
        contact_gap = beta_gap(gamma_contact, upstream_shell.gamma, beta_contact, upstream_beta)
    speed_den = 1.0 + shock_direction * beta_contact * shock_speed_cd
    one_minus = 0.5 * (gamma_contact**-2 + upstream_shell.gamma**-2 + contact_gap**2)
    shock_gap = (contact_gap + shock_speed_cd * one_minus) / speed_den
    contact_axis = gamma_contact**-2 / (1.0 + beta_contact)
    axis_rate = contact_axis * (1.0 - shock_direction * shock_speed_cd) / speed_den
    beta_shock_lab = (beta_contact + shock_direction * shock_speed_cd) / speed_den
    crossing = shell_lab_width_cm(upstream_shell) / (constants.para_c * shock_gap)
    return BranchJump(
        name=name,
        valid_shock=True,
        gamma_rel=gamma_rel,
        beta_upstream=upstream_beta,
        beta_shock_lab=beta_shock_lab,
        axis_rate=axis_rate,
        crossing_time_lab_s=crossing,
        compression=compression,
        specific_internal_energy=specific_internal,
        magnetic_fraction=magnetic_fraction,
        pressure_dyn_cm2=pressure,
    )


def _build_branch_history(
    name: str,
    shell: InternalShockShell,
    jump: BranchJump,
    gamma_contact: float,
    radius_collision_cm: float,
    redshift: float,
    epsilon_e: float,
    epsilon_b: float,
    num_steps: int,
) -> BranchHistory:
    if not jump.valid_shock:
        empty = np.empty(0, dtype=float)
        return BranchHistory(
            name=name,
            valid_shock=False,
            characteristic_time_s=empty,
            gamma=empty,
            radius_cm=empty,
            shell_density_cm3=empty,
            thermal_energy_density_erg_cm3=empty,
            thermal_luminosity_comoving_erg_s=empty,
            electron_luminosity_comoving_erg_s=empty,
            ordered_b_g=empty,
            turbulent_b_g=empty,
            total_b_g=empty,
            swept_mass_g=empty,
            internal_energy_erg=empty,
            comoving_volume_cm3=empty,
            upstream_gamma=shell.gamma,
            upstream_baryonic_mass_g=baryonic_mass_g(shell),
            upstream_lab_width_cm=shell_lab_width_cm(shell),
            jump=jump,
        )
    crossing = cast(float, jump.crossing_time_lab_s)
    beta_shock = cast(float, jump.beta_shock_lab)
    axis_rate = cast(float, jump.axis_rate)
    xi = np.linspace(0.0, 1.0, num_steps, dtype=float)
    t_lab = xi * crossing
    radius = radius_collision_cm + beta_shock * constants.para_c * t_lab
    characteristic_time = (1.0 + redshift) * t_lab * axis_rate
    gamma = np.full(num_steps, gamma_contact, dtype=float)
    upstream_density = shell_proper_number_density(shell, radius)
    rho_up = upstream_density * constants.para_m_p
    thermal_density = upstream_density * jump.compression * jump.specific_internal_energy
    thermal_density = thermal_density * constants.para_m_p * constants.para_c**2
    ordered_b = _upstream_ordered_b(rho_up, shell.sigma) * jump.compression
    turbulent_b = np.sqrt(8.0 * pi * epsilon_b * thermal_density)
    total_b = np.sqrt(ordered_b**2 + turbulent_b**2)
    relative_beta = shell_lab_width_cm(shell) / (constants.para_c * crossing)
    dmdt_lab = 4.0 * pi * radius**2 * shell.gamma * upstream_density * constants.para_m_p
    dmdt_lab = dmdt_lab * constants.para_c * relative_beta
    dmdt_comoving = gamma_contact * dmdt_lab
    thermal_luminosity = jump.specific_internal_energy * constants.para_c**2 * dmdt_comoving
    thermal_luminosity[[0, -1]] = 0.0
    swept_mass = dmdt_lab * t_lab
    comoving_volume = swept_mass / (jump.compression * rho_up)
    internal_energy = thermal_density * comoving_volume
    return BranchHistory(
        name=name,
        valid_shock=jump.valid_shock,
        characteristic_time_s=characteristic_time,
        gamma=gamma,
        radius_cm=radius,
        shell_density_cm3=upstream_density,
        thermal_energy_density_erg_cm3=thermal_density,
        thermal_luminosity_comoving_erg_s=thermal_luminosity,
        electron_luminosity_comoving_erg_s=epsilon_e * thermal_luminosity,
        ordered_b_g=ordered_b,
        turbulent_b_g=turbulent_b,
        total_b_g=total_b,
        swept_mass_g=swept_mass,
        internal_energy_erg=internal_energy,
        comoving_volume_cm3=comoving_volume,
        upstream_gamma=shell.gamma,
        upstream_baryonic_mass_g=baryonic_mass_g(shell),
        upstream_lab_width_cm=shell_lab_width_cm(shell),
        jump=jump,
    )


def _upstream_ordered_b(rho_up_g_cm3: np.ndarray, sigma: float) -> np.ndarray:
    return np.sqrt(4.0 * np.pi * sigma * rho_up_g_cm3 * constants.para_c**2)
