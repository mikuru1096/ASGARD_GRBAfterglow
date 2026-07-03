from __future__ import annotations

from dataclasses import dataclass
from math import pi

import numpy as np

from asgard_core.asgard_runtime import _rs_fast_shock_allowed, _rs_shock_upstream_u, _rs_vegas_downstream_u
from src import constants


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
    beta_shock_lab: float
    crossing_time_lab_s: float
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
    radius_collision = constants.para_c * engine_gap_s * beta_slow * beta_fast / (beta_fast - beta_slow)
    lab_collision_time = radius_collision / (beta_slow * constants.para_c)
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
    t0_axis = (lab_collision_time - radius_collision / constants.para_c) * (1.0 + redshift)

    fs = _build_branch_history(
        "fs",
        slow,
        fs_jump,
        gamma_contact,
        beta_contact,
        radius_collision,
        lab_collision_time,
        t0_axis,
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
        beta_contact,
        radius_collision,
        lab_collision_time,
        t0_axis,
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
    if gamma_rel <= 1.0 + 1.0e-10:
        return False
    return bool(_rs_fast_shock_allowed(gamma_rel, sigma))


def baryonic_mass_g(shell: InternalShockShell) -> float:
    matter_fraction = 1.0 / (1.0 + shell.sigma)
    return shell.total_energy_iso_erg * matter_fraction / (shell.gamma * constants.para_c**2)


def beta_from_gamma(gamma: float) -> float:
    return float(np.sqrt(1.0 - gamma**-2))


def shell_lab_width_cm(shell: InternalShockShell) -> float:
    return constants.para_c * shell.duration_s


def shell_proper_number_density(shell: InternalShockShell, radius_cm: float) -> float:
    return baryonic_mass_g(shell) / (
        4.0 * pi * radius_cm**2 * shell.gamma * shell_lab_width_cm(shell) * constants.para_m_p
    )


def relative_gamma(gamma_a: float, gamma_b: float) -> float:
    return gamma_a * gamma_b * (1.0 - beta_from_gamma(gamma_a) * beta_from_gamma(gamma_b))


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
    lo = gamma_slow * (1.0 + 1.0e-12)
    hi = gamma_fast * (1.0 - 1.0e-12)

    def balance(gamma_contact: float) -> float:
        g_fs = relative_gamma(gamma_contact, gamma_slow)
        g_rs = relative_gamma(gamma_contact, gamma_fast)
        p_fs = _postshock_pressure(g_fs, n_slow, sigma_slow)
        p_rs = _postshock_pressure(g_rs, n_fast, sigma_fast)
        return p_fs - p_rs

    f_lo = balance(lo)
    f_hi = balance(hi)
    if f_lo * f_hi > 0.0:
        scan = np.linspace(lo, hi, 129)
        values = np.array([balance(float(gamma)) for gamma in scan], dtype=float)
        changes = np.flatnonzero(values[1:] * values[:-1] <= 0.0)
        if changes.size == 0:
            raise ValueError("contact pressure balance has no bracket between the 2 shell Lorentz factors.")
        lo = float(scan[changes[0]])
        hi = float(scan[changes[0] + 1])
        f_lo = float(values[changes[0]])
        f_hi = float(values[changes[0] + 1])
    for _ in range(96):
        mid = 0.5 * (lo + hi)
        f_mid = balance(mid)
        if f_mid == 0.0:
            return mid
        if f_lo * f_mid <= 0.0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    return 0.5 * (lo + hi)


def _postshock_pressure(gamma_rel: float, upstream_density_cm3: float, sigma: float) -> float:
    compression = _rs_mag_comp(gamma_rel, sigma)
    specific_internal = _rs_mag_specific_internal(gamma_rel, sigma)
    rest_energy_density = upstream_density_cm3 * constants.para_m_p * constants.para_c**2
    thermal_pressure = compression * specific_internal * rest_energy_density / 3.0
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
    compression = _rs_mag_comp(gamma_rel, upstream_shell.sigma)
    specific_internal = _rs_mag_specific_internal(gamma_rel, upstream_shell.sigma)
    pressure = _postshock_pressure(gamma_rel, upstream_density_cm3, upstream_shell.sigma)
    valid = fast_shock_allowed(gamma_rel, upstream_shell.sigma)
    upstream_beta = beta_from_gamma(upstream_shell.gamma)
    shock_speed_cd = _shock_speed_contact_frame(gamma_rel, upstream_shell.sigma)
    if shock_direction > 0.0:
        beta_shock_lab = (beta_contact + shock_speed_cd) / (1.0 + beta_contact * shock_speed_cd)
        crossing = shell_lab_width_cm(upstream_shell) / (constants.para_c * (beta_shock_lab - upstream_beta))
    else:
        beta_shock_lab = (beta_contact - shock_speed_cd) / (1.0 - beta_contact * shock_speed_cd)
        crossing = shell_lab_width_cm(upstream_shell) / (constants.para_c * (upstream_beta - beta_shock_lab))
    magnetic_fraction = upstream_shell.sigma / (1.0 + upstream_shell.sigma)
    return BranchJump(
        name=name,
        valid_shock=valid and crossing > 0.0 and compression > 0.0 and specific_internal > 0.0,
        gamma_rel=gamma_rel,
        beta_upstream=upstream_beta,
        beta_shock_lab=beta_shock_lab,
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
    beta_contact: float,
    radius_collision_cm: float,
    lab_collision_time_s: float,
    t0_axis_s: float,
    redshift: float,
    epsilon_e: float,
    epsilon_b: float,
    num_steps: int,
) -> BranchHistory:
    xi = np.linspace(0.0, 1.0, num_steps, dtype=float)
    t_lab = xi * jump.crossing_time_lab_s
    radius = radius_collision_cm + jump.beta_shock_lab * constants.para_c * t_lab
    characteristic_time = ((lab_collision_time_s + t_lab) - radius / constants.para_c) * (1.0 + redshift)
    characteristic_time = characteristic_time - t0_axis_s
    gamma = np.full(num_steps, gamma_contact, dtype=float)
    upstream_density = shell_proper_number_density(shell, radius)
    if not jump.valid_shock:
        zeros = np.zeros(num_steps, dtype=float)
        return BranchHistory(
            name=name,
            valid_shock=False,
            characteristic_time_s=characteristic_time,
            gamma=gamma,
            radius_cm=radius,
            shell_density_cm3=upstream_density,
            thermal_energy_density_erg_cm3=zeros,
            thermal_luminosity_comoving_erg_s=zeros,
            electron_luminosity_comoving_erg_s=zeros,
            ordered_b_g=zeros,
            turbulent_b_g=zeros,
            total_b_g=zeros,
            swept_mass_g=zeros,
            internal_energy_erg=zeros,
            comoving_volume_cm3=zeros,
            upstream_gamma=shell.gamma,
            upstream_baryonic_mass_g=baryonic_mass_g(shell),
            upstream_lab_width_cm=shell_lab_width_cm(shell),
            jump=jump,
        )
    rho_up = upstream_density * constants.para_m_p
    thermal_density = upstream_density * jump.compression * jump.specific_internal_energy
    thermal_density = thermal_density * constants.para_m_p * constants.para_c**2
    ordered_b = _upstream_ordered_b(rho_up, shell.sigma) * jump.compression
    turbulent_b = np.sqrt(8.0 * pi * epsilon_b * thermal_density)
    total_b = np.sqrt(ordered_b**2 + turbulent_b**2)
    relative_beta = np.abs(jump.beta_shock_lab - jump.beta_upstream)
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


def _shock_speed_contact_frame(gamma_rel: float, sigma: float) -> float:
    downstream_u = float(_rs_vegas_downstream_u(gamma_rel, sigma))
    return downstream_u / np.sqrt(1.0 + downstream_u**2)


def _upstream_ordered_b(rho_up_g_cm3: np.ndarray, sigma: float) -> np.ndarray:
    return np.sqrt(4.0 * np.pi * sigma * rho_up_g_cm3 * constants.para_c**2)


def _rs_mag_comp(gamma_rel: float, sigma: float) -> float:
    gamma_eff = max(1.0, float(gamma_rel))
    if gamma_eff <= 1.0 or sigma <= 1.0e-6:
        return 4.0 * gamma_eff
    return _rs_shock_upstream_u(gamma_eff, sigma) / _rs_vegas_downstream_u(gamma_eff, sigma)


def _rs_mag_specific_internal(gamma_rel: float, sigma: float) -> float:
    gamma_eff = max(1.0, float(gamma_rel))
    if gamma_eff <= 1.0:
        return 0.0
    if sigma <= 0.0:
        return gamma_eff - 1.0
    adiabatic_index = 4.0 / 3.0 + 1.0 / (3.0 * gamma_eff)
    u_down = _rs_vegas_downstream_u(gamma_eff, sigma)
    u_up = _rs_shock_upstream_u(gamma_eff, sigma)
    gamma_down = np.sqrt(1.0 + u_down * u_down)
    gamma_up = np.sqrt(1.0 + u_up * u_up)
    compression = u_up / u_down
    h_down = (1.0 + sigma) * gamma_up / gamma_down - compression * sigma
    return float((h_down - 1.0) / adiabatic_index)
