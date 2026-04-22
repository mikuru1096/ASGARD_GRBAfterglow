from __future__ import annotations

from dataclasses import dataclass
import numpy as np

import src.Electron.FS_electron_weno5_1d as electron_weno5_module
import src.Electron.FS_electron_t2g1_1d as electron_t2g1_module
import src.Electron.FS_electron_slc1_1d as electron_slc1_module
import src.Electron.FS_electron_charint_1d as electron_charint_module
import src.Electron.FS_electron_charint_2d as electron_charint_2d_module
import src.Electron.FS_electron_fullhide_1d as electron_fullhide_1d_module
import src.Electron.FS_electron_fullhide_2d as electron_fullhide_2d_module
import src.Electron.electron_get_y as electron_get_y_module
from asgard_core.asgard_config import FitConfig
from asgard_core.asgard_numpy import trapezoid
from asgard_core.asgard_types import (
    ReverseShockParameters,
    ReverseShockDynamics,
    DynamicsSolution,
    ElectronSolution,
    ReverseShockEmission,
)
from asgard_core.asgard_physics_utils import ambient_density, doppler_denominator, compute_magnetic_field
from src import Dynamics, Electron, constants


_ELECTRON_SOLVER_ALIASES = {
    "fullhide": "fullhide_1d",
    "fullhide_1d": "fullhide_1d",
    "fullhide_2d": "fullhide_2d",
    "slc1": "slc1_1d",
    "slc1_1d": "slc1_1d",
    "charint": "charint_1d",
    "charint_1d": "charint_1d",
    "charint_2d": "charint_2d",
    "t2g1": "t2g1_1d",
    "t2g1_1d": "t2g1_1d",
    "weno5": "weno5_1d",
    "weno5_1d": "weno5_1d",
}

_COOLING_KERNEL_ALIASES = {
    "legacy": "legacy",
}


def solve_dynamics(boundary: np.ndarray, config: FitConfig) -> DynamicsSolution:
    reverse_params = _resolve_reverse_shock_parameters(config)
    if reverse_params is None:
        r_tobs, r_gamma, radius, swept_mass_g = Dynamics.dynamics_forward(boundary, config.num_r, config.index_dyn)
        return DynamicsSolution(r_tobs, r_gamma, radius, swept_mass_g)

    (
        t_cross,
        r_cross,
        e3_cross,
        gam20,
        r_tobs,
        r_gamma,
        radius,
        swept_mass_g,
        swept_reverse_mass_g,
        magnetic_field_g,
        gam_e,
        d_n_gam_e,
    ) = Dynamics.dynamics_reverse(
        reverse_params.delta_t_s,
        reverse_params.epsilon_e,
        reverse_params.epsilon_b,
        reverse_params.p,
        reverse_params.f_e,
        boundary,
        config.num_r,
        config.num_gam_e,
    )
    reverse_shock = ReverseShockDynamics(
        t_cross,
        r_cross,
        e3_cross,
        gam20,
        swept_reverse_mass_g,
        magnetic_field_g,
        gam_e,
        _renormalize_reverse_shock_distribution(gam_e, d_n_gam_e, swept_reverse_mass_g, reverse_params.f_e),
    )
    return DynamicsSolution(r_tobs, r_gamma, radius, swept_mass_g, reverse_shock=reverse_shock)


def solve_electron(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: FitConfig,
) -> ElectronSolution:
    solver_name = _resolve_electron_solver(config)
    _resolve_cooling_kernel(config)
    if solver_name == "weno5_1d":
        gam_e, d_n_gam_e, l_syn_spec, seed_syn = electron_weno5_module.fs_electron_weno5_1d(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.num_threads,
        )
        nu_m, nu_c, nu_a = _compute_characteristic_frequencies_weno5(
            config,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            gam_e,
            d_n_gam_e,
        )
        return _build_electron_solution(
            config,
            dynamics,
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu_m,
            nu_c,
            nu_a,
        )

    if solver_name == "t2g1_1d":
        gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_t2g1_module.fs_electron_t2g1_1d(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
        )
        return _build_electron_solution(
            config,
            dynamics,
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu_m,
            nu_c,
            nu_a,
        )

    if solver_name == "slc1_1d":
        gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_slc1_module.fs_electron_slc1_1d(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
        )
        return _build_electron_solution(
            config,
            dynamics,
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu_m,
            nu_c,
            nu_a,
        )

    if solver_name == "charint_1d":
        gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_charint_module.fs_electron_charint_1d(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
            1 if config.electron_adaptive_substeps else 0,
            config.electron_substep_rtol,
            config.electron_substep_min,
            config.electron_substep_max,
        )
        return _build_electron_solution(
            config,
            dynamics,
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu_m,
            nu_c,
            nu_a,
        )

    if solver_name == "charint_2d":
        num_chi = _resolve_num_chi(config, solver_name)
        gam_e, d_n_gam_e_chi, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_charint_2d_module.fs_electron_charint_2d(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            config.num_gam_e,
            num_chi,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
        )
        chi_grid = _build_log_chi_grid(dynamics.r_gamma, num_chi)
        return _build_electron_solution(
            config,
            dynamics,
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu_m,
            nu_c,
            nu_a,
            d_n_gam_e_chi=d_n_gam_e_chi,
            chi_grid=chi_grid,
        )

    if solver_name == "fullhide_2d":
        num_chi = _resolve_num_chi(config, solver_name)
        gam_e, d_n_gam_e_chi, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_fullhide_2d_module.fs_electron_fullhide_2d(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            config.num_gam_e,
            num_chi,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
        )
        chi_grid = _build_log_chi_grid(dynamics.r_gamma, num_chi)
        return _build_electron_solution(
            config,
            dynamics,
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu_m,
            nu_c,
            nu_a,
            d_n_gam_e_chi=d_n_gam_e_chi,
            chi_grid=chi_grid,
        )

    gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_fullhide_1d_module.fs_electron_fullhide_1d(
        boundary,
        dynamics.r_tobs,
        dynamics.r_gamma,
        dynamics.radius,
        v_seed,
        config.num_gam_e,
        config.index_y,
        config.index_syn_integr,
        config.num_threads,
        1 if config.electron_adaptive_substeps else 0,
        config.electron_substep_rtol,
        config.electron_substep_min,
        config.electron_substep_max,
    )
    return _build_electron_solution(
        config,
        dynamics,
        gam_e,
        d_n_gam_e,
        l_syn_spec,
        seed_syn,
        nu_m,
        nu_c,
        nu_a,
    )


def _resolve_electron_solver(config: FitConfig) -> str:
    if config.weno5:
        return "weno5_1d"
    solver_name = _ELECTRON_SOLVER_ALIASES.get(config.electron_solver.lower())
    if solver_name is None:
        raise ValueError(f"Unsupported electron solver: {config.electron_solver}")
    return solver_name


def _resolve_cooling_kernel(config: FitConfig) -> str:
    cooling_kernel = _COOLING_KERNEL_ALIASES.get(config.cooling_kernel.lower())
    if cooling_kernel is None:
        raise ValueError(f"Unsupported cooling kernel: {config.cooling_kernel}")
    return cooling_kernel


def _resolve_num_chi(config: FitConfig, solver_name: str | None = None) -> int:
    resolved_solver = _resolve_electron_solver(config) if solver_name is None else solver_name
    user_value = config.num_chi
    if resolved_solver.endswith("_1d"):
        return 1 if user_value is None else 1
    if user_value is None:
        return 64
    if int(user_value) < 2:
        raise ValueError("num_chi must be >= 2 for 2d electron solvers.")
    return int(user_value)


def _build_log_chi_grid(r_gamma: np.ndarray, num_chi: int) -> np.ndarray:
    gamma_arr = np.asarray(r_gamma, dtype=float)
    chi_max = 1.0 + 8.0 * np.max(gamma_arr * gamma_arr)
    deta = np.log10(chi_max) / float(num_chi)
    eta_grid = (np.arange(num_chi, dtype=float) + 0.5) * deta
    return np.power(10.0, eta_grid)


def _build_electron_solution(
    config: FitConfig,
    dynamics: DynamicsSolution,
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
    l_syn_spec: np.ndarray,
    seed_syn: np.ndarray,
    nu_m: np.ndarray,
    nu_c: np.ndarray,
    nu_a: np.ndarray,
    *,
    d_n_gam_e_chi: np.ndarray | None = None,
    chi_grid: np.ndarray | None = None,
) -> ElectronSolution:
    cooling_timescale_s, dynamical_timescale_s = _compute_forward_timescales(
        dynamics.r_gamma,
        dynamics.radius,
        nu_c,
        config,
    )
    return ElectronSolution(
        gam_e=np.asarray(gam_e, dtype=float),
        d_n_gam_e=np.asarray(d_n_gam_e, dtype=float),
        l_syn_spec=np.asarray(l_syn_spec, dtype=float),
        seed_syn=np.asarray(seed_syn, dtype=float),
        nu_m=np.asarray(nu_m, dtype=float),
        nu_c=np.asarray(nu_c, dtype=float),
        nu_a=np.asarray(nu_a, dtype=float),
        d_n_gam_e_chi=None if d_n_gam_e_chi is None else np.asarray(d_n_gam_e_chi, dtype=float),
        chi_grid=None if chi_grid is None else np.asarray(chi_grid, dtype=float),
        cooling_timescale_s=cooling_timescale_s,
        dynamical_timescale_s=dynamical_timescale_s,
    )


def _compute_forward_timescales(
    r_gamma: np.ndarray,
    radius_cm: np.ndarray,
    nu_c: np.ndarray,
    config: FitConfig,
) -> tuple[np.ndarray, np.ndarray]:
    gamma = np.asarray(r_gamma, dtype=float)
    radius = np.asarray(radius_cm, dtype=float)
    nu_c_arr = np.asarray(nu_c, dtype=float)
    magnetic_field_g = np.asarray(compute_magnetic_field(gamma, radius, config), dtype=float)
    doppler_den = np.asarray(doppler_denominator(gamma, config.z), dtype=float)
    beta = np.zeros_like(gamma)
    valid_gamma = gamma > 1.0
    beta[valid_gamma] = np.sqrt(1.0 - gamma[valid_gamma] ** (-2.0))

    gamma_c = np.zeros_like(nu_c_arr)
    valid = (magnetic_field_g > 0.0) & (doppler_den > 0.0) & (nu_c_arr > 0.0)
    gamma_c[valid] = np.sqrt(nu_c_arr[valid] * doppler_den[valid] / (4.2e6 * magnetic_field_g[valid]))

    cooling_timescale_s = np.zeros_like(nu_c_arr)
    valid_cooling = valid & (gamma > 0.0)
    cooling_timescale_s[valid_cooling] = (
        7.7e8
        * (1.0 + float(config.z))
        / (gamma[valid_cooling] * magnetic_field_g[valid_cooling] ** 2 * gamma_c[valid_cooling])
    )

    dynamical_timescale_s = np.zeros_like(radius)
    valid_dyn = beta > 0.0
    dynamical_timescale_s[valid_dyn] = radius[valid_dyn] / (gamma[valid_dyn] * beta[valid_dyn] * constants.para_c)
    return cooling_timescale_s, dynamical_timescale_s


def solve_reverse_shock_emission(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: FitConfig,
) -> ReverseShockEmission | None:
    reverse_params = _resolve_reverse_shock_parameters(config)
    if reverse_params is None or dynamics.reverse_shock is None:
        return None

    (
        nu_m,
        nu_c,
        nu_a,
        magnetic_field_g,
        nu_M,
    ) = _compute_reverse_shock_characteristic_frequencies(
        config,
        reverse_params,
        dynamics,
    )
    l_syn_spec, seed_syn = _compute_reverse_shock_synchrotron_emission(dynamics, v_seed, config)
    return ReverseShockEmission(
        l_syn_spec=l_syn_spec,
        seed_syn=seed_syn,
        magnetic_field_g=magnetic_field_g,
        nu_m=nu_m,
        nu_c=nu_c,
        nu_a=nu_a,
        nu_M=nu_M,
    )


def _compute_reverse_shock_synchrotron_emission(
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: FitConfig,
) -> tuple[np.ndarray, np.ndarray]:
    if dynamics.reverse_shock is None:
        raise ValueError("Reverse shock dynamics are required to compute reverse emission.")
    num_nu = v_seed.shape[0]
    num_r = dynamics.radius.shape[0]
    l_syn_spec = np.zeros((num_nu, num_r), dtype=float)
    seed_syn = np.zeros((num_nu, num_r), dtype=float)
    magnetic_field_g = np.asarray(dynamics.reverse_shock.magnetic_field_g, dtype=float)

    for i in range(num_r):
        db = float(magnetic_field_g[i])
        radius_loc = float(dynamics.radius[i])
        if not np.isfinite(db) or not np.isfinite(radius_loc) or db <= 0.0 or radius_loc <= 0.0:
            continue
        p_syn_i, seed_syn_i = electron_get_y_module.get_syn_selected(
            config.index_syn_integr,
            radius_loc,
            db,
            config.num_threads,
            dynamics.reverse_shock.gam_e,
            dynamics.reverse_shock.d_n_gam_e[:, i],
            v_seed,
        )
        l_syn_spec[:, i] = np.asarray(p_syn_i, dtype=float)
        seed_syn[:, i] = np.asarray(seed_syn_i, dtype=float)
    return l_syn_spec, seed_syn


def _resolve_reverse_shock_parameters(config: FitConfig) -> ReverseShockParameters | None:
    reverse_enabled = config.reverse or config.reverse_shock.enabled
    if not reverse_enabled:
        return None

    if config.reverse_shock.delta_t_s is None:
        raise ValueError("ReverseShockConfig.delta_t_s must be set when reverse shock is enabled.")

    return ReverseShockParameters(
        delta_t_s=config.reverse_shock.delta_t_s,
        epsilon_e=config.epsilon_e if config.reverse_shock.epsilon_e is None else config.reverse_shock.epsilon_e,
        epsilon_b=config.epsilon_b if config.reverse_shock.epsilon_b is None else config.reverse_shock.epsilon_b,
        p=config.p if config.reverse_shock.p is None else config.reverse_shock.p,
        f_e=config.f_e if config.reverse_shock.f_e is None else config.reverse_shock.f_e,
    )


def _compute_characteristic_frequencies_weno5(
    config: FitConfig,
    r_tobs: np.ndarray,
    r_gamma: np.ndarray,
    radius: np.ndarray,
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    num_r = radius.shape[0]
    nu_m = np.zeros(num_r, dtype=float)
    nu_c = np.zeros(num_r, dtype=float)
    nu_a = np.zeros(num_r, dtype=float)

    for i in range(1, num_r):
        radius_loc = radius[i - 1]
        gamma_loc = 0.5 * (r_gamma[i - 1] + r_gamma[i])
        d_ne = _ambient_density(radius_loc, config)
        db = 0.39 * np.sqrt(config.epsilon_b * d_ne * (gamma_loc * (gamma_loc - 1.0)))
        gam_e_max = 3.0 * constants.para_m_energy / np.sqrt(8.0 * db * constants.para_e**3)
        gam_e_m = _minimum_electron_lorentz_factor(config, gamma_loc, gam_e_max)
        gam_e_c = 7.7e8 * (1.0 + config.z) / gamma_loc / db**2 / r_tobs[i]
        doppler_den = _doppler_denominator(gamma_loc, config.z)

        nu_m[i - 1] = _synchrotron_frequency(db, gam_e_m, doppler_den)
        nu_c[i - 1] = _synchrotron_frequency(db, gam_e_c, doppler_den)

        nu_a_comoving = electron_weno5_module.get_y.get_nu_a(radius_loc, db, gam_e, d_n_gam_e[:, i - 1])
        nu_a[i - 1] = nu_a_comoving / doppler_den

    return nu_m, nu_c, nu_a


def _compute_reverse_shock_characteristic_frequencies(
    config: FitConfig,
    reverse_params: ReverseShockParameters,
    dynamics: DynamicsSolution,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    num_r = dynamics.radius.shape[0]
    nu_m = np.zeros(num_r, dtype=float)
    nu_c = np.zeros(num_r, dtype=float)
    nu_a = np.zeros(num_r, dtype=float)
    magnetic_field_g = np.zeros(num_r, dtype=float)
    nu_M = np.zeros(num_r, dtype=float)

    eta_0 = config.eta_0
    beta4 = np.sqrt(1.0 - eta_0**-2)
    gamma_floor = float(dynamics.reverse_shock.gam_e[0])
    for i in range(1, num_r):
        radius_loc = dynamics.radius[i - 1]
        gamma2 = 0.5 * (dynamics.r_gamma[i - 1] + dynamics.r_gamma[i])
        beta2 = np.sqrt(1.0 - gamma2**-2)
        gamma34 = (1.0 - beta2 * beta4) * eta_0 * gamma2

        db = float(dynamics.reverse_shock.magnetic_field_g[i - 1])
        magnetic_field_g[i - 1] = db
        gam_e_max = 3.0 * constants.para_m_energy / np.sqrt(8.0 * db * constants.para_e**3)
        gam_e_m = _minimum_reverse_shock_electron_lorentz_factor(reverse_params, gamma34, gam_e_max)
        gam_e_m = max(gam_e_m, gamma_floor)
        gam_e_c = 7.7e8 * (1.0 + config.z) / gamma2 / db**2 / dynamics.r_tobs[i]
        doppler_den = _doppler_denominator(gamma2, config.z)

        nu_m[i - 1] = _synchrotron_frequency(db, gam_e_m, doppler_den)
        nu_c[i - 1] = _synchrotron_frequency(db, gam_e_c, doppler_den)
        nu_M[i - 1] = _synchrotron_frequency(db, gam_e_max, doppler_den)

        nu_a_comoving = electron_get_y_module.get_nu_a(
            radius_loc,
            db,
            dynamics.reverse_shock.gam_e,
            dynamics.reverse_shock.d_n_gam_e[:, i - 1],
        )
        nu_a[i - 1] = nu_a_comoving / doppler_den

    return nu_m, nu_c, nu_a, magnetic_field_g, nu_M


def _compute_synchrotron_emission_from_distribution(
    radius_cm: np.ndarray,
    magnetic_field_g: np.ndarray,
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
    v_seed: np.ndarray,
    config: FitConfig,
) -> tuple[np.ndarray, np.ndarray]:
    num_nu = v_seed.shape[0]
    num_r = radius_cm.shape[0]
    l_syn_spec = np.zeros((num_nu, num_r), dtype=float)
    seed_syn = np.zeros((num_nu, num_r), dtype=float)

    for i in range(num_r):
        db = float(magnetic_field_g[i])
        radius_loc = float(radius_cm[i])
        if not np.isfinite(db) or not np.isfinite(radius_loc) or db <= 0.0 or radius_loc <= 0.0:
            continue
        p_syn_i, seed_syn_i = electron_weno5_module.get_y.get_syn_selected(
            config.index_syn_integr,
            radius_loc,
            db,
            config.num_threads,
            gam_e,
            d_n_gam_e[:, i],
            v_seed,
        )
        l_syn_spec[:, i] = np.asarray(p_syn_i, dtype=float)
        seed_syn[:, i] = np.asarray(seed_syn_i, dtype=float)

    return l_syn_spec, seed_syn




def _renormalize_reverse_shock_distribution(
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
    swept_mass_g: np.ndarray,
    f_e: float,
) -> np.ndarray:
    gam = np.asarray(gam_e, dtype=float)
    dist = np.asarray(d_n_gam_e, dtype=float).copy()
    targets = np.asarray(swept_mass_g, dtype=float) / constants.para_m_p * float(f_e)
    for i in range(dist.shape[1]):
        total = float(trapezoid(dist[:, i], gam))
        target = float(targets[i])
        if not np.isfinite(total) or total <= 0.0:
            continue
        if not np.isfinite(target) or target <= 0.0:
            continue
        dist[:, i] *= target / total
    return dist


def _ambient_density(radius_cm: float, config: FitConfig) -> float:
    """DEPRECATED: Use asgard_physics_utils.ambient_density instead."""
    return ambient_density(radius_cm, config)


def _doppler_denominator(gamma_bulk: float, redshift: float) -> float:
    """DEPRECATED: Use asgard_physics_utils.doppler_denominator instead."""
    return doppler_denominator(gamma_bulk, redshift)


def _synchrotron_frequency(magnetic_field_g: float, electron_lorentz_factor: float, doppler_den: float) -> float:
    return 4.2e6 * magnetic_field_g * electron_lorentz_factor * electron_lorentz_factor / doppler_den


def _reverse_ambient_density(radius_cm: float, config: FitConfig) -> float:
    if config.a_star >= 0.0:
        d_ne_wind = config.a_star * 3.0e35 / radius_cm**2
        if d_ne_wind <= config.d_ne / 4.0:
            return config.d_ne
        return d_ne_wind
    return config.d_ne


def _minimum_electron_lorentz_factor(config: FitConfig, gamma_bulk: float, gam_e_max: float) -> float:
    temp_gam = config.epsilon_e / config.f_e * constants.para_m_p_div_m_e * (gamma_bulk - 1.0)

    if config.p > 2.0:
        return (config.p - 2.0) / (config.p - 1.0) * temp_gam + 1.0

    if 1.0 < config.p < 2.0:
        return ((2.0 - config.p) / (config.p - 1.0) * temp_gam * gam_e_max ** (config.p - 2.0)) ** (
            1.0 / (config.p - 1.0)
        ) + 1.0

    if np.isclose(config.p, 2.0):
        gam_e_m = 1.0
        temp = temp_gam / np.log(gam_e_max / gam_e_m)
        while abs(1.0 - gam_e_m / temp) > 1.0e-5:
            temp = temp_gam / np.log(gam_e_max / gam_e_m)
            if gam_e_m > temp:
                gam_e_m = 0.5 * (gam_e_m + temp)
            else:
                gam_e_m = 0.5 * (gam_e_m + gam_e_max)
        return gam_e_m

    raise ValueError(f"Unsupported electron index p={config.p}.")


def _minimum_reverse_shock_electron_lorentz_factor(
    reverse_params: ReverseShockParameters,
    gamma34: float,
    gam_e_max: float,
) -> float:
    temp_gam = reverse_params.epsilon_e / reverse_params.f_e * constants.para_m_p_div_m_e * (gamma34 - 1.0)

    if reverse_params.p > 2.05:
        return (reverse_params.p - 2.0) / (reverse_params.p - 1.0) * temp_gam + 1.0

    if 1.0 < reverse_params.p < 2.0:
        return ((2.0 - reverse_params.p) / (reverse_params.p - 1.0) * temp_gam * gam_e_max ** (reverse_params.p - 2.0)) ** (
            1.0 / (reverse_params.p - 1.0)
        ) + 1.0

    if reverse_params.p >= 2.0:
        return 0.05 / 1.05 * temp_gam + 1.0

    raise ValueError(f"Unsupported reverse-shock electron index p={reverse_params.p}.")
