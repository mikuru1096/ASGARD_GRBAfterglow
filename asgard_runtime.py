from __future__ import annotations

from dataclasses import dataclass
import numpy as np

import src.Electron.FS_electron_weno5 as electron_weno5_module
try:
    import src.Electron.FS_electron_t2g1_v3 as electron_t2g1_module
except ImportError:
    try:
        import src.Electron.FS_electron_t2g1_v2 as electron_t2g1_module
    except ImportError:
        import src.Electron.FS_electron_t2g1 as electron_t2g1_module

try:
    import src.Electron.FS_electron_slc1 as electron_slc1_module
except ImportError:
    try:
        import src.Electron.FS_electron_slc1_v3 as electron_slc1_module
    except ImportError:
        try:
            import src.Electron.FS_electron_slc1_v2 as electron_slc1_module
        except ImportError:
            import src.Electron.FS_electron_slc1 as electron_slc1_module
try:
    import src.Electron.FS_electron_charint as electron_charint_module
except ImportError:
    import src.Electron.FS_electron_charint as electron_charint_module
from asgard_models import FitConfig
from src import Dynamics, Electron, Radiation, constants


@dataclass
class ReverseShockParameters:
    delta_t_s: float
    epsilon_e: float
    epsilon_b: float
    p: float
    f_e: float


@dataclass
class ReverseShockDynamics:
    t_cross: float
    r_cross: float
    e3_cross: float
    gam20: float
    swept_mass_g: np.ndarray
    gam_e: np.ndarray
    d_n_gam_e: np.ndarray


@dataclass
class DynamicsSolution:
    r_tobs: np.ndarray
    r_gamma: np.ndarray
    radius: np.ndarray
    swept_mass_g: np.ndarray
    reverse_shock: ReverseShockDynamics | None = None


@dataclass
class ElectronSolution:
    gam_e: np.ndarray
    d_n_gam_e: np.ndarray
    l_syn_spec: np.ndarray
    seed_syn: np.ndarray
    nu_m: np.ndarray
    nu_c: np.ndarray
    nu_a: np.ndarray
    work_x_edge_log10: np.ndarray | None = None
    work_d_n_x: np.ndarray | None = None


@dataclass
class ReverseShockEmission:
    l_syn_spec: np.ndarray
    seed_syn: np.ndarray
    magnetic_field_g: np.ndarray
    nu_m: np.ndarray
    nu_c: np.ndarray
    nu_a: np.ndarray
    nu_M: np.ndarray


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
        gam_e,
        d_n_gam_e,
    )
    return DynamicsSolution(r_tobs, r_gamma, radius, swept_mass_g, reverse_shock=reverse_shock)


def solve_electron(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: FitConfig,
) -> ElectronSolution:
    solver_name = _resolve_electron_solver(config)
    if solver_name == "weno5":
        gam_e, d_n_gam_e, l_syn_spec, seed_syn = Electron.fs_electron_weno5(
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
        return ElectronSolution(gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a)

    if solver_name == "t2g1":
        gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_t2g1_module.fs_electron_t2g1(
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
        return ElectronSolution(gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a)

    if solver_name == "slc1":
        gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_slc1_module.fs_electron_slc1(
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
        return ElectronSolution(gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a)

    if solver_name == "charint":
        gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_charint_module.fs_electron_charint(
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
        return ElectronSolution(gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a)

    gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = Electron.fs_electron_fullhide(
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
    return ElectronSolution(gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a)


def _resolve_electron_solver(config: FitConfig) -> str:
    if config.weno5:
        return "weno5"
    solver_name = config.electron_solver.lower()
    if solver_name not in {"fullhide", "t2g1", "weno5", "slc1", "charint"}:
        raise ValueError(f"Unsupported electron solver: {config.electron_solver}")
    return solver_name


def solve_reverse_shock_emission(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: FitConfig,
) -> ReverseShockEmission | None:
    reverse_params = _resolve_reverse_shock_parameters(config)
    if reverse_params is None or dynamics.reverse_shock is None:
        return None

    l_syn_spec, seed_syn = Radiation.seed_reverse(
        dynamics.reverse_shock.t_cross,
        dynamics.reverse_shock.r_cross,
        dynamics.reverse_shock.e3_cross,
        dynamics.reverse_shock.gam20,
        reverse_params.delta_t_s,
        reverse_params.epsilon_b,
        boundary,
        dynamics.r_tobs,
        dynamics.r_gamma,
        dynamics.radius,
        dynamics.reverse_shock.gam_e,
        dynamics.reverse_shock.d_n_gam_e,
        v_seed,
        config.num_threads,
    )
    nu_m, nu_c, nu_a, magnetic_field_g, nu_M = _compute_reverse_shock_characteristic_frequencies(
        config,
        reverse_params,
        dynamics,
    )
    return ReverseShockEmission(
        l_syn_spec=l_syn_spec,
        seed_syn=seed_syn,
        magnetic_field_g=magnetic_field_g,
        nu_m=nu_m,
        nu_c=nu_c,
        nu_a=nu_a,
        nu_M=nu_M,
    )


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

        d_ne = _reverse_ambient_density(radius_loc, config)
        e2 = 4.0 * gamma2 * gamma2 * d_ne * constants.para_m_p_e
        if radius_loc < dynamics.reverse_shock.r_cross:
            e3 = e2
        else:
            e3 = dynamics.reverse_shock.e3_cross * (dynamics.reverse_shock.r_cross / radius_loc) ** 3
            e3 = e3 * gamma2 / dynamics.reverse_shock.gam20

        db = np.sqrt(8.0 * np.pi * reverse_params.epsilon_b * e3)
        magnetic_field_g[i - 1] = db
        gam_e_max = 3.0 * constants.para_m_energy / np.sqrt(8.0 * db * constants.para_e**3)
        gam_e_m = _minimum_reverse_shock_electron_lorentz_factor(reverse_params, gamma34, gam_e_max)
        gam_e_m = max(gam_e_m, gamma_floor)
        gam_e_c = 7.7e8 * (1.0 + config.z) / gamma2 / db**2 / dynamics.r_tobs[i]
        doppler_den = _doppler_denominator(gamma2, config.z)

        nu_m[i - 1] = _synchrotron_frequency(db, gam_e_m, doppler_den)
        nu_c[i - 1] = _synchrotron_frequency(db, gam_e_c, doppler_den)
        nu_M[i - 1] = _synchrotron_frequency(db, gam_e_max, doppler_den)

        nu_a_comoving = electron_weno5_module.get_y.get_nu_a(
            radius_loc,
            db,
            dynamics.reverse_shock.gam_e,
            dynamics.reverse_shock.d_n_gam_e[:, i - 1],
        )
        nu_a[i - 1] = nu_a_comoving / doppler_den

    return nu_m, nu_c, nu_a, magnetic_field_g, nu_M


def _ambient_density(radius_cm: float, config: FitConfig) -> float:
    if config.a_star > 0.0:
        d_ne_wind = config.a_star * 3.0e35 / radius_cm**2
        if d_ne_wind <= config.d_ne / 4.0:
            d_ne = config.d_ne
        else:
            d_ne = d_ne_wind
    else:
        d_ne = config.d_ne * (
            1.0
            + (config.f_jump - 1.0)
            * np.exp(-(np.log10(radius_cm) - np.log10(config.r_tr)) ** 2 / (2.0 * config.f_wide**2))
        )

    if radius_cm < config.r0:
        d_ne = config.a_star * 3.0e35 / config.r0**2

    return d_ne


def _doppler_denominator(gamma_bulk: float, redshift: float) -> float:
    beta_bulk = np.sqrt(1.0 - gamma_bulk**-2)
    return gamma_bulk * (1.0 - beta_bulk) * (1.0 + redshift)


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
