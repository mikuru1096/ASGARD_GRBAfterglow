from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_physics_utils import ambient_density
from src import constants


@dataclass
class CoupledShockGeometry:
    proper_time_s: np.ndarray
    fs_width_cm: np.ndarray
    rs_width_cm: np.ndarray
    center_delay_s: np.ndarray


def build_coupled_shock_geometry(dynamics, config: RuntimeConfig) -> CoupledShockGeometry:
    if dynamics.reverse_shock is None:
        raise ValueError("Reverse-shock dynamics are required to build coupled-shock geometry.")
    if config.reverse_shock.delta_t_s is None:
        raise ValueError("ReverseShockConfig.delta_t_s must be set when reverse shock is enabled.")

    radius_cm = dynamics.radius
    gamma = dynamics.r_gamma
    proper_time_s = _integrate_proper_time(radius_cm, gamma)

    fs_width_cm = np.zeros_like(radius_cm)
    rs_width_cm = np.zeros_like(radius_cm)

    eta_0 = config.eta_0
    shell_mass_g = config.e_iso / eta_0 / constants.para_c**2
    delta_0_cm = config.reverse_shock.delta_t_s * constants.para_c

    for i, radius_loc in enumerate(radius_cm):
        gamma_loc = gamma[i]
        n1 = ambient_density(radius_loc, config)
        n2 = 4.0 * gamma_loc * n1
        fs_width_cm[i] = dynamics.swept_mass_g[i] / (4.0 * np.pi * radius_loc**2 * n2 * constants.para_m_p)

        delta_shell_cm = max(delta_0_cm, radius_loc / eta_0**2)
        n4 = shell_mass_g / (4.0 * np.pi * constants.para_m_p * radius_loc**2 * eta_0 * delta_shell_cm)
        u2 = np.sqrt(gamma_loc * gamma_loc - 1.0)
        u4 = np.sqrt(eta_0 * eta_0 - 1.0)
        gamma34 = (gamma_loc * gamma_loc + eta_0 * eta_0 - 1.0) / (eta_0 * gamma_loc + u2 * u4)
        n3 = (4.0 * gamma34 + 3.0) * n4
        rs_width_cm[i] = (
            dynamics.reverse_shock.swept_mass_g[i] / (4.0 * np.pi * radius_loc**2 * n3 * constants.para_m_p)
        )

    center_delay_s = 0.5 * (fs_width_cm + rs_width_cm) / constants.para_c
    return CoupledShockGeometry(
        proper_time_s=proper_time_s,
        fs_width_cm=fs_width_cm,
        rs_width_cm=rs_width_cm,
        center_delay_s=center_delay_s,
    )


def build_cross_zone_seed_fields(
    fs_seed_syn: np.ndarray,
    rs_seed_syn: np.ndarray,
    geometry: CoupledShockGeometry,
    angular_factor: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    tau = geometry.proper_time_s
    tau_ret = tau - geometry.center_delay_s
    seed_fs_to_rs = _retarded_seed_interpolation(fs_seed_syn, tau, tau_ret, angular_factor)
    seed_rs_to_fs = _retarded_seed_interpolation(rs_seed_syn, tau, tau_ret, angular_factor)
    return seed_fs_to_rs, seed_rs_to_fs


def _integrate_proper_time(radius_cm: np.ndarray, gamma: np.ndarray) -> np.ndarray:
    proper_time_s = np.zeros_like(radius_cm)
    for i in range(1, radius_cm.shape[0]):
        gamma_mean = 0.5 * (gamma[i - 1] + gamma[i])
        beta_mean = np.sqrt(1.0 - gamma_mean**-2)
        d_radius = radius_cm[i] - radius_cm[i - 1]
        proper_time_s[i] = proper_time_s[i - 1] + d_radius / (beta_mean * gamma_mean * constants.para_c)
    return proper_time_s


def _retarded_seed_interpolation(
    seed_syn: np.ndarray,
    source_time_s: np.ndarray,
    retarded_time_s: np.ndarray,
    angular_factor: float,
) -> np.ndarray:
    shifted_seed = np.zeros_like(seed_syn)
    for i_nu in range(seed_syn.shape[0]):
        shifted_seed[i_nu] = angular_factor * np.interp(
            retarded_time_s,
            source_time_s,
            seed_syn[i_nu],
            left=0.0,
            right=seed_syn[i_nu, -1],
        )
    return shifted_seed
