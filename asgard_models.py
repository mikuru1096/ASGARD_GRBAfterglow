from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from src import constants


def default_num_threads() -> int:
    env_value = os.environ.get("ASGARD_NUM_THREADS")
    if env_value is not None:
        return max(1, int(env_value))
    cpu_count = os.cpu_count()
    if cpu_count is None:
        return 8
    return max(1, cpu_count)


@dataclass
class SpectrumOutputConfig:
    enabled: bool = False
    num_nu_obs: int = 180
    nu_min_hz: float = 1.0e-6 * constants.para_ev2hz
    nu_max_hz: float = 1.0e-3 * constants.para_tev2hz


@dataclass
class ReverseShockConfig:
    enabled: bool = False
    delta_t_s: Optional[float] = None
    epsilon_e: Optional[float] = None
    epsilon_b: Optional[float] = None
    p: Optional[float] = None
    f_e: Optional[float] = None
    include_ssc: bool = False
    include_cross_zone_ic: bool = False


@dataclass(frozen=True)
class ExecutionPolicy:
    num_threads: int = field(default_factory=default_num_threads)
    serial_thresholds: dict[str, int] = field(
        default_factory=lambda: {
            "ssa_work_items": 16384,
            "y_nakar_work_items": 1024,
        }
    )
    patch_batch_size: int = 1
    omp_nested: bool = False


@dataclass
class FitConfig:
    num_threads: int = field(default_factory=default_num_threads)
    index_dyn: int = 3
    index_y: int = 2
    index_syn_integr: int = 2
    electron_solver: str = "fullhide"
    electron_adaptive_substeps: bool = False
    electron_substep_rtol: float = 2.0e-2
    electron_substep_min: int = 25
    electron_substep_max: int = 150
    include_forward_ssc: bool = True
    num_gam_e: int = 201
    num_nu: int = 201
    num_r: int = 300
    num_theta: int = 300
    num_phi: int = 1

    z: float = 0.4
    eta_0: float = 1.0e2
    epsilon_e: float = 1.0e-1
    epsilon_b: float = 1.0e-3
    p: float = 2.5
    opening_angle_jet: float = 1.0e-1
    theta_v: float = 0.0
    f_e: float = 1.0e-1
    e_iso: float = 1.0e53
    d_ne: float = 1.0e-1
    a_star: float = -1.0
    r0: float = 1.0e9
    initial_radius_cm: float = 1.0e14

    ebv: float = 0.0
    rv: float = 2.93
    lyman_ar: float = 0.0
    f_sys: float = -1.0

    e_inj_t1: float = 1.0
    e_inj_t2: float = 100.0
    l_inj_0: float = 0.0
    q_inj: float = -1.0

    r_tr: float = 1.0e18
    f_jump: float = 1.0
    f_wide: float = 0.1

    num_tobs: int = 200
    t_obs_min_log10: float = 2.0
    t_obs_max_log10: float = 8.0
    luminosity_distance_cm_override: Optional[float] = None

    weno5: bool = False
    reverse: bool = False
    plot_lc: bool = False
    show_plots: bool = False

    spectrum_output: SpectrumOutputConfig = field(default_factory=SpectrumOutputConfig)
    reverse_shock: ReverseShockConfig = field(default_factory=ReverseShockConfig)


@dataclass
class PhysicalSolution:
    characteristic_time_s: np.ndarray
    gamma: np.ndarray
    radius_cm: np.ndarray
    absorbed_spectral_flux: np.ndarray
    nu_m: np.ndarray
    nu_c: np.ndarray
    nu_a: np.ndarray
    rs_nu_m: Optional[np.ndarray] = None
    rs_nu_c: Optional[np.ndarray] = None
    rs_nu_a: Optional[np.ndarray] = None


@dataclass
class SimulationSetup:
    luminosity_distance_cm: float
    boundary: np.ndarray
    seed_frequency_hz: np.ndarray
    observer_time_s: np.ndarray


@dataclass
class FitResult:
    t_obs_s: np.ndarray
    characteristic_time_s: np.ndarray
    bands: tuple[str, ...]
    bands_flux: np.ndarray
    redchi: float
    nu_m: np.ndarray
    nu_c: np.ndarray
    nu_a: np.ndarray
    rs_nu_m: Optional[np.ndarray] = None
    rs_nu_c: Optional[np.ndarray] = None
    rs_nu_a: Optional[np.ndarray] = None
    spectrum_freq_hz: Optional[np.ndarray] = None
    spectrum_fnu: Optional[np.ndarray] = None
