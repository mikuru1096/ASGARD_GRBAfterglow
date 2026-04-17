"""
Configuration module for ASGARD simulations.

This module defines the configuration hierarchy:
- PhysicsConfig: Physical parameters (epsilon_e, epsilon_b, p, etc.)
- NumericalConfig: Grid resolution and solver settings
- OutputConfig: Output control (spectrum, reverse shock, etc.)
- SimulationConfig: Top-level config composing the above three
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from src import constants


def default_num_threads() -> int:
    """Get default number of threads from environment or CPU count."""
    env_value = os.environ.get("ASGARD_NUM_THREADS")
    if env_value is not None:
        return max(1, int(env_value))
    cpu_count = os.cpu_count()
    if cpu_count is None:
        return 8
    return max(1, cpu_count)


@dataclass
class SpectrumOutputConfig:
    """Configuration for spectrum output."""
    enabled: bool = False
    num_nu_obs: int = 180
    nu_min_hz: float = 1.0e-6 * constants.para_ev2hz
    nu_max_hz: float = 1.0e-3 * constants.para_tev2hz


@dataclass
class ReverseShockConfig:
    """Configuration for reverse shock physics."""
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
    """Execution policy for parallel computation."""
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
class PhysicsConfig:
    """Physical parameters for GRB afterglow simulation."""
    # Microphysics
    epsilon_e: float = 1.0e-1
    epsilon_b: float = 1.0e-3
    p: float = 2.5
    f_e: float = 1.0e-1

    # Jet properties
    eta_0: float = 1.0e2
    e_iso: float = 1.0e53
    opening_angle_jet: float = 1.0e-1
    initial_radius_cm: float = 1.0e14

    # Medium properties
    d_ne: float = 1.0e-1
    a_star: float = -1.0
    r0: float = 1.0e9

    # Observer properties
    z: float = 0.4
    theta_v: float = 0.0
    luminosity_distance_cm_override: Optional[float] = None

    # Extinction
    ebv: float = 0.0
    rv: float = 2.93
    lyman_ar: float = 0.0
    f_sys: float = -1.0

    # Energy injection
    e_inj_t1: float = 1.0
    e_inj_t2: float = 100.0
    l_inj_0: float = 0.0
    q_inj: float = -1.0

    # Density jump
    r_tr: float = 1.0e18
    f_jump: float = 1.0
    f_wide: float = 0.1

    # Reverse shock
    reverse_shock: ReverseShockConfig = field(default_factory=ReverseShockConfig)


@dataclass
class NumericalConfig:
    """Numerical parameters for grid resolution and solvers."""
    # Grid resolution
    num_r: int = 300
    num_theta: int = 300
    num_phi: int = 1
    num_gam_e: int = 201
    num_nu: int = 201
    num_tobs: int = 200

    # Time grid
    t_obs_min_log10: float = 2.0
    t_obs_max_log10: float = 8.0

    # Solver selection
    electron_solver: str = "fullhide"
    index_dyn: int = 3
    index_y: int = 2
    index_syn_integr: int = 2

    # Adaptive stepping
    electron_adaptive_substeps: bool = False
    electron_substep_rtol: float = 2.0e-2
    electron_substep_min: int = 25
    electron_substep_max: int = 150

    # Physics toggles
    include_forward_ssc: bool = True

    # Execution policy
    execution: ExecutionPolicy = field(default_factory=ExecutionPolicy)

    @property
    def num_threads(self) -> int:
        """Convenience property for num_threads."""
        return self.execution.num_threads


@dataclass
class OutputConfig:
    """Configuration for output control."""
    spectrum: SpectrumOutputConfig = field(default_factory=SpectrumOutputConfig)
    plot_lc: bool = False
    show_plots: bool = False


@dataclass
class SimulationConfig:
    """Top-level simulation configuration."""
    physics: PhysicsConfig = field(default_factory=PhysicsConfig)
    numerical: NumericalConfig = field(default_factory=NumericalConfig)
    output: OutputConfig = field(default_factory=OutputConfig)


# Legacy support - FitConfig for backward compatibility during migration
@dataclass
class FitConfig:
    """
    DEPRECATED: Use SimulationConfig instead.

    This class is kept for backward compatibility during migration.
    It will be removed in a future version.
    """
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

    def to_simulation_config(self) -> SimulationConfig:
        """Convert FitConfig to new SimulationConfig format."""
        physics = PhysicsConfig(
            epsilon_e=self.epsilon_e,
            epsilon_b=self.epsilon_b,
            p=self.p,
            f_e=self.f_e,
            eta_0=self.eta_0,
            e_iso=self.e_iso,
            opening_angle_jet=self.opening_angle_jet,
            initial_radius_cm=self.initial_radius_cm,
            d_ne=self.d_ne,
            a_star=self.a_star,
            r0=self.r0,
            z=self.z,
            theta_v=self.theta_v,
            luminosity_distance_cm_override=self.luminosity_distance_cm_override,
            ebv=self.ebv,
            rv=self.rv,
            lyman_ar=self.lyman_ar,
            f_sys=self.f_sys,
            e_inj_t1=self.e_inj_t1,
            e_inj_t2=self.e_inj_t2,
            l_inj_0=self.l_inj_0,
            q_inj=self.q_inj,
            r_tr=self.r_tr,
            f_jump=self.f_jump,
            f_wide=self.f_wide,
            reverse_shock=self.reverse_shock,
        )

        execution = ExecutionPolicy(
            num_threads=self.num_threads,
        )

        numerical = NumericalConfig(
            num_r=self.num_r,
            num_theta=self.num_theta,
            num_phi=self.num_phi,
            num_gam_e=self.num_gam_e,
            num_nu=self.num_nu,
            num_tobs=self.num_tobs,
            t_obs_min_log10=self.t_obs_min_log10,
            t_obs_max_log10=self.t_obs_max_log10,
            electron_solver=self.electron_solver,
            index_dyn=self.index_dyn,
            index_y=self.index_y,
            index_syn_integr=self.index_syn_integr,
            electron_adaptive_substeps=self.electron_adaptive_substeps,
            electron_substep_rtol=self.electron_substep_rtol,
            electron_substep_min=self.electron_substep_min,
            electron_substep_max=self.electron_substep_max,
            include_forward_ssc=self.include_forward_ssc,
            execution=execution,
        )

        output = OutputConfig(
            spectrum=self.spectrum_output,
            plot_lc=self.plot_lc,
            show_plots=self.show_plots,
        )

        return SimulationConfig(
            physics=physics,
            numerical=numerical,
            output=output,
        )


# Result dataclasses
@dataclass
class PhysicalSolution:
    """Physical solution data."""
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
    """Simulation setup data."""
    luminosity_distance_cm: float
    boundary: np.ndarray
    seed_frequency_hz: np.ndarray
    observer_time_s: np.ndarray


@dataclass
class FitResult:
    """Fitting result data."""
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
