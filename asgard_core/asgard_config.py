"""Internal runtime configuration for ASGARD simulations."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from src import constants

MAX_DENSITY_JUMPS = 8
MAX_DENSITY_PROFILE_POINTS = 96


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
    time_s: float | None = None
    dataset_names: tuple[str, ...] = field(default_factory=tuple)


@dataclass
class ReverseShockConfig:
    """Configuration for reverse shock physics."""
    enabled: bool = False
    delta_t_s: float | None = None
    sigma: float = 0.0
    epsilon_e: float | None = None
    epsilon_b: float | None = None
    p: float | None = None
    f_e: float | None = None
    include_ssc: bool = False
    include_cross_zone_ic: bool = False


@dataclass
class HadronicConfig:
    """Configuration for forward-shock hadronic emission."""
    enabled: bool = False
    solver: str = "legacy_1d"
    epsilon_p: float = 0.0
    p_p: float = 2.2
    eta_acc: float = 1.0
    num_gam_p: int = 161
    include_proton_synch: bool = True
    include_pg: bool = False
    include_neutrino: bool = False
    include_bethe_heitler: bool = False
    include_hadronic_inverse_compton: bool = False
    include_pair_production: bool = False
    include_pp: bool = False
    pgamma_scheme: str = "disabled"
    num_nu_nu: int = 121
    reverse_enabled: bool = False
    reverse_epsilon_p: float = 0.0
    pair_cascade_iterations: int = 1
    quantum_syn: bool = False
    pp_model: int = -1


_HADRONIC_DELEGATES = (
    ("hadronic_enabled", "enabled"),
    ("hadronic_solver", "solver"),
    ("epsilon_p", "epsilon_p"),
    ("p_p", "p_p"),
    ("eta_acc", "eta_acc"),
    ("num_gam_p", "num_gam_p"),
    ("include_proton_synch", "include_proton_synch"),
    ("include_pg", "include_pg"),
    ("include_neutrino", "include_neutrino"),
    ("pgamma_scheme", "pgamma_scheme"),
    ("num_nu_nu", "num_nu_nu"),
    ("include_bethe_heitler", "include_bethe_heitler"),
    ("include_hadronic_inverse_compton", "include_hadronic_inverse_compton"),
    ("include_pair_production", "include_pair_production"),
    ("include_pp", "include_pp"),
    ("pair_cascade_iterations", "pair_cascade_iterations"),
    ("reverse_enabled", "reverse_enabled"),
    ("reverse_epsilon_p", "reverse_epsilon_p"),
)

_REVERSE_SHOCK_DELEGATES = (
    ("rvs_shock", "enabled"),
    ("rvs_ssc", "include_ssc"),
    ("reverse_delta_t_s", "delta_t_s"),
    ("reverse_sigma", "sigma"),
    ("include_cross_zone_ic", "include_cross_zone_ic"),
)


@dataclass(frozen=True)
class ExecutionPolicy:
    """Execution policy for parallel computation."""
    num_threads: int = field(default_factory=default_num_threads)


@dataclass
class _RuntimeConfig:
    """
    Internal normalized runtime configuration consumed by the Fortran/kernel chain.
    """
    num_threads: int = field(default_factory=default_num_threads)
    index_dyn: int = 3
    index_y: int = 2
    index_syn_integr: int = 2
    electron_solver: str = "fullhide_1d"
    cooling_kernel: str = "legacy"
    radiation_kernel: str = "legacy"
    dynamics_kernel: str = "forward_legacy"
    geometry_kernel: str = "sed_legacy"
    electron_photon_coupling: str = "separated"
    structured_backend: str = "fortran_1d"
    patch_sampling: str = "uniform"
    patch_projection: str = "auto"
    patch_sampling_pilot_theta: int = 0
    patch_sampling_num_times: int = 12
    patch_sampling_beaming_factor: float = 3.0
    patch_sampling_beaming_resolution: float = 8.0
    structured_parallel_mode: str = "outer"
    structured_outer_threads: int | None = None
    structured_inner_threads: int | None = None
    projection_adaptive_rtol: float = 2.0e-2
    projection_adaptive_max_depth: int = 4
    structured_adaptive_rtol: float = 0.0  # >0 enables adaptive theta grid for chi_2d
    structured_adaptive_max_depth: int = 4
    fullhide2d_transport_model: str = "legacy"
    fullhide2d_stochastic_accel_norm: float = 0.0
    fullhide2d_escape_mode: str = "closed"
    electron_adaptive_substeps: bool = False
    electron_substep_rtol: float = 2.0e-2
    electron_substep_min: int = 100
    electron_substep_max: int = 1000
    include_forward_ssc: bool = True
    thermal_electrons: bool = False
    num_gam_e: int = 201
    num_nu: int = 201
    num_r: int = 300
    eats_num_theta: int = 300
    eats_num_phi: int = 1
    structured_num_theta: int = 12
    structured_num_phi: int = 24
    downstream_num_chi: int | None = None

    z: float = 0.4
    eta_0: float = 1.0e2
    epsilon_e: float = 1.0e-1
    epsilon_b: float = 1.0e-3
    epsilon_b_floor: float | None = None
    magnetic_decay_alpha_t: float = 0.0
    magnetic_decay_t0_s: float = 1.0
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
    jump_r_cm: tuple[float, ...] = field(default_factory=tuple)
    jump_factor: tuple[float, ...] = field(default_factory=tuple)
    jump_width_log10: tuple[float, ...] = field(default_factory=tuple)
    density_profile_radius_cm: tuple[float, ...] = field(default_factory=tuple)
    density_profile_n_cm3: tuple[float, ...] = field(default_factory=tuple)

    num_tobs: int = 200
    t_obs_min_log10: float = 2.0
    t_obs_max_log10: float = 8.0
    luminosity_distance_cm_override: float | None = None

    reverse: bool = False

    spectrum_output: SpectrumOutputConfig = field(default_factory=SpectrumOutputConfig)
    reverse_shock: ReverseShockConfig = field(default_factory=ReverseShockConfig)
    hadronic: HadronicConfig = field(default_factory=HadronicConfig)
    nu_callback: Callable[[str, np.ndarray, np.ndarray, np.ndarray], None] | None = None


def _make_hadronic_delegate(field_name: str) -> property:
    def getter(self):
        return getattr(self.hadronic, field_name)

    def setter(self, value) -> None:
        current = getattr(self.hadronic, field_name)
        setattr(self.hadronic, field_name, type(current)(value))

    return property(getter, setter)


for _public_name, _field_name in _HADRONIC_DELEGATES:
    setattr(_RuntimeConfig, _public_name, _make_hadronic_delegate(_field_name))


def _make_reverse_shock_delegate(field_name: str) -> property:
    def getter(self):
        return getattr(self.reverse_shock, field_name)

    def setter(self, value) -> None:
        current = getattr(self.reverse_shock, field_name)
        setattr(self.reverse_shock, field_name, type(current)(value))

    return property(getter, setter)


for _public_name, _field_name in _REVERSE_SHOCK_DELEGATES:
    setattr(_RuntimeConfig, _public_name, _make_reverse_shock_delegate(_field_name))

del _public_name, _field_name


RuntimeConfig = _RuntimeConfig


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
    spectrum_time_s: float | None = None
    spectrum_freq_hz: np.ndarray | None = None
    spectrum_fnu: np.ndarray | None = None
    spectrum_redchi: float | None = None
