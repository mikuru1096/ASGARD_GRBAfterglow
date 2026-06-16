from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, NamedTuple, Optional

import numpy as np

from asgard_core.angular_sampling import angular_separation, build_patch_grid, is_axisymmetric_jet
from asgard_core.asgard_config import (
    RuntimeConfig,
    HadronicConfig,
    MAX_DENSITY_PROFILE_POINTS,
    ReverseShockConfig,
    SpectrumOutputConfig,
    default_num_threads,
)
from asgard_core.asgard_postprocess import build_multiband_observer_frequencies, combine_multiband_flux
from asgard_core.asgard_state import (
    FluxComponents,
    SolveState,
    make_query_cfg,
    make_query_setup,
    make_tgrid,
    solve_state_from_setup,
)
from src import Interpolation, constants

class Scale(str, Enum):
    LINEAR = "linear"
    LOG10 = "log10"
    FIXED = "fixed"


@dataclass
class Medium:
    """Callable density profile. Use UniformMedium, WindMedium, or TabulatedMedium."""
    rho: Callable[[float, float, float], float]
    kind: str = "custom"
    label: str = "custom"
    n_ism: float = 1.0
    A_star: float = 0.0
    n0: float | None = None
    k: float = 2.0
    density_profile_radius_cm: tuple[float, ...] = field(default_factory=tuple)
    density_profile_n_cm3: tuple[float, ...] = field(default_factory=tuple)
    kernel_params: dict | None = None

    def density(self, phi: float, theta: float, radius_cm: float) -> float:
        values = np.asarray(self.rho(phi, theta, radius_cm), dtype=float)
        if values.ndim == 0:
            return float(values)
        return values

    def __call__(self, phi: float, theta: float, radius_cm: float) -> float:
        return self.density(phi, theta, radius_cm)

    def to_kernel_params(self) -> dict[str, float]:
        if self.kernel_params is not None:
            return self.kernel_params
        raise NotImplementedError("User-defined Medium is not supported by the current ASGARD kernel.")

    @property
    def number_density_cm3(self) -> float:
        return self.n_ism

    @number_density_cm3.setter
    def number_density_cm3(self, value: float) -> None:
        self.n_ism = float(value)
        if self.kind == "ism":
            self.kernel_params = {"d_ne": self.n_ism, "a_star": -1.0}

    @property
    def density_floor_cm3(self) -> float:
        return self.n_ism

    @density_floor_cm3.setter
    def density_floor_cm3(self, value: float) -> None:
        self.n_ism = float(value)
        if self.kind == "wind":
            self.kernel_params = _wind_kernel_params(self.A_star, self.n_ism, self.n0, self.k)

    @property
    def density_cap_cm3(self) -> float | None:
        return self.n0

    @density_cap_cm3.setter
    def density_cap_cm3(self, value: float | None) -> None:
        self.n0 = None if value is None else float(value)
        if self.kind == "wind":
            self.kernel_params = _wind_kernel_params(self.A_star, self.n_ism, self.n0, self.k)


def UniformMedium(number_density_cm3: float) -> Medium:
    _n_ism = float(number_density_cm3)
    medium = Medium(
        rho=lambda phi, theta, rc: float(np.full(np.asarray(rc).shape, medium.n_ism)),
        kind="ism",
        label="ism",
        n_ism=_n_ism,
        kernel_params={"d_ne": _n_ism, "a_star": -1.0},
    )
    return medium


def WindMedium(
    a_star: float,
    density_floor_cm3: float,
    density_cap_cm3: float | None,
) -> Medium:
    _A_star = float(a_star)
    _n_ism = float(density_floor_cm3)
    _n0 = None if density_cap_cm3 is None else float(density_cap_cm3)
    medium = Medium(
        rho=lambda phi, theta, rc: _wind_rho(medium, rc),
        kind="wind",
        label="wind",
        n_ism=_n_ism,
        A_star=_A_star,
        n0=_n0,
        k=2.0,
        kernel_params=_wind_kernel_params(_A_star, _n_ism, _n0, 2.0),
    )
    return medium


def _wind_rho(medium: Medium, radius_cm) -> float:
    if medium.k != 2.0:
        raise NotImplementedError("The current ASGARD wind kernel only supports k=2.")
    radius = np.asarray(radius_cm, dtype=float)
    d_ne_wind = medium.A_star * 3.0e35 / radius**2
    values = np.where(d_ne_wind <= medium.n_ism / 4.0, medium.n_ism, d_ne_wind)
    if medium.n0 is not None and np.isfinite(medium.n0) and medium.n0 > 0.0:
        values = np.minimum(values, float(medium.n0))
    if values.ndim == 0:
        return float(values)
    return values


def _wind_kernel_params(A_star: float, n_ism: float, n0: float | None, k: float) -> dict[str, float]:
    if k != 2.0:
        raise NotImplementedError("The current ASGARD wind kernel only supports k=2.")
    r0 = 0.0
    if n0 is not None and np.isfinite(n0) and n0 > 0.0:
        r0 = float(np.sqrt(A_star * 3.0e35 / float(n0)))
    return {"d_ne": n_ism, "a_star": A_star, "r0": r0}


def TabulatedMedium(radius_cm, density_cm3, label: str) -> Medium:
    radius = tuple(float(value) for value in radius_cm)
    density = tuple(float(value) for value in density_cm3)
    _validate_density_profile(radius, density)
    radius_arr = np.asarray(radius, dtype=float)
    density_arr = np.asarray(density, dtype=float)

    def _rho(_phi, _theta, rc):
        values = np.exp(np.interp(np.log(np.asarray(rc, dtype=float)), np.log(radius_arr), np.log(density_arr)))
        return float(values) if values.ndim == 0 else values

    return Medium(
        rho=_rho,
        kind="density_profile",
        label=label,
        n_ism=float(density[-1]),
        density_profile_radius_cm=radius,
        density_profile_n_cm3=density,
        kernel_params={"d_ne": float(density[-1]), "a_star": -1.0},
    )


def _validate_density_profile(radius: tuple[float, ...], density: tuple[float, ...]) -> None:
    radius_arr = np.asarray(radius, dtype=float)
    density_arr = np.asarray(density, dtype=float)
    if radius_arr.shape != density_arr.shape:
        raise ValueError("density profile radius and density arrays must have the same length.")
    if radius_arr.size < 2:
        raise ValueError("density profile requires at least two points.")
    if radius_arr.size > MAX_DENSITY_PROFILE_POINTS:
        raise ValueError(f"At most {MAX_DENSITY_PROFILE_POINTS} density profile points are supported.")
    if not np.all(np.isfinite(radius_arr)) or not np.all(np.isfinite(density_arr)):
        raise ValueError("density profile arrays must contain finite values.")
    if np.any(radius_arr <= 0.0) or np.any(density_arr <= 0.0):
        raise ValueError("density profile radii and densities must be positive.")
    if np.any(np.diff(radius_arr) <= 0.0):
        raise ValueError("density profile radii must be strictly increasing.")


class Magnetar(NamedTuple):
    L0: float
    t0: float
    q: float


@dataclass
class JetProfile:
    """Unified jet profile. Use top_hat_jet, gaussian_jet, or power_law_jet."""
    kind: str
    theta_max: float
    # Tophat / Gaussian / PowerLaw params
    E_iso: float = 0.0
    lf: float = 1.0
    theta_j: float = 0.1
    theta_c: float = 0.1
    k_e: float = 2.0
    k_g: float = 2.0
    # TwoComponent params
    E_iso_n: float = 0.0
    lf_n: float = 1.0
    theta_n: float = 0.0
    E_iso_w: float = 0.0
    lf_w: float = 1.0
    theta_w: float = 0.0
    # StepPowerLaw params
    E_iso_c: float = 0.0
    lf_c: float = 1.0
    # Function-defined jet params
    e_iso_fn: Callable[[float, float], float] | None = None
    gamma0_fn: Callable[[float, float], float] | None = None
    # Common
    spreading: bool = False
    duration: float | None = None
    magnetar: Magnetar | None = None

    def is_active(self, phi: float, theta: float) -> bool:
        return self.energy_iso(phi, theta) > 0.0 and self.gamma0(phi, theta) > 1.0

    @property
    def energy_iso_erg(self) -> float:
        return self.E_iso

    @energy_iso_erg.setter
    def energy_iso_erg(self, value: float) -> None:
        self.E_iso = float(value)

    @property
    def initial_lorentz_factor(self) -> float:
        return self.lf

    @initial_lorentz_factor.setter
    def initial_lorentz_factor(self, value: float) -> None:
        self.lf = float(value)

    @property
    def opening_angle_rad(self) -> float:
        return self.theta_j

    @opening_angle_rad.setter
    def opening_angle_rad(self, value: float) -> None:
        self.theta_j = float(value)
        if self.kind == "tophat":
            self.theta_max = float(value)

    @property
    def core_angle_rad(self) -> float:
        return self.theta_c

    @core_angle_rad.setter
    def core_angle_rad(self, value: float) -> None:
        self.theta_c = float(value)

    @property
    def outer_angle_rad(self) -> float:
        return self.theta_max

    @outer_angle_rad.setter
    def outer_angle_rad(self, value: float) -> None:
        self.theta_max = float(value)

    @property
    def shell_duration_s(self) -> float | None:
        return self.duration

    @shell_duration_s.setter
    def shell_duration_s(self, value: float | None) -> None:
        self.duration = None if value is None else float(value)

    @property
    def energy_index(self) -> float:
        return self.k_e

    @energy_index.setter
    def energy_index(self, value: float) -> None:
        self.k_e = float(value)

    @property
    def lorentz_index(self) -> float:
        return self.k_g

    @lorentz_index.setter
    def lorentz_index(self, value: float) -> None:
        self.k_g = float(value)


def top_hat_jet(
    energy_iso_erg: float,
    initial_lorentz_factor: float,
    opening_angle_rad: float,
    shell_duration_s: float | None,
    magnetar: Optional[Magnetar],
    spreading: bool,
) -> JetProfile:
    _theta_j = float(opening_angle_rad)
    jet = JetProfile(
        kind="tophat",
        theta_max=_theta_j,
        E_iso=float(energy_iso_erg),
        lf=float(initial_lorentz_factor),
        theta_j=_theta_j,
        duration=None if shell_duration_s is None else float(shell_duration_s),
        magnetar=magnetar,
        spreading=bool(spreading),
    )
    jet.energy_iso = lambda phi, theta: jet.E_iso if theta < jet.theta_j else 0.0
    jet.gamma0 = lambda phi, theta: jet.lf if theta < jet.theta_j else 1.0
    return jet


def gaussian_jet(
    energy_iso_erg: float,
    initial_lorentz_factor: float,
    core_angle_rad: float,
    outer_angle_rad: float,
    shell_duration_s: float | None,
    magnetar: Optional[Magnetar],
    spreading: bool,
) -> JetProfile:
    jet = JetProfile(
        kind="gaussian",
        theta_max=float(outer_angle_rad),
        E_iso=float(energy_iso_erg),
        lf=float(initial_lorentz_factor),
        theta_c=float(core_angle_rad),
        duration=None if shell_duration_s is None else float(shell_duration_s),
        magnetar=magnetar,
        spreading=bool(spreading),
    )
    jet.energy_iso = lambda phi, theta: jet.E_iso * np.exp(-0.5 * (theta / jet.theta_c) ** 2)
    jet.gamma0 = lambda phi, theta: 1.0 + (jet.lf - 1.0) * np.exp(-0.5 * (theta / jet.theta_c) ** 2)
    return jet


def power_law_jet(
    energy_iso_erg: float,
    initial_lorentz_factor: float,
    core_angle_rad: float,
    energy_index: float,
    lorentz_index: float | None,
    outer_angle_rad: float,
    shell_duration_s: float | None,
    magnetar: Optional[Magnetar],
    spreading: bool,
) -> JetProfile:
    _k_e = float(energy_index)
    _k_g = float(_k_e if lorentz_index is None else lorentz_index)
    jet = JetProfile(
        kind="powerlaw",
        theta_max=float(outer_angle_rad),
        E_iso=float(energy_iso_erg),
        lf=float(initial_lorentz_factor),
        theta_c=float(core_angle_rad),
        k_e=_k_e,
        k_g=_k_g,
        duration=None if shell_duration_s is None else float(shell_duration_s),
        magnetar=magnetar,
        spreading=bool(spreading),
    )
    jet.energy_iso = lambda phi, theta: (
        jet.E_iso if theta <= jet.theta_c
        else jet.E_iso * (theta / jet.theta_c) ** (-jet.k_e)
    )
    jet.gamma0 = lambda phi, theta: (
        jet.lf if theta <= jet.theta_c
        else 1.0 + (jet.lf - 1.0) * (theta / jet.theta_c) ** (-jet.k_g)
    )
    return jet


Jet = JetProfile


class Observer:
    def __init__(
        self,
        z: float,
        viewing_angle_rad: float,
        viewing_azimuth_rad: float,
        luminosity_distance_cm: Optional[float],
    ) -> None:
        self.z = float(z)
        self.theta_obs = float(viewing_angle_rad)
        self.phi_obs = float(viewing_azimuth_rad)
        self.lumi_dist_cm = None if luminosity_distance_cm is None else float(luminosity_distance_cm)

    @property
    def viewing_angle_rad(self) -> float:
        return self.theta_obs

    @viewing_angle_rad.setter
    def viewing_angle_rad(self, value: float) -> None:
        self.theta_obs = float(value)

    @property
    def viewing_azimuth_rad(self) -> float:
        return self.phi_obs

    @viewing_azimuth_rad.setter
    def viewing_azimuth_rad(self, value: float) -> None:
        self.phi_obs = float(value)

    @property
    def luminosity_distance_cm(self) -> float | None:
        return self.lumi_dist_cm

    @luminosity_distance_cm.setter
    def luminosity_distance_cm(self, value: float | None) -> None:
        self.lumi_dist_cm = None if value is None else float(value)


class Radiation:
    def __init__(
        self,
        epsilon_e: float,
        epsilon_B: float,
        p: float,
        *,
        proton_energy_fraction: float,
        epsilon_b_floor: float | None = None,
        magnetic_decay_alpha_t: float,
        magnetic_decay_t0_s: float,
        accelerated_electron_fraction: float,
        thermal_electrons: bool,
        include_ssc: bool,
        include_kn_correction: bool,
        proton_synch: bool,
        include_pgamma: bool,
        bethe_heitler: bool,
        hadronic_inverse_compton: bool,
        pp: bool,
        neutrino: bool,
        acceleration_efficiency: float,
        reverse_proton_energy_fraction: float,
        pgamma_scheme: str,
        pair_production: bool,
    ) -> None:
        self.eps_e = float(epsilon_e)
        self.eps_B = float(epsilon_B)
        self.p = float(p)
        self.epsilon_p = float(proton_energy_fraction)
        self.epsilon_b_floor = epsilon_b_floor
        self.magnetic_decay_alpha_t = float(magnetic_decay_alpha_t)
        self.magnetic_decay_t0_s = float(magnetic_decay_t0_s)
        self.xi_N = float(accelerated_electron_fraction)
        self.thermal_electrons = bool(thermal_electrons)
        self.ssc = bool(include_ssc)
        self.kn = bool(include_kn_correction)
        self.proton_synch = bool(proton_synch)
        self.pg = bool(include_pgamma)
        self.bethe_heitler = bool(bethe_heitler)
        self.hadronic_inverse_compton = bool(hadronic_inverse_compton)
        self.pp = bool(pp)
        self.neutrino = bool(neutrino)
        self.eta_acc = float(acceleration_efficiency)
        self.reverse_epsilon_p = float(reverse_proton_energy_fraction)
        self.pgamma_scheme = str(pgamma_scheme)
        self.pair_production = bool(pair_production)

    @property
    def epsilon_e(self) -> float:
        return self.eps_e

    @epsilon_e.setter
    def epsilon_e(self, value: float) -> None:
        self.eps_e = float(value)

    @property
    def epsilon_B(self) -> float:
        return self.eps_B

    @epsilon_B.setter
    def epsilon_B(self, value: float) -> None:
        self.eps_B = float(value)

    @property
    def accelerated_electron_fraction(self) -> float:
        return self.xi_N

    @accelerated_electron_fraction.setter
    def accelerated_electron_fraction(self, value: float) -> None:
        self.xi_N = float(value)

    @property
    def include_ssc(self) -> bool:
        return self.ssc

    @include_ssc.setter
    def include_ssc(self, value: bool) -> None:
        self.ssc = bool(value)

    @property
    def include_kn_correction(self) -> bool:
        return self.kn

    @include_kn_correction.setter
    def include_kn_correction(self, value: bool) -> None:
        self.kn = bool(value)

    @property
    def include_pgamma(self) -> bool:
        return self.pg

    @include_pgamma.setter
    def include_pgamma(self, value: bool) -> None:
        self.pg = bool(value)

    @property
    def proton_energy_fraction(self) -> float:
        return self.epsilon_p

    @proton_energy_fraction.setter
    def proton_energy_fraction(self, value: float) -> None:
        self.epsilon_p = float(value)

    @property
    def acceleration_efficiency(self) -> float:
        return self.eta_acc

    @acceleration_efficiency.setter
    def acceleration_efficiency(self, value: float) -> None:
        self.eta_acc = float(value)


@dataclass
class Numerics:
    num_radius: int
    num_theta: int
    num_phi: int
    num_observer_time: int
    num_electron_gamma: int
    num_photon_frequency: int
    num_chi: int | None
    num_threads: int
    electron_adaptive_substeps: bool
    electron_substep_rtol: float
    electron_substep_min: int
    electron_substep_max: int
    initial_radius_cm: float


@dataclass
class ObserverGrid:
    time_min_s: float
    time_max_s: float


@dataclass
class SolverOptions:
    electron_solver: str
    dynamics_solver: str
    geometry_projection: str
    electron_photon_coupling: str
    ssc_cooling_mode: str
    synchrotron_integration: str
    cooling_kernel: str
    radiation_kernel: str
    structured_backend: str
    patch_sampling: str
    patch_projection: str
    patch_sampling_pilot_theta: int
    patch_sampling_num_times: int
    patch_sampling_beaming_factor: float
    patch_sampling_beaming_resolution: float
    structured_parallel_mode: str
    structured_outer_threads: int | None
    structured_inner_threads: int | None
    fullhide2d_transport_model: str
    fullhide2d_stochastic_accel_norm: float
    fullhide2d_escape_mode: str


@dataclass
class ReverseShock:
    enabled: bool
    shell_duration_s: float
    upstream_sigma: float
    include_cross_zone_ic: bool
    include_ssc: bool


@dataclass
class Hadronic:
    enabled: bool
    solver: str
    num_proton_gamma: int
    num_neutrino_frequency: int
    pgamma_scheme: str
    pair_cascade_iterations: int


class _ActivePatch(NamedTuple):
    phi_center: float
    theta_center: float
    half_angle: float
    domega: float
    theta_v: float
    e_iso: float
    gamma0: float


@dataclass
class _RuntimeSetups:
    medium: Optional[str] = None
    jet: Optional[str] = None
    z: float = 0.4
    theta_obs: float = 0.0
    phi_obs: float = 0.0
    lumi_dist: Optional[float] = None
    rvs_shock: bool = False
    fwd_ssc: bool = False
    rvs_ssc: bool = False
    ssc_cooling: bool = True
    kn: bool = False
    num_threads: int = field(default_factory=default_num_threads)
    num_gam_e: int = 201
    num_nu: int = 201
    num_r: int = 300
    num_theta: int = 300
    num_phi: int = 1
    num_tobs: int = 200
    patch_theta: int = 12
    patch_phi: int = 24
    observer_time_min_s: float = 1.0e2
    observer_time_max_s: float = 1.0e8
    initial_radius_cm: float = 1.0e14
    r_tr: float = 1.0e18
    f_jump: float = 1.0
    f_wide: float = 0.1
    jump_r_cm: tuple[float, ...] = field(default_factory=tuple)
    jump_factor: tuple[float, ...] = field(default_factory=tuple)
    jump_width_log10: tuple[float, ...] = field(default_factory=tuple)
    density_profile_radius_cm: tuple[float, ...] = field(default_factory=tuple)
    density_profile_n_cm3: tuple[float, ...] = field(default_factory=tuple)
    reverse_delta_t_s: float = 10.0
    reverse_sigma: float = 0.0
    include_cross_zone_ic: bool = False
    weno5: bool = False
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
    structured_outer_threads: Optional[int] = None
    structured_inner_threads: Optional[int] = None
    num_chi: Optional[int] = None
    fullhide2d_transport_model: str = "legacy"
    fullhide2d_stochastic_accel_norm: float = 0.0
    fullhide2d_escape_mode: str = "closed"
    electron_adaptive_substeps: bool = False
    electron_substep_rtol: float = 2.0e-2
    electron_substep_min: int = 100
    electron_substep_max: int = 1000
    hadronic_enabled: bool = False
    hadronic_solver: str = "legacy_1d"
    num_gam_p: int = 161
    num_nu_nu: int = 121
    pgamma_scheme: str = "disabled"
    pair_cascade_iterations: int = 1
    index_dyn: int = 3
    index_syn_integr: int = 2
    clean: bool = True


@dataclass
class FluxPair:
    sync: np.ndarray
    ssc: np.ndarray

@dataclass
class CharTrack:
    t_obs: np.ndarray
    radius: np.ndarray
    Gamma: np.ndarray
    N_p: np.ndarray
    Doppler: np.ndarray
    B_comv: np.ndarray
    nu_m: np.ndarray
    nu_c: np.ndarray
    nu_a: np.ndarray
    nu_M: np.ndarray
    cooling_timescale_s: Optional[np.ndarray] = None
    dynamical_timescale_s: Optional[np.ndarray] = None
    gamma_e: Optional[np.ndarray] = None
    dN_dgamma_e: Optional[np.ndarray] = None
    dN_dgamma_e_bh: Optional[np.ndarray] = None
    dN_dgamma_e_chi: Optional[np.ndarray] = None
    chi_grid: Optional[np.ndarray] = None
    l_syn_spec_chi: Optional[np.ndarray] = None
    seed_syn_chi: Optional[np.ndarray] = None
    tau_syn_chi: Optional[np.ndarray] = None
    chi_radius_cm: Optional[np.ndarray] = None
    chi_gamma_bulk: Optional[np.ndarray] = None
    chi_dvolume_weight: Optional[np.ndarray] = None
    gamma_p: Optional[np.ndarray] = None
    dN_dgamma_p: Optional[np.ndarray] = None
    gamma_secondary: Optional[np.ndarray] = None
    dN_dgamma_n: Optional[np.ndarray] = None
    dN_dgamma_pi_plus: Optional[np.ndarray] = None
    dN_dgamma_pi_minus: Optional[np.ndarray] = None
    dN_dgamma_mu_minus_left: Optional[np.ndarray] = None
    dN_dgamma_mu_minus_right: Optional[np.ndarray] = None
    dN_dgamma_mu_plus_left: Optional[np.ndarray] = None
    dN_dgamma_mu_plus_right: Optional[np.ndarray] = None
    hadronic_gamma: Optional[np.ndarray] = None
    hadronic_pion_synch: Optional[np.ndarray] = None
    hadronic_muon_synch: Optional[np.ndarray] = None
    hadronic_pion_inverse_compton: Optional[np.ndarray] = None
    hadronic_muon_inverse_compton: Optional[np.ndarray] = None
    neutrino_frequency_hz: Optional[np.ndarray] = None
    neutrino_luminosity: Optional[np.ndarray] = None
    seed_frequency_hz: Optional[np.ndarray] = None
    l_had_syn_spec: Optional[np.ndarray] = None
    l_had_pg_gamma_spec: Optional[np.ndarray] = None
    l_had_bethe_heitler_spec: Optional[np.ndarray] = None
    l_had_hadronic_ic_spec: Optional[np.ndarray] = None
    am3_process_power: Optional[np.ndarray] = None
    tau_pg: Optional[np.ndarray] = None
    tau_bh: Optional[np.ndarray] = None
    pg_photon_survival: Optional[np.ndarray] = None
    secondary_rs_event_active: Optional[np.ndarray] = None
    secondary_rs_start_radius: Optional[np.ndarray] = None
    secondary_rs_shock_end_radius: Optional[np.ndarray] = None
    secondary_rs_start_tobs_axis: Optional[np.ndarray] = None
    secondary_rs_shock_end_tobs_axis: Optional[np.ndarray] = None
    secondary_rs_gamma_contact: Optional[np.ndarray] = None
    secondary_rs_pressure_3: Optional[np.ndarray] = None
    secondary_rs_gamma_43: Optional[np.ndarray] = None
    secondary_rs_beta_rs: Optional[np.ndarray] = None
    secondary_rs_u_diss: Optional[np.ndarray] = None
    secondary_rs_dissipated_energy_erg: Optional[np.ndarray] = None
    secondary_rs_electron_injected_energy_erg: Optional[np.ndarray] = None
    secondary_rs_swept_mass_g: Optional[np.ndarray] = None
    secondary_rs_internal_energy_erg: Optional[np.ndarray] = None
    secondary_rs_comoving_volume_cm3: Optional[np.ndarray] = None
    secondary_rs_pressure_total: Optional[np.ndarray] = None
    secondary_rs_enthalpy_density_total: Optional[np.ndarray] = None
    secondary_rs_branch_swept_mass_g: Optional[np.ndarray] = None
    secondary_rs_branch_internal_energy_erg: Optional[np.ndarray] = None
    secondary_rs_branch_comoving_volume_cm3: Optional[np.ndarray] = None
    secondary_rs_branch_B: Optional[np.ndarray] = None
    secondary_rs_B: Optional[np.ndarray] = None
    secondary_rs_gamma_e: Optional[np.ndarray] = None
    secondary_rs_dN_dgamma_e: Optional[np.ndarray] = None
    secondary_rs_branch_gamma_m: Optional[np.ndarray] = None
    secondary_rs_branch_gamma_contact: Optional[np.ndarray] = None
    secondary_rs_branch_gamma_43: Optional[np.ndarray] = None
    secondary_rs_branch_compression: Optional[np.ndarray] = None
    secondary_rs_branch_beta_rs: Optional[np.ndarray] = None
    secondary_rs_branch_u_diss: Optional[np.ndarray] = None
    secondary_rs_branch_reacceleration_seed_energy_erg: Optional[np.ndarray] = None
    secondary_rs_branch_reaccelerated_energy_erg: Optional[np.ndarray] = None
    secondary_rs_branch_luminosity_syn: Optional[np.ndarray] = None
    secondary_rs_nu_m: Optional[np.ndarray] = None
    secondary_rs_nu_c: Optional[np.ndarray] = None
    secondary_rs_nu_a: Optional[np.ndarray] = None
    timings: Optional[dict] = None


@dataclass
class TrackBundle:
    fwd: CharTrack
    rev: Optional[CharTrack]
    patches: list[dict[str, float]]

@dataclass
class FluxResult:
    total: np.ndarray
    fwd: FluxPair
    rev: FluxPair
    cross_ic: Optional[np.ndarray]

    @property
    def shape(self):
        return self.total.shape

    def __array__(self, dtype=None):
        return np.asarray(self.total, dtype=dtype)

    def __getitem__(self, item):
        return self.total[item]


@dataclass
class AdaptiveFluxResult:
    time_s: np.ndarray
    flux: FluxResult

    @property
    def total(self) -> np.ndarray:
        return self.flux.total

    @property
    def fwd(self) -> FluxPair:
        return self.flux.fwd

    @property
    def rev(self) -> FluxPair:
        return self.flux.rev

    @property
    def cross_ic(self) -> Optional[np.ndarray]:
        return self.flux.cross_ic


@dataclass
class SkyImage:
    image: np.ndarray
    extent: np.ndarray
    pixel_solid_angle: float
    pixel_size: float
    direct_flux: np.ndarray
    rendered_flux: np.ndarray
    normalization_scale: np.ndarray
    x_centroid: np.ndarray
    y_centroid: np.ndarray

    @property
    def shape(self):
        return self.image.shape

    @property
    def flux_ratio(self) -> np.ndarray:
        return self.rendered_flux / np.maximum(self.direct_flux, np.finfo(float).tiny)

    @property
    def centroid(self) -> np.ndarray:
        return np.column_stack((self.x_centroid, self.y_centroid))


@dataclass
class PolarizationResult:
    I_sync: np.ndarray
    Q: np.ndarray
    U: np.ndarray
    linear_polarization: np.ndarray
    polarization_angle_rad: np.ndarray
    components: dict[str, dict[str, np.ndarray]]

    @property
    def shape(self):
        return self.I_sync.shape


def make_flux_density_entry(
    times_s: np.ndarray,
    frequencies_hz: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray,
) -> dict:
    return dict(
        times_s=np.asarray(times_s, dtype=float),
        frequencies_hz=np.asarray(frequencies_hz, dtype=float),
        flux=np.asarray(flux, dtype=float),
        flux_err=np.asarray(flux_err, dtype=float),
    )


def make_spectrum_entry(
    time_s: float,
    frequencies_hz: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray,
) -> dict:
    return dict(
        time_s=float(time_s),
        frequencies_hz=np.asarray(frequencies_hz, dtype=float),
        flux=np.asarray(flux, dtype=float),
        flux_err=np.asarray(flux_err, dtype=float),
    )


def make_flux_entry(
    nu_min_hz: float,
    nu_max_hz: float,
    time_s: float,
    flux: float,
    flux_err: float,
    num_points: int = 64,
) -> dict:
    return dict(
        nu_min_hz=float(nu_min_hz),
        nu_max_hz=float(nu_max_hz),
        time_s=float(time_s),
        flux=float(flux),
        flux_err=float(flux_err),
        num_points=int(num_points),
    )


def make_empty_obs() -> dict:
    return dict(flux_density=[], spectrum=[], flux=[])


def _compose_runtime_setups(
    *,
    numerics: Numerics | None,
    observer_grid: ObserverGrid | None,
    solver_options: SolverOptions | None,
    reverse_shock: ReverseShock | None,
    hadronic: Hadronic | None,
) -> _RuntimeSetups:
    if any(value is None for value in (numerics, observer_grid, solver_options, reverse_shock, hadronic)):
        raise TypeError("Model requires explicit numerics, observer_grid, solver_options, reverse_shock, and hadronic.")
    result = deepcopy(_RuntimeSetups())
    result.num_r = int(numerics.num_radius)
    result.num_theta = int(numerics.num_theta)
    result.num_phi = int(numerics.num_phi)
    result.num_tobs = int(numerics.num_observer_time)
    result.num_gam_e = int(numerics.num_electron_gamma)
    result.num_nu = int(numerics.num_photon_frequency)
    result.num_chi = numerics.num_chi
    result.num_threads = int(numerics.num_threads)
    result.electron_adaptive_substeps = bool(numerics.electron_adaptive_substeps)
    result.electron_substep_rtol = float(numerics.electron_substep_rtol)
    result.electron_substep_min = int(numerics.electron_substep_min)
    result.electron_substep_max = int(numerics.electron_substep_max)
    result.initial_radius_cm = float(numerics.initial_radius_cm)
    result.observer_time_min_s = float(observer_grid.time_min_s)
    result.observer_time_max_s = float(observer_grid.time_max_s)
    ssc_mode = str(solver_options.ssc_cooling_mode).lower()
    if ssc_mode == "none":
        index_y, result.ssc_cooling, result.kn = 0, False, False
    elif ssc_mode == "numeric_ic_kn":
        index_y, result.ssc_cooling, result.kn = 1, True, True
    elif ssc_mode == "nakar_y_thomson":
        index_y, result.ssc_cooling, result.kn = 2, True, False
    else:
        raise ValueError("ssc_cooling_mode must be 'none', 'numeric_ic_kn', or 'nakar_y_thomson'.")
    result.electron_solver = str(solver_options.electron_solver)
    result.dynamics_kernel = str(solver_options.dynamics_solver)
    result.geometry_kernel = str(solver_options.geometry_projection)
    result.electron_photon_coupling = str(solver_options.electron_photon_coupling)
    if str(solver_options.synchrotron_integration).lower() != "fixed_grid":
        raise ValueError("synchrotron_integration currently supports only 'fixed_grid'.")
    result.index_syn_integr = 2
    result.cooling_kernel = str(solver_options.cooling_kernel)
    result.radiation_kernel = str(solver_options.radiation_kernel)
    result.structured_backend = str(solver_options.structured_backend)
    result.patch_sampling = str(solver_options.patch_sampling)
    result.patch_projection = str(solver_options.patch_projection)
    result.patch_sampling_pilot_theta = int(solver_options.patch_sampling_pilot_theta)
    result.patch_sampling_num_times = int(solver_options.patch_sampling_num_times)
    result.patch_sampling_beaming_factor = float(solver_options.patch_sampling_beaming_factor)
    result.patch_sampling_beaming_resolution = float(solver_options.patch_sampling_beaming_resolution)
    result.structured_parallel_mode = str(solver_options.structured_parallel_mode)
    result.structured_outer_threads = solver_options.structured_outer_threads
    result.structured_inner_threads = solver_options.structured_inner_threads
    result.fullhide2d_transport_model = str(solver_options.fullhide2d_transport_model)
    result.fullhide2d_stochastic_accel_norm = float(solver_options.fullhide2d_stochastic_accel_norm)
    result.fullhide2d_escape_mode = str(solver_options.fullhide2d_escape_mode)
    result.weno5 = False
    result._index_y_override = index_y
    result.rvs_shock = bool(reverse_shock.enabled)
    result.reverse_delta_t_s = float(reverse_shock.shell_duration_s)
    result.reverse_sigma = float(reverse_shock.upstream_sigma)
    result.include_cross_zone_ic = bool(reverse_shock.include_cross_zone_ic)
    result.rvs_ssc = bool(reverse_shock.include_ssc)
    result.hadronic_enabled = bool(hadronic.enabled)
    result.hadronic_solver = str(hadronic.solver)
    result.num_gam_p = int(hadronic.num_proton_gamma)
    result.num_nu_nu = int(hadronic.num_neutrino_frequency)
    result.pgamma_scheme = str(hadronic.pgamma_scheme)
    result.pair_cascade_iterations = int(hadronic.pair_cascade_iterations)
    return result



class Model:
    def __init__(
        self,
        *,
        jet: Jet,
        medium: Medium,
        observer: Observer,
        fwd_rad: Radiation,
        rvs_rad: Optional[Radiation] = None,
        numerics: Numerics,
        observer_grid: ObserverGrid,
        solver_options: SolverOptions,
        reverse_shock: ReverseShock,
        hadronic: Hadronic,
    ) -> None:
        self.medium = medium
        self.jet = jet
        self.observer = observer
        self.fwd_rad = fwd_rad
        self.rvs_rad = rvs_rad
        self.setups = _compose_runtime_setups(
            numerics=numerics,
            observer_grid=observer_grid,
            solver_options=solver_options,
            reverse_shock=reverse_shock,
            hadronic=hadronic,
        )
        self._last_details: Optional[TrackBundle] = None
        self._raw_cache: dict[tuple[str, str], tuple[FluxResult, TrackBundle]] = {}
        self._details_cache: dict[str, TrackBundle] = {}

    def flux_density_grid(
        self,
        times_s: np.ndarray,
        nu_hz: np.ndarray,
        *,
        projection_kind: str = "lightcurve",
    ) -> FluxResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        return self._compute_raw(times_s, nu_hz, projection_kind=projection_kind)

    def flux_density_grid_adaptive(
        self,
        times_s: np.ndarray,
        nu_hz: np.ndarray,
        *,
        projection_kind: str = "lightcurve",
    ) -> AdaptiveFluxResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        if self.jet.kind != "tophat" or not self._supports_direct_kernel():
            raise NotImplementedError("adaptive observer-time grid currently supports direct top-hat kernel-backed models.")
        from .api_adaptive import _adaptive_observer_time_grid

        self._ensure_direct_adaptive_eats_resolution()
        self._compute_raw(times_s, nu_hz, projection_kind=projection_kind)
        adaptive_times = _adaptive_observer_time_grid(self, times_s)
        flux = self._compute_raw(
            adaptive_times,
            nu_hz,
            solve_reference_times_s=adaptive_times,
            projection_kind=projection_kind,
        )
        return AdaptiveFluxResult(time_s=adaptive_times, flux=flux)

    def flux_density(
        self,
        times_s: np.ndarray,
        nu_hz: np.ndarray,
        *,
        projection_kind: str = "lightcurve",
    ) -> FluxResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        def pack(result: FluxResult, extract: Callable[[np.ndarray], np.ndarray]) -> FluxResult:
            pair_total = extract(result.total)
            pair_fwd_sync = extract(result.fwd.sync)
            pair_fwd_ssc = extract(result.fwd.ssc)
            pair_rev_sync = extract(result.rev.sync)
            pair_rev_ssc = extract(result.rev.ssc)
            pair_cross_ic = None if result.cross_ic is None else extract(result.cross_ic)
            return FluxResult(
                total=pair_total,
                fwd=FluxPair(sync=pair_fwd_sync, ssc=pair_fwd_ssc),
                rev=FluxPair(sync=pair_rev_sync, ssc=pair_rev_ssc),
                cross_ic=pair_cross_ic,
            )

        if times_s.ndim == 1 and nu_hz.ndim == 1 and times_s.shape == nu_hz.shape:
            unique_freqs, inverse = np.unique(nu_hz, return_inverse=True)
            result = self.flux_density_grid(times_s, unique_freqs, projection_kind=projection_kind)
            pair_index = np.arange(times_s.shape[0], dtype=int)
            return pack(result, lambda values: values[inverse, pair_index])

        result = self._compute_raw(times_s, nu_hz, projection_kind=projection_kind)
        return pack(result, lambda values: _extract_pair_flux(values, times_s, nu_hz))

    def flux_density_exposures(
        self,
        times_s: np.ndarray,
        nu_hz: np.ndarray,
        exposures_s: np.ndarray,
        num_subsamples: int = 4,
    ) -> FluxResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        exposures_s = np.asarray(exposures_s, dtype=float)
        if times_s.shape != nu_hz.shape or times_s.shape != exposures_s.shape:
            raise ValueError("times_s, nu_hz, and exposures_s must have the same shape.")
        if int(num_subsamples) <= 0:
            raise ValueError("num_subsamples must be positive.")
        from .api_adaptive import _adaptive_exposure_average

        return _adaptive_exposure_average(self, times_s, nu_hz, exposures_s, int(num_subsamples))

    def spectrum(self, time_s: float, nu_hz: np.ndarray, *, projection_kind: str = "sed") -> np.ndarray:
        return self.flux_density_grid(np.array([time_s], dtype=float), nu_hz, projection_kind=projection_kind).total[:, 0]

    def sky_image(self, t_obs: np.ndarray | float, nu_obs: float, fov: float, npixel: int = 128) -> SkyImage:
        times_s = np.atleast_1d(np.asarray(t_obs, dtype=float))
        if times_s.size == 0:
            raise ValueError("t_obs array must be non-empty.")
        if np.any(times_s <= 0.0):
            raise ValueError("t_obs must be positive.")
        if not np.isfinite(nu_obs) or nu_obs <= 0.0:
            raise ValueError("nu_obs must be positive.")
        if not np.isfinite(fov) or fov <= 0.0:
            raise ValueError("fov must be positive.")
        if int(npixel) <= 0:
            raise ValueError("npixel must be positive.")
        from .api_observe import _render_sky_image

        return _render_sky_image(self, times_s, float(nu_obs), float(fov), int(npixel))

    def polarization(
        self,
        times_s: np.ndarray,
        nu_hz: np.ndarray,
        *,
        magnetic_geometry: str = "shock_random",
        local_emissivity: str = "analytic_then_kernel",
    ) -> PolarizationResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        from .api_observe import _compute_polarization

        return _compute_polarization(
            self,
            times_s,
            nu_hz,
            magnetic_geometry=magnetic_geometry,
            local_emissivity=local_emissivity,
        )

    def flux(
        self,
        time_s: np.ndarray | float,
        nu_min_hz: float,
        nu_max_hz: float,
        num_points: int = 64,
        *,
        projection_kind: str = "sed",
    ):
        times_s = np.atleast_1d(np.asarray(time_s, dtype=float))
        nu_hz = np.logspace(np.log10(nu_min_hz), np.log10(nu_max_hz), num_points)
        result = self.flux_density_grid(times_s, nu_hz, projection_kind=projection_kind)
        total = np.trapezoid(result.total, nu_hz, axis=0)
        fwd_sync = np.trapezoid(result.fwd.sync, nu_hz, axis=0)
        fwd_ssc = np.trapezoid(result.fwd.ssc, nu_hz, axis=0)
        rev_sync = np.trapezoid(result.rev.sync, nu_hz, axis=0)
        rev_ssc = np.trapezoid(result.rev.ssc, nu_hz, axis=0)
        cross_ic = None if result.cross_ic is None else np.trapezoid(result.cross_ic, nu_hz, axis=0)
        band_result = FluxResult(
            total=total,
            fwd=FluxPair(sync=fwd_sync, ssc=fwd_ssc),
            rev=FluxPair(sync=rev_sync, ssc=rev_ssc),
            cross_ic=cross_ic,
        )
        if np.asarray(time_s).ndim == 0:
            return float(band_result.total[0])
        return band_result

    def details(self, t_min: Optional[float] = None, t_max: Optional[float] = None) -> TrackBundle:
        if t_min is not None or t_max is not None:
            t1 = self.setups.observer_time_min_s if t_min is None else float(t_min)
            t2 = self.setups.observer_time_max_s if t_max is None else float(t_max)
            num_tobs = self._detail_time_count(t1, t2)
            times = np.logspace(np.log10(t1), np.log10(t2), num_tobs)
            self._last_details = self._compute_details_only(times)
        elif self._last_details is None:
            self._last_details = self._compute_details_only(self.default_detail_times())
        return self._last_details

    def jet_E_iso(self, phi: float, theta: float) -> float:
        return self.jet.energy_iso(phi, theta)

    def jet_Gamma0(self, phi: float, theta: float) -> float:
        return self.jet.gamma0(phi, theta)

    def medium_rho(self, phi: float, theta: float, r_cm: float) -> float:
        return self.medium.density(phi, theta, r_cm)

    def component_fluxes(
        self,
        times_s: np.ndarray,
        nu_hz: np.ndarray,
        *,
        projection_kind: str = "lightcurve",
    ) -> FluxResult:
        return self._compute_raw(times_s, nu_hz, projection_kind=projection_kind)

    def flux_density_bands(self, times_s: np.ndarray, *, projection_kind: str = "lightcurve") -> np.ndarray:
        frequencies_hz = build_multiband_observer_frequencies()[1]
        total_matrix = self.flux_density_grid(times_s, frequencies_hz, projection_kind=projection_kind).total
        return combine_multiband_flux(total_matrix, frequencies_hz, 8)

    def default_times(self) -> np.ndarray:
        return np.logspace(
            np.log10(self.setups.observer_time_min_s),
            np.log10(self.setups.observer_time_max_s),
            self.setups.num_tobs,
        )

    def default_detail_times(self) -> np.ndarray:
        return np.logspace(
            np.log10(self.setups.observer_time_min_s),
            np.log10(self.setups.observer_time_max_s),
            self._detail_time_count(self.setups.observer_time_min_s, self.setups.observer_time_max_s),
        )

    def default_frequencies(self) -> np.ndarray:
        return np.logspace(np.log10(1.0e9), np.log10(1.0e24), 64)

    def _compute_raw(
        self,
        times_s: np.ndarray,
        nu_hz: np.ndarray,
        *,
        solve_reference_times_s: np.ndarray | None = None,
        projection_kind: str = "lightcurve",
    ) -> FluxResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        from asgard_core.asgard_state import _normalize_projection_kind
        from .api_adaptive import _array_signature, _observe_parts, _remember_cache_entry

        projection_kind = _normalize_projection_kind(projection_kind)
        reference_signature = None
        if solve_reference_times_s is not None:
            reference_signature = _array_signature(np.asarray(solve_reference_times_s, dtype=float))
        cache_key = (_array_signature(times_s), _array_signature(nu_hz), reference_signature, projection_kind)
        cached = self._raw_cache.get(cache_key)
        if cached is not None:
            self._last_details = cached[1]
            return cached[0]
        if self.jet.kind == "tophat" and self._supports_direct_kernel():
            config = _direct_tophat_patch_config(self)
            state = _solve_patch_state(
                self,
                config,
                times_s,
                nu_hz,
                solve_reference_times_s=solve_reference_times_s,
            )
            observed = _observe_parts(state, times_s, nu_hz, projection_kind=projection_kind)
            details = _make_details(state.components, patches=[{"phi": 0.0, "theta": 0.0, "weight": 1.0}], state=state)
            model_result = (observed, details)
        else:
            model_result = _solve_patch_model(
                self,
                times_s,
                nu_hz,
                solve_reference_times_s=solve_reference_times_s,
                projection_kind=projection_kind,
            )
        self._last_details = model_result[1]
        _remember_cache_entry(self._raw_cache, cache_key, model_result)
        return model_result[0]

    def _compute_details_only(self, times_s: np.ndarray) -> TrackBundle:
        times_s = np.asarray(times_s, dtype=float)
        from .api_adaptive import _array_signature, _remember_cache_entry

        cache_key = _array_signature(times_s)
        cached = self._details_cache.get(cache_key)
        if cached is not None:
            self._last_details = cached
            return cached
        if self.jet.kind == "tophat" and self._supports_direct_kernel():
            config = _direct_tophat_patch_config(self)
            state = _solve_patch_state(self, config, times_s, None)
            details = _make_details(state.components, patches=[{"phi": 0.0, "theta": 0.0, "weight": 1.0}], state=state)
        else:
            details = _patch_details(self, times_s)
        self._last_details = details
        _remember_cache_entry(self._details_cache, cache_key, details)
        return details

    def _supports_direct_kernel(self) -> bool:
        return self.medium.kind in ("ism", "wind", "density_profile")

    def _ensure_direct_adaptive_eats_resolution(self) -> None:
        if self.jet.kind != "tophat" or not self._supports_direct_kernel():
            return
        if float(self.observer.theta_obs) == 0.0:
            return
        gamma = max(float(self.jet.lf), 1.0)
        theta_j = float(self.jet.theta_j)
        factor = float(self.setups.patch_sampling_beaming_factor)
        resolution = float(self.setups.patch_sampling_beaming_resolution)
        required_theta = int(np.ceil(theta_j * resolution * gamma / factor)) + 1
        required_phi = int(np.ceil(np.pi * resolution * gamma * np.sin(theta_j) / factor)) + 1
        required_phi = max(required_phi, 2)
        if required_theta <= int(self.setups.num_theta) and required_phi <= int(self.setups.num_phi):
            return
        self.setups.num_theta = max(int(self.setups.num_theta), required_theta)
        self.setups.num_phi = max(int(self.setups.num_phi), required_phi)
        self._last_details = None
        self._raw_cache.clear()
        self._details_cache.clear()

    def _detail_time_count(self, t_min: float, t_max: float) -> int:
        if t_min <= 0.0 or t_max <= 0.0 or t_max <= t_min:
            return int(self.setups.num_tobs)
        decades = np.log10(float(t_max) / float(t_min))
        density_jumps = self.setups.jump_r_cm and self.setups.jump_factor and self.setups.jump_width_log10
        density_profile = self.setups.density_profile_radius_cm and self.setups.density_profile_n_cm3
        samples_per_decade = 96.0 if self.setups.rvs_shock and (density_jumps or density_profile) else 8.0
        return max(int(self.setups.num_tobs), int(np.ceil(samples_per_decade * decades)))


def _solve_patch_state(
    model: Model,
    config: RuntimeConfig,
    times_s: np.ndarray,
    requested_frequencies_hz: np.ndarray | None,
    timings: Optional[dict[str, float]] = None,
    solve_reference_times_s: np.ndarray | None = None,
) -> SolveState:
    solve_reference = times_s if solve_reference_times_s is None else np.asarray(solve_reference_times_s, dtype=float)
    base_count = max(int(model.setups.num_tobs), int(np.unique(solve_reference).size))
    solve_times_s = make_tgrid(solve_reference, base_count)
    if solve_reference.size > 1:
        solve_t_min = min(float(model.setups.observer_time_min_s), float(np.min(solve_reference)))
        solve_t_max_requested = float(np.max(solve_reference))
        solve_count = max(base_count, model._detail_time_count(solve_t_min, solve_t_max_requested))
        if solve_t_max_requested <= solve_t_min:
            solve_t_max = max(float(model.setups.observer_time_max_s), solve_t_min * 10.0)
        else:
            log_t_min = np.log10(solve_t_min)
            log_t_max_requested = np.log10(solve_t_max_requested)
            if solve_count <= 2:
                log_t_max = log_t_max_requested
            else:
                log_step = (log_t_max_requested - log_t_min) / float(solve_count - 2)
                log_t_max = log_t_max_requested + log_step
            solve_t_max = 10.0**log_t_max
        solve_times_s = np.logspace(np.log10(solve_t_min), np.log10(solve_t_max), solve_count)
    query_config = make_query_cfg(config, solve_times_s)
    query_config.num_r = max(int(query_config.num_r), int(solve_times_s.size))
    setup = make_query_setup(query_config, solve_times_s, requested_frequencies_hz)
    return solve_state_from_setup(
        query_config,
        setup,
        timings=timings,
        requested_frequencies_hz=requested_frequencies_hz,
    )


def _direct_tophat_patch_config(model: Model) -> RuntimeConfig:
    return _build_fit_config_for_patch(
        model,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )


def _solve_patch_model(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    solve_reference_times_s: np.ndarray | None = None,
    projection_kind: str = "lightcurve",
) -> tuple[FluxResult, TrackBundle]:
    structured_backend = str(model.setups.structured_backend).lower()
    patch_sampling = str(model.setups.patch_sampling).lower()
    from asgard_core.angular_sampling import SUPPORTED_PATCH_SAMPLING

    if patch_sampling not in SUPPORTED_PATCH_SAMPLING:
        raise ValueError(
            f"Unknown patch_sampling={patch_sampling!r}; expected one of {SUPPORTED_PATCH_SAMPLING}."
        )
    if structured_backend != "python_patch":
        if patch_sampling != "uniform":
            raise NotImplementedError(
                "patch_sampling='dominant_region_ioka_v1' is only supported by "
                "structured_backend='python_patch'."
            )
        if solve_reference_times_s is not None:
            raise NotImplementedError("structured_backend='fortran_1d' does not yet support external solve_reference_times_s.")
        from asgard_core.structured_jet_kernel import solve_structured_jet_fortran

        return solve_structured_jet_fortran(model, times_s, nu_hz, _build_fit_config_for_patch)

    return _solve_patch_model_python(
        model,
        times_s,
        nu_hz,
        solve_reference_times_s=solve_reference_times_s,
        projection_kind=projection_kind,
    )


def _empty_patch_flux_accumulator(nu_hz: np.ndarray, times_s: np.ndarray) -> FluxResult:
    shape = (nu_hz.shape[0], times_s.shape[0])
    return FluxResult(
        total=np.zeros(shape, dtype=float),
        fwd=FluxPair(sync=np.zeros(shape, dtype=float), ssc=np.zeros(shape, dtype=float)),
        rev=FluxPair(sync=np.zeros(shape, dtype=float), ssc=np.zeros(shape, dtype=float)),
        cross_ic=np.zeros(shape, dtype=float),
    )


def _accumulate_patch_flux(accumulator: FluxResult, observed: FluxResult) -> None:
    accumulator.total += observed.total
    accumulator.fwd.sync += observed.fwd.sync
    accumulator.fwd.ssc += observed.fwd.ssc
    accumulator.rev.sync += observed.rev.sync
    accumulator.rev.ssc += observed.rev.ssc
    if observed.cross_ic is not None:
        accumulator.cross_ic += observed.cross_ic


def _solve_patch_model_python(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    solve_reference_times_s: np.ndarray | None = None,
    projection_kind: str = "lightcurve",
) -> tuple[FluxResult, TrackBundle]:
    from .api_adaptive import _observe_parts

    flux_accumulator = _empty_patch_flux_accumulator(nu_hz, times_s)
    patches_meta: list[dict[str, float]] = []
    details_fwd: Optional[CharTrack] = None
    details_rev: Optional[CharTrack] = None

    patch_sampling = str(model.setups.patch_sampling).lower()
    patch_projection = _resolve_patch_projection(model, patch_sampling)
    if patch_projection == "surface_element" and is_axisymmetric_jet(model.jet):
        return _solve_axisymmetric_surface_patch_model_python(
            model,
            times_s,
            nu_hz,
            solve_reference_times_s=solve_reference_times_s,
        )
    for patch, state in _iter_solved_patch_elements(
        model,
        times_s,
        nu_hz,
        _iter_patch_elements(model, times_s),
        solve_reference_times_s=solve_reference_times_s,
    ):
        if patch_projection == "tophat_cell":
            observed = _observe_parts(state, times_s, nu_hz, projection_kind=projection_kind)
        else:
            observed = _observe_surface_element_parts(state, times_s, nu_hz, patch.domega)
        _accumulate_patch_flux(flux_accumulator, observed)
        patches_meta.append(_patch_metadata(patch, patch_sampling, patch_projection))
        if details_fwd is None:
            details = _make_details(state.components, patches_meta, state=state)
            details_fwd = details.fwd
            details_rev = details.rev

    if details_fwd is None:
        raise ValueError("No active jet patches were found for the requested structured jet.")
    return (
        flux_accumulator,
        TrackBundle(fwd=details_fwd, rev=details_rev, patches=patches_meta),
    )


def _solve_axisymmetric_surface_patch_model_python(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    solve_reference_times_s: np.ndarray | None = None,
) -> tuple[FluxResult, TrackBundle]:
    grid = build_patch_grid(model, times_s)
    flux_accumulator = _empty_patch_flux_accumulator(nu_hz, times_s)
    patches_meta: list[dict[str, float]] = []
    details_fwd: Optional[CharTrack] = None
    details_rev: Optional[CharTrack] = None
    patch_sampling = str(model.setups.patch_sampling).lower()

    for i_theta, theta_center_value in enumerate(grid.theta_centers):
        theta_center = float(theta_center_value)
        e_iso = model.jet.energy_iso(0.0, theta_center)
        gamma0 = model.jet.gamma0(0.0, theta_center)
        if e_iso <= 0.0 or gamma0 <= 1.0:
            continue
        config = _build_fit_config_for_patch(
            model,
            theta_v=0.0,
            opening_angle_jet=float(grid.patch_half_angle[i_theta, 0]),
            e_iso=e_iso,
            gamma0=gamma0,
            theta_center=theta_center,
        )
        state = _solve_patch_state(
            model,
            config,
            times_s,
            nu_hz,
            solve_reference_times_s=solve_reference_times_s,
        )
        boundary = state.setup.boundary
        for i_phi, phi_center_value in enumerate(grid.phi_centers):
            phi_center = float(phi_center_value)
            theta_v = float(angular_separation(theta_center, phi_center, model.observer.theta_obs, model.observer.phi_obs))
            patch = _ActivePatch(
                phi_center,
                theta_center,
                float(grid.patch_half_angle[i_theta, i_phi]),
                float(grid.domega[i_theta, i_phi]),
                theta_v,
                float(e_iso),
                float(gamma0),
            )
            boundary[9] = patch.theta_v
            observed = _observe_surface_element_parts(state, times_s, nu_hz, patch.domega)
            _accumulate_patch_flux(flux_accumulator, observed)
            patches_meta.append(_patch_metadata(patch, patch_sampling, "surface_element"))
        if details_fwd is None:
            details = _make_details(state.components, patches_meta, state=state)
            details_fwd = details.fwd
            details_rev = details.rev

    if details_fwd is None:
        raise ValueError("No active jet patches were found for the requested structured jet.")
    return (
        flux_accumulator,
        TrackBundle(fwd=details_fwd, rev=details_rev, patches=patches_meta),
    )


def _patch_details(model: Model, times_s: np.ndarray) -> TrackBundle:
    patches_meta: list[dict[str, float]] = []
    first_component: Optional[FluxComponents] = None
    first_details: Optional[TrackBundle] = None

    patch_sampling = str(model.setups.patch_sampling).lower()
    patch_projection = _resolve_patch_projection(model, patch_sampling)
    for patch in _active_patch_elements(model, _iter_patch_elements(model, times_s)):
        patches_meta.append(_patch_metadata(patch, patch_sampling, patch_projection))
        if first_component is None:
            config = _patch_runtime_config(model, patch)
            first_state = _solve_patch_state(model, config, times_s, None)
            first_component = first_state.components
            first_details = _make_details(first_component, patches_meta, state=first_state)

    if first_component is None or first_details is None:
        raise ValueError("No active jet patches were found for the requested structured jet.")
    return first_details


def _resolve_patch_projection(model: Model, patch_sampling: str) -> str:
    projection = str(model.setups.patch_projection).lower()
    if projection == "auto":
        return "tophat_cell" if patch_sampling == "uniform" else "surface_element"
    if projection in {"tophat_cell", "surface_element"}:
        return projection
    raise ValueError("patch_projection must be 'auto', 'tophat_cell', or 'surface_element'.")


def _observe_surface_element_parts(
    state: SolveState,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    domega: float,
) -> FluxResult:
    from asgard_core.asgard_state import _build_observer_setup_from_state

    setup = _build_observer_setup_from_state(state, times_s)
    components = state.components
    fwd_sync = _project_surface_element(
        setup,
        components.fwd.characteristic_time_s,
        components.fwd.gamma,
        components.fwd.radius_cm,
        setup.seed_frequency_hz,
        components.fwd_sync,
        nu_hz,
        domega,
    )
    fwd_ssc = _project_optional_surface_element(setup, components, components.fwd_ssc, nu_hz, domega)
    rev_sync = _project_optional_surface_element(setup, components, components.rev_sync, nu_hz, domega)
    rev_ssc = _project_optional_surface_element(setup, components, components.rev_ssc, nu_hz, domega)
    cross_ic = _project_optional_surface_element(setup, components, components.cross_ic, nu_hz, domega)
    total = fwd_sync + fwd_ssc + rev_sync + rev_ssc + cross_ic
    return FluxResult(
        total=total,
        fwd=FluxPair(sync=fwd_sync, ssc=fwd_ssc),
        rev=FluxPair(sync=rev_sync, ssc=rev_ssc),
        cross_ic=cross_ic,
    )


def _project_optional_surface_element(
    setup,
    components: FluxComponents,
    source: np.ndarray | None,
    nu_hz: np.ndarray,
    domega: float,
) -> np.ndarray:
    if source is None:
        return np.zeros((nu_hz.shape[0], setup.observer_time_s.shape[0]), dtype=float)
    return _project_surface_element(
        setup,
        components.fwd.characteristic_time_s,
        components.fwd.gamma,
        components.fwd.radius_cm,
        setup.seed_frequency_hz,
        source,
        nu_hz,
        domega,
    )


def _project_surface_element(
    setup,
    characteristic_time_s: np.ndarray,
    gamma: np.ndarray,
    radius_cm: np.ndarray,
    seed_frequency_hz: np.ndarray,
    absorbed_spectral_flux: np.ndarray,
    frequencies_hz: np.ndarray,
    domega: float,
) -> np.ndarray:
    if not np.any(absorbed_spectral_flux):
        return np.zeros((frequencies_hz.shape[0], setup.observer_time_s.shape[0]), dtype=float)
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    order = np.argsort(frequencies_hz)
    sorted_frequencies = frequencies_hz[order]
    flux_sorted = Interpolation.sed_interpolation_surface_element(
        setup.boundary,
        characteristic_time_s,
        gamma,
        radius_cm,
        absorbed_spectral_flux,
        seed_frequency_hz,
        sorted_frequencies,
        setup.observer_time_s,
        float(domega),
    )
    if np.array_equal(order, np.arange(order.shape[0])):
        return flux_sorted
    flux_matrix = np.empty_like(flux_sorted)
    flux_matrix[order] = flux_sorted
    return flux_matrix


def _extract_pair_flux(grid: np.ndarray, times_s: np.ndarray, frequencies_hz: np.ndarray) -> np.ndarray:
    if grid.ndim != 2:
        return np.asarray(grid, dtype=float)
    if frequencies_hz.shape == times_s.shape:
        unique_freqs, inverse = np.unique(frequencies_hz, return_inverse=True)
        if unique_freqs.shape[0] == grid.shape[0]:
            return grid[inverse, np.arange(times_s.shape[0])]
    if grid.shape == (frequencies_hz.shape[0], times_s.shape[0]):
        return np.diag(grid)
    raise ValueError("Cannot extract pairwise flux from the provided grid.")


def _build_fit_config_for_patch(
    model: Model,
    *,
    theta_v: float,
    opening_angle_jet: float,
    e_iso: float,
    gamma0: float,
    theta_center: Optional[float] = None,
) -> RuntimeConfig:
    if model.jet.spreading:
        raise NotImplementedError("Jet spreading is not implemented in the current ASGARD backend.")
    reverse_delta_t = model.setups.reverse_delta_t_s if model.jet.duration is None else float(model.jet.duration)
    index_y = model.setups._index_y_override
    config = RuntimeConfig(
        num_threads=model.setups.num_threads,
        num_gam_e=model.setups.num_gam_e,
        num_nu=model.setups.num_nu,
        num_r=model.setups.num_r,
        num_theta=model.setups.num_theta,
        num_phi=model.setups.num_phi,
        num_tobs=model.setups.num_tobs,
        t_obs_min_log10=float(np.log10(model.setups.observer_time_min_s)),
        t_obs_max_log10=float(np.log10(model.setups.observer_time_max_s)),
        index_dyn=model.setups.index_dyn,
        index_syn_integr=model.setups.index_syn_integr,
        electron_solver=model.setups.electron_solver,
        cooling_kernel=model.setups.cooling_kernel,
        radiation_kernel=model.setups.radiation_kernel,
        dynamics_kernel=model.setups.dynamics_kernel,
        geometry_kernel=model.setups.geometry_kernel,
        electron_photon_coupling=model.setups.electron_photon_coupling,
        structured_backend=model.setups.structured_backend,
        patch_sampling=model.setups.patch_sampling,
        patch_projection=model.setups.patch_projection,
        patch_sampling_pilot_theta=model.setups.patch_sampling_pilot_theta,
        patch_sampling_num_times=model.setups.patch_sampling_num_times,
        patch_sampling_beaming_factor=model.setups.patch_sampling_beaming_factor,
        patch_sampling_beaming_resolution=model.setups.patch_sampling_beaming_resolution,
        structured_parallel_mode=model.setups.structured_parallel_mode,
        structured_outer_threads=model.setups.structured_outer_threads,
        structured_inner_threads=model.setups.structured_inner_threads,
        num_chi=model.setups.num_chi,
        fullhide2d_transport_model=model.setups.fullhide2d_transport_model,
        fullhide2d_stochastic_accel_norm=model.setups.fullhide2d_stochastic_accel_norm,
        fullhide2d_escape_mode=model.setups.fullhide2d_escape_mode,
        electron_adaptive_substeps=model.setups.electron_adaptive_substeps,
        electron_substep_rtol=model.setups.electron_substep_rtol,
        electron_substep_min=model.setups.electron_substep_min,
        electron_substep_max=model.setups.electron_substep_max,
        weno5=model.setups.weno5,
        z=model.observer.z,
        theta_v=theta_v,
        opening_angle_jet=opening_angle_jet,
        e_iso=e_iso,
        eta_0=gamma0,
        epsilon_e=model.fwd_rad.eps_e,
        epsilon_b=model.fwd_rad.eps_B,
        epsilon_b_floor=model.fwd_rad.epsilon_b_floor,
        magnetic_decay_alpha_t=model.fwd_rad.magnetic_decay_alpha_t,
        magnetic_decay_t0_s=model.fwd_rad.magnetic_decay_t0_s,
        p=model.fwd_rad.p,
        f_e=model.fwd_rad.xi_N,
        thermal_electrons=bool(model.fwd_rad.thermal_electrons),
        index_y=index_y,
        include_forward_ssc=model.fwd_rad.ssc,
        luminosity_distance_cm_override=model.observer.lumi_dist_cm,
        initial_radius_cm=model.setups.initial_radius_cm,
        spectrum_output=SpectrumOutputConfig(enabled=False),
    )
    config.hadronic = HadronicConfig(
        enabled=bool(model.setups.hadronic_enabled and model.fwd_rad.epsilon_p > 0.0),
        solver=str(model.setups.hadronic_solver),
        epsilon_p=float(model.fwd_rad.epsilon_p),
        p_p=float(model.fwd_rad.p),
        eta_acc=float(model.fwd_rad.eta_acc),
        num_gam_p=int(model.setups.num_gam_p),
        include_proton_synch=bool(model.fwd_rad.proton_synch),
        include_pg=bool(model.fwd_rad.pg),
        include_bethe_heitler=bool(model.fwd_rad.bethe_heitler),
        include_hadronic_inverse_compton=bool(model.fwd_rad.hadronic_inverse_compton),
        include_pp=bool(model.fwd_rad.pp),
        include_neutrino=bool(model.fwd_rad.neutrino),
        include_pair_production=bool(model.fwd_rad.pair_production),
        pgamma_scheme=str(model.setups.pgamma_scheme if model.setups.pgamma_scheme != "disabled" else model.fwd_rad.pgamma_scheme),
        num_nu_nu=int(model.setups.num_nu_nu),
        reverse_enabled=bool(model.setups.rvs_shock and model.fwd_rad.reverse_epsilon_p > 0.0),
        reverse_epsilon_p=float(model.fwd_rad.reverse_epsilon_p),
        pair_cascade_iterations=int(model.setups.pair_cascade_iterations),
        pp_model=int(getattr(model.fwd_rad, "pp_model", -1)),
    )
    magnetar = model.jet.magnetar
    if magnetar is not None and _jet_magnetar_active(model.jet, 0.0 if theta_center is None else theta_center):
        config.l_inj_0 = float(magnetar.L0)
        config.e_inj_t2 = float(magnetar.t0)
        config.q_inj = float(magnetar.q)
    config.reverse = model.rvs_rad is not None
    if model.rvs_rad is not None:
        config.reverse_shock = ReverseShockConfig(
            enabled=True,
            delta_t_s=reverse_delta_t,
            epsilon_e=model.rvs_rad.eps_e,
            epsilon_b=model.rvs_rad.eps_B,
            p=model.rvs_rad.p,
            f_e=model.rvs_rad.xi_N,
            sigma=model.setups.reverse_sigma,
            include_ssc=model.rvs_rad.ssc,
            include_cross_zone_ic=model.setups.include_cross_zone_ic,
        )
    if model.medium.kind not in ("ism", "wind", "density_profile"):
        raise NotImplementedError("User-defined Medium is not supported by the current ASGARD kernel.")
    for key, value in model.medium.to_kernel_params().items():
        setattr(config, key, value)
    if model.medium.density_profile_radius_cm or model.medium.density_profile_n_cm3:
        config.density_profile_radius_cm = tuple(float(value) for value in model.medium.density_profile_radius_cm)
        config.density_profile_n_cm3 = tuple(float(value) for value in model.medium.density_profile_n_cm3)
    if model.setups.density_profile_radius_cm or model.setups.density_profile_n_cm3:
        config.density_profile_radius_cm = tuple(float(value) for value in model.setups.density_profile_radius_cm)
        config.density_profile_n_cm3 = tuple(float(value) for value in model.setups.density_profile_n_cm3)
    if model.medium.kind == "ism" and float(model.setups.f_jump) != 1.0:
        config.r_tr = float(model.setups.r_tr)
        config.f_jump = float(model.setups.f_jump)
        config.f_wide = float(model.setups.f_wide)
    if (
        model.setups.jump_r_cm
        or model.setups.jump_factor
        or model.setups.jump_width_log10
    ):
        config.jump_r_cm = tuple(float(value) for value in model.setups.jump_r_cm)
        config.jump_factor = tuple(float(value) for value in model.setups.jump_factor)
        config.jump_width_log10 = tuple(float(value) for value in model.setups.jump_width_log10)
    return config


def _as_float_array_or_none(value):
    return None if value is None else np.asarray(value, dtype=float)


def _secondary_rs_detail_kwargs(secondary_rs) -> dict:
    if secondary_rs is None:
        return {}
    return dict(
        secondary_rs_event_active=secondary_rs.event_active,
        secondary_rs_start_radius=secondary_rs.start_radius_cm,
        secondary_rs_shock_end_radius=secondary_rs.shock_end_radius_cm,
        secondary_rs_start_tobs_axis=secondary_rs.start_tobs_axis_s,
        secondary_rs_shock_end_tobs_axis=secondary_rs.shock_end_tobs_axis_s,
        secondary_rs_gamma_contact=secondary_rs.gamma_contact,
        secondary_rs_pressure_3=secondary_rs.pressure_3,
        secondary_rs_gamma_43=secondary_rs.gamma_43,
        secondary_rs_beta_rs=secondary_rs.beta_rs,
        secondary_rs_u_diss=secondary_rs.dissipated_energy_density,
        secondary_rs_dissipated_energy_erg=secondary_rs.dissipated_energy_erg,
        secondary_rs_electron_injected_energy_erg=secondary_rs.electron_injected_energy_erg,
        secondary_rs_swept_mass_g=secondary_rs.swept_mass_g,
        secondary_rs_internal_energy_erg=secondary_rs.internal_energy_erg,
        secondary_rs_comoving_volume_cm3=secondary_rs.comoving_volume_cm3,
        secondary_rs_pressure_total=secondary_rs.pressure_total,
        secondary_rs_enthalpy_density_total=secondary_rs.enthalpy_density_total,
        secondary_rs_branch_swept_mass_g=secondary_rs.branch_swept_mass_g,
        secondary_rs_branch_internal_energy_erg=secondary_rs.branch_internal_energy_erg,
        secondary_rs_branch_comoving_volume_cm3=secondary_rs.branch_comoving_volume_cm3,
        secondary_rs_branch_B=secondary_rs.branch_magnetic_field_g,
        secondary_rs_B=secondary_rs.magnetic_field_g,
        secondary_rs_gamma_e=secondary_rs.gam_e,
        secondary_rs_dN_dgamma_e=secondary_rs.d_n_gam_e,
        secondary_rs_branch_gamma_m=secondary_rs.branch_gamma_m,
        secondary_rs_branch_gamma_contact=secondary_rs.branch_gamma_contact,
        secondary_rs_branch_gamma_43=secondary_rs.branch_gamma_43,
        secondary_rs_branch_compression=secondary_rs.branch_compression,
        secondary_rs_branch_beta_rs=secondary_rs.branch_beta_rs,
        secondary_rs_branch_u_diss=secondary_rs.branch_dissipated_energy_density,
        secondary_rs_branch_reacceleration_seed_energy_erg=secondary_rs.branch_reacceleration_seed_energy_erg,
        secondary_rs_branch_reaccelerated_energy_erg=secondary_rs.branch_reaccelerated_energy_erg,
        secondary_rs_branch_luminosity_syn=secondary_rs.branch_luminosity_syn,
        secondary_rs_nu_m=secondary_rs.nu_m,
        secondary_rs_nu_c=secondary_rs.nu_c,
        secondary_rs_nu_a=secondary_rs.nu_a,
    )


def _make_details(
    components: FluxComponents,
    patches: list[dict[str, float]],
    state: Optional[SolveState] = None,
) -> TrackBundle:
    electron = None if state is None else state.electron
    hadronic = None if state is None else state.hadronic
    reverse_shock = None if state is None else state.dynamics.reverse_shock
    reverse_emission = None if state is None else state.reverse_emission

    fwd_gamma_e = None if electron is None else np.asarray(electron.gam_e, dtype=float)
    fwd_dnde = None if electron is None else np.asarray(electron.d_n_gam_e, dtype=float)
    fwd_dnde_bh = None if electron is None else _as_float_array_or_none(electron.d_n_gam_e_bh)
    if hadronic is not None and hadronic.d_n_gam_e_bh is not None:
        fwd_dnde_bh = np.asarray(hadronic.d_n_gam_e_bh, dtype=float)
    fwd_dnde_chi = None if electron is None else _as_float_array_or_none(electron.d_n_gam_e_chi)
    fwd_chi_grid = None if electron is None else _as_float_array_or_none(electron.chi_grid)
    fwd_lsyn_chi = None if electron is None else _as_float_array_or_none(electron.l_syn_spec_chi)
    fwd_seed_chi = None if electron is None else _as_float_array_or_none(electron.seed_syn_chi)
    fwd_tau_chi = None if electron is None else _as_float_array_or_none(electron.tau_syn_chi)
    fwd_chi_radius = None if electron is None else _as_float_array_or_none(electron.chi_radius_cm)
    fwd_chi_gamma = None if electron is None else _as_float_array_or_none(electron.chi_gamma_bulk)
    fwd_chi_weight = None if electron is None else _as_float_array_or_none(electron.chi_dvolume_weight)
    fwd_gamma_p = None if hadronic is None else np.asarray(hadronic.gam_p, dtype=float)
    fwd_dndp = None if hadronic is None else np.asarray(hadronic.d_n_gam_p, dtype=float)
    fwd_gamma_secondary = None if hadronic is None else _as_float_array_or_none(hadronic.gam_secondary)
    fwd_dndn = None if hadronic is None else _as_float_array_or_none(hadronic.d_n_gam_n)
    fwd_dndpi_plus = None if hadronic is None else _as_float_array_or_none(hadronic.d_n_gam_pi_plus)
    fwd_dndpi_minus = None if hadronic is None else _as_float_array_or_none(hadronic.d_n_gam_pi_minus)
    fwd_dndmu_ml = None if hadronic is None else _as_float_array_or_none(hadronic.d_n_gam_mu_minus_left)
    fwd_dndmu_mr = None if hadronic is None else _as_float_array_or_none(hadronic.d_n_gam_mu_minus_right)
    fwd_dndmu_pl = None if hadronic is None else _as_float_array_or_none(hadronic.d_n_gam_mu_plus_left)
    fwd_dndmu_pr = None if hadronic is None else _as_float_array_or_none(hadronic.d_n_gam_mu_plus_right)
    fwd_had_gamma = None
    if hadronic is not None:
        fwd_had_gamma = np.asarray(hadronic.l_had_syn_spec + hadronic.l_had_pg_gamma, dtype=float)
        if hadronic.l_had_pion_synch is not None:
            fwd_had_gamma = fwd_had_gamma + np.asarray(hadronic.l_had_pion_synch, dtype=float)
        if hadronic.l_had_muon_synch is not None:
            fwd_had_gamma = fwd_had_gamma + np.asarray(hadronic.l_had_muon_synch, dtype=float)
        if hadronic.l_had_pion_inverse_compton is not None:
            fwd_had_gamma = fwd_had_gamma + np.asarray(hadronic.l_had_pion_inverse_compton, dtype=float)
        if hadronic.l_had_muon_inverse_compton is not None:
            fwd_had_gamma = fwd_had_gamma + np.asarray(hadronic.l_had_muon_inverse_compton, dtype=float)
    fwd_had_pi_syn = None if hadronic is None else _as_float_array_or_none(hadronic.l_had_pion_synch)
    fwd_had_mu_syn = None if hadronic is None else _as_float_array_or_none(hadronic.l_had_muon_synch)
    fwd_had_pi_ic = None if hadronic is None else _as_float_array_or_none(hadronic.l_had_pion_inverse_compton)
    fwd_had_mu_ic = None if hadronic is None else _as_float_array_or_none(hadronic.l_had_muon_inverse_compton)
    fwd_nu_freq = None
    fwd_nu_lum = None
    if hadronic is not None and state.config.hadronic.include_neutrino:
        fwd_nu_freq = np.asarray(hadronic.neutrino_frequency_hz, dtype=float)
        fwd_nu_lum = np.asarray(hadronic.neutrino_luminosity, dtype=float)
    fwd_had_syn = None
    fwd_had_pg_gamma = None
    fwd_had_bh = None
    fwd_had_hic = None
    fwd_am3_power = None
    fwd_tau_pg = None
    fwd_tau_bh = None
    fwd_pg_survival = None
    fwd_timings = None
    fwd_seed_freq = None
    if hadronic is not None:
        fwd_seed_freq = np.asarray(state.photon_field.seed_frequency_hz, dtype=float)
        fwd_had_syn = np.asarray(hadronic.l_had_syn_spec, dtype=float)
        fwd_had_pg_gamma = np.asarray(hadronic.l_had_pg_gamma, dtype=float)
        fwd_had_bh = _as_float_array_or_none(hadronic.l_had_bethe_heitler)
        fwd_had_hic = _as_float_array_or_none(hadronic.l_had_hadronic_inverse_compton)
        fwd_am3_power = _as_float_array_or_none(hadronic.am3_process_power)
        fwd_tau_pg = _as_float_array_or_none(hadronic.tau_pg)
        fwd_tau_bh = _as_float_array_or_none(hadronic.tau_bh)
        fwd_pg_survival = _as_float_array_or_none(hadronic.pg_photon_survival)
        fwd_timings = dict(hadronic.timings) if hadronic.timings else {}
    rev_gamma_e = None
    rev_dnde = None
    if reverse_shock is not None:
        rev_gamma_e = np.asarray(reverse_shock.gam_e, dtype=float)
        rev_dnde = np.asarray(reverse_shock.d_n_gam_e, dtype=float)
    secondary_rs = None if reverse_emission is None else reverse_emission.secondary_rs
    return TrackBundle(
        fwd=CharTrack(
            t_obs=components.fwd.characteristic_time_s,
            radius=components.fwd.radius_cm,
            Gamma=components.fwd.gamma,
            N_p=np.asarray(components.fwd.swept_mass_g, dtype=float) / constants.para_m_p,
            Doppler=components.fwd.doppler,
            B_comv=components.fwd.magnetic_field_g,
            nu_m=components.fwd.nu_m,
            nu_c=components.fwd.nu_c,
            nu_a=components.fwd.nu_a,
            nu_M=components.fwd.nu_M,
            cooling_timescale_s=components.fwd.cooling_timescale_s,
            dynamical_timescale_s=components.fwd.dynamical_timescale_s,
            gamma_e=fwd_gamma_e,
            dN_dgamma_e=fwd_dnde,
            dN_dgamma_e_bh=fwd_dnde_bh,
            dN_dgamma_e_chi=fwd_dnde_chi,
            chi_grid=fwd_chi_grid,
            l_syn_spec_chi=fwd_lsyn_chi,
            seed_syn_chi=fwd_seed_chi,
            tau_syn_chi=fwd_tau_chi,
            chi_radius_cm=fwd_chi_radius,
            chi_gamma_bulk=fwd_chi_gamma,
            chi_dvolume_weight=fwd_chi_weight,
            gamma_p=fwd_gamma_p,
            dN_dgamma_p=fwd_dndp,
            gamma_secondary=fwd_gamma_secondary,
            dN_dgamma_n=fwd_dndn,
            dN_dgamma_pi_plus=fwd_dndpi_plus,
            dN_dgamma_pi_minus=fwd_dndpi_minus,
            dN_dgamma_mu_minus_left=fwd_dndmu_ml,
            dN_dgamma_mu_minus_right=fwd_dndmu_mr,
            dN_dgamma_mu_plus_left=fwd_dndmu_pl,
            dN_dgamma_mu_plus_right=fwd_dndmu_pr,
            hadronic_gamma=fwd_had_gamma,
            hadronic_pion_synch=fwd_had_pi_syn,
            hadronic_muon_synch=fwd_had_mu_syn,
            hadronic_pion_inverse_compton=fwd_had_pi_ic,
            hadronic_muon_inverse_compton=fwd_had_mu_ic,
            neutrino_frequency_hz=fwd_nu_freq,
            neutrino_luminosity=fwd_nu_lum,
            seed_frequency_hz=fwd_seed_freq,
            l_had_syn_spec=fwd_had_syn,
            l_had_pg_gamma_spec=fwd_had_pg_gamma,
            l_had_bethe_heitler_spec=fwd_had_bh,
            l_had_hadronic_ic_spec=fwd_had_hic,
            am3_process_power=fwd_am3_power,
            tau_pg=fwd_tau_pg,
            tau_bh=fwd_tau_bh,
            pg_photon_survival=fwd_pg_survival,
            timings=fwd_timings,
        ),
        rev=None
        if components.rev is None
        else CharTrack(
            t_obs=components.rev.characteristic_time_s,
            radius=components.rev.radius_cm,
            Gamma=components.rev.gamma,
            N_p=np.asarray(components.rev.swept_mass_g, dtype=float) / constants.para_m_p,
            Doppler=components.rev.doppler,
            B_comv=components.rev.magnetic_field_g,
            nu_m=components.rev.nu_m,
            nu_c=components.rev.nu_c,
            nu_a=components.rev.nu_a,
            nu_M=components.rev.nu_M,
            cooling_timescale_s=components.rev.cooling_timescale_s,
            dynamical_timescale_s=components.rev.dynamical_timescale_s,
            gamma_e=rev_gamma_e,
            dN_dgamma_e=rev_dnde,
            **_secondary_rs_detail_kwargs(secondary_rs),
        ),
        patches=patches,
    )


def _iter_patch_elements(model: Model, observer_time_s: np.ndarray | None = None):
    grid = build_patch_grid(model, observer_time_s)
    for i_theta, theta_center in enumerate(grid.theta_centers):
        for i_phi, phi_center in enumerate(grid.phi_centers):
            yield (
                float(phi_center),
                float(theta_center),
                float(grid.patch_half_angle[i_theta, i_phi]),
                float(grid.domega[i_theta, i_phi]),
            )


def _active_patch_elements(model: Model, patch_elements):
    for phi_center, theta_center, patch_half_angle, domega in patch_elements:
        e_iso = model.jet.energy_iso(phi_center, theta_center)
        gamma0 = model.jet.gamma0(phi_center, theta_center)
        if e_iso <= 0.0 or gamma0 <= 1.0:
            continue
        theta_v = angular_separation(theta_center, phi_center, model.observer.theta_obs, model.observer.phi_obs)
        yield _ActivePatch(
            float(phi_center),
            float(theta_center),
            float(patch_half_angle),
            float(domega),
            float(theta_v),
            float(e_iso),
            float(gamma0),
        )


def _patch_runtime_config(
    model: Model,
    patch: _ActivePatch,
    opening_angle_jet: float | None = None,
) -> RuntimeConfig:
    return _build_fit_config_for_patch(
        model,
        theta_v=patch.theta_v,
        opening_angle_jet=patch.half_angle if opening_angle_jet is None else float(opening_angle_jet),
        e_iso=patch.e_iso,
        gamma0=patch.gamma0,
        theta_center=patch.theta_center,
    )


def _iter_solved_patch_elements(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray | None,
    patch_elements,
    *,
    opening_angle_jet: float | None = None,
    timings: Optional[dict[str, float]] = None,
    solve_reference_times_s: np.ndarray | None = None,
):
    for patch in _active_patch_elements(model, patch_elements):
        state = _solve_patch_state(
            model,
            _patch_runtime_config(model, patch, opening_angle_jet=opening_angle_jet),
            times_s,
            nu_hz,
            timings=timings,
            solve_reference_times_s=solve_reference_times_s,
        )
        yield patch, state


def _patch_metadata(patch: _ActivePatch, patch_sampling: str, patch_projection: str) -> dict[str, float | str]:
    return {
        "phi": patch.phi_center,
        "theta": patch.theta_center,
        "theta_v": patch.theta_v,
        "half_angle": patch.half_angle,
        "domega": patch.domega,
        "patch_sampling": patch_sampling,
        "patch_projection": patch_projection,
        "E_iso": patch.e_iso,
        "Gamma0": patch.gamma0,
    }


def _jet_magnetar_active(jet: JetProfile, theta_center: float) -> bool:
    if jet.kind in ("tophat", "gaussian", "powerlaw", "steppowerlaw"):
        return theta_center <= jet.theta_c
    if jet.kind == "twocomponent":
        return theta_center <= jet.theta_n
    return True
