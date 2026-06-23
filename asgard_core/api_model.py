from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, NamedTuple, Optional

import numpy as np

from asgard_core.asgard_config import (
    RuntimeConfig,
    HadronicConfig,
    MAX_DENSITY_PROFILE_POINTS,
    ReverseShockConfig,
    SpectrumOutputConfig,
    default_num_threads,
)
from asgard_core.asgard_state import (
    FluxComponents,
    SolveState,
    make_query_cfg,
    make_query_setup,
    make_tgrid,
    project_flux_grid,
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
    """Unified jet profile."""
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
    theta_table: tuple[float, ...] = field(default_factory=tuple)
    e_iso_table: tuple[float, ...] = field(default_factory=tuple)
    gamma0_table: tuple[float, ...] = field(default_factory=tuple)
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


def tabulated_angular_jet(
    theta_rad,
    energy_iso_erg,
    lorentz_factor,
    shell_duration_s: float | None,
    magnetar: Optional[Magnetar],
    spreading: bool,
) -> JetProfile:
    theta = tuple(float(value) for value in theta_rad)
    e_iso = tuple(float(value) for value in energy_iso_erg)
    gamma0 = tuple(float(value) for value in lorentz_factor)
    _validate_tabulated_angular_jet(theta, e_iso, gamma0)
    theta_arr = np.asarray(theta, dtype=float)
    log_e_iso = np.log(np.asarray(e_iso, dtype=float))
    gamma0_arr = np.asarray(gamma0, dtype=float)
    jet = JetProfile(
        kind="tabulated",
        theta_max=float(theta[-1]),
        E_iso=float(np.max(np.asarray(e_iso, dtype=float))),
        lf=float(np.max(gamma0_arr)),
        theta_table=theta,
        e_iso_table=e_iso,
        gamma0_table=gamma0,
        duration=None if shell_duration_s is None else float(shell_duration_s),
        magnetar=magnetar,
        spreading=bool(spreading),
    )

    def _energy_iso(_phi, theta_value):
        theta_eval = np.asarray(theta_value, dtype=float)
        log_value = np.interp(theta_eval, theta_arr, log_e_iso, left=-np.inf, right=-np.inf)
        value = np.exp(log_value)
        return float(value) if value.ndim == 0 else value

    def _gamma0(_phi, theta_value):
        theta_eval = np.asarray(theta_value, dtype=float)
        value = np.interp(theta_eval, theta_arr, gamma0_arr, left=1.0, right=1.0)
        return float(value) if value.ndim == 0 else value

    jet.energy_iso = _energy_iso
    jet.gamma0 = _gamma0
    return jet


def _validate_tabulated_angular_jet(
    theta_rad: tuple[float, ...],
    energy_iso_erg: tuple[float, ...],
    lorentz_factor: tuple[float, ...],
) -> None:
    theta = np.asarray(theta_rad, dtype=float)
    e_iso = np.asarray(energy_iso_erg, dtype=float)
    gamma0 = np.asarray(lorentz_factor, dtype=float)
    if theta.shape != e_iso.shape or theta.shape != gamma0.shape:
        raise ValueError("theta_rad, energy_iso_erg, and lorentz_factor must have the same shape.")
    if theta.size < 2:
        raise ValueError("tabulated angular jet requires at least two theta samples.")
    if not np.all(np.isfinite(theta)) or not np.all(np.isfinite(e_iso)) or not np.all(np.isfinite(gamma0)):
        raise ValueError("tabulated angular jet arrays must contain finite values.")
    if np.any(theta < 0.0) or np.any(np.diff(theta) <= 0.0):
        raise ValueError("theta_rad must be non-negative and strictly increasing.")
    if np.any(e_iso <= 0.0):
        raise ValueError("energy_iso_erg must be positive inside the tabulated active angular domain.")
    if np.any(gamma0 <= 1.0):
        raise ValueError("lorentz_factor must be greater than 1 inside the tabulated active angular domain.")


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
    structured_num_theta: int
    structured_num_phi: int
    eats_num_theta: int
    eats_num_phi: int
    downstream_num_chi: int | None
    num_observer_time: int
    num_electron_gamma: int
    num_photon_frequency: int
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
    projection_adaptive_rtol: float = 2.0e-2
    projection_adaptive_max_depth: int = 4
    nu_callback: Optional[Callable[[str, np.ndarray, np.ndarray, np.ndarray], None]] = None


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
    b_chi_g: Optional[np.ndarray] = None
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
        return self.rendered_flux / self.direct_flux

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


def _build_base_runtime_config(
    *,
    numerics: Numerics,
    observer_grid: ObserverGrid,
    solver_options: SolverOptions,
    reverse_shock: ReverseShock,
    hadronic: Hadronic,
    observer: Observer,
    fwd_rad: Radiation,
    medium: Medium,
    jet: Jet,
    rvs_rad: Radiation | None,
) -> RuntimeConfig:
    ssc_mode = str(solver_options.ssc_cooling_mode).lower()
    if ssc_mode == "none":
        index_y = 0
    elif ssc_mode == "numeric_ic_kn":
        index_y = 1
    elif ssc_mode == "nakar_y_thomson":
        index_y = 2
    else:
        raise ValueError("ssc_cooling_mode must be 'none', 'numeric_ic_kn', or 'nakar_y_thomson'.")
    synch_mode = str(solver_options.synchrotron_integration).lower()
    if synch_mode == "fixed_grid":
        index_syn_integr = 2
    elif synch_mode == "cyclotron":
        index_syn_integr = 4
    else:
        raise ValueError("synchrotron_integration must be 'fixed_grid' or 'cyclotron'.")
    n, so = numerics, solver_options
    config = RuntimeConfig(
        num_threads=int(n.num_threads), num_gam_e=int(n.num_electron_gamma),
        num_nu=int(n.num_photon_frequency), num_r=int(n.num_radius),
        eats_num_theta=int(n.eats_num_theta), eats_num_phi=int(n.eats_num_phi),
        num_tobs=int(n.num_observer_time),
        t_obs_min_log10=float(np.log10(observer_grid.time_min_s)),
        t_obs_max_log10=float(np.log10(observer_grid.time_max_s)),
        index_dyn=3, index_syn_integr=index_syn_integr,
        electron_solver=str(so.electron_solver), cooling_kernel=str(so.cooling_kernel),
        radiation_kernel=str(so.radiation_kernel), dynamics_kernel=str(so.dynamics_solver),
        geometry_kernel=str(so.geometry_projection),
        electron_photon_coupling=str(so.electron_photon_coupling),
        structured_backend=str(so.structured_backend), patch_sampling=str(so.patch_sampling),
        patch_projection=str(so.patch_projection),
        patch_sampling_pilot_theta=int(so.patch_sampling_pilot_theta),
        patch_sampling_num_times=int(so.patch_sampling_num_times),
        patch_sampling_beaming_factor=float(so.patch_sampling_beaming_factor),
        patch_sampling_beaming_resolution=float(so.patch_sampling_beaming_resolution),
        structured_parallel_mode=str(so.structured_parallel_mode),
        structured_outer_threads=so.structured_outer_threads,
        structured_inner_threads=so.structured_inner_threads,
        projection_adaptive_rtol=float(so.projection_adaptive_rtol),
        projection_adaptive_max_depth=int(so.projection_adaptive_max_depth),
        downstream_num_chi=n.downstream_num_chi,
        fullhide2d_transport_model=str(so.fullhide2d_transport_model),
        fullhide2d_stochastic_accel_norm=float(so.fullhide2d_stochastic_accel_norm),
        fullhide2d_escape_mode=str(so.fullhide2d_escape_mode),
        nu_callback=so.nu_callback,
        electron_adaptive_substeps=bool(n.electron_adaptive_substeps),
        electron_substep_rtol=float(n.electron_substep_rtol),
        electron_substep_min=int(n.electron_substep_min),
        electron_substep_max=int(n.electron_substep_max),
        z=observer.z, theta_v=observer.theta_obs, opening_angle_jet=jet.theta_max,
        e_iso=jet.energy_iso_erg, eta_0=jet.initial_lorentz_factor,
        epsilon_e=fwd_rad.eps_e, epsilon_b=fwd_rad.eps_B,
        epsilon_b_floor=fwd_rad.epsilon_b_floor,
        magnetic_decay_alpha_t=fwd_rad.magnetic_decay_alpha_t,
        magnetic_decay_t0_s=fwd_rad.magnetic_decay_t0_s,
        p=fwd_rad.p, f_e=fwd_rad.xi_N,
        thermal_electrons=bool(fwd_rad.thermal_electrons),
        index_y=index_y, include_forward_ssc=fwd_rad.ssc,
        luminosity_distance_cm_override=observer.lumi_dist_cm,
        initial_radius_cm=float(n.initial_radius_cm),
        spectrum_output=SpectrumOutputConfig(enabled=False),
    )
    config.hadronic = HadronicConfig(
        enabled=bool(hadronic.enabled and fwd_rad.epsilon_p > 0.0),
        solver=str(hadronic.solver), epsilon_p=float(fwd_rad.epsilon_p),
        p_p=float(fwd_rad.p), eta_acc=float(fwd_rad.eta_acc),
        num_gam_p=int(hadronic.num_proton_gamma),
        include_proton_synch=bool(fwd_rad.proton_synch),
        include_pg=bool(fwd_rad.pg),
        include_bethe_heitler=bool(fwd_rad.bethe_heitler),
        include_hadronic_inverse_compton=bool(fwd_rad.hadronic_inverse_compton),
        include_pp=bool(fwd_rad.pp), include_neutrino=bool(fwd_rad.neutrino),
        include_pair_production=bool(fwd_rad.pair_production),
        pgamma_scheme=str(hadronic.pgamma_scheme if hadronic.pgamma_scheme != "disabled" else fwd_rad.pgamma_scheme),
        num_nu_nu=int(hadronic.num_neutrino_frequency),
        reverse_enabled=bool(reverse_shock.enabled and fwd_rad.reverse_epsilon_p > 0.0),
        reverse_epsilon_p=float(fwd_rad.reverse_epsilon_p),
        pair_cascade_iterations=int(hadronic.pair_cascade_iterations),
        pp_model=int(getattr(fwd_rad, "pp_model", -1)),
    )
    config.reverse = rvs_rad is not None
    config.reverse_shock = ReverseShockConfig(
        enabled=bool(reverse_shock.enabled),
        delta_t_s=float(reverse_shock.shell_duration_s) if jet.duration is None else float(jet.duration),
        epsilon_e=rvs_rad.eps_e if rvs_rad is not None else fwd_rad.eps_e,
        epsilon_b=rvs_rad.eps_B if rvs_rad is not None else fwd_rad.eps_B,
        p=rvs_rad.p if rvs_rad is not None else fwd_rad.p,
        f_e=rvs_rad.xi_N if rvs_rad is not None else fwd_rad.xi_N,
        sigma=float(reverse_shock.upstream_sigma),
        include_ssc=bool(reverse_shock.include_ssc) if rvs_rad is not None else False,
        include_cross_zone_ic=bool(reverse_shock.include_cross_zone_ic),
    )
    if medium.kind not in ("ism", "wind", "density_profile"):
        raise NotImplementedError("User-defined Medium is not supported by the current ASGARD kernel.")
    for key, value in medium.to_kernel_params().items():
        setattr(config, key, value)
    if medium.density_profile_radius_cm or medium.density_profile_n_cm3:
        config.density_profile_radius_cm = tuple(float(v) for v in medium.density_profile_radius_cm)
        config.density_profile_n_cm3 = tuple(float(v) for v in medium.density_profile_n_cm3)
    return config



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
        self.setups = _build_base_runtime_config(
            numerics=numerics, observer_grid=observer_grid, solver_options=solver_options,
            reverse_shock=reverse_shock, hadronic=hadronic,
            observer=observer, fwd_rad=fwd_rad, medium=medium, jet=jet, rvs_rad=rvs_rad,
        )
        self._last_details: Optional[TrackBundle] = None
        self._raw_cache: dict[tuple[object, ...], tuple[FluxResult, TrackBundle]] = {}
        self._details_cache: dict[tuple[object, ...], TrackBundle] = {}

    def _observer_cache_signature(self) -> tuple[float, float, float, float | None]:
        distance = self.observer.lumi_dist_cm
        return (
            float(self.observer.z),
            float(self.observer.theta_obs),
            float(self.observer.phi_obs),
            None if distance is None else float(distance),
        )

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
            t1 = 10**self.setups.t_obs_min_log10 if t_min is None else float(t_min)
            t2 = 10**self.setups.t_obs_max_log10 if t_max is None else float(t_max)
            num_tobs = self._detail_time_count(t1, t2)
            times = np.logspace(np.log10(t1), np.log10(t2), num_tobs)
            self._last_details = self._compute_details_only(times)
        elif self._last_details is None:
            self._last_details = self._compute_details_only(self.default_detail_times())
        return self._last_details

    def default_times(self) -> np.ndarray:
        return np.logspace(
            self.setups.t_obs_min_log10,
            self.setups.t_obs_max_log10,
            self.setups.num_tobs,
        )

    def default_detail_times(self) -> np.ndarray:
        return np.logspace(
            self.setups.t_obs_min_log10,
            self.setups.t_obs_max_log10,
            self._detail_time_count(10**self.setups.t_obs_min_log10, 10**self.setups.t_obs_max_log10),
        )

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
        cache_key = (
            _array_signature(times_s),
            _array_signature(nu_hz),
            reference_signature,
            projection_kind,
            self._observer_cache_signature(),
        )
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

    def _total_matrix(
        self,
        times_s: np.ndarray,
        nu_hz: np.ndarray,
        *,
        timings: dict[str, float] | None = None,
        projection_kind: str = "lightcurve",
    ) -> np.ndarray:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        if self.jet.kind == "tophat" and self._supports_direct_kernel():
            config = _direct_tophat_patch_config(self)
            state = _solve_patch_state(self, config, times_s, nu_hz, timings=timings)
            observed = project_flux_grid(
                state,
                times_s,
                nu_hz,
                timings=timings,
                mode="total_only",
                projection_kind=projection_kind,
            )
            return np.asarray(observed.components["total"], dtype=float)

        result, details = _solve_patch_model(
            self,
            times_s,
            nu_hz,
            projection_kind=projection_kind,
        )
        self._last_details = details
        return np.asarray(result.total, dtype=float)

    def _compute_details_only(self, times_s: np.ndarray) -> TrackBundle:
        times_s = np.asarray(times_s, dtype=float)
        from .api_adaptive import _array_signature, _remember_cache_entry

        cache_key = (_array_signature(times_s), self._observer_cache_signature())
        cached = self._details_cache.get(cache_key)
        if cached is not None:
            self._last_details = cached
            return cached
        if self.jet.kind == "tophat" and self._supports_direct_kernel():
            config = _direct_tophat_patch_config(self)
            state = _solve_patch_state(self, config, times_s, None)
            details = _make_details(state.components, patches=[{"phi": 0.0, "theta": 0.0, "weight": 1.0}], state=state)
        else:
            _result, details = _solve_patch_model(
                self,
                times_s,
                np.array([1.0e9], dtype=float),
                projection_kind="lightcurve",
            )
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
        if required_theta <= int(self.setups.eats_num_theta) and required_phi <= int(self.setups.eats_num_phi):
            return
        self.setups.eats_num_theta = max(int(self.setups.eats_num_theta), required_theta)
        self.setups.eats_num_phi = max(int(self.setups.eats_num_phi), required_phi)
        self._last_details = None
        self._raw_cache.clear()
        self._details_cache.clear()

    def _detail_time_count(self, t_min: float, t_max: float) -> int:
        if t_min <= 0.0 or t_max <= 0.0 or t_max <= t_min:
            return int(self.setups.num_tobs)
        decades = np.log10(float(t_max) / float(t_min))
        density_jumps = self.setups.jump_r_cm and self.setups.jump_factor and self.setups.jump_width_log10
        density_profile = self.setups.density_profile_radius_cm and self.setups.density_profile_n_cm3
        samples_per_decade = 96.0 if self.setups.reverse_shock.enabled and (density_jumps or density_profile) else 8.0
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
        solve_t_min = min(float(10**model.setups.t_obs_min_log10), float(np.min(solve_reference)))
        solve_t_max_requested = float(np.max(solve_reference))
        solve_count = max(base_count, model._detail_time_count(solve_t_min, solve_t_max_requested))
        if solve_t_max_requested <= solve_t_min:
            solve_t_max = max(float(model10**self.setups.t_obs_max_log10), solve_t_min * 10.0)
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
    if solve_reference_times_s is not None:
        raise NotImplementedError(
            "structured_backend does not yet support external solve_reference_times_s."
        )
    from asgard_core.structured_jet_kernel import solve_structured_jet_fortran

    return solve_structured_jet_fortran(model, times_s, nu_hz, _build_fit_config_for_patch)


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
    config = deepcopy(model.setups)
    config.theta_v = theta_v
    config.opening_angle_jet = opening_angle_jet
    config.e_iso = e_iso
    config.eta_0 = gamma0
    if model.jet.duration is not None:
        config.reverse_shock.delta_t_s = float(model.jet.duration)
    magnetar = model.jet.magnetar
    if magnetar is not None and _jet_magnetar_active(model.jet, 0.0 if theta_center is None else theta_center):
        config.l_inj_0, config.e_inj_t2, config.q_inj = float(magnetar.L0), float(magnetar.t0), float(magnetar.q)
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

    _opt = lambda src, attr: None if src is None else _as_float_array_or_none(getattr(src, attr, None))
    fwd_gamma_e, fwd_dnde = (None, None) if electron is None else (electron.gam_e, electron.d_n_gam_e)
    fwd_dnde_bh = _as_float_array_or_none(hadronic.d_n_gam_e_bh) if hadronic is not None else _opt(electron, "d_n_gam_e_bh")
    _echi = ["d_n_gam_e_chi", "chi_grid", "l_syn_spec_chi", "seed_syn_chi", "tau_syn_chi",
             "chi_radius_cm", "chi_gamma_bulk", "chi_dvolume_weight", "b_chi_g"]
    (fwd_dnde_chi, fwd_chi_grid, fwd_lsyn_chi, fwd_seed_chi, fwd_tau_chi,
     fwd_chi_radius, fwd_chi_gamma, fwd_chi_weight, fwd_b_chi) = (_opt(electron, a) for a in _echi)
    fwd_gamma_p, fwd_dndp = (None, None) if hadronic is None else (hadronic.gam_p, hadronic.d_n_gam_p)
    _hdist = ["gam_secondary", "d_n_gam_n", "d_n_gam_pi_plus", "d_n_gam_pi_minus",
              "d_n_gam_mu_minus_left", "d_n_gam_mu_minus_right", "d_n_gam_mu_plus_left", "d_n_gam_mu_plus_right"]
    (fwd_gamma_secondary, fwd_dndn, fwd_dndpi_plus, fwd_dndpi_minus,
     fwd_dndmu_ml, fwd_dndmu_mr, fwd_dndmu_pl, fwd_dndmu_pr) = (_opt(hadronic, a) for a in _hdist)
    _hcomp = ["l_had_pion_synch", "l_had_muon_synch", "l_had_pion_inverse_compton", "l_had_muon_inverse_compton"]
    fwd_had_gamma = None
    if hadronic is not None:
        fwd_had_gamma = hadronic.l_had_syn_spec + hadronic.l_had_pg_gamma
        for a in _hcomp:
            v = getattr(hadronic, a, None)
            if v is not None:
                fwd_had_gamma = fwd_had_gamma + v
    fwd_had_pi_syn, fwd_had_mu_syn, fwd_had_pi_ic, fwd_had_mu_ic = (_opt(hadronic, a) for a in _hcomp)
    fwd_nu_freq, fwd_nu_lum = (None, None)
    if hadronic is not None and state.config.hadronic.include_neutrino:
        fwd_nu_freq, fwd_nu_lum = hadronic.neutrino_frequency_hz, hadronic.neutrino_luminosity
    fwd_seed_freq = None if hadronic is None else np.asarray(state.photon_field.seed_frequency_hz, dtype=float)
    _hsed = ["l_had_syn_spec", "l_had_pg_gamma", "l_had_bethe_heitler", "l_had_hadronic_inverse_compton",
             "am3_process_power", "tau_pg", "tau_bh", "pg_photon_survival"]
    (fwd_had_syn, fwd_had_pg_gamma, fwd_had_bh, fwd_had_hic,
     fwd_am3_power, fwd_tau_pg, fwd_tau_bh, fwd_pg_survival) = (_opt(hadronic, a) for a in _hsed)
    fwd_timings = dict(hadronic.timings) if (hadronic is not None and hadronic.timings) else None
    rev_gamma_e, rev_dnde = (None, None)
    if reverse_shock is not None:
        rev_gamma_e, rev_dnde = reverse_shock.gam_e, reverse_shock.d_n_gam_e
    secondary_rs = None if reverse_emission is None else reverse_emission.secondary_rs
    return TrackBundle(
        fwd=CharTrack(
            t_obs=components.fwd.characteristic_time_s,
            radius=components.fwd.radius_cm,
            Gamma=components.fwd.gamma,
            N_p=np.asarray(components.fwd.swept_mass_g, dtype=float) / constants.para_m_p,
            Doppler=components.fwd.doppler,
            B_comv=components.fwd.magnetic_field_g,
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
            b_chi_g=fwd_b_chi,
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
            gamma_e=rev_gamma_e,
            dN_dgamma_e=rev_dnde,
            **_secondary_rs_detail_kwargs(secondary_rs),
        ),
        patches=patches,
    )


def _jet_magnetar_active(jet: JetProfile, theta_center: float) -> bool:
    if jet.kind in ("tophat", "gaussian", "powerlaw", "steppowerlaw"):
        return theta_center <= jet.theta_c
    if jet.kind == "twocomponent":
        return theta_center <= jet.theta_n
    return True
