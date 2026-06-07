from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, NamedTuple, Optional

import numpy as np

from asgard_core.asgard_config import default_num_threads
from asgard_core.asgard_numpy import trapezoid
from asgard_core.asgard_observables import build_multiband_observer_frequencies, combine_multiband_flux

class Scale(str, Enum):
    LINEAR = "linear"
    LOG = "log"
    LOG10 = "log10"
    FIXED = "fixed"


@dataclass
class Medium:
    """Callable density profile. Use factory functions (make_ism_medium, make_wind_medium)."""
    rho: Callable[[float, float, float], float]
    kind: str = "custom"
    label: str = "custom"
    n_ism: float = 1.0
    A_star: float = 0.0
    n0: float | None = None
    k: float = 2.0
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


def make_ism_medium(n_ism: Optional[float] = None, n0: Optional[float] = None) -> Medium:
    _n_ism = 1.0 if n_ism is None and n0 is None else float(n_ism if n_ism is not None else n0)
    medium = Medium(
        rho=lambda phi, theta, rc: float(np.full(np.asarray(rc).shape, medium.n_ism)),
        kind="ism",
        label="ism",
        n_ism=_n_ism,
        kernel_params={"d_ne": _n_ism, "a_star": -1.0},
    )
    return medium


def make_wind_medium(
    A_star: Optional[float] = None,
    n_ism: Optional[float] = None,
    n0: Optional[float] = None,
    k: float = 2.0,
    Astar: Optional[float] = None,
) -> Medium:
    _A_star = float(1.0 if A_star is None and Astar is None else A_star if A_star is not None else Astar)
    _n_ism = 0.1 if n_ism is None else float(n_ism)
    _n0 = None if n0 is None else float(n0)
    _k = float(k)
    medium = Medium(
        rho=lambda phi, theta, rc: _wind_rho(medium, rc),
        kind="wind",
        label="wind",
        n_ism=_n_ism,
        A_star=_A_star,
        n0=_n0,
        k=_k,
        kernel_params=_wind_kernel_params(_A_star, _n_ism, _n0, _k),
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


# Backward-compatible aliases
ISM = make_ism_medium
Wind = make_wind_medium


class Magnetar(NamedTuple):
    L0: float
    t0: float
    q: float


@dataclass
class JetProfile:
    """Unified jet profile. Use factory functions (make_tophat_jet, etc.) to construct."""
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
    # Ejecta params
    e_iso_fn: Callable[[float, float], float] | None = None
    gamma0_fn: Callable[[float, float], float] | None = None
    # Common
    spreading: bool = False
    duration: float | None = None
    magnetar: Magnetar | None = None

    def is_active(self, phi: float, theta: float) -> bool:
        return self.energy_iso(phi, theta) > 0.0 and self.gamma0(phi, theta) > 1.0


def make_tophat_jet(
    *args,
    E_iso: Optional[float] = None,
    lf: Optional[float] = None,
    theta_j: Optional[float] = None,
    Gamma0: Optional[float] = None,
    theta_c: Optional[float] = None,
    duration: Optional[float] = None,
    magnetar: Optional[Magnetar] = None,
    spreading: bool = False,
) -> JetProfile:
    if len(args) == 3 and E_iso is None and lf is None and theta_j is None and Gamma0 is None and theta_c is None:
        if float(args[0]) < 10.0 and float(args[1]) > 1.0e20:
            theta_c, E_iso, Gamma0 = args
        else:
            E_iso, lf, theta_j = args
    elif len(args) != 0:
        raise TypeError("make_tophat_jet accepts either keyword arguments or three positional arguments.")
    _E_iso = float(E_iso)
    _lf = float(lf if lf is not None else Gamma0)
    _theta_j = float(theta_j if theta_j is not None else theta_c)
    jet = JetProfile(
        kind="tophat",
        theta_max=_theta_j,
        E_iso=_E_iso,
        lf=_lf,
        theta_j=_theta_j,
        duration=None if duration is None else float(duration),
        magnetar=magnetar,
        spreading=bool(spreading),
    )
    jet.energy_iso = lambda phi, theta: jet.E_iso if theta < jet.theta_j else 0.0
    jet.gamma0 = lambda phi, theta: jet.lf if theta < jet.theta_j else 1.0
    return jet


def make_gaussian_jet(
    E_iso: float,
    lf: Optional[float] = None,
    theta_c: float = 0.1,
    theta_max: float = 0.6,
    Gamma0: Optional[float] = None,
    duration: Optional[float] = None,
    magnetar: Optional[Magnetar] = None,
    spreading: bool = False,
) -> JetProfile:
    jet = JetProfile(
        kind="gaussian",
        theta_max=theta_max,
        E_iso=float(E_iso),
        lf=float(lf if lf is not None else Gamma0),
        theta_c=float(theta_c),
        duration=None if duration is None else float(duration),
        magnetar=magnetar,
        spreading=bool(spreading),
    )
    jet.energy_iso = lambda phi, theta: jet.E_iso * np.exp(-0.5 * (theta / jet.theta_c) ** 2)
    jet.gamma0 = lambda phi, theta: 1.0 + (jet.lf - 1.0) * np.exp(-0.5 * (theta / jet.theta_c) ** 2)
    return jet


def make_powerlaw_jet(
    E_iso: float,
    lf: Optional[float] = None,
    theta_c: float = 0.1,
    k: Optional[float] = None,
    theta_max: float = np.pi / 2.0,
    Gamma0: Optional[float] = None,
    k_e: Optional[float] = None,
    k_g: Optional[float] = None,
    duration: Optional[float] = None,
    magnetar: Optional[Magnetar] = None,
    spreading: bool = False,
) -> JetProfile:
    _k_e = float(2.0 if k is None and k_e is None else k if k_e is None else k_e)
    _k_g = float(_k_e if k_g is None else k_g)
    jet = JetProfile(
        kind="powerlaw",
        theta_max=theta_max,
        E_iso=float(E_iso),
        lf=float(lf if lf is not None else Gamma0),
        theta_c=float(theta_c),
        k_e=_k_e,
        k_g=_k_g,
        duration=None if duration is None else float(duration),
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


def make_twocomponent_jet(
    E_iso_n: Optional[float] = None,
    lf_n: Optional[float] = None,
    theta_n: Optional[float] = None,
    E_iso_w: Optional[float] = None,
    lf_w: Optional[float] = None,
    theta_w: Optional[float] = None,
    E_iso_c: Optional[float] = None,
    Gamma0_c: Optional[float] = None,
    E_iso_outer: Optional[float] = None,
    Gamma0_outer: Optional[float] = None,
    theta_c: Optional[float] = None,
    theta_o: Optional[float] = None,
    duration: Optional[float] = None,
    magnetar: Optional[Magnetar] = None,
    spreading: bool = False,
) -> JetProfile:
    jet = JetProfile(
        kind="twocomponent",
        theta_max=float(theta_w if theta_w is not None else theta_o),
        E_iso_n=float(E_iso_n if E_iso_n is not None else E_iso_c),
        lf_n=float(lf_n if lf_n is not None else Gamma0_c),
        theta_n=float(theta_n if theta_n is not None else theta_c),
        E_iso_w=float(E_iso_w if E_iso_w is not None else E_iso_outer),
        lf_w=float(lf_w if lf_w is not None else Gamma0_outer),
        theta_w=float(theta_w if theta_w is not None else theta_o),
        duration=None if duration is None else float(duration),
        magnetar=magnetar,
        spreading=bool(spreading),
    )
    jet.energy_iso = lambda phi, theta: (
        jet.E_iso_n if theta <= jet.theta_n
        else jet.E_iso_w if theta <= jet.theta_w
        else 0.0
    )
    jet.gamma0 = lambda phi, theta: (
        jet.lf_n if theta <= jet.theta_n
        else jet.lf_w if theta <= jet.theta_w
        else 1.0
    )
    return jet


def make_steppowerlaw_jet(
    E_iso_c: float,
    lf_c: Optional[float] = None,
    theta_c: float = 0.1,
    E_iso_w: float = 1.0e51,
    lf_w: Optional[float] = None,
    theta_w: float = 0.3,
    k: float = 2.0,
    Gamma0_c: Optional[float] = None,
    Gamma0_w: Optional[float] = None,
    k_e: Optional[float] = None,
    k_g: Optional[float] = None,
    duration: Optional[float] = None,
    magnetar: Optional[Magnetar] = None,
    spreading: bool = False,
) -> JetProfile:
    _k_e = float(k if k_e is None else k_e)
    _k_g = float(_k_e if k_g is None else k_g)
    jet = JetProfile(
        kind="steppowerlaw",
        theta_max=float(theta_w),
        E_iso_c=float(E_iso_c),
        lf_c=float(lf_c if lf_c is not None else Gamma0_c),
        theta_c=float(theta_c),
        E_iso_w=float(E_iso_w),
        lf_w=float(lf_w if lf_w is not None else Gamma0_w),
        theta_w=float(theta_w),
        k_e=_k_e,
        k_g=_k_g,
        duration=None if duration is None else float(duration),
        magnetar=magnetar,
        spreading=bool(spreading),
    )
    jet.energy_iso = lambda phi, theta: (
        jet.E_iso_c if theta <= jet.theta_c
        else jet.E_iso_w * (theta / jet.theta_c) ** (-jet.k_e) if theta <= jet.theta_w
        else 0.0
    )
    jet.gamma0 = lambda phi, theta: (
        jet.lf_c if theta <= jet.theta_c
        else 1.0 + (jet.lf_w - 1.0) * (theta / jet.theta_c) ** (-jet.k_g) if theta <= jet.theta_w
        else 1.0
    )
    return jet


def make_ejecta_jet(
    e_iso_fn: Optional[Callable[[float, float], float]] = None,
    gamma0_fn: Optional[Callable[[float, float], float]] = None,
    theta_max: float = np.pi / 2.0,
    E_iso: Optional[Callable[[float, float], float]] = None,
    Gamma0: Optional[Callable[[float, float], float]] = None,
    duration: Optional[float] = None,
    magnetar: Optional[Magnetar] = None,
    spreading: bool = False,
) -> JetProfile:
    _e_iso_fn = e_iso_fn if e_iso_fn is not None else E_iso
    _gamma0_fn = gamma0_fn if gamma0_fn is not None else Gamma0
    jet = JetProfile(
        kind="ejecta",
        theta_max=theta_max,
        e_iso_fn=_e_iso_fn,
        gamma0_fn=_gamma0_fn,
        duration=None if duration is None else float(duration),
        magnetar=magnetar,
        spreading=bool(spreading),
    )
    jet.energy_iso = lambda phi, theta: float(jet.e_iso_fn(phi, theta))
    jet.gamma0 = lambda phi, theta: float(jet.gamma0_fn(phi, theta))
    return jet


# Backward-compatible aliases
TophatJet = make_tophat_jet
GaussianJet = make_gaussian_jet
PowerLawJet = make_powerlaw_jet
TwoComponentJet = make_twocomponent_jet
StepPowerLawJet = make_steppowerlaw_jet
Ejecta = make_ejecta_jet
Jet = JetProfile


class Observer:
    def __init__(
        self,
        *args,
        z: Optional[float] = None,
        theta_obs: float = 0.0,
        phi_obs: float = 0.0,
        lumi_dist_cm: Optional[float] = None,
        lumi_dist: Optional[float] = None,
    ) -> None:
        if len(args) == 3 and z is None and lumi_dist_cm is None and lumi_dist is None:
            lumi_dist, z, theta_obs = args
        elif len(args) == 4 and z is None and lumi_dist_cm is None and lumi_dist is None:
            lumi_dist, z, theta_obs, phi_obs = args
        elif len(args) != 0:
            raise TypeError("Observer accepts either keyword arguments or positional (lumi_dist, z, theta_obs[, phi_obs]).")
        self.z = float(z)
        self.theta_obs = float(theta_obs)
        self.phi_obs = float(phi_obs)
        self.lumi_dist_cm = None if lumi_dist_cm is None and lumi_dist is None else float(
            lumi_dist_cm if lumi_dist_cm is not None else lumi_dist
        )


@dataclass
class Radiation:
    eps_e: float
    eps_B: float
    p: float
    epsilon_p: float = 0.0
    epsilon_b_floor: float | None = None
    magnetic_decay_alpha_t: float = 0.0
    magnetic_decay_t0_s: float = 1.0
    xi_N: float = 0.1
    thermal_electrons: bool = False
    ssc: bool = False
    kn: bool = False
    proton_synch: bool = True
    pg: bool = False
    bethe_heitler: bool = False
    hadronic_inverse_compton: bool = False
    pp: bool = False
    neutrino: bool = False
    eta_acc: float = 1.0
    reverse_epsilon_p: float = 0.0
    pgamma_scheme: str = "disabled"
    pair_production: bool = False


@dataclass
class Setups:
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
    reverse_delta_t_s: float = 10.0
    reverse_sigma: float = 0.0
    include_cross_zone_ic: bool = False
    weno5: bool = False
    electron_solver: str = "fullhide_1d"
    cooling_kernel: str = "legacy"
    radiation_kernel: str = "legacy"
    dynamics_kernel: str = "forward_legacy"
    geometry_kernel: str = "sed_legacy"
    structured_backend: str = "fortran_1d"
    structured_parallel_mode: str = "outer"
    structured_outer_threads: Optional[int] = None
    structured_inner_threads: Optional[int] = None
    num_chi: Optional[int] = None
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
    pg_photon_survival: Optional[np.ndarray] = None
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



class Model:
    def __init__(
        self,
        medium: Optional[Medium] = None,
        jet: Optional[Jet] = None,
        observer: Optional[Observer] = None,
        fwd_rad: Optional[Radiation] = None,
        rvs_rad: Optional[Radiation] = None,
        setups: Optional[Setups] = None,
        resolutions: Optional[tuple[float, float, int]] = None,
        *args,
    ) -> None:
        if len(args) == 4 and isinstance(args[0], JetProfile) and isinstance(args[1], Medium) and isinstance(args[2], Observer) and isinstance(args[3], Radiation):
            jet, medium, observer, fwd_rad = args
        elif len(args) == 4 and isinstance(args[0], Medium) and isinstance(args[1], JetProfile) and isinstance(args[2], Observer) and isinstance(args[3], Radiation):
            medium, jet, observer, fwd_rad = args
        elif len(args) != 0:
            raise TypeError("Model accepts either keyword arguments or positional (jet, medium, observer, radiation).")
        elif isinstance(medium, JetProfile) and isinstance(jet, Medium) and isinstance(observer, Observer) and isinstance(fwd_rad, Radiation):
            medium, jet = jet, medium
        self.medium = medium
        self.jet = jet
        self.observer = observer
        self.fwd_rad = fwd_rad
        self.rvs_rad = rvs_rad
        self.setups = deepcopy(Setups() if setups is None else setups)
        if resolutions is not None:
            self._apply_resolutions(resolutions)
        self._last_details: Optional[TrackBundle] = None
        self._raw_cache: dict[tuple[str, str], tuple[FluxResult, TrackBundle]] = {}
        self._details_cache: dict[str, TrackBundle] = {}

    def flux_density_grid(self, times_s: np.ndarray, nu_hz: np.ndarray) -> FluxResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        return self._compute_raw(times_s, nu_hz)

    def flux_density(self, times_s: np.ndarray, nu_hz: np.ndarray) -> FluxResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        from .api_observe import _extract_pair_flux

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
            result = self.flux_density_grid(times_s, unique_freqs)
            pair_index = np.arange(times_s.shape[0], dtype=int)
            return pack(result, lambda values: values[inverse, pair_index])

        result = self._compute(times_s, nu_hz)
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

    def spectrum(self, time_s: float, nu_hz: np.ndarray) -> np.ndarray:
        return self.flux_density_grid(np.array([time_s], dtype=float), nu_hz).total[:, 0]

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

    def flux(self, time_s: np.ndarray | float, nu_min_hz: float, nu_max_hz: float, num_points: int = 64):
        times_s = np.atleast_1d(np.asarray(time_s, dtype=float))
        nu_hz = np.logspace(np.log10(nu_min_hz), np.log10(nu_max_hz), num_points)
        result = self.flux_density_grid(times_s, nu_hz)
        total = trapezoid(result.total, nu_hz, axis=0)
        fwd_sync = trapezoid(result.fwd.sync, nu_hz, axis=0)
        fwd_ssc = trapezoid(result.fwd.ssc, nu_hz, axis=0)
        rev_sync = trapezoid(result.rev.sync, nu_hz, axis=0)
        rev_ssc = trapezoid(result.rev.ssc, nu_hz, axis=0)
        cross_ic = None if result.cross_ic is None else trapezoid(result.cross_ic, nu_hz, axis=0)
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

    def component_fluxes(self, times_s: np.ndarray, nu_hz: np.ndarray) -> FluxResult:
        return self._compute(times_s, nu_hz)

    def flux_density_bands(self, times_s: np.ndarray) -> np.ndarray:
        frequencies_hz = build_multiband_observer_frequencies()[1]
        total_matrix = self.flux_density_grid(times_s, frequencies_hz).total
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

    def _compute(self, times_s: np.ndarray, nu_hz: np.ndarray) -> FluxResult:
        return self._compute_raw(times_s, nu_hz)

    def _compute_raw(
        self,
        times_s: np.ndarray,
        nu_hz: np.ndarray,
        *,
        solve_reference_times_s: np.ndarray | None = None,
    ) -> FluxResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        from .api_adaptive import _array_signature, _remember_cache_entry
        from .api_observe import _solve_direct_model, _solve_patch_model

        reference_signature = None
        if solve_reference_times_s is not None:
            reference_signature = _array_signature(np.asarray(solve_reference_times_s, dtype=float))
        cache_key = (_array_signature(times_s), _array_signature(nu_hz), reference_signature)
        cached = self._raw_cache.get(cache_key)
        if cached is not None:
            self._last_details = cached[1]
            return cached[0]
        if self.jet.kind == "tophat" and self._supports_direct_kernel():
            model_result = _solve_direct_model(
                self,
                times_s,
                nu_hz,
                solve_reference_times_s=solve_reference_times_s,
            )
        else:
            model_result = _solve_patch_model(
                self,
                times_s,
                nu_hz,
                solve_reference_times_s=solve_reference_times_s,
            )
        self._last_details = model_result[1]
        _remember_cache_entry(self._raw_cache, cache_key, model_result)
        return model_result[0]

    def _compute_details_only(self, times_s: np.ndarray) -> TrackBundle:
        times_s = np.asarray(times_s, dtype=float)
        from .api_adaptive import _array_signature, _remember_cache_entry
        from .api_observe import _evaluate_direct_details, _patch_details

        cache_key = _array_signature(times_s)
        cached = self._details_cache.get(cache_key)
        if cached is not None:
            self._last_details = cached
            return cached
        if self.jet.kind == "tophat" and self._supports_direct_kernel():
            details = _evaluate_direct_details(self, times_s)
        else:
            details = _patch_details(self, times_s)
        self._last_details = details
        _remember_cache_entry(self._details_cache, cache_key, details)
        return details

    def _supports_direct_kernel(self) -> bool:
        return self.medium.kind in ("ism", "wind")

    def _apply_resolutions(self, resolutions: tuple[float, float, int]) -> None:
        theta_ppd, phi_ppd, t_ppd = resolutions
        theta_max_deg = np.degrees(self.jet.theta_max)
        self.setups.patch_theta = max(1, int(np.ceil(theta_max_deg * float(theta_ppd))))
        self.setups.patch_phi = max(1, int(np.ceil(360.0 * float(phi_ppd))))
        self.setups.num_theta = self.setups.patch_theta
        self.setups.num_phi = self.setups.patch_phi
        decades = np.log10(self.setups.observer_time_max_s / self.setups.observer_time_min_s)
        self.setups.num_tobs = max(8, int(np.ceil(decades * float(t_ppd))))
        self._last_details = None
        self._raw_cache.clear()
        self._details_cache.clear()

    def _detail_time_count(self, t_min: float, t_max: float) -> int:
        if t_min <= 0.0 or t_max <= 0.0 or t_max <= t_min:
            return int(self.setups.num_tobs)
        decades = np.log10(float(t_max) / float(t_min))
        return max(int(self.setups.num_tobs), int(np.ceil(8.0 * decades)))
