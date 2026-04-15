from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, field
from enum import Enum
import hashlib
import json
from typing import Any, Callable, Optional

import numpy as np

from asgard_component_backend import (
    ComponentSpectra,
    build_solve_time_grid,
    build_query_setup,
    observe_spectra_from_setup,
    solve_component_spectra_from_setup,
    ModelState,
    observe_flux_grid_from_state,
    solve_model_state_from_setup,
    state_covers_request,
)
from asgard_models import FitConfig, ReverseShockConfig, SpectrumOutputConfig, default_num_threads
from asgard_observables import OUTPUT_BANDS, build_multiband_observer_frequencies, combine_multiband_flux
from asgard_postprocess import compute_light_curve_redchi
from asgard_presets import build_baseline_config
from src import constants


class Scale(str, Enum):
    LINEAR = "linear"
    linear = "linear"
    LOG = "log"
    log = "log"
    LOG10 = "log10"
    log10 = "log10"
    FIXED = "fixed"
    fixed = "fixed"


@dataclass
class Medium:
    rho: Callable[[float, float, float], float]
    label: str = "custom"

    def density(self, phi: float, theta: float, radius_cm: float) -> float:
        values = np.asarray(self.rho(phi, theta, radius_cm), dtype=float)
        if values.ndim == 0:
            return float(values)
        return values

    def __call__(self, phi: float, theta: float, radius_cm: float) -> float:
        return self.density(phi, theta, radius_cm)

    def to_kernel_params(self) -> dict[str, float]:
        raise NotImplementedError("User-defined Medium is not supported by the current ASGARD kernel.")


@dataclass
class ISM(Medium):
    n_ism: float = 1.0

    def __init__(self, n_ism: Optional[float] = None, n0: Optional[float] = None) -> None:
        density = 1.0 if n_ism is None and n0 is None else float(n_ism if n_ism is not None else n0)
        self.n_ism = density
        super().__init__(rho=self._rho, label="ism")

    def _rho(self, _phi: float, _theta: float, radius_cm: float):
        radius = np.asarray(radius_cm, dtype=float)
        values = np.full(radius.shape, self.n_ism, dtype=float)
        if values.ndim == 0:
            return float(values)
        return values

    def to_kernel_params(self) -> dict[str, float]:
        return {"d_ne": self.n_ism, "a_star": -1.0}


@dataclass
class Wind(Medium):
    A_star: float = 1.0
    n_ism: float = 0.1
    n0: Optional[float] = None
    k: float = 2.0

    def __init__(
        self,
        A_star: Optional[float] = None,
        n_ism: Optional[float] = None,
        n0: Optional[float] = None,
        k: float = 2.0,
        Astar: Optional[float] = None,
    ) -> None:
        self.A_star = float(A_star if A_star is not None else Astar)
        self.n_ism = 0.1 if n_ism is None else float(n_ism)
        self.n0 = None if n0 is None else float(n0)
        self.k = float(k)
        super().__init__(rho=self._rho, label="wind")

    def _rho(self, _phi: float, _theta: float, radius_cm: float) -> float:
        if self.k != 2.0:
            raise NotImplementedError("The current ASGARD wind kernel only supports k=2.")
        radius = np.asarray(radius_cm, dtype=float)
        d_ne_wind = self.A_star * 3.0e35 / radius**2
        values = d_ne_wind
        if self.n0 is not None and np.isfinite(self.n0) and self.n0 > 0.0:
            values = np.minimum(values, float(self.n0))
        values = np.maximum(self.n_ism, values)
        if values.ndim == 0:
            return float(values)
        return values

    def to_kernel_params(self) -> dict[str, float]:
        if self.k != 2.0:
            raise NotImplementedError("The current ASGARD wind kernel only supports k=2.")
        r0 = 0.0
        if self.n0 is not None and np.isfinite(self.n0) and self.n0 > 0.0:
            r0 = float(np.sqrt(self.A_star * 3.0e35 / float(self.n0)))
        return {"d_ne": self.n_ism, "a_star": self.A_star, "r0": r0}


@dataclass
class Magnetar:
    L0: float
    t0: float
    q: float


class Jet:
    def __init__(self, theta_max: float) -> None:
        self.theta_max = float(theta_max)

    def energy_iso(self, phi: float, theta: float) -> float:
        raise NotImplementedError

    def gamma0(self, phi: float, theta: float) -> float:
        raise NotImplementedError

    def is_active(self, phi: float, theta: float) -> bool:
        return self.energy_iso(phi, theta) > 0.0 and self.gamma0(phi, theta) > 1.0


class TophatJet(Jet):
    def __init__(
        self,
        *args,
        E_iso: Optional[float] = None,
        lf: Optional[float] = None,
        theta_j: Optional[float] = None,
        Gamma0: Optional[float] = None,
        theta_c: Optional[float] = None,
        duration: Optional[float] = None,
        magnetar: Optional[Magnetar] = None,
        spreading: bool = False,
    ) -> None:
        if len(args) == 3 and E_iso is None and lf is None and theta_j is None and Gamma0 is None and theta_c is None:
            if float(args[0]) < 10.0 and float(args[1]) > 1.0e20:
                theta_c, E_iso, Gamma0 = args
            else:
                E_iso, lf, theta_j = args
        elif len(args) != 0:
            raise TypeError("TophatJet accepts either keyword arguments or three positional arguments.")
        self.E_iso = float(E_iso)
        self.lf = float(lf if lf is not None else Gamma0)
        self.theta_j = float(theta_j if theta_j is not None else theta_c)
        self.duration = None if duration is None else float(duration)
        self.magnetar = magnetar
        self.spreading = bool(spreading)
        super().__init__(theta_max=self.theta_j)

    def energy_iso(self, phi: float, theta: float) -> float:
        return self.E_iso if theta <= self.theta_j else 0.0

    def gamma0(self, phi: float, theta: float) -> float:
        return self.lf if theta <= self.theta_j else 1.0


class GaussianJet(Jet):
    def __init__(
        self,
        E_iso: float,
        lf: Optional[float] = None,
        theta_c: float = 0.1,
        theta_max: float = 0.6,
        Gamma0: Optional[float] = None,
        duration: Optional[float] = None,
        magnetar: Optional[Magnetar] = None,
        spreading: bool = False,
    ) -> None:
        self.E_iso = float(E_iso)
        self.lf = float(lf if lf is not None else Gamma0)
        self.theta_c = float(theta_c)
        self.duration = None if duration is None else float(duration)
        self.magnetar = magnetar
        self.spreading = bool(spreading)
        super().__init__(theta_max=theta_max)

    def energy_iso(self, phi: float, theta: float) -> float:
        return self.E_iso * np.exp(-0.5 * (theta / self.theta_c) ** 2)

    def gamma0(self, phi: float, theta: float) -> float:
        return 1.0 + (self.lf - 1.0) * np.exp(-0.5 * (theta / self.theta_c) ** 2)


class PowerLawJet(Jet):
    def __init__(
        self,
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
    ) -> None:
        self.E_iso = float(E_iso)
        self.lf = float(lf if lf is not None else Gamma0)
        self.theta_c = float(theta_c)
        self.k_e = float(k if k_e is None else k_e)
        self.k_g = float(self.k_e if k_g is None else k_g)
        self.duration = None if duration is None else float(duration)
        self.magnetar = magnetar
        self.spreading = bool(spreading)
        super().__init__(theta_max=theta_max)

    def energy_iso(self, phi: float, theta: float) -> float:
        if theta <= self.theta_c:
            return self.E_iso
        return self.E_iso * (theta / self.theta_c) ** (-self.k_e)

    def gamma0(self, phi: float, theta: float) -> float:
        if theta <= self.theta_c:
            return self.lf
        return 1.0 + (self.lf - 1.0) * (theta / self.theta_c) ** (-self.k_g)


class TwoComponentJet(Jet):
    def __init__(
        self,
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
    ) -> None:
        self.E_iso_n = float(E_iso_n if E_iso_n is not None else E_iso_c)
        self.lf_n = float(lf_n if lf_n is not None else Gamma0_c)
        self.theta_n = float(theta_n if theta_n is not None else theta_c)
        self.E_iso_w = float(E_iso_w if E_iso_w is not None else E_iso_outer)
        self.lf_w = float(lf_w if lf_w is not None else Gamma0_outer)
        self.theta_w = float(theta_w if theta_w is not None else theta_o)
        self.duration = None if duration is None else float(duration)
        self.magnetar = magnetar
        self.spreading = bool(spreading)
        super().__init__(theta_max=self.theta_w)

    def energy_iso(self, phi: float, theta: float) -> float:
        if theta <= self.theta_n:
            return self.E_iso_n
        if theta <= self.theta_w:
            return self.E_iso_w
        return 0.0

    def gamma0(self, phi: float, theta: float) -> float:
        if theta <= self.theta_n:
            return self.lf_n
        if theta <= self.theta_w:
            return self.lf_w
        return 1.0


class StepPowerLawJet(Jet):
    def __init__(
        self,
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
    ) -> None:
        self.E_iso_c = float(E_iso_c)
        self.lf_c = float(lf_c if lf_c is not None else Gamma0_c)
        self.theta_c = float(theta_c)
        self.E_iso_w = float(E_iso_w)
        self.lf_w = float(lf_w if lf_w is not None else Gamma0_w)
        self.theta_w = float(theta_w)
        self.k_e = float(k if k_e is None else k_e)
        self.k_g = float(self.k_e if k_g is None else k_g)
        self.duration = None if duration is None else float(duration)
        self.magnetar = magnetar
        self.spreading = bool(spreading)
        super().__init__(theta_max=self.theta_w)

    def energy_iso(self, phi: float, theta: float) -> float:
        if theta <= self.theta_c:
            return self.E_iso_c
        if theta <= self.theta_w:
            return self.E_iso_w * (theta / self.theta_c) ** (-self.k_e)
        return 0.0

    def gamma0(self, phi: float, theta: float) -> float:
        if theta <= self.theta_c:
            return self.lf_c
        if theta <= self.theta_w:
            return 1.0 + (self.lf_w - 1.0) * (theta / self.theta_c) ** (-self.k_g)
        return 1.0


class Ejecta(Jet):
    def __init__(
        self,
        e_iso_fn: Optional[Callable[[float, float], float]] = None,
        gamma0_fn: Optional[Callable[[float, float], float]] = None,
        theta_max: float = np.pi / 2.0,
        E_iso: Optional[Callable[[float, float], float]] = None,
        Gamma0: Optional[Callable[[float, float], float]] = None,
        duration: Optional[float] = None,
        magnetar: Optional[Magnetar] = None,
        spreading: bool = False,
    ) -> None:
        self.e_iso_fn = e_iso_fn if e_iso_fn is not None else E_iso
        self.gamma0_fn = gamma0_fn if gamma0_fn is not None else Gamma0
        self.duration = None if duration is None else float(duration)
        self.magnetar = magnetar
        self.spreading = bool(spreading)
        super().__init__(theta_max=theta_max)

    def energy_iso(self, phi: float, theta: float) -> float:
        return float(self.e_iso_fn(phi, theta))

    def gamma0(self, phi: float, theta: float) -> float:
        return float(self.gamma0_fn(phi, theta))


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
    xi_N: float = 0.1
    ssc: bool = False
    kn: bool = False


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
    include_cross_zone_ic: bool = False
    weno5: bool = False
    electron_solver: str = "fullhide"
    electron_adaptive_substeps: bool = False
    electron_substep_rtol: float = 2.0e-2
    electron_substep_min: int = 25
    electron_substep_max: int = 150
    index_dyn: int = 3
    index_syn_integr: int = 2
    clean: bool = True


@dataclass
class BranchView:
    sync: np.ndarray
    ssc: np.ndarray

    @property
    def synch(self) -> np.ndarray:
        return self.sync


@dataclass
class DetailView:
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


@dataclass
class ModelDetails:
    fwd: DetailView
    rev: Optional[DetailView]
    patches: list[dict[str, float]]

    @property
    def rvs(self) -> Optional[DetailView]:
        return self.rev


@dataclass
class ModelFluxResult:
    total: np.ndarray
    fwd: BranchView
    rev: BranchView
    cross_ic: Optional[np.ndarray]

    @property
    def rvs(self) -> BranchView:
        return self.rev

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

    @property
    def shape(self):
        return self.image.shape


@dataclass
class ObsData:
    flux_density: list[dict[str, Any]] = field(default_factory=list)
    spectrum: list[dict[str, Any]] = field(default_factory=list)
    flux: list[dict[str, Any]] = field(default_factory=list)

    def add_flux_density(
        self,
        times_s: Optional[np.ndarray] = None,
        frequencies_hz: Optional[np.ndarray] = None,
        flux: Optional[np.ndarray] = None,
        flux_err: Optional[np.ndarray] = None,
        *,
        t: Optional[np.ndarray] = None,
        nu: Optional[np.ndarray] = None,
        f_nu: Optional[np.ndarray] = None,
        err: Optional[np.ndarray] = None,
    ) -> None:
        self.flux_density.append(
            {
                "times_s": np.asarray(times_s if times_s is not None else t, dtype=float),
                "frequencies_hz": np.asarray(frequencies_hz if frequencies_hz is not None else nu, dtype=float),
                "flux": np.asarray(flux if flux is not None else f_nu, dtype=float),
                "flux_err": np.asarray(flux_err if flux_err is not None else err, dtype=float),
            }
        )

    def add_spectrum(
        self,
        time_s: Optional[float] = None,
        frequencies_hz: Optional[np.ndarray] = None,
        flux: Optional[np.ndarray] = None,
        flux_err: Optional[np.ndarray] = None,
        *,
        t: Optional[float] = None,
        nu: Optional[np.ndarray] = None,
        f_nu: Optional[np.ndarray] = None,
        err: Optional[np.ndarray] = None,
    ) -> None:
        self.spectrum.append(
            {
                "time_s": float(time_s if time_s is not None else t),
                "frequencies_hz": np.asarray(frequencies_hz if frequencies_hz is not None else nu, dtype=float),
                "flux": np.asarray(flux if flux is not None else f_nu, dtype=float),
                "flux_err": np.asarray(flux_err if flux_err is not None else err, dtype=float),
            }
        )

    def add_flux(
        self,
        nu_min_hz: Optional[float] = None,
        nu_max_hz: Optional[float] = None,
        time_s: Optional[float] = None,
        flux: Optional[float] = None,
        flux_err: Optional[float] = None,
        *,
        nu_min: Optional[float] = None,
        nu_max: Optional[float] = None,
        t: Optional[float] = None,
        value: Optional[float] = None,
        err: Optional[float] = None,
        num_points: int = 64,
    ) -> None:
        self.flux.append(
            {
                "nu_min_hz": float(nu_min_hz if nu_min_hz is not None else nu_min),
                "nu_max_hz": float(nu_max_hz if nu_max_hz is not None else nu_max),
                "time_s": float(time_s if time_s is not None else t),
                "flux": float(flux if flux is not None else value),
                "flux_err": float(flux_err if flux_err is not None else err),
                "num_points": int(num_points),
            }
        )

    @staticmethod
    def logscale_screen(times_s: np.ndarray, data_density: float = 5.0) -> np.ndarray:
        times_s = np.asarray(times_s, dtype=float)
        if times_s.ndim != 1:
            raise ValueError("times_s must be one-dimensional.")
        if times_s.size == 0:
            return np.array([], dtype=int)
        log_times = np.log10(times_s)
        decades = max(log_times.max() - log_times.min(), 0.0)
        num_bins = max(1, int(np.ceil(decades * data_density)))
        edges = np.linspace(log_times.min(), log_times.max(), num_bins + 1)
        indices: list[int] = []
        for i in range(num_bins):
            left = edges[i]
            right = edges[i + 1]
            mask = (log_times >= left) & (log_times < right if i < num_bins - 1 else log_times <= right)
            if np.any(mask):
                candidates = np.where(mask)[0]
                center = 0.5 * (left + right)
                idx = int(candidates[np.argmin(np.abs(log_times[candidates] - center))])
                indices.append(idx)
        return np.asarray(sorted(set(indices)), dtype=int)

class ParamDef:
    def __init__(self, name: str, *args, scale: Scale = Scale.LINEAR) -> None:
        self.name = name
        self.path: Optional[str]
        if len(args) >= 3 and isinstance(args[0], str):
            self.path = args[0]
            self.lower = float(args[1])
            self.upper = float(args[2])
            if len(args) >= 4:
                scale = args[3]
        elif len(args) >= 2:
            self.path = None
            self.lower = float(args[0])
            self.upper = float(args[1])
            if len(args) >= 3:
                scale = args[2]
        else:
            raise TypeError("ParamDef requires either (name, path, lower, upper, scale) or (name, lower, upper, scale).")
        self.scale = Scale(scale)

    def transform(self, value: float) -> float:
        if self.scale in (Scale.LOG, Scale.LOG10):
            return 10.0**value
        if self.scale == Scale.FIXED:
            return float(self.lower)
        return value

    def is_fixed(self) -> bool:
        return self.scale == Scale.FIXED or np.isclose(self.lower, self.upper)


@dataclass
class FitResult:
    best_params: dict[str, float]
    best_loglike: float
    samples: Optional[np.ndarray] = None
    log_prob: Optional[np.ndarray] = None
    logz: Optional[float] = None
    logzerr: Optional[float] = None
    labels: Optional[list[str]] = None
    top_k_params: Optional[list[dict[str, float]]] = None


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
        if len(args) == 4 and isinstance(args[0], Jet) and isinstance(args[1], Medium) and isinstance(args[2], Observer) and isinstance(args[3], Radiation):
            jet, medium, observer, fwd_rad = args
        elif len(args) == 4 and isinstance(args[0], Medium) and isinstance(args[1], Jet) and isinstance(args[2], Observer) and isinstance(args[3], Radiation):
            medium, jet, observer, fwd_rad = args
        elif len(args) != 0:
            raise TypeError("Model accepts either keyword arguments or positional (jet, medium, observer, radiation).")
        elif isinstance(medium, Jet) and isinstance(jet, Medium) and isinstance(observer, Observer) and isinstance(fwd_rad, Radiation):
            medium, jet = jet, medium
        self.medium = medium
        self.jet = jet
        self.observer = observer
        self.fwd_rad = fwd_rad
        self.rvs_rad = rvs_rad
        self.setups = deepcopy(Setups() if setups is None else setups)
        if resolutions is not None:
            self._apply_resolutions(resolutions)
        self._last_details: Optional[ModelDetails] = None
        self._details_cache: dict[tuple[float, float, int, str], ModelDetails] = {}
        self._state_cache: dict[str, list[ModelState]] = {}
        self._observed_cache: dict[tuple[str, str, str], ModelFluxResult] = {}

    def flux_density_grid(self, times_s: np.ndarray, nu_hz: np.ndarray) -> ModelFluxResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        if _should_use_adaptive_spectrum_grid(times_s, nu_hz):
            return _adaptive_spectrum_grid(self, times_s, nu_hz)
        if _should_use_adaptive_flux_grid(times_s, nu_hz):
            return _adaptive_flux_density_grid(self, times_s, nu_hz)
        return self._compute_raw(times_s, nu_hz)

    def flux_density(self, times_s: np.ndarray, nu_hz: np.ndarray) -> ModelFluxResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        if times_s.ndim == 1 and nu_hz.ndim == 1 and times_s.shape == nu_hz.shape:
            unique_freqs, inverse = np.unique(nu_hz, return_inverse=True)
            result = self.flux_density_grid(times_s, unique_freqs)
            pair_index = np.arange(times_s.shape[0], dtype=int)
            pair_total = result.total[inverse, pair_index]
            pair_fwd_sync = result.fwd.sync[inverse, pair_index]
            pair_fwd_ssc = result.fwd.ssc[inverse, pair_index]
            pair_rev_sync = result.rev.sync[inverse, pair_index]
            pair_rev_ssc = result.rev.ssc[inverse, pair_index]
            pair_cross_ic = None if result.cross_ic is None else result.cross_ic[inverse, pair_index]
            return ModelFluxResult(
                total=pair_total,
                fwd=BranchView(sync=pair_fwd_sync, ssc=pair_fwd_ssc),
                rev=BranchView(sync=pair_rev_sync, ssc=pair_rev_ssc),
                cross_ic=pair_cross_ic,
            )

        result = self._compute(times_s, nu_hz)
        pair_total = _extract_pair_flux(result.total, times_s, nu_hz)
        pair_fwd_sync = _extract_pair_flux(result.fwd.sync, times_s, nu_hz)
        pair_fwd_ssc = _extract_pair_flux(result.fwd.ssc, times_s, nu_hz)
        pair_rev_sync = _extract_pair_flux(result.rev.sync, times_s, nu_hz)
        pair_rev_ssc = _extract_pair_flux(result.rev.ssc, times_s, nu_hz)
        pair_cross_ic = None if result.cross_ic is None else _extract_pair_flux(result.cross_ic, times_s, nu_hz)
        return ModelFluxResult(
            total=pair_total,
            fwd=BranchView(sync=pair_fwd_sync, ssc=pair_fwd_ssc),
            rev=BranchView(sync=pair_rev_sync, ssc=pair_rev_ssc),
            cross_ic=pair_cross_ic,
        )

    def flux_density_exposures(
        self,
        times_s: np.ndarray,
        nu_hz: np.ndarray,
        exposures_s: np.ndarray,
        num_subsamples: int = 8,
    ) -> ModelFluxResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        exposures_s = np.asarray(exposures_s, dtype=float)
        if times_s.shape != nu_hz.shape or times_s.shape != exposures_s.shape:
            raise ValueError("times_s, nu_hz, and exposures_s must have the same shape.")
        if int(num_subsamples) <= 0:
            raise ValueError("num_subsamples must be positive.")
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
        return _render_sky_image(self, times_s, float(nu_obs), float(fov), int(npixel))

    def flux(self, time_s: np.ndarray | float, nu_min_hz: float, nu_max_hz: float, num_points: int = 64):
        times_s = np.atleast_1d(np.asarray(time_s, dtype=float))
        nu_hz = np.logspace(np.log10(nu_min_hz), np.log10(nu_max_hz), num_points)
        result = self.flux_density_grid(times_s, nu_hz)
        total = np.trapezoid(result.total, nu_hz, axis=0)
        fwd_sync = np.trapezoid(result.fwd.sync, nu_hz, axis=0)
        fwd_ssc = np.trapezoid(result.fwd.ssc, nu_hz, axis=0)
        rev_sync = np.trapezoid(result.rev.sync, nu_hz, axis=0)
        rev_ssc = np.trapezoid(result.rev.ssc, nu_hz, axis=0)
        cross_ic = None if result.cross_ic is None else np.trapezoid(result.cross_ic, nu_hz, axis=0)
        band_result = ModelFluxResult(
            total=total,
            fwd=BranchView(sync=fwd_sync, ssc=fwd_ssc),
            rev=BranchView(sync=rev_sync, ssc=rev_ssc),
            cross_ic=cross_ic,
        )
        if np.asarray(time_s).ndim == 0:
            return float(band_result.total[0])
        return band_result

    def details(self, t_min: Optional[float] = None, t_max: Optional[float] = None) -> ModelDetails:
        if t_min is not None or t_max is not None:
            t1 = self.setups.observer_time_min_s if t_min is None else float(t_min)
            t2 = self.setups.observer_time_max_s if t_max is None else float(t_max)
            key = (t1, t2, self.setups.num_tobs, _model_signature(self))
            if key not in self._details_cache:
                times = np.logspace(np.log10(t1), np.log10(t2), self.setups.num_tobs)
                self._details_cache[key] = self._compute_details_only(times)
            self._last_details = self._details_cache[key]
        elif self._last_details is None:
            self._last_details = self._compute_details_only(self.default_times())
        return self._last_details

    def jet_E_iso(self, phi: float, theta: float) -> float:
        return self.jet.energy_iso(phi, theta)

    def jet_Gamma0(self, phi: float, theta: float) -> float:
        return self.jet.gamma0(phi, theta)

    def medium_rho(self, phi: float, theta: float, r_cm: float) -> float:
        return self.medium.density(phi, theta, r_cm)

    def component_fluxes(self, times_s: np.ndarray, nu_hz: np.ndarray) -> ModelFluxResult:
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

    def default_frequencies(self) -> np.ndarray:
        return np.logspace(np.log10(1.0e9), np.log10(1.0e24), 64)

    def _compute(self, times_s: np.ndarray, nu_hz: np.ndarray) -> ModelFluxResult:
        return self._compute_raw(times_s, nu_hz)

    def _compute_raw(self, times_s: np.ndarray, nu_hz: np.ndarray) -> ModelFluxResult:
        times_s = np.asarray(times_s, dtype=float)
        nu_hz = np.asarray(nu_hz, dtype=float)
        if isinstance(self.jet, TophatJet) and self._supports_direct_kernel():
            model_result = _evaluate_direct_model(self, times_s, nu_hz)
        else:
            model_result = _evaluate_patch_model(self, times_s, nu_hz)
        self._last_details = model_result[1]
        return model_result[0]

    def _compute_details_only(self, times_s: np.ndarray) -> ModelDetails:
        times_s = np.asarray(times_s, dtype=float)
        if isinstance(self.jet, TophatJet) and self._supports_direct_kernel():
            details = _evaluate_direct_details(self, times_s)
        else:
            details = _evaluate_patch_details(self, times_s)
        self._last_details = details
        return details

    def _supports_direct_kernel(self) -> bool:
        return isinstance(self.medium, (ISM, Wind))

    def _apply_resolutions(self, resolutions: tuple[float, float, int]) -> None:
        theta_ppd, phi_ppd, t_ppd = resolutions
        theta_max_deg = np.degrees(self.jet.theta_max)
        self.setups.patch_theta = max(1, int(np.ceil(theta_max_deg * float(theta_ppd))))
        self.setups.patch_phi = max(1, int(np.ceil(360.0 * float(phi_ppd))))
        self.setups.num_theta = self.setups.patch_theta
        self.setups.num_phi = self.setups.patch_phi
        decades = np.log10(self.setups.observer_time_max_s / self.setups.observer_time_min_s)
        self.setups.num_tobs = max(8, int(np.ceil(decades * float(t_ppd))))


class Fitter:
    def __init__(
        self,
        model: Optional[Model] = None,
        data: Optional[ObsData] = None,
        params: Optional[list[ParamDef]] = None,
        cfg: Optional[Any] = None,
        num_workers: int = 1,
        **config_kwargs,
    ) -> None:
        if isinstance(model, ObsData) and data is None:
            data = model
            model = None
        if cfg is None and config_kwargs:
            cfg = config_kwargs
        self.model = model if model is not None else _coerce_model(cfg)
        self.data = ObsData() if data is None else data
        self.params = [] if params is None else list(params)
        self.num_workers = int(num_workers)

    def build_model(self, values: dict[str, float]) -> Model:
        model = deepcopy(self.model)
        for param in self.params:
            if param.name in values or param.is_fixed():
                raw_value = param.lower if param.is_fixed() else values[param.name]
                path = _resolve_param_path(model, param)
                _set_dotted_attr(model, path, param.transform(raw_value))
        return model

    def loglike(self, values: dict[str, float]) -> float:
        model = self.build_model(values)
        loglike = 0.0
        for dataset in self.data.flux_density:
            pred = _evaluate_flux_observations(model, dataset["times_s"], dataset["frequencies_hz"])
            resid = (pred - dataset["flux"]) / dataset["flux_err"]
            loglike -= 0.5 * float(np.sum(resid * resid))
        for dataset in self.data.spectrum:
            pred = model.spectrum(dataset["time_s"], dataset["frequencies_hz"])
            resid = (pred - dataset["flux"]) / dataset["flux_err"]
            loglike -= 0.5 * float(np.sum(resid * resid))
        for dataset in self.data.flux:
            pred = model.flux(
                dataset["time_s"],
                dataset["nu_min_hz"],
                dataset["nu_max_hz"],
                num_points=dataset["num_points"],
            )
            resid = (pred - dataset["flux"]) / dataset["flux_err"]
            loglike -= 0.5 * float(resid * resid)
        return loglike

    def flux_density_grid(self, values: dict[str, float], times_s: np.ndarray, nu_hz: np.ndarray) -> ModelFluxResult:
        return self.build_model(values).flux_density_grid(times_s, nu_hz)

    def add_flux_density(self, *args, **kwargs) -> None:
        self.data.add_flux_density(*args, **kwargs)

    def add_spectrum(self, *args, **kwargs) -> None:
        self.data.add_spectrum(*args, **kwargs)

    def add_flux(self, *args, **kwargs) -> None:
        self.data.add_flux(*args, **kwargs)

    def run_emcee(self, initial: np.ndarray, nwalkers: int, nsteps: int) -> FitResult:
        import emcee

        ndim = len(self.params)

        def _lnprob(theta):
            trial = {param.name: theta[i] for i, param in enumerate(self.params)}
            return self.loglike(trial)

        sampler = emcee.EnsembleSampler(nwalkers, ndim, _lnprob)
        sampler.run_mcmc(initial, nsteps, progress=False)
        flat = sampler.get_chain(flat=True)
        log_prob = sampler.get_log_prob(flat=True)
        best_idx = int(np.argmax(log_prob))
        best = {param.name: flat[best_idx, i] for i, param in enumerate(self.params)}
        return FitResult(
            best_params=best,
            best_loglike=float(log_prob[best_idx]),
            samples=flat,
            log_prob=log_prob,
            labels=[param.name for param in self.params],
            top_k_params=[best],
        )

    def run_multinest(self, outputfiles_basename: str, verbose: bool = True) -> FitResult:
        from pymultinest.solve import solve

        n_dims = len(self.params)

        def _prior_transform(cube):
            params = np.asarray(cube, dtype=float).copy()
            for i, param in enumerate(self.params):
                params[i] = param.lower + params[i] * (param.upper - param.lower)
            return params

        def _loglike(theta):
            trial = {param.name: theta[i] for i, param in enumerate(self.params)}
            return self.loglike(trial)

        result = solve(
            LogLikelihood=_loglike,
            Prior=_prior_transform,
            n_dims=n_dims,
            outputfiles_basename=outputfiles_basename,
            verbose=verbose,
        )
        samples = np.asarray(result["samples"], dtype=float)
        best_idx = int(np.argmax(result["log_likelihood"]))
        best = {param.name: float(samples[best_idx, i]) for i, param in enumerate(self.params)}
        return FitResult(
            best_params=best,
            best_loglike=float(result["log_likelihood"][best_idx]),
            samples=samples,
            log_prob=np.asarray(result["log_likelihood"], dtype=float),
            logz=float(result["logZ"]),
            logzerr=float(result["logZerr"]),
            labels=[param.name for param in self.params],
            top_k_params=[best],
        )

    def fit(
        self,
        param_defs: Optional[list[ParamDef]] = None,
        *,
        sampler: str = "emcee",
        total_steps: int = 128,
        burn_frac: float = 0.5,
        thin: int = 1,
        nwalkers: Optional[int] = None,
        outputfiles_basename: str = "chains/vegas-",
        verbose: bool = False,
        nsteps: Optional[int] = None,
        nburn: Optional[int] = None,
        npool: Optional[int] = None,
        resolution: Optional[tuple[float, float, int]] = None,
    ) -> FitResult:
        if param_defs is not None:
            self.params = list(param_defs)
        if not self.params:
            raise ValueError("No parameter definitions were provided to Fitter.fit().")
        if nsteps is not None:
            total_steps = int(nsteps)
        if nburn is not None:
            burn_frac = float(nburn) / float(total_steps)
        if npool is not None:
            self.num_workers = int(npool)
        if resolution is not None:
            self.model._apply_resolutions(resolution)

        active_params = [param for param in self.params if not param.is_fixed()]
        if sampler.lower() in ("multinest", "pymultinest"):
            if active_params != self.params:
                raise NotImplementedError("run_multinest currently requires all fitted parameters to be active.")
            return self.run_multinest(outputfiles_basename=outputfiles_basename, verbose=verbose)

        if sampler.lower() not in ("emcee", "mcmc"):
            raise ValueError(f"Unsupported sampler: {sampler}")
        if not active_params:
            best = {param.name: param.lower for param in self.params}
            best_loglike = self.loglike(best)
            return FitResult(best_params=best, best_loglike=best_loglike, labels=[param.name for param in self.params], top_k_params=[best])

        import emcee

        ndim = len(active_params)
        walker_num = max(2 * ndim + 2, 8) if nwalkers is None else int(nwalkers)
        initial = np.zeros((walker_num, ndim), dtype=float)
        for i, param in enumerate(active_params):
            span = param.upper - param.lower
            if span <= 0.0:
                initial[:, i] = param.lower
            else:
                initial[:, i] = param.lower + span * np.random.default_rng(1234 + i).uniform(0.2, 0.8, size=walker_num)

        def _lnprob(theta):
            trial = {param.name: theta[i] for i, param in enumerate(active_params)}
            for param in active_params:
                if param.name not in trial:
                    continue
                if trial[param.name] < param.lower or trial[param.name] > param.upper:
                    return -np.inf
            return self.loglike(trial)

        sampler_obj = emcee.EnsembleSampler(walker_num, ndim, _lnprob)
        sampler_obj.run_mcmc(initial, int(total_steps), progress=False)
        chain = sampler_obj.get_chain(flat=True)
        log_prob = sampler_obj.get_log_prob(flat=True)
        burn = int(burn_frac * chain.shape[0])
        flat = chain[burn::thin]
        flat_log_prob = log_prob[burn::thin]
        best_idx = int(np.argmax(flat_log_prob))
        best = {param.name: float(flat[best_idx, i]) for i, param in enumerate(active_params)}
        for param in self.params:
            if param.is_fixed():
                best[param.name] = param.lower
        return FitResult(
            best_params=best,
            best_loglike=float(flat_log_prob[best_idx]),
            samples=flat,
            log_prob=flat_log_prob,
            labels=[param.name for param in active_params],
            top_k_params=[best],
        )


def _array_signature(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values, dtype=float))
    digest = hashlib.blake2b(digest_size=16)
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
    digest.update(array.view(np.uint8))
    return digest.hexdigest()


def _fit_config_signature(config: FitConfig) -> str:
    payload = asdict(config)
    payload.pop("num_tobs", None)
    payload.pop("t_obs_min_log10", None)
    payload.pop("t_obs_max_log10", None)
    digest = hashlib.blake2b(digest_size=16)
    digest.update(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8"))
    return digest.hexdigest()


def _model_signature(model: Model) -> str:
    payload = {
        "medium": {
            "type": type(model.medium).__name__,
            "label": getattr(model.medium, "label", None),
            "data": {k: v for k, v in vars(model.medium).items() if k != "rho"},
        },
        "jet": {"type": type(model.jet).__name__, "data": vars(model.jet)},
        "observer": vars(model.observer),
        "fwd_rad": vars(model.fwd_rad),
        "rvs_rad": None if model.rvs_rad is None else vars(model.rvs_rad),
        "setups": asdict(model.setups),
    }
    digest = hashlib.blake2b(digest_size=16)
    digest.update(json.dumps(payload, sort_keys=True, default=float, separators=(",", ":")).encode("utf-8"))
    return digest.hexdigest()


def _find_cached_state(
    model: Model,
    state_signature: str,
    times_s: np.ndarray,
    frequencies_hz: np.ndarray | None,
) -> Optional[ModelState]:
    for state in model._state_cache.get(state_signature, []):
        if state_covers_request(state, times_s, frequencies_hz):
            return state
    return None


def _remember_state(model: Model, state_signature: str, state: ModelState) -> ModelState:
    bucket = model._state_cache.setdefault(state_signature, [])
    bucket.append(state)
    bucket.sort(key=lambda item: float(np.max(item.setup.observer_time_s)))
    return state


def _observe_result_from_components(observed: dict[str, np.ndarray | None]) -> ModelFluxResult:
    total = np.asarray(observed["total"], dtype=float)
    rev_sync = np.zeros_like(total) if observed["rev_sync"] is None else np.asarray(observed["rev_sync"], dtype=float)
    rev_ssc = np.zeros_like(total) if observed["rev_ssc"] is None else np.asarray(observed["rev_ssc"], dtype=float)
    cross_ic = None if observed["cross_ic"] is None else np.asarray(observed["cross_ic"], dtype=float)
    return ModelFluxResult(
        total=total,
        fwd=BranchView(sync=np.asarray(observed["fwd_sync"], dtype=float), ssc=np.asarray(observed["fwd_ssc"], dtype=float)),
        rev=BranchView(sync=rev_sync, ssc=rev_ssc),
        cross_ic=cross_ic,
    )


def _observe_cached_state(
    model: Model,
    state: ModelState,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
) -> ModelFluxResult:
    state_signature = _fit_config_signature(state.config)
    query_key = (state_signature, _array_signature(times_s), _array_signature(nu_hz))
    cached = model._observed_cache.get(query_key)
    if cached is not None:
        return cached
    observed_state = observe_flux_grid_from_state(state, times_s, nu_hz)
    result = _observe_result_from_components(observed_state.components)
    model._observed_cache[query_key] = result
    return result


def _log_time_midpoint(t_left: float, t_right: float) -> float:
    if t_left <= 0.0 or t_right <= 0.0:
        return 0.5 * (t_left + t_right)
    return float(np.sqrt(t_left * t_right))


def _interpolate_segment_value(
    t_left: float,
    f_left: float,
    t_right: float,
    f_right: float,
    t_mid: float,
) -> float:
    if f_left > 0.0 and f_right > 0.0 and t_left > 0.0 and t_right > 0.0 and t_mid > 0.0:
        slope = (np.log(f_right) - np.log(f_left)) / (np.log(t_right) - np.log(t_left))
        return float(np.exp(np.log(f_left) + slope * (np.log(t_mid) - np.log(t_left))))
    weight = (t_mid - t_left) / (t_right - t_left)
    return float(f_left + weight * (f_right - f_left))


def _relative_segment_error(actual: float, estimate: float) -> float:
    scale = max(abs(actual), 1.0e-99)
    return abs(actual - estimate) / scale


def _should_use_adaptive_flux_grid(times_s: np.ndarray, nu_hz: np.ndarray) -> bool:
    times_s = np.asarray(times_s, dtype=float)
    nu_hz = np.asarray(nu_hz, dtype=float)
    if times_s.ndim != 1 or nu_hz.ndim != 1:
        return False
    if times_s.size < 48 or nu_hz.size == 0 or nu_hz.size > 8:
        return False
    if np.any(~np.isfinite(times_s)) or np.any(times_s <= 0.0):
        return False
    return float(np.max(times_s)) / float(np.min(times_s)) >= 2.0


def _should_use_adaptive_spectrum_grid(times_s: np.ndarray, nu_hz: np.ndarray) -> bool:
    times_s = np.asarray(times_s, dtype=float)
    nu_hz = np.asarray(nu_hz, dtype=float)
    if times_s.ndim != 1 or nu_hz.ndim != 1:
        return False
    if times_s.size != 1 or nu_hz.size < 64:
        return False
    if np.any(~np.isfinite(nu_hz)) or np.any(nu_hz <= 0.0):
        return False
    return float(np.max(nu_hz)) / float(np.min(nu_hz)) >= 10.0


def _adaptive_segment_error(
    t_left: float,
    t_mid: float,
    t_right: float,
    left_values: dict[str, Optional[np.ndarray]],
    mid_values: dict[str, Optional[np.ndarray]],
    right_values: dict[str, Optional[np.ndarray]],
) -> float:
    max_error = 0.0
    for key in ("total", "fwd_sync", "fwd_ssc", "rev_sync", "rev_ssc", "cross_ic"):
        left = left_values[key]
        mid = mid_values[key]
        right = right_values[key]
        if left is None or mid is None or right is None:
            continue
        estimate = np.empty_like(mid)
        for i_freq in range(mid.shape[0]):
            estimate[i_freq] = _interpolate_segment_value(
                t_left,
                float(left[i_freq]),
                t_right,
                float(right[i_freq]),
                t_mid,
            )
        scale = np.maximum(np.abs(mid), 1.0e-99)
        error = float(np.max(np.abs(mid - estimate) / scale))
        if error > max_error:
            max_error = error
    return max_error


def _adaptive_frequency_error(
    nu_left: float,
    nu_mid: float,
    nu_right: float,
    left_values: dict[str, Optional[float]],
    mid_values: dict[str, Optional[float]],
    right_values: dict[str, Optional[float]],
) -> float:
    max_error = 0.0
    for key in ("total", "fwd_sync", "fwd_ssc", "rev_sync", "rev_ssc", "cross_ic"):
        left = left_values[key]
        mid = mid_values[key]
        right = right_values[key]
        if left is None or mid is None or right is None:
            continue
        estimate = _interpolate_segment_value(nu_left, float(left), nu_right, float(right), nu_mid)
        error = _relative_segment_error(float(mid), estimate)
        if error > max_error:
            max_error = error
    return max_error


def _interpolate_time_series(
    source_times_s: np.ndarray,
    source_values: np.ndarray,
    target_times_s: np.ndarray,
) -> np.ndarray:
    source_times_s = np.asarray(source_times_s, dtype=float)
    target_times_s = np.asarray(target_times_s, dtype=float)
    result = np.empty((source_values.shape[0], target_times_s.size), dtype=float)
    for i_freq in range(source_values.shape[0]):
        y = np.asarray(source_values[i_freq], dtype=float)
        if np.all(y > 0.0):
            result[i_freq] = np.exp(
                np.interp(
                    np.log(target_times_s),
                    np.log(source_times_s),
                    np.log(y),
                )
            )
        else:
            result[i_freq] = np.interp(target_times_s, source_times_s, y)
    return result


def _interpolate_frequency_series(
    source_freqs_hz: np.ndarray,
    source_values: np.ndarray,
    target_freqs_hz: np.ndarray,
) -> np.ndarray:
    source_freqs_hz = np.asarray(source_freqs_hz, dtype=float)
    source_values = np.asarray(source_values, dtype=float)
    target_freqs_hz = np.asarray(target_freqs_hz, dtype=float)
    if np.all(source_values > 0.0):
        return np.exp(
            np.interp(
                np.log(target_freqs_hz),
                np.log(source_freqs_hz),
                np.log(source_values),
            )
        )
    return np.interp(target_freqs_hz, source_freqs_hz, source_values)


def _adaptive_flux_density_grid(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    *,
    tolerance: float = 1.0e-2,
    max_depth: int = 6,
    min_ratio: float = 1.02,
    coarse_nodes: int = 16,
) -> ModelFluxResult:
    requested_times = np.asarray(times_s, dtype=float)
    order = np.argsort(requested_times)
    sorted_times = requested_times[order]
    unique_times = np.unique(sorted_times)
    if unique_times.size <= coarse_nodes:
        return model._compute_raw(requested_times, nu_hz)

    node_times = np.logspace(
        np.log10(float(unique_times[0])),
        np.log10(float(unique_times[-1])),
        min(int(coarse_nodes), unique_times.size),
    )
    node_result = model._compute_raw(node_times, nu_hz)
    component_cache = {
        float(t): {
            "total": node_result.total[:, idx],
            "fwd_sync": node_result.fwd.sync[:, idx],
            "fwd_ssc": node_result.fwd.ssc[:, idx],
            "rev_sync": node_result.rev.sync[:, idx],
            "rev_ssc": node_result.rev.ssc[:, idx],
            "cross_ic": None if node_result.cross_ic is None else node_result.cross_ic[:, idx],
        }
        for idx, t in enumerate(node_times)
    }
    segments = [(float(node_times[i]), float(node_times[i + 1])) for i in range(node_times.size - 1)]

    for _ in range(max_depth):
        midpoint_times: list[float] = []
        midpoint_meta: list[tuple[float, float, float]] = []
        next_segments: list[tuple[float, float]] = []
        for t_left, t_right in segments:
            if t_right / t_left < min_ratio:
                next_segments.append((t_left, t_right))
                continue
            t_mid = _log_time_midpoint(t_left, t_right)
            if t_mid in component_cache:
                next_segments.append((t_left, t_right))
                continue
            midpoint_times.append(t_mid)
            midpoint_meta.append((t_left, t_mid, t_right))
        if not midpoint_times:
            break
        midpoint_result = model._compute_raw(np.array(midpoint_times, dtype=float), nu_hz)
        for idx, t_mid in enumerate(midpoint_times):
            component_cache[float(t_mid)] = {
                "total": midpoint_result.total[:, idx],
                "fwd_sync": midpoint_result.fwd.sync[:, idx],
                "fwd_ssc": midpoint_result.fwd.ssc[:, idx],
                "rev_sync": midpoint_result.rev.sync[:, idx],
                "rev_ssc": midpoint_result.rev.ssc[:, idx],
                "cross_ic": None if midpoint_result.cross_ic is None else midpoint_result.cross_ic[:, idx],
            }
        refined = False
        for t_left, t_mid, t_right in midpoint_meta:
            error = _adaptive_segment_error(
                t_left,
                t_mid,
                t_right,
                component_cache[t_left],
                component_cache[t_mid],
                component_cache[t_right],
            )
            if error > tolerance and t_right / t_left >= min_ratio:
                next_segments.append((t_left, t_mid))
                next_segments.append((t_mid, t_right))
                refined = True
            else:
                next_segments.append((t_left, t_right))
        segments = next_segments
        if not refined:
            break

    source_times = np.array(sorted(component_cache.keys()), dtype=float)
    total_grid = np.column_stack([component_cache[float(t)]["total"] for t in source_times])
    fwd_sync_grid = np.column_stack([component_cache[float(t)]["fwd_sync"] for t in source_times])
    fwd_ssc_grid = np.column_stack([component_cache[float(t)]["fwd_ssc"] for t in source_times])
    rev_sync_grid = np.column_stack([component_cache[float(t)]["rev_sync"] for t in source_times])
    rev_ssc_grid = np.column_stack([component_cache[float(t)]["rev_ssc"] for t in source_times])
    has_cross_ic = any(component_cache[float(t)]["cross_ic"] is not None for t in source_times)
    cross_ic_grid = None
    if has_cross_ic:
        cross_ic_grid = np.column_stack([
            np.zeros_like(component_cache[float(source_times[0])]["total"])
            if component_cache[float(t)]["cross_ic"] is None
            else component_cache[float(t)]["cross_ic"]
            for t in source_times
        ])

    total_sorted = _interpolate_time_series(source_times, total_grid, sorted_times)
    fwd_sync_sorted = _interpolate_time_series(source_times, fwd_sync_grid, sorted_times)
    fwd_ssc_sorted = _interpolate_time_series(source_times, fwd_ssc_grid, sorted_times)
    rev_sync_sorted = _interpolate_time_series(source_times, rev_sync_grid, sorted_times)
    rev_ssc_sorted = _interpolate_time_series(source_times, rev_ssc_grid, sorted_times)
    cross_ic_sorted = None if cross_ic_grid is None else _interpolate_time_series(source_times, cross_ic_grid, sorted_times)

    inverse_order = np.argsort(order)
    total = total_sorted[:, inverse_order]
    fwd_sync = fwd_sync_sorted[:, inverse_order]
    fwd_ssc = fwd_ssc_sorted[:, inverse_order]
    rev_sync = rev_sync_sorted[:, inverse_order]
    rev_ssc = rev_ssc_sorted[:, inverse_order]
    cross_ic = None if cross_ic_sorted is None else cross_ic_sorted[:, inverse_order]
    return ModelFluxResult(
        total=total,
        fwd=BranchView(sync=fwd_sync, ssc=fwd_ssc),
        rev=BranchView(sync=rev_sync, ssc=rev_ssc),
        cross_ic=cross_ic,
    )


def _adaptive_spectrum_grid(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    *,
    tolerance: float = 1.0e-2,
    max_depth: int = 6,
    min_ratio: float = 1.02,
    coarse_nodes: int = 16,
) -> ModelFluxResult:
    requested_freqs = np.asarray(nu_hz, dtype=float)
    order = np.argsort(requested_freqs)
    sorted_freqs = requested_freqs[order]
    unique_freqs = np.unique(sorted_freqs)
    if unique_freqs.size <= coarse_nodes:
        return model._compute_raw(times_s, nu_hz)

    node_freqs = np.logspace(
        np.log10(float(unique_freqs[0])),
        np.log10(float(unique_freqs[-1])),
        min(int(coarse_nodes), unique_freqs.size),
    )
    node_result = model._compute_raw(times_s, node_freqs)
    component_cache = {
        float(freq): {
            "total": float(node_result.total[idx, 0]),
            "fwd_sync": float(node_result.fwd.sync[idx, 0]),
            "fwd_ssc": float(node_result.fwd.ssc[idx, 0]),
            "rev_sync": float(node_result.rev.sync[idx, 0]),
            "rev_ssc": float(node_result.rev.ssc[idx, 0]),
            "cross_ic": None if node_result.cross_ic is None else float(node_result.cross_ic[idx, 0]),
        }
        for idx, freq in enumerate(node_freqs)
    }
    segments = [(float(node_freqs[i]), float(node_freqs[i + 1])) for i in range(node_freqs.size - 1)]

    for _ in range(max_depth):
        midpoint_freqs: list[float] = []
        midpoint_meta: list[tuple[float, float, float]] = []
        next_segments: list[tuple[float, float]] = []
        for nu_left, nu_right in segments:
            if nu_right / nu_left < min_ratio:
                next_segments.append((nu_left, nu_right))
                continue
            nu_mid = _log_time_midpoint(nu_left, nu_right)
            if nu_mid in component_cache:
                next_segments.append((nu_left, nu_right))
                continue
            midpoint_freqs.append(nu_mid)
            midpoint_meta.append((nu_left, nu_mid, nu_right))
        if not midpoint_freqs:
            break
        midpoint_result = model._compute_raw(times_s, np.array(midpoint_freqs, dtype=float))
        for idx, nu_mid in enumerate(midpoint_freqs):
            component_cache[float(nu_mid)] = {
                "total": float(midpoint_result.total[idx, 0]),
                "fwd_sync": float(midpoint_result.fwd.sync[idx, 0]),
                "fwd_ssc": float(midpoint_result.fwd.ssc[idx, 0]),
                "rev_sync": float(midpoint_result.rev.sync[idx, 0]),
                "rev_ssc": float(midpoint_result.rev.ssc[idx, 0]),
                "cross_ic": None if midpoint_result.cross_ic is None else float(midpoint_result.cross_ic[idx, 0]),
            }
        refined = False
        for nu_left, nu_mid, nu_right in midpoint_meta:
            error = _adaptive_frequency_error(
                nu_left,
                nu_mid,
                nu_right,
                component_cache[nu_left],
                component_cache[nu_mid],
                component_cache[nu_right],
            )
            if error > tolerance and nu_right / nu_left >= min_ratio:
                next_segments.append((nu_left, nu_mid))
                next_segments.append((nu_mid, nu_right))
                refined = True
            else:
                next_segments.append((nu_left, nu_right))
        segments = next_segments
        if not refined:
            break

    source_freqs = np.array(sorted(component_cache.keys()), dtype=float)
    total_source = np.array([component_cache[float(freq)]["total"] for freq in source_freqs], dtype=float)
    fwd_sync_source = np.array([component_cache[float(freq)]["fwd_sync"] for freq in source_freqs], dtype=float)
    fwd_ssc_source = np.array([component_cache[float(freq)]["fwd_ssc"] for freq in source_freqs], dtype=float)
    rev_sync_source = np.array([component_cache[float(freq)]["rev_sync"] for freq in source_freqs], dtype=float)
    rev_ssc_source = np.array([component_cache[float(freq)]["rev_ssc"] for freq in source_freqs], dtype=float)
    has_cross_ic = any(component_cache[float(freq)]["cross_ic"] is not None for freq in source_freqs)
    cross_ic_source = None
    if has_cross_ic:
        cross_ic_source = np.array(
            [
                0.0 if component_cache[float(freq)]["cross_ic"] is None else float(component_cache[float(freq)]["cross_ic"])
                for freq in source_freqs
            ],
            dtype=float,
        )

    total_sorted = _interpolate_frequency_series(source_freqs, total_source, sorted_freqs)
    fwd_sync_sorted = _interpolate_frequency_series(source_freqs, fwd_sync_source, sorted_freqs)
    fwd_ssc_sorted = _interpolate_frequency_series(source_freqs, fwd_ssc_source, sorted_freqs)
    rev_sync_sorted = _interpolate_frequency_series(source_freqs, rev_sync_source, sorted_freqs)
    rev_ssc_sorted = _interpolate_frequency_series(source_freqs, rev_ssc_source, sorted_freqs)
    cross_ic_sorted = None if cross_ic_source is None else _interpolate_frequency_series(source_freqs, cross_ic_source, sorted_freqs)

    inverse_order = np.argsort(order)
    total = total_sorted[inverse_order][:, None]
    fwd_sync = fwd_sync_sorted[inverse_order][:, None]
    fwd_ssc = fwd_ssc_sorted[inverse_order][:, None]
    rev_sync = rev_sync_sorted[inverse_order][:, None]
    rev_ssc = rev_ssc_sorted[inverse_order][:, None]
    cross_ic = None if cross_ic_sorted is None else cross_ic_sorted[inverse_order][:, None]
    return ModelFluxResult(
        total=total,
        fwd=BranchView(sync=fwd_sync, ssc=fwd_ssc),
        rev=BranchView(sync=rev_sync, ssc=rev_ssc),
        cross_ic=cross_ic,
    )


def _batch_fetch_pair_result(
    model: Model,
    cache: dict[tuple[float, float], tuple[float, float, float, float, float, Optional[float]]],
    query_pairs: list[tuple[float, float]],
) -> None:
    missing: list[tuple[float, float]] = []
    seen: set[tuple[float, float]] = set()
    for pair in query_pairs:
        if pair not in cache and pair not in seen:
            missing.append(pair)
            seen.add(pair)
    if not missing:
        return
    times_s = np.array([pair[0] for pair in missing], dtype=float)
    frequencies_hz = np.array([pair[1] for pair in missing], dtype=float)
    result = model.flux_density(times_s, frequencies_hz)
    for idx, pair in enumerate(missing):
        cross_ic = None if result.cross_ic is None else float(result.cross_ic[idx])
        cache[pair] = (
            float(result.total[idx]),
            float(result.fwd.sync[idx]),
            float(result.fwd.ssc[idx]),
            float(result.rev.sync[idx]),
            float(result.rev.ssc[idx]),
            cross_ic,
        )


def _adaptive_exposure_average(
    model: Model,
    times_s: np.ndarray,
    frequencies_hz: np.ndarray,
    exposures_s: np.ndarray,
    num_subsamples: int,
    *,
    tolerance: float = 1.0e-2,
    max_depth: int = 6,
    min_ratio: float = 1.02,
) -> ModelFluxResult:
    pair_cache: dict[tuple[float, float], tuple[float, float, float, float, float, Optional[float]]] = {}
    exposure_nodes: list[list[float]] = []
    segments_by_exposure: list[list[tuple[float, float]]] = []
    initial_pairs: list[tuple[float, float]] = []

    for time_s, freq_hz, exposure_s in zip(times_s, frequencies_hz, exposures_s):
        t_start = max(float(time_s) - 0.5 * float(exposure_s), 1.0e-30)
        t_stop = float(time_s) + 0.5 * float(exposure_s)
        if np.isclose(t_start, t_stop):
            nodes = [float(time_s)]
            segments = []
        else:
            nodes = np.linspace(t_start, t_stop, num_subsamples, dtype=float).tolist()
            segments = [(nodes[i], nodes[i + 1]) for i in range(len(nodes) - 1)]
        exposure_nodes.append(nodes)
        segments_by_exposure.append(segments)
        initial_pairs.extend((node, float(freq_hz)) for node in nodes)

    _batch_fetch_pair_result(model, pair_cache, initial_pairs)

    for _ in range(max_depth):
        midpoint_pairs: list[tuple[float, float]] = []
        midpoint_meta: list[tuple[int, float, float, float, float]] = []
        next_segments_by_exposure: list[list[tuple[float, float]]] = [[] for _ in exposure_nodes]

        for idx, segments in enumerate(segments_by_exposure):
            freq_hz = float(frequencies_hz[idx])
            for t_left, t_right in segments:
                if t_right / t_left < min_ratio:
                    next_segments_by_exposure[idx].append((t_left, t_right))
                    continue
                t_mid = _log_time_midpoint(t_left, t_right)
                midpoint_pairs.append((t_mid, freq_hz))
                midpoint_meta.append((idx, t_left, t_mid, t_right, freq_hz))

        if not midpoint_pairs:
            break

        _batch_fetch_pair_result(model, pair_cache, midpoint_pairs)
        refined = False
        for idx, t_left, t_mid, t_right, freq_hz in midpoint_meta:
            f_left = pair_cache[(t_left, freq_hz)][0]
            f_mid = pair_cache[(t_mid, freq_hz)][0]
            f_right = pair_cache[(t_right, freq_hz)][0]
            f_interp = _interpolate_segment_value(t_left, f_left, t_right, f_right, t_mid)
            err = _relative_segment_error(f_mid, f_interp)
            if err > tolerance and t_right / t_left >= min_ratio:
                exposure_nodes[idx].append(t_mid)
                next_segments_by_exposure[idx].append((t_left, t_mid))
                next_segments_by_exposure[idx].append((t_mid, t_right))
                refined = True
            else:
                next_segments_by_exposure[idx].append((t_left, t_right))
        segments_by_exposure = next_segments_by_exposure
        if not refined:
            break

    total = np.zeros(times_s.shape[0], dtype=float)
    fwd_sync = np.zeros_like(total)
    fwd_ssc = np.zeros_like(total)
    rev_sync = np.zeros_like(total)
    rev_ssc = np.zeros_like(total)
    cross_ic = np.zeros_like(total)
    has_cross_ic = False

    for idx, nodes in enumerate(exposure_nodes):
        freq_hz = float(frequencies_hz[idx])
        node_array = np.array(sorted(set(nodes)), dtype=float)
        duration = float(node_array[-1] - node_array[0]) if node_array.size > 1 else 0.0
        values = np.array([pair_cache[(float(node), freq_hz)] for node in node_array], dtype=object)
        total_values = np.array([entry[0] for entry in values], dtype=float)
        fwd_sync_values = np.array([entry[1] for entry in values], dtype=float)
        fwd_ssc_values = np.array([entry[2] for entry in values], dtype=float)
        rev_sync_values = np.array([entry[3] for entry in values], dtype=float)
        rev_ssc_values = np.array([entry[4] for entry in values], dtype=float)
        cross_values = np.array([0.0 if entry[5] is None else entry[5] for entry in values], dtype=float)
        has_cross_ic = has_cross_ic or np.any(cross_values != 0.0)

        if duration == 0.0:
            total[idx] = total_values[0]
            fwd_sync[idx] = fwd_sync_values[0]
            fwd_ssc[idx] = fwd_ssc_values[0]
            rev_sync[idx] = rev_sync_values[0]
            rev_ssc[idx] = rev_ssc_values[0]
            cross_ic[idx] = cross_values[0]
            continue

        total[idx] = float(np.trapezoid(total_values, node_array) / duration)
        fwd_sync[idx] = float(np.trapezoid(fwd_sync_values, node_array) / duration)
        fwd_ssc[idx] = float(np.trapezoid(fwd_ssc_values, node_array) / duration)
        rev_sync[idx] = float(np.trapezoid(rev_sync_values, node_array) / duration)
        rev_ssc[idx] = float(np.trapezoid(rev_ssc_values, node_array) / duration)
        cross_ic[idx] = float(np.trapezoid(cross_values, node_array) / duration)

    return ModelFluxResult(
        total=total,
        fwd=BranchView(sync=fwd_sync, ssc=fwd_ssc),
        rev=BranchView(sync=rev_sync, ssc=rev_ssc),
        cross_ic=None if not has_cross_ic else cross_ic,
    )


def _resolve_patch_state(
    model: Model,
    config: FitConfig,
    times_s: np.ndarray,
    requested_frequencies_hz: np.ndarray | None,
) -> ModelState:
    state_signature = _fit_config_signature(config)
    cached = _find_cached_state(model, state_signature, times_s, requested_frequencies_hz)
    if cached is not None:
        return cached
    solve_times_s = build_solve_time_grid(times_s, model.setups.num_tobs)
    if requested_frequencies_hz is not None:
        solve_t_min = min(float(model.setups.observer_time_min_s), float(np.min(times_s)))
        solve_count = max(int(model.setups.num_tobs), int(np.unique(times_s).size))
        solve_t_max_requested = float(np.max(times_s))
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
    setup = build_query_setup(config, solve_times_s, requested_frequencies_hz)
    state = solve_model_state_from_setup(
        config,
        setup,
        requested_frequencies_hz=requested_frequencies_hz,
    )
    return _remember_state(model, state_signature, state)


def _evaluate_direct_model(model: Model, times_s: np.ndarray, nu_hz: np.ndarray) -> tuple[ModelFluxResult, ModelDetails]:
    config = _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )
    state = _resolve_patch_state(model, config, times_s, nu_hz)
    observed = _observe_cached_state(model, state, times_s, nu_hz)
    details = _build_details(state.component_spectra, patches=[{"phi": 0.0, "theta": 0.0, "weight": 1.0}])
    return observed, details


def _evaluate_direct_details(model: Model, times_s: np.ndarray) -> ModelDetails:
    config = _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )
    state = _resolve_patch_state(model, config, times_s, None)
    return _build_details(state.component_spectra, patches=[{"phi": 0.0, "theta": 0.0, "weight": 1.0}])


def _evaluate_patch_model(model: Model, times_s: np.ndarray, nu_hz: np.ndarray) -> tuple[ModelFluxResult, ModelDetails]:
    total = np.zeros((nu_hz.shape[0], times_s.shape[0]), dtype=float)
    fwd_sync_total = np.zeros_like(total)
    fwd_ssc_total = np.zeros_like(total)
    rev_sync_total = np.zeros_like(total)
    rev_ssc_total = np.zeros_like(total)
    cross_ic_total = np.zeros_like(total)
    patches_meta: list[dict[str, float]] = []
    details_fwd: Optional[DetailView] = None
    details_rev: Optional[DetailView] = None

    for phi_center, theta_center, patch_half_angle in _iter_jet_patches(model):
        e_iso = model.jet.energy_iso(phi_center, theta_center)
        gamma0 = model.jet.gamma0(phi_center, theta_center)
        if e_iso <= 0.0 or gamma0 <= 1.0:
            continue
        theta_v = _angular_separation(theta_center, phi_center, model.observer.theta_obs, model.observer.phi_obs)
        config = _build_fit_config_for_patch(
            model,
            phi_center=phi_center,
            theta_v=theta_v,
            opening_angle_jet=patch_half_angle,
            e_iso=e_iso,
            gamma0=gamma0,
            theta_center=theta_center,
        )
        state = _resolve_patch_state(model, config, times_s, nu_hz)
        observed = _observe_cached_state(model, state, times_s, nu_hz)
        total += observed.total
        fwd_sync_total += observed.fwd.sync
        fwd_ssc_total += observed.fwd.ssc
        rev_sync_total += observed.rev.sync
        rev_ssc_total += observed.rev.ssc
        if observed.cross_ic is not None:
            cross_ic_total += observed.cross_ic
        patches_meta.append(
            {
                "phi": float(phi_center),
                "theta": float(theta_center),
                "theta_v": float(theta_v),
                "half_angle": float(patch_half_angle),
                "E_iso": float(e_iso),
                "Gamma0": float(gamma0),
            }
        )
        if details_fwd is None:
            details = _build_details(state.component_spectra, patches_meta)
            details_fwd = details.fwd
            details_rev = details.rev

    if details_fwd is None:
        raise ValueError("No active jet patches were found for the requested structured jet.")
    return (
        ModelFluxResult(
            total=total,
            fwd=BranchView(sync=fwd_sync_total, ssc=fwd_ssc_total),
            rev=BranchView(sync=rev_sync_total, ssc=rev_ssc_total),
            cross_ic=cross_ic_total,
        ),
        ModelDetails(fwd=details_fwd, rev=details_rev, patches=patches_meta),
    )


def _evaluate_patch_details(model: Model, times_s: np.ndarray) -> ModelDetails:
    patches_meta: list[dict[str, float]] = []
    first_component: Optional[ComponentSpectra] = None

    for phi_center, theta_center, patch_half_angle in _iter_jet_patches(model):
        e_iso = model.jet.energy_iso(phi_center, theta_center)
        gamma0 = model.jet.gamma0(phi_center, theta_center)
        if e_iso <= 0.0 or gamma0 <= 1.0:
            continue
        theta_v = _angular_separation(theta_center, phi_center, model.observer.theta_obs, model.observer.phi_obs)
        patches_meta.append(
            {
                "phi": float(phi_center),
                "theta": float(theta_center),
                "theta_v": float(theta_v),
                "half_angle": float(patch_half_angle),
                "E_iso": float(e_iso),
                "Gamma0": float(gamma0),
            }
        )
        if first_component is None:
            config = _build_fit_config_for_patch(
                model,
                phi_center=phi_center,
                theta_v=theta_v,
                opening_angle_jet=patch_half_angle,
                e_iso=e_iso,
                gamma0=gamma0,
                theta_center=theta_center,
            )
            first_component = _resolve_patch_state(model, config, times_s, None).component_spectra

    if first_component is None:
        raise ValueError("No active jet patches were found for the requested structured jet.")
    return _build_details(first_component, patches_meta)


def _render_sky_image(model: Model, times_s: np.ndarray, nu_obs: float, fov: float, npixel: int) -> SkyImage:
    pixel_size = fov / float(npixel)
    extent = np.array([-0.5 * fov, 0.5 * fov, -0.5 * fov, 0.5 * fov], dtype=float)
    pixel_axis = np.linspace(extent[0] + 0.5 * pixel_size, extent[1] - 0.5 * pixel_size, npixel)
    x_grid, y_grid = np.meshgrid(pixel_axis, pixel_axis, indexing="ij")
    image = np.zeros((times_s.shape[0], npixel, npixel), dtype=float)
    angular_diameter_distance_cm = _angular_diameter_distance_cm(model.observer)
    sightline, sky_x_axis, sky_y_axis = _sky_basis(model.observer)
    frequencies_hz = np.array([nu_obs], dtype=float)
    use_ring_cache = _sky_image_can_cache_azimuthally(model)
    ring_cache: dict[tuple[float, float, float, float], tuple[np.ndarray, np.ndarray]] = {}

    for phi_center, theta_center, patch_half_angle in _iter_sky_image_patches(model, npixel):
        e_iso = model.jet.energy_iso(phi_center, theta_center)
        gamma0 = model.jet.gamma0(phi_center, theta_center)
        if e_iso <= 0.0 or gamma0 <= 1.0:
            continue
        cache_key = (
            round(theta_center, 12),
            round(patch_half_angle, 12),
            round(e_iso, 6),
            round(gamma0, 6),
        )
        cached = ring_cache.get(cache_key) if use_ring_cache else None
        if cached is None:
            theta_v = _angular_separation(theta_center, phi_center, model.observer.theta_obs, model.observer.phi_obs)
            config = _build_fit_config_for_patch(
                model,
                phi_center=phi_center,
                theta_v=theta_v,
                opening_angle_jet=patch_half_angle,
                e_iso=e_iso,
                gamma0=gamma0,
                theta_center=theta_center,
            )
            state = _resolve_patch_state(model, config, times_s, frequencies_hz)
            observed = _observe_cached_state(model, state, times_s, frequencies_hz)
            patch_flux = np.asarray(observed.total[0, :], dtype=float)
            radius_cm = _interpolate_positive_series(
                state.component_spectra.fwd.characteristic_time_s,
                state.component_spectra.fwd.radius_cm,
                times_s,
            )
            if use_ring_cache:
                ring_cache[cache_key] = (patch_flux, radius_cm)
        else:
            patch_flux, radius_cm = cached
        if not np.any(np.isfinite(patch_flux) & (patch_flux > 0.0)):
            continue

        x_center, y_center = _project_patch_to_sky(
            radius_cm,
            theta_center,
            phi_center,
            sightline,
            sky_x_axis,
            sky_y_axis,
            angular_diameter_distance_cm,
        )
        sigma = np.maximum(
            radius_cm * np.sin(max(patch_half_angle, 1.0e-12)) / angular_diameter_distance_cm / 2.0,
            0.5 * pixel_size,
        )
        image += patch_flux[:, None, None] * _gaussian_splat_stack(x_grid, y_grid, x_center, y_center, sigma)

    return SkyImage(image=image, extent=extent, pixel_solid_angle=pixel_size * pixel_size)


def _build_fit_config_for_patch(
    model: Model,
    *,
    phi_center: float,
    theta_v: float,
    opening_angle_jet: float,
    e_iso: float,
    gamma0: float,
    theta_center: Optional[float] = None,
) -> FitConfig:
    if getattr(model.jet, "spreading", False):
        raise NotImplementedError("Jet spreading is not implemented in the current ASGARD backend.")
    reverse_delta_t = model.setups.reverse_delta_t_s
    if getattr(model.jet, "duration", None) is not None:
        reverse_delta_t = float(model.jet.duration)
    index_y = 0
    # `ssc` controls whether the SSC radiation component is emitted.
    # `ssc_cooling` controls whether IC cooling is included in the electron solver.
    if model.setups.ssc_cooling:
        index_y = 1 if model.fwd_rad.kn else 2
    config = build_baseline_config(
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
        p=model.fwd_rad.p,
        f_e=model.fwd_rad.xi_N,
        index_y=index_y,
        include_forward_ssc=model.fwd_rad.ssc,
        luminosity_distance_cm_override=model.observer.lumi_dist_cm,
        initial_radius_cm=model.setups.initial_radius_cm,
        spectrum_output=SpectrumOutputConfig(enabled=False),
    )
    magnetar = getattr(model.jet, "magnetar", None)
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
            include_ssc=model.rvs_rad.ssc,
            include_cross_zone_ic=model.setups.include_cross_zone_ic,
        )
    kernel_medium = _project_medium_to_kernel(
        model.medium,
        phi_center=phi_center,
        theta_center=0.0 if theta_center is None else theta_center,
    )
    for key, value in kernel_medium.items():
        setattr(config, key, value)
    return config


def _build_details(component_spectra: ComponentSpectra, patches: list[dict[str, float]]) -> ModelDetails:
    return ModelDetails(
        fwd=DetailView(
            t_obs=component_spectra.fwd.characteristic_time_s,
            radius=component_spectra.fwd.radius_cm,
            Gamma=component_spectra.fwd.gamma,
            N_p=component_spectra.fwd.swept_mass_g / constants.para_m_p,
            Doppler=component_spectra.fwd.doppler,
            B_comv=component_spectra.fwd.magnetic_field_g,
            nu_m=component_spectra.fwd.nu_m,
            nu_c=component_spectra.fwd.nu_c,
            nu_a=component_spectra.fwd.nu_a,
            nu_M=component_spectra.fwd.nu_M,
        ),
        rev=None
        if component_spectra.rev is None
        else DetailView(
            t_obs=component_spectra.rev.characteristic_time_s,
            radius=component_spectra.rev.radius_cm,
            Gamma=component_spectra.rev.gamma,
            N_p=component_spectra.rev.swept_mass_g / constants.para_m_p,
            Doppler=component_spectra.rev.doppler,
            B_comv=component_spectra.rev.magnetic_field_g,
            nu_m=component_spectra.rev.nu_m,
            nu_c=component_spectra.rev.nu_c,
            nu_a=component_spectra.rev.nu_a,
            nu_M=component_spectra.rev.nu_M,
        ),
        patches=patches,
    )


def _iter_jet_patches(model: Model):
    theta_edges = np.linspace(0.0, model.jet.theta_max, model.setups.patch_theta + 1)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, model.setups.patch_phi + 1)
    for i_theta in range(model.setups.patch_theta):
        theta1 = theta_edges[i_theta]
        theta2 = theta_edges[i_theta + 1]
        theta_center = 0.5 * (theta1 + theta2)
        for i_phi in range(model.setups.patch_phi):
            phi1 = phi_edges[i_phi]
            phi2 = phi_edges[i_phi + 1]
            phi_center = 0.5 * (phi1 + phi2)
            domega = (np.cos(theta1) - np.cos(theta2)) * (phi2 - phi1)
            patch_half_angle = np.sqrt(max(domega, 1.0e-12) / np.pi)
            yield phi_center, theta_center, patch_half_angle


def _iter_sky_image_patches(model: Model, npixel: int):
    theta_bins = min(model.setups.patch_theta, max(2, int(np.ceil(np.sqrt(npixel) / 6.0))))
    phi_bins = min(model.setups.patch_phi, max(12, 6 * theta_bins))
    theta_edges = np.linspace(0.0, model.jet.theta_max, theta_bins + 1)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, phi_bins + 1)
    for i_theta in range(theta_bins):
        theta1 = theta_edges[i_theta]
        theta2 = theta_edges[i_theta + 1]
        theta_center = 0.5 * (theta1 + theta2)
        for i_phi in range(phi_bins):
            phi1 = phi_edges[i_phi]
            phi2 = phi_edges[i_phi + 1]
            phi_center = 0.5 * (phi1 + phi2)
            domega = (np.cos(theta1) - np.cos(theta2)) * (phi2 - phi1)
            patch_half_angle = np.sqrt(max(domega, 1.0e-12) / np.pi)
            yield phi_center, theta_center, patch_half_angle


def _angular_separation(theta1: float, phi1: float, theta2: float, phi2: float) -> float:
    cos_alpha = (
        np.cos(theta1) * np.cos(theta2)
        + np.sin(theta1) * np.sin(theta2) * np.cos(phi1 - phi2)
    )
    return float(np.arccos(np.clip(cos_alpha, -1.0, 1.0)))


def _sky_image_can_cache_azimuthally(model: Model) -> bool:
    if abs(model.observer.theta_obs) > 1.0e-12:
        return False
    if not isinstance(model.medium, (ISM, Wind)):
        return False
    return _jet_is_axisymmetric(model.jet)


def _jet_is_axisymmetric(jet: Jet) -> bool:
    return isinstance(jet, (TophatJet, GaussianJet, PowerLawJet, TwoComponentJet, StepPowerLawJet))


def _sky_basis(observer: Observer) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sightline = _direction_vector(observer.theta_obs, observer.phi_obs)
    trial = np.array([0.0, 0.0, 1.0], dtype=float)
    sky_x = np.cross(trial, sightline)
    if np.linalg.norm(sky_x) < 1.0e-12:
        sky_x = np.array([1.0, 0.0, 0.0], dtype=float)
    else:
        sky_x = sky_x / np.linalg.norm(sky_x)
    sky_y = np.cross(sightline, sky_x)
    sky_y = sky_y / np.linalg.norm(sky_y)
    return sightline, sky_x, sky_y


def _direction_vector(theta: float, phi: float) -> np.ndarray:
    return np.array(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ],
        dtype=float,
    )


def _angular_diameter_distance_cm(observer: Observer) -> float:
    if observer.lumi_dist_cm is None or observer.lumi_dist_cm <= 0.0:
        raise ValueError("Observer.lumi_dist_cm must be set for sky_image().")
    return observer.lumi_dist_cm / (1.0 + observer.z) ** 2


def _project_patch_to_sky(
    radius_cm: np.ndarray,
    theta_center: float,
    phi_center: float,
    sightline: np.ndarray,
    sky_x_axis: np.ndarray,
    sky_y_axis: np.ndarray,
    angular_diameter_distance_cm: float,
) -> tuple[np.ndarray, np.ndarray]:
    direction = _direction_vector(theta_center, phi_center)
    position = radius_cm[:, None] * direction[None, :]
    line_of_sight = np.sum(position * sightline[None, :], axis=1)
    transverse = position - line_of_sight[:, None] * sightline[None, :]
    x_center = transverse @ sky_x_axis / angular_diameter_distance_cm
    y_center = transverse @ sky_y_axis / angular_diameter_distance_cm
    return x_center, y_center


def _gaussian_splat_stack(
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    x_center: np.ndarray,
    y_center: np.ndarray,
    sigma: np.ndarray,
) -> np.ndarray:
    sigma = np.maximum(np.asarray(sigma, dtype=float), 1.0e-30)
    dx = x_grid[None, :, :] - np.asarray(x_center, dtype=float)[:, None, None]
    dy = y_grid[None, :, :] - np.asarray(y_center, dtype=float)[:, None, None]
    sigma_view = sigma[:, None, None]
    exponent = -0.5 * (dx * dx + dy * dy) / (sigma_view * sigma_view)
    return np.exp(exponent) / (2.0 * np.pi * sigma_view * sigma_view)


def _interpolate_positive_series(source_t: np.ndarray, source_y: np.ndarray, target_t: np.ndarray) -> np.ndarray:
    source_t = np.asarray(source_t, dtype=float)
    source_y = np.asarray(source_y, dtype=float)
    target_t = np.asarray(target_t, dtype=float)
    if np.all(source_t > 0.0) and np.all(source_y > 0.0) and np.all(target_t > 0.0):
        return np.exp(
            np.interp(
                np.log(target_t),
                np.log(source_t),
                np.log(source_y),
                left=np.log(source_y[0]),
                right=np.log(source_y[-1]),
            )
        )
    return np.interp(target_t, source_t, source_y, left=source_y[0], right=source_y[-1])


def _set_dotted_attr(obj: Any, path: str, value: Any) -> None:
    target = obj
    parts = path.split(".")
    for name in parts[:-1]:
        target = getattr(target, name)
    setattr(target, parts[-1], value)


def _evaluate_flux_observations(model: Model, times_s: np.ndarray, frequencies_hz: np.ndarray) -> np.ndarray:
    times_s = np.asarray(times_s, dtype=float)
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    if frequencies_hz.ndim != 1 or times_s.ndim != 1:
        raise ValueError("Flux-density observations must be one-dimensional arrays.")
    if frequencies_hz.shape == times_s.shape:
        unique_freqs, inverse = np.unique(frequencies_hz, return_inverse=True)
        grid = model.flux_density_grid(times_s, unique_freqs).total
        return grid[inverse, np.arange(times_s.shape[0])]
    return model.flux_density_grid(times_s, frequencies_hz).total


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


def _jet_magnetar_active(jet: Jet, theta_center: float) -> bool:
    if isinstance(jet, (TophatJet, GaussianJet, PowerLawJet, StepPowerLawJet)):
        return theta_center <= getattr(jet, "theta_c", getattr(jet, "theta_j", jet.theta_max))
    if isinstance(jet, TwoComponentJet):
        return theta_center <= jet.theta_n
    return True


def _project_medium_to_kernel(medium: Medium, *, phi_center: float, theta_center: float) -> dict[str, float]:
    if isinstance(medium, (ISM, Wind)):
        return medium.to_kernel_params()

    radius = np.logspace(9.0, 19.0, 256)
    density = np.asarray(medium.density(phi_center, theta_center, radius), dtype=float)
    if density.shape != radius.shape:
        raise ValueError("Custom medium callable must return a density array with the same shape as the radius grid.")
    if not np.all(np.isfinite(density)) or np.any(density <= 0.0):
        raise ValueError("Custom medium callable must return positive finite densities.")

    candidates = [
        _fit_ism_medium(radius, density),
        _fit_wind_medium(radius, density),
        _fit_jump_medium(radius, density),
    ]
    errors = [_medium_fit_error(radius, density, candidate) for candidate in candidates]
    return candidates[int(np.argmin(errors))]


def _fit_ism_medium(radius: np.ndarray, density: np.ndarray) -> dict[str, float]:
    return {
        "d_ne": float(np.exp(np.mean(np.log(density)))),
        "a_star": -1.0,
        "r0": 1.0e9,
        "r_tr": 1.0e18,
        "f_jump": 1.0,
        "f_wide": 0.1,
    }


def _fit_wind_medium(radius: np.ndarray, density: np.ndarray) -> dict[str, float]:
    logr = np.log(radius)
    logd = np.log(density)
    slopes = np.gradient(logd, logr)
    wind_mask = np.abs(slopes + 2.0) < 0.35
    if np.any(wind_mask):
        a_star = float(np.median(density[wind_mask] * radius[wind_mask] ** 2 / 3.0e35))
    else:
        a_star = 0.0

    tail_slope = float(np.median(slopes[-16:]))
    floor = 0.0 if abs(tail_slope + 2.0) < 0.35 else float(np.min(density))
    r0 = 0.0
    if a_star > 0.0:
        wind = a_star * 3.0e35 / radius**2
        head_slope = float(np.median(slopes[:16]))
        if abs(head_slope) < 0.35 and density[0] < 0.95 * wind[0]:
            plateau_mask = np.abs(slopes) < 0.35
            if np.any(plateau_mask):
                n0 = float(np.median(density[plateau_mask]))
                r0 = float(np.sqrt(a_star * 3.0e35 / n0))
    return {
        "d_ne": floor,
        "a_star": a_star,
        "r0": r0,
        "r_tr": 1.0e18,
        "f_jump": 1.0,
        "f_wide": 0.1,
    }


def _fit_jump_medium(radius: np.ndarray, density: np.ndarray) -> dict[str, float]:
    floor = float(np.min(density))
    peak_idx = int(np.argmax(density))
    peak = float(density[peak_idx])
    logr = np.log10(radius)
    excess = np.maximum(density - floor, 0.0)
    if peak <= floor or np.sum(excess) <= 0.0:
        return _fit_ism_medium(radius, density)
    center = logr[peak_idx]
    variance = np.sum(excess * (logr - center) ** 2) / np.sum(excess)
    width = float(np.sqrt(max(variance, 0.03**2)))
    return {
        "d_ne": floor,
        "a_star": -1.0,
        "r0": 1.0e9,
        "r_tr": float(radius[peak_idx]),
        "f_jump": float(max(peak / floor, 1.0)),
        "f_wide": width,
    }


def _medium_fit_error(radius: np.ndarray, density: np.ndarray, params: dict[str, float]) -> float:
    model = _evaluate_kernel_density(radius, params)
    return float(np.mean((np.log10(model) - np.log10(density)) ** 2))


def _evaluate_kernel_density(radius: np.ndarray, params: dict[str, float]) -> np.ndarray:
    radius = np.asarray(radius, dtype=float)
    if params["a_star"] > 0.0:
        wind = params["a_star"] * 3.0e35 / radius**2
        density = np.maximum(params["d_ne"], wind)
    else:
        density = params["d_ne"] * (
            1.0
            + (params["f_jump"] - 1.0)
            * np.exp(-(np.log10(radius) - np.log10(params["r_tr"])) ** 2 / (2.0 * params["f_wide"] ** 2))
        )
    if np.any(radius < params["r0"]) and params["a_star"] > 0.0:
        density = density.copy()
        density[radius < params["r0"]] = params["a_star"] * 3.0e35 / params["r0"] ** 2
    return density


def _resolve_param_path(model: Model, param: ParamDef) -> str:
    if param.path is not None:
        return param.path
    name = param.name.lower()
    alias_map = {
        "e_iso": "jet.E_iso",
        "log10_eiso": "jet.E_iso",
        "log10_e_iso": "jet.E_iso",
        "gamma0": "jet.lf",
        "log10_gamma0": "jet.lf",
        "theta_c": "jet.theta_j" if isinstance(model.jet, TophatJet) else "jet.theta_c",
        "theta_j": "jet.theta_j",
        "theta_obs": "observer.theta_obs",
        "z": "observer.z",
        "lumi_dist": "observer.lumi_dist_cm",
        "lumi_dist_cm": "observer.lumi_dist_cm",
        "eps_e": "fwd_rad.eps_e",
        "epsilon_e": "fwd_rad.eps_e",
        "eps_b": "fwd_rad.eps_B",
        "epsilon_b": "fwd_rad.eps_B",
        "p": "fwd_rad.p",
        "xi_n": "fwd_rad.xi_N",
        "f_e": "fwd_rad.xi_N",
        "n0": "medium.n_ism",
        "n_ism": "medium.n_ism",
        "d_ne": "medium.n_ism",
        "a_star": "medium.A_star",
        "astar": "medium.A_star",
        "e_iso_c": "jet.E_iso_c",
        "e_iso_n": "jet.E_iso_n",
        "e_iso_outer": "jet.E_iso_w",
        "e_iso_w": "jet.E_iso_w",
        "gamma0_c": "jet.lf_c",
        "gamma0_n": "jet.lf_n",
        "gamma0_outer": "jet.lf_w",
        "gamma0_w": "jet.lf_w",
        "theta_n": "jet.theta_n",
        "theta_o": "jet.theta_w",
        "theta_w": "jet.theta_w",
        "k": "jet.k_e",
        "k_e": "jet.k_e",
        "k_g": "jet.k_g",
        "duration": "jet.duration",
        "l0": "jet.magnetar.L0",
        "t0": "jet.magnetar.t0",
        "q": "jet.magnetar.q",
        "eps_e_r": "rvs_rad.eps_e",
        "epsilon_e_r": "rvs_rad.eps_e",
        "eps_b_r": "rvs_rad.eps_B",
        "epsilon_b_r": "rvs_rad.eps_B",
        "p_r": "rvs_rad.p",
        "xi_n_r": "rvs_rad.xi_N",
        "f_e_r": "rvs_rad.xi_N",
    }
    if name not in alias_map:
        raise KeyError(f"Cannot infer parameter path for {param.name}.")
    return alias_map[name]


def _coerce_model(cfg: Any) -> Model:
    if isinstance(cfg, Model):
        return cfg
    if cfg is None:
        raise ValueError("Either a Model or cfg must be provided.")
    if isinstance(cfg, Setups):
        cfg = cfg.__dict__.copy()
    if not isinstance(cfg, dict):
        raise TypeError("cfg must be a Model or a dictionary of model options.")
    setups_source = cfg.get("setups", Setups())
    setups = deepcopy(setups_source if isinstance(setups_source, Setups) else Setups(**setups_source))
    observer = cfg.get(
        "observer",
        Observer(
            z=cfg.get("z", setups.z),
            theta_obs=cfg.get("theta_obs", setups.theta_obs),
            phi_obs=cfg.get("phi_obs", setups.phi_obs),
            lumi_dist=cfg.get("lumi_dist", setups.lumi_dist),
            lumi_dist_cm=cfg.get("lumi_dist_cm"),
        ),
    )
    medium = cfg.get("medium")
    medium_kind = str(cfg.get("medium_type", cfg.get("medium_name", setups.medium or "ism"))).lower()
    if isinstance(medium, str):
        medium_kind = medium.lower()
        medium = None
    if medium is None:
        if medium_kind == "wind":
            medium = Wind(A_star=cfg.get("A_star", cfg.get("Astar", 1.0)), n0=cfg.get("n0", cfg.get("n_ism", 0.1)))
        else:
            medium = ISM(n0=cfg.get("n0", cfg.get("n_ism", 0.1)))

    jet = cfg.get("jet")
    jet_kind = str(cfg.get("jet_type", cfg.get("jet_name", setups.jet or "tophat"))).lower()
    if isinstance(jet, str):
        jet_kind = jet.lower()
        jet = None
    if jet is None:
        magnetar = None
        if "magnetar" in cfg and cfg["magnetar"] is not None:
            source = cfg["magnetar"]
            if isinstance(source, Magnetar):
                magnetar = source
            else:
                magnetar = Magnetar(L0=source["L0"], t0=source["t0"], q=source["q"])
        elif {"L0", "t0", "q"} <= set(cfg.keys()):
            magnetar = Magnetar(L0=cfg["L0"], t0=cfg["t0"], q=cfg["q"])
        if jet_kind == "gaussian":
            jet = GaussianJet(
                E_iso=cfg["E_iso"],
                Gamma0=cfg["Gamma0"],
                theta_c=cfg["theta_c"],
                theta_max=cfg.get("theta_max", 0.6),
                duration=cfg.get("duration"),
                magnetar=magnetar,
                spreading=cfg.get("spreading", False),
            )
        elif jet_kind == "powerlaw":
            jet = PowerLawJet(
                E_iso=cfg["E_iso"],
                Gamma0=cfg["Gamma0"],
                theta_c=cfg["theta_c"],
                k=cfg.get("k"),
                k_e=cfg.get("k_e"),
                k_g=cfg.get("k_g"),
                theta_max=cfg.get("theta_max", np.pi / 2.0),
                duration=cfg.get("duration"),
                magnetar=magnetar,
                spreading=cfg.get("spreading", False),
            )
        elif jet_kind == "twocomponent":
            jet = TwoComponentJet(
                E_iso_c=cfg["E_iso_c"],
                Gamma0_c=cfg["Gamma0_c"],
                theta_c=cfg["theta_c"],
                E_iso_outer=cfg["E_iso_outer"],
                Gamma0_outer=cfg["Gamma0_outer"],
                theta_o=cfg["theta_o"],
                duration=cfg.get("duration"),
                magnetar=magnetar,
                spreading=cfg.get("spreading", False),
            )
        elif jet_kind == "steppowerlaw":
            jet = StepPowerLawJet(
                E_iso_c=cfg["E_iso_c"],
                Gamma0_c=cfg["Gamma0_c"],
                theta_c=cfg["theta_c"],
                E_iso_w=cfg["E_iso_w"],
                Gamma0_w=cfg["Gamma0_w"],
                theta_w=cfg["theta_w"],
                k=cfg.get("k", 2.0),
                k_e=cfg.get("k_e"),
                k_g=cfg.get("k_g"),
                duration=cfg.get("duration"),
                magnetar=magnetar,
                spreading=cfg.get("spreading", False),
            )
        else:
            jet = TophatJet(
                E_iso=cfg["E_iso"],
                Gamma0=cfg["Gamma0"],
                theta_c=cfg["theta_c"],
                duration=cfg.get("duration"),
                magnetar=magnetar,
                spreading=cfg.get("spreading", False),
            )

    fwd_rad = cfg.get(
        "fwd_rad",
        Radiation(
            eps_e=cfg.get("eps_e", cfg.get("epsilon_e", 1.0e-1)),
            eps_B=cfg.get("eps_B", cfg.get("epsilon_B", 1.0e-3)),
            p=cfg.get("p", 2.5),
            xi_N=cfg.get("xi_N", cfg.get("f_e", 1.0e-1)),
            ssc=cfg.get("ssc", setups.fwd_ssc),
            kn=cfg.get("kn", setups.kn),
        ),
    )
    rvs_rad = cfg.get("rvs_rad")
    if rvs_rad is None and cfg.get("reverse", setups.rvs_shock):
        rvs_rad = Radiation(
            eps_e=cfg.get("eps_e_r", cfg.get("eps_e", 1.0e-1)),
            eps_B=cfg.get("eps_B_r", cfg.get("eps_B", 1.0e-2)),
            p=cfg.get("p_r", cfg.get("p", 2.4)),
            xi_N=cfg.get("xi_N_r", cfg.get("f_e_r", 1.0)),
            ssc=cfg.get("rvs_ssc", setups.rvs_ssc),
            kn=cfg.get("kn", setups.kn),
        )

    resolutions = cfg.get("resolutions")
    return Model(medium=medium, jet=jet, observer=observer, fwd_rad=fwd_rad, rvs_rad=rvs_rad, setups=setups, resolutions=resolutions)
