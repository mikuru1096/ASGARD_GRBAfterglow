from __future__ import annotations

from typing import Any

import numpy as np

from .api_fit import Param
from .api_model import (
    Hadronic,
    Model,
    Numerics,
    Observer,
    ObserverGrid,
    PolarizationResult,
    Radiation,
    ReverseShock,
    SolverOptions,
)


def _polarization(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    *,
    magnetic_geometry: str,
    local_emissivity: str,
) -> PolarizationResult:
    times_s = np.asarray(times_s, dtype=float)
    nu_hz = np.asarray(nu_hz, dtype=float)
    if times_s.ndim != 1 or nu_hz.ndim != 1:
        raise ValueError("polarization() requires 1D times_s and nu_hz grids.")
    if times_s.size == 0 or nu_hz.size == 0:
        raise ValueError("polarization() grids must be non-empty.")
    if np.any(times_s <= 0.0) or np.any(nu_hz <= 0.0):
        raise ValueError("polarization() times and frequencies must be positive.")
    if magnetic_geometry not in {"shock_random", "toroidal"}:
        raise ValueError("magnetic_geometry must be 'shock_random' or 'toroidal'.")
    if local_emissivity not in {"analytic", "analytic_then_kernel"}:
        raise ValueError("local_emissivity must be 'analytic' or 'analytic_then_kernel'.")
    raise NotImplementedError(
        "polarization currently has no supported backend after the Python patch projection path was removed."
    )


def _skyimage(model: Model, times_s: np.ndarray, nu_obs: float, fov: float, npixel: int):
    if model.observer.lumi_dist_cm is None or model.observer.lumi_dist_cm <= 0.0:
        raise ValueError("Observer.luminosity_distance_cm must be set for sky_image().")
    raise NotImplementedError(
        "sky_image currently has no supported backend after the Python patch projection path was removed."
    )


def _parampath(model: Model, param: Param) -> str:
    if param.path is not None:
        return param.path
    name = param.name
    path_map = {
        "energy_iso_erg": "jet.energy_iso_erg",
        "log10_energy_iso_erg": "jet.energy_iso_erg",
        "initial_lorentz_factor": "jet.initial_lorentz_factor",
        "log10_initial_lorentz_factor": "jet.initial_lorentz_factor",
        "core_angle_rad": "jet.opening_angle_rad" if model.jet.kind == "tophat" else "jet.core_angle_rad",
        "opening_angle_rad": "jet.opening_angle_rad",
        "outer_angle_rad": "jet.outer_angle_rad",
        "energy_index": "jet.energy_index",
        "lorentz_index": "jet.lorentz_index",
        "shell_duration_s": "jet.shell_duration_s",
        "viewing_angle_rad": "observer.viewing_angle_rad",
        "viewing_azimuth_rad": "observer.viewing_azimuth_rad",
        "z": "observer.z",
        "luminosity_distance_cm": "observer.luminosity_distance_cm",
        "epsilon_e": "fwd_rad.epsilon_e",
        "epsilon_B": "fwd_rad.epsilon_B",
        "p": "fwd_rad.p",
        "proton_energy_fraction": "fwd_rad.proton_energy_fraction",
        "epsilon_b_floor": "fwd_rad.epsilon_b_floor",
        "magnetic_decay_alpha_t": "fwd_rad.magnetic_decay_alpha_t",
        "magnetic_decay_t0_s": "fwd_rad.magnetic_decay_t0_s",
        "accelerated_electron_fraction": "fwd_rad.accelerated_electron_fraction",
        "acceleration_efficiency": "fwd_rad.acceleration_efficiency",
        "reverse_proton_energy_fraction": "fwd_rad.reverse_proton_energy_fraction",
        "number_density_cm3": "medium.number_density_cm3",
        "a_star": "medium.A_star",
        "density_floor_cm3": "medium.density_floor_cm3",
        "density_cap_cm3": "medium.density_cap_cm3",
        "reverse_shock_enabled": "reverse_shock.enabled",
        "reverse_shock_shell_duration_s": "reverse_shock.shell_duration_s",
        "reverse_shock_upstream_sigma": "reverse_shock.upstream_sigma",
        "hadronic_num_proton_gamma": "hadronic.num_proton_gamma",
        "hadronic_num_neutrino_frequency": "hadronic.num_neutrino_frequency",
        "hadronic_pair_cascade_iterations": "hadronic.pair_cascade_iterations",
    }
    if name not in path_map:
        raise KeyError(f"Cannot infer canonical parameter path for {param.name}; pass Param(name, path, lower, upper, scale).")
    return path_map[name]


def _asmodel(cfg: Any) -> Model:
    if isinstance(cfg, Model):
        return cfg
    if cfg is None:
        raise ValueError("Either a Model or cfg must be provided.")
    if not isinstance(cfg, dict):
        raise TypeError("cfg must be a Model or a dictionary containing explicit public objects.")

    def build(key: str, cls):
        value = cfg[key]
        return value if isinstance(value, cls) else cls(**value)

    return Model(
        medium=cfg["medium"],
        jet=cfg["jet"],
        observer=build("observer", Observer),
        fwd_rad=build("fwd_rad", Radiation),
        rvs_rad=None if cfg.get("rvs_rad") is None else build("rvs_rad", Radiation),
        numerics=build("numerics", Numerics),
        observer_grid=build("observer_grid", ObserverGrid),
        solver_options=build("solver_options", SolverOptions),
        reverse_shock=build("reverse_shock", ReverseShock),
        hadronic=build("hadronic", Hadronic),
    )
