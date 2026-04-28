from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import constants
from asgard_core.hadronic_pair_production import solve_pair_production


@dataclass(frozen=True)
class PairCascadeKernels:
    frequency_hz: np.ndarray
    synchrotron_hz_per_gev: np.ndarray
    inverse_compton_hz_per_gev: np.ndarray


@dataclass(frozen=True)
class PairCascadeOutput:
    photon_energy_gev: np.ndarray
    electron_energy_gev: np.ndarray
    frequency_hz: np.ndarray
    photon_loss_rate_local_s_inv: np.ndarray
    tau_pair_path: np.ndarray
    shell_path_time_s: float
    pair_injection_rate_per_gev_total_local: np.ndarray
    pair_syn_luminosity_hz: np.ndarray
    pair_ic_luminosity_hz: np.ndarray
    pair_total_luminosity_hz: np.ndarray
    absorbed_power_gev_per_cm3_s: float
    injected_power_gev_per_cm3_s: float


def shell_path_time_seconds(radius_cm: float, gamma_bulk: float) -> float:
    radius = float(radius_cm)
    gamma = float(gamma_bulk)
    if radius <= 0.0:
        raise ValueError("radius_cm must be > 0.")
    if gamma <= 0.0:
        raise ValueError("gamma_bulk must be > 0.")
    return radius / (12.0 * gamma * constants.para_c)


def compute_pair_cascade_bookkeeping(
    *,
    photon_energy_gev: np.ndarray,
    photon_density_per_gev: np.ndarray,
    electron_energy_gev: np.ndarray,
    radius_cm: float,
    gamma_bulk: float,
    kernels: PairCascadeKernels,
    max_com_energy_factor: int = 138,
) -> PairCascadeOutput:
    pair = solve_pair_production(
        photon_energy_gev=photon_energy_gev,
        photon_density_per_gev=photon_density_per_gev,
        electron_energy_gev=electron_energy_gev,
        max_com_energy_factor=max_com_energy_factor,
    )

    pair_injection = np.asarray(pair.pair_injection_rate_per_gev_total, dtype=float)
    photon_loss_rate_local = np.asarray(pair.photon_loss_rate, dtype=float)
    path_time_s = shell_path_time_seconds(radius_cm=radius_cm, gamma_bulk=gamma_bulk)
    tau_pair_path = photon_loss_rate_local * path_time_s

    frequency_hz = _as_strictly_increasing(kernels.frequency_hz, "kernels.frequency_hz")
    syn_kernel = _as_projection_kernel(
        kernels.synchrotron_hz_per_gev,
        "kernels.synchrotron_hz_per_gev",
        out_size=frequency_hz.size,
        in_size=pair_injection.size,
    )
    ic_kernel = _as_projection_kernel(
        kernels.inverse_compton_hz_per_gev,
        "kernels.inverse_compton_hz_per_gev",
        out_size=frequency_hz.size,
        in_size=pair_injection.size,
    )

    pair_syn = syn_kernel @ pair_injection
    pair_ic = ic_kernel @ pair_injection

    return PairCascadeOutput(
        photon_energy_gev=np.asarray(pair.photon_energy_gev, dtype=float),
        electron_energy_gev=np.asarray(pair.electron_energy_gev, dtype=float),
        frequency_hz=frequency_hz,
        photon_loss_rate_local_s_inv=photon_loss_rate_local,
        tau_pair_path=tau_pair_path,
        shell_path_time_s=path_time_s,
        pair_injection_rate_per_gev_total_local=pair_injection,
        pair_syn_luminosity_hz=pair_syn,
        pair_ic_luminosity_hz=pair_ic,
        pair_total_luminosity_hz=pair_syn + pair_ic,
        absorbed_power_gev_per_cm3_s=float(pair.absorbed_power_gev_per_cm3_s),
        injected_power_gev_per_cm3_s=float(pair.injected_power_gev_per_cm3_s),
    )


def _as_projection_kernel(
    kernel: np.ndarray,
    name: str,
    *,
    out_size: int,
    in_size: int,
) -> np.ndarray:
    arr = np.asarray(kernel, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be a 2D array.")
    if arr.shape != (out_size, in_size):
        raise ValueError(
            f"{name} must have shape {(out_size, in_size)}, got {arr.shape}."
        )
    return arr


def _as_strictly_increasing(values: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1D array.")
    if arr.size < 2:
        raise ValueError(f"{name} must have at least two points.")
    if np.any(arr[1:] <= arr[:-1]):
        raise ValueError(f"{name} must be strictly increasing.")
    return arr
