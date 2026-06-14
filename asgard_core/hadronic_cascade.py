"""γγ pair/synch cascade helpers and shell-sequence time evolution."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from src import constants


@dataclass(frozen=True)
class TimeDependentPairCascadeOutput:
    photon_energy_gev: np.ndarray
    electron_energy_gev: np.ndarray
    frequency_hz: np.ndarray
    photon_loss_rate_local_s_inv: np.ndarray
    tau_pair_path: np.ndarray
    pair_density_per_gamma: np.ndarray
    pair_syn_luminosity_hz: np.ndarray
    pair_syn_seed_per_hz: np.ndarray
    cascade_photon_density_per_gev: np.ndarray
    absorbed_power_gev_per_cm3_s: np.ndarray
    injected_power_gev_per_cm3_s: np.ndarray


def shell_path_time_seconds(radius_cm: float, gamma_bulk: float) -> float:
    if radius_cm <= 0.0:
        raise ValueError("pair cascade shell radius must be positive.")
    if gamma_bulk < 1.0:
        raise ValueError("pair cascade bulk Lorentz factor must be >= 1.")
    return float(radius_cm) / (12.0 * float(gamma_bulk) * constants.para_c)


def compute_time_dependent_pair_cascade_sequence(
    photon_energy_gev: np.ndarray,
    primary_photon_density_per_gev: np.ndarray,
    electron_energy_gev: np.ndarray,
    frequency_hz: np.ndarray,
    radius_cm: np.ndarray,
    gamma_bulk: np.ndarray,
    observer_time_s: np.ndarray,
    b_field_g: np.ndarray,
    *,
    num_threads: int,
    index_syn_integr: int,
    substeps_per_shell: int,
) -> TimeDependentPairCascadeOutput:
    """Advance the gamma-gamma pair/synch cascade as a Fortran shell sequence."""
    from src.Hadronic import hadronic_forward_1d as hadronic_module

    e_ph = _as_positive_grid(photon_energy_gev, "photon_energy_gev")
    e_e = _as_positive_grid(electron_energy_gev, "electron_energy_gev")
    nu = _as_positive_grid(frequency_hz, "frequency_hz")
    primary = np.asarray(primary_photon_density_per_gev, dtype=float)
    radius = _as_positive_grid(radius_cm, "radius_cm")
    gamma = np.asarray(gamma_bulk, dtype=float)
    tobs = np.asarray(observer_time_s, dtype=float)
    b_field = np.asarray(b_field_g, dtype=float)
    if primary.shape != (e_ph.size, radius.size):
        raise ValueError("primary photon density must have shape (num_photon, num_shell).")
    if e_ph.shape != nu.shape:
        raise ValueError("photon_energy_gev and frequency_hz must have matching shapes.")
    if gamma.shape != radius.shape or tobs.shape != radius.shape or b_field.shape != radius.shape:
        raise ValueError("gamma_bulk, observer_time_s, and b_field_g must match radius_cm.")
    if np.any(primary < 0.0) or np.any(gamma < 1.0) or np.any(b_field < 0.0):
        raise ValueError("pair cascade sequence received non-physical inputs.")
    if int(substeps_per_shell) < 1:
        raise ValueError("substeps_per_shell must be positive.")

    (
        photon_loss,
        tau_pair,
        pair_density,
        pair_lum,
        pair_seed,
        cascade_photon_density,
        absorbed_power,
        injected_power,
    ) = hadronic_module.fs_hadronic_pair_cascade_sequence(
        e_ph,
        primary,
        e_e,
        nu,
        radius,
        gamma,
        tobs,
        b_field,
        int(num_threads),
        int(index_syn_integr),
        int(substeps_per_shell),
    )

    return TimeDependentPairCascadeOutput(
        photon_energy_gev=e_ph,
        electron_energy_gev=e_e,
        frequency_hz=nu,
        photon_loss_rate_local_s_inv=photon_loss,
        tau_pair_path=tau_pair,
        pair_density_per_gamma=pair_density,
        pair_syn_luminosity_hz=pair_lum,
        pair_syn_seed_per_hz=pair_seed,
        cascade_photon_density_per_gev=cascade_photon_density,
        absorbed_power_gev_per_cm3_s=absorbed_power,
        injected_power_gev_per_cm3_s=injected_power,
    )


def _as_positive_grid(values: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1D grid.")
    if arr.size < 2:
        raise ValueError(f"{name} must contain at least two points.")
    if np.any(arr <= 0.0) or np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be positive and strictly increasing.")
    return arr
