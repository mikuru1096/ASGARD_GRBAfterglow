from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import constants

try:
    import src.Hadronic.FS_hadronic_1d as hadronic_fortran_module
except ImportError:
    hadronic_fortran_module = None


ELECTRON_MASS_GEV = constants.para_m_e_gev
_HAS_FORTRAN_PAIR_PRODUCTION = hadronic_fortran_module is not None and hasattr(
    hadronic_fortran_module, "fs_hadronic_pair_production_shell"
)
PAIR_PRODUCTION_BACKEND = "fortran_pair_production" if _HAS_FORTRAN_PAIR_PRODUCTION else "unavailable"


@dataclass(frozen=True)
class PairProductionOutput:
    photon_energy_gev: np.ndarray
    electron_energy_gev: np.ndarray
    photon_loss_rate: np.ndarray
    pair_injection_rate_per_gev_per_species: np.ndarray
    pair_injection_rate_per_gev_total: np.ndarray
    absorbed_power_gev_per_cm3_s: float
    injected_power_gev_per_cm3_s: float


def solve_pair_production(
    photon_energy_gev: np.ndarray,
    photon_density_per_gev: np.ndarray,
    electron_energy_gev: np.ndarray,
    max_com_energy_factor: int = 138,
) -> PairProductionOutput:
    e_ph = _as_strictly_increasing(photon_energy_gev, "photon_energy_gev")
    n_ph = _as_matching_and_nonnegative(photon_density_per_gev, e_ph, "photon_density_per_gev")
    e_e = _as_strictly_increasing(electron_energy_gev, "electron_energy_gev")
    if int(max_com_energy_factor) < 1:
        raise ValueError("max_com_energy_factor must be >= 1.")

    if not _HAS_FORTRAN_PAIR_PRODUCTION:
        raise RuntimeError("Pair production core must be provided by the Fortran backend.")

    (
        photon_loss_rate,
        pair_rate_per_gev_per_species,
        pair_rate_per_gev_total,
        absorbed_power,
        injected_power,
    ) = hadronic_fortran_module.fs_hadronic_pair_production_shell(
        e_ph,
        n_ph,
        e_e,
        int(max_com_energy_factor),
    )
    return PairProductionOutput(
        photon_energy_gev=e_ph,
        electron_energy_gev=e_e,
        photon_loss_rate=np.asarray(photon_loss_rate, dtype=float),
        pair_injection_rate_per_gev_per_species=np.asarray(pair_rate_per_gev_per_species, dtype=float),
        pair_injection_rate_per_gev_total=np.asarray(pair_rate_per_gev_total, dtype=float),
        absorbed_power_gev_per_cm3_s=float(absorbed_power),
        injected_power_gev_per_cm3_s=float(injected_power),
    )


def _as_strictly_increasing(values: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1d array.")
    if arr.size < 2:
        raise ValueError(f"{name} must contain at least two points.")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive.")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing.")
    return arr


def _as_matching_and_nonnegative(values: np.ndarray, grid: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != grid.shape:
        raise ValueError(f"{name} must match grid shape.")
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must be non-negative.")
    return arr
