from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from src import constants

try:
    import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module
except ImportError:
    hadronic_fortran_module = None


FINE_STRUCTURE_ALPHA = 1.0 / 137.0
ELECTRON_MASS_GEV = constants.para_m_e_gev
PROTON_MASS_GEV = constants.para_m_p_gev

_HAS_FORTRAN_BH = hadronic_fortran_module is not None and hasattr(
    hadronic_fortran_module, "fs_hadronic_bethe_heitler_shell"
)
BETHE_HEITLER_BACKEND = "fortran_bethe_heitler" if _HAS_FORTRAN_BH else "unavailable"


@dataclass(frozen=True)
class BetheHeitlerOutput:
    electron_energy_gev: np.ndarray
    pair_rate_per_gev: np.ndarray
    proton_loss_rate: np.ndarray


def solve_bethe_heitler(
    proton_energy_gev: np.ndarray,
    proton_density_per_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    photon_density_per_gev: np.ndarray,
    electron_energy_gev: np.ndarray,
) -> BetheHeitlerOutput:
    e_p = _as_strictly_increasing(proton_energy_gev, "proton_energy_gev")
    n_p = _as_matching(proton_density_per_gev, e_p, "proton_density_per_gev")
    e_ph = _as_strictly_increasing(photon_energy_gev, "photon_energy_gev")
    n_ph = _as_matching(photon_density_per_gev, e_ph, "photon_density_per_gev")
    e_e = _as_strictly_increasing(electron_energy_gev, "electron_energy_gev")

    if not _HAS_FORTRAN_BH:
        raise RuntimeError("Bethe-Heitler core must be provided by the Fortran backend.")

    pair_rate_per_gev, proton_loss_rate = hadronic_fortran_module.fs_hadronic_bethe_heitler_shell(
        e_p,
        n_p,
        e_ph,
        n_ph,
        e_e,
    )
    pair_rate_per_gev = np.asarray(pair_rate_per_gev, dtype=float)
    proton_loss_rate = np.asarray(proton_loss_rate, dtype=float)

    return BetheHeitlerOutput(
        electron_energy_gev=e_e,
        pair_rate_per_gev=pair_rate_per_gev,
        proton_loss_rate=proton_loss_rate,
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


def _as_matching(values: np.ndarray, grid: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != grid.shape:
        raise ValueError(f"{name} must match grid shape.")
    return arr


def _log_spacing(grid: np.ndarray, name: str) -> float:
    dln = np.diff(np.log(grid))
    if not np.allclose(dln, dln[0], rtol=1.0e-6, atol=1.0e-12):
        raise ValueError(f"{name} must be logarithmically uniform.")
    return float(dln[0])
