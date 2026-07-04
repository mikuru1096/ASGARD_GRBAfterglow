"""γγ pair/synch cascade helpers and shell-sequence time evolution."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from asgard_core.hadronic_processes import loginterp
from src import constants


ME_GEV = constants.para_m_e_gev


@dataclass(frozen=True, slots=True)
class PairCascade:
    """γγ pair/synch fields passed back to photon-field assembly."""

    syn_lum: np.ndarray
    syn_seed: np.ndarray
    tau_pair: np.ndarray
    density_grid: np.ndarray


def solve_paircascade(
    seed_hz: np.ndarray,
    seed_field: np.ndarray,
    electron_gamma: np.ndarray,
    radius_cm: np.ndarray,
    gamma_bulk: np.ndarray,
    tobs_s: np.ndarray,
    bfield_g: np.ndarray,
    *,
    threads: int,
    syn_index: int,
    substeps: int,
) -> PairCascade:
    """Advance the gamma-gamma pair/synch cascade as a Fortran shell sequence."""
    from src.Hadronic import hadronic_forward_1d as hadronic_module

    nu = _posgrid(seed_hz, "seed_hz")
    seed = np.asarray(seed_field, dtype=float)
    gamma_e = _posgrid(electron_gamma, "electron_gamma")
    radius = _posgrid(radius_cm, "radius_cm")
    gamma = np.asarray(gamma_bulk, dtype=float)
    tobs = np.asarray(tobs_s, dtype=float)
    bfield = np.asarray(bfield_g, dtype=float)
    if seed.shape != (nu.size, radius.size):
        raise ValueError("seed_field must have shape (num_photon, num_shell).")
    if gamma.shape != radius.shape or tobs.shape != radius.shape or bfield.shape != radius.shape:
        raise ValueError("gamma_bulk, tobs_s, and bfield_g must match radius_cm.")
    if np.any(seed < 0.0) or np.any(gamma < 1.0) or np.any(bfield < 0.0):
        raise ValueError("pair cascade sequence received non-physical inputs.")
    if int(substeps) < 1:
        raise ValueError("substeps must be positive.")

    photon_gev = constants.para_h_gev * nu
    dln = float(np.log(photon_gev[1] / photon_gev[0]))
    pair_offset = max(0, int(np.ceil(np.log(ME_GEV / photon_gev[0]) / dln)))
    pair_gev = photon_gev[0] * np.exp((pair_offset + np.arange(photon_gev.size, dtype=float)) * dln)
    pair_gamma = pair_gev / ME_GEV

    (
        _,
        tau,
        density,
        lum,
        seedout,
        *_,
    ) = hadronic_module.cascade_sequence(
        photon_gev,
        seed / constants.para_h_gev,
        pair_gev,
        nu,
        radius,
        gamma,
        tobs,
        bfield,
        int(threads),
        int(syn_index),
        int(substeps),
    )

    density_grid = np.zeros((gamma_e.size, density.shape[1]), dtype=float)
    for i_shell in range(density.shape[1]):
        density_grid[:, i_shell] = loginterp(pair_gamma, density[:, i_shell], gamma_e)
    return PairCascade(
        syn_lum=lum,
        syn_seed=seedout,
        tau_pair=tau,
        density_grid=density_grid,
    )


def _posgrid(values: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1D grid.")
    if arr.size < 2:
        raise ValueError(f"{name} must contain at least 2 points.")
    if np.any(arr <= 0.0) or np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be positive and strictly increasing.")
    return arr
