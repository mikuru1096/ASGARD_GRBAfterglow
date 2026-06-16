from __future__ import annotations

import numpy as np

import src.Hadronic.hadronic_forward_1d as hadronic_legacy_module
from src import constants


PROTON_MASS_GEV = constants.para_m_p_e * constants.para_erg2ev * 1.0e-9


def photon_density_hz_to_gev(photon_nu_hz: np.ndarray, photon_density_per_hz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Convert `n_nu` to `n_E` using the Fortran hadronic utility."""
    nu = np.asarray(photon_nu_hz, dtype=float)
    density = np.asarray(photon_density_per_hz, dtype=float)
    if nu.ndim != 1 or density.ndim != 1 or nu.shape != density.shape:
        raise ValueError("photon_nu_hz and photon_density_per_hz must be 1d arrays with the same shape.")
    if np.any(nu <= 0.0):
        raise ValueError("photon_nu_hz must be positive.")
    photon_energy_gev, photon_density_per_gev = hadronic_legacy_module.fs_hadronic_photon_density_hz_to_gev(
        nu,
        density,
    )
    return np.asarray(photon_energy_gev, dtype=float), np.asarray(photon_density_per_gev, dtype=float)
