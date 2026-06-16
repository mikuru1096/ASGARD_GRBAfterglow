from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import constants
from asgard_core.hadronic_validation import as_matching, as_strictly_increasing_grid

import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module


ELECTRON_MASS_GEV = constants.para_m_e_gev
PROTON_MASS_GEV = constants.para_m_p_gev

BETHE_HEITLER_BACKEND = "fortran_bethe_heitler"


@dataclass(frozen=True)
class BetheHeitlerOutput:
    electron_energy_gev: np.ndarray
    pair_rate_per_gev: np.ndarray
    proton_loss_rate: np.ndarray
    photon_loss_rate: np.ndarray


def solve_bethe_heitler(
    proton_energy_gev: np.ndarray,
    proton_density_per_gev: np.ndarray,
    photon_energy_gev: np.ndarray,
    photon_density_per_gev: np.ndarray,
    electron_energy_gev: np.ndarray,
) -> BetheHeitlerOutput:
    e_p = as_strictly_increasing_grid(proton_energy_gev, "proton_energy_gev")
    n_p = as_matching(proton_density_per_gev, e_p, "proton_density_per_gev")
    e_ph = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev")
    n_ph = as_matching(photon_density_per_gev, e_ph, "photon_density_per_gev")
    e_e = as_strictly_increasing_grid(electron_energy_gev, "electron_energy_gev")

    pair_rate_per_gev, proton_loss_rate, photon_loss_rate = hadronic_fortran_module.fs_hadronic_bethe_heitler_shell(
        e_p,
        n_p,
        e_ph,
        n_ph,
        e_e,
    )
    pair_rate_per_gev = np.asarray(pair_rate_per_gev, dtype=float)
    proton_loss_rate = np.asarray(proton_loss_rate, dtype=float)
    photon_loss_rate = np.asarray(photon_loss_rate, dtype=float)

    return BetheHeitlerOutput(
        electron_energy_gev=e_e,
        pair_rate_per_gev=pair_rate_per_gev,
        proton_loss_rate=proton_loss_rate,
        photon_loss_rate=photon_loss_rate,
    )
