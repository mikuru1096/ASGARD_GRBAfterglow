from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import constants
from asgard_core.hadronic_validation import as_matching_nonnegative, as_strictly_increasing_grid

import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module


ELECTRON_MASS_GEV = constants.para_m_e_gev


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
    e_ph = as_strictly_increasing_grid(photon_energy_gev, "photon_energy_gev")
    n_ph = as_matching_nonnegative(photon_density_per_gev, e_ph, "photon_density_per_gev")
    e_e = as_strictly_increasing_grid(electron_energy_gev, "electron_energy_gev")
    if int(max_com_energy_factor) < 1:
        raise ValueError("max_com_energy_factor must be >= 1.")

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
