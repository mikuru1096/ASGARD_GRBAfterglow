from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from asgard_core.hadronic_validation import as_matching, as_strictly_increasing_grid

import src.Hadronic.hadronic_forward_1d as hadronic_fortran_module


SECONDARY_RADIATION_BACKEND = "fortran_secondary_radiation"


@dataclass(frozen=True)
class SecondarySpeciesDistribution:
    pion_plus_per_gev: np.ndarray
    pion_minus_per_gev: np.ndarray
    muon_minus_left_per_gev: np.ndarray
    muon_minus_right_per_gev: np.ndarray
    muon_plus_left_per_gev: np.ndarray
    muon_plus_right_per_gev: np.ndarray


@dataclass(frozen=True)
class SecondaryTargetPhotonField:
    photon_energy_gev: np.ndarray
    photons_on_had_grid_per_gev: np.ndarray
    ind_min_energy_pho_hadgrid: int = 0


@dataclass(frozen=True)
class SecondaryRadiationSpectrum:
    photon_energy_gev: np.ndarray
    pion_synch_rate_per_gev: np.ndarray
    muon_synch_rate_per_gev: np.ndarray
    pion_ic_rate_per_gev: np.ndarray
    muon_ic_rate_per_gev: np.ndarray


def solve_secondary_radiation_spectrum(
    hadron_energy_gev: np.ndarray,
    species: SecondarySpeciesDistribution,
    target: SecondaryTargetPhotonField,
    magnetic_field_g: float,
) -> SecondaryRadiationSpectrum:
    e_had = as_strictly_increasing_grid(hadron_energy_gev, "hadron_energy_gev", require_finite=True)
    e_ph = as_strictly_increasing_grid(target.photon_energy_gev, "photon_energy_gev", require_finite=True)
    (
        pion_syn,
        muon_syn,
        pion_ic,
        muon_ic,
    ) = hadronic_fortran_module.fs_hadronic_secondary_radiation_shell(
        e_had,
        e_ph,
        as_matching(species.pion_plus_per_gev, e_had, "pion_plus_per_gev", require_finite=True),
        as_matching(species.pion_minus_per_gev, e_had, "pion_minus_per_gev", require_finite=True),
        as_matching(species.muon_minus_left_per_gev, e_had, "muon_minus_left_per_gev", require_finite=True),
        as_matching(species.muon_minus_right_per_gev, e_had, "muon_minus_right_per_gev", require_finite=True),
        as_matching(species.muon_plus_left_per_gev, e_had, "muon_plus_left_per_gev", require_finite=True),
        as_matching(species.muon_plus_right_per_gev, e_had, "muon_plus_right_per_gev", require_finite=True),
        as_matching(target.photons_on_had_grid_per_gev, e_ph, "photons_on_had_grid_per_gev", require_finite=True),
        int(target.ind_min_energy_pho_hadgrid),
        float(magnetic_field_g),
    )
    return SecondaryRadiationSpectrum(
        photon_energy_gev=e_ph,
        pion_synch_rate_per_gev=np.asarray(pion_syn, dtype=float),
        muon_synch_rate_per_gev=np.asarray(muon_syn, dtype=float),
        pion_ic_rate_per_gev=np.asarray(pion_ic, dtype=float),
        muon_ic_rate_per_gev=np.asarray(muon_ic, dtype=float),
    )
