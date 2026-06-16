from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.hadronic_processes import (  # noqa: E402
    SecondarySpeciesDistribution,
    SecondaryTargetPhotonField,
    solve_secondary_radiation_spectrum,
)


def _sample_inputs() -> tuple[np.ndarray, SecondarySpeciesDistribution, SecondaryTargetPhotonField, float]:
    dln = 0.2
    hadron_energy_gev = np.exp(np.arange(8, dtype=float) * dln)
    photon_energy_gev = np.exp((np.arange(12, dtype=float) - 4.0) * dln)
    species = SecondarySpeciesDistribution(
        pion_plus_per_gev=np.linspace(2.0e-11, 5.0e-11, hadron_energy_gev.size),
        pion_minus_per_gev=np.linspace(1.0e-11, 3.0e-11, hadron_energy_gev.size),
        muon_minus_left_per_gev=np.linspace(2.0e-12, 4.0e-12, hadron_energy_gev.size),
        muon_minus_right_per_gev=np.linspace(1.5e-12, 3.5e-12, hadron_energy_gev.size),
        muon_plus_left_per_gev=np.linspace(2.2e-12, 4.2e-12, hadron_energy_gev.size),
        muon_plus_right_per_gev=np.linspace(1.8e-12, 3.8e-12, hadron_energy_gev.size),
    )
    target = SecondaryTargetPhotonField(
        photon_energy_gev=photon_energy_gev,
        photons_on_had_grid_per_gev=np.linspace(2.0e-5, 8.0e-5, photon_energy_gev.size),
        ind_min_energy_pho_hadgrid=3,
    )
    return hadron_energy_gev, species, target, 25.0


def test_secondary_radiation_smoke() -> None:
    hadron_energy_gev, species, target, magnetic_field_g = _sample_inputs()
    output = solve_secondary_radiation_spectrum(
        hadron_energy_gev=hadron_energy_gev,
        species=species,
        target=target,
        magnetic_field_g=magnetic_field_g,
    )

    assert output.photon_energy_gev.shape == target.photon_energy_gev.shape
    for arr in (
        output.pion_synch_rate_per_gev,
        output.muon_synch_rate_per_gev,
        output.pion_ic_rate_per_gev,
        output.muon_ic_rate_per_gev,
    ):
        assert arr.shape == target.photon_energy_gev.shape
        assert np.all(np.isfinite(arr))
        assert np.all(arr >= 0.0)


def test_secondary_radiation_zero_inputs() -> None:
    hadron_energy_gev, species, target, magnetic_field_g = _sample_inputs()
    zeros = np.zeros_like(hadron_energy_gev)
    species = SecondarySpeciesDistribution(
        pion_plus_per_gev=zeros,
        pion_minus_per_gev=zeros,
        muon_minus_left_per_gev=zeros,
        muon_minus_right_per_gev=zeros,
        muon_plus_left_per_gev=zeros,
        muon_plus_right_per_gev=zeros,
    )
    target = SecondaryTargetPhotonField(
        photon_energy_gev=target.photon_energy_gev,
        photons_on_had_grid_per_gev=np.zeros_like(target.photon_energy_gev),
        ind_min_energy_pho_hadgrid=target.ind_min_energy_pho_hadgrid,
    )
    output = solve_secondary_radiation_spectrum(
        hadron_energy_gev=hadron_energy_gev,
        species=species,
        target=target,
        magnetic_field_g=magnetic_field_g,
    )
    assert np.all(output.pion_synch_rate_per_gev == 0.0)
    assert np.all(output.muon_synch_rate_per_gev == 0.0)
    assert np.all(output.pion_ic_rate_per_gev == 0.0)
    assert np.all(output.muon_ic_rate_per_gev == 0.0)


def main() -> None:
    test_secondary_radiation_smoke()
    test_secondary_radiation_zero_inputs()


if __name__ == "__main__":
    main()
