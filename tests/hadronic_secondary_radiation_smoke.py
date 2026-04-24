from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.hadronic_secondary_radiation import (  # noqa: E402
    AM3_MASS_ELECTRON_GEV,
    AM3_MASS_MUON_GEV,
    AM3_MASS_PION_CHARGED_GEV,
    AM3_SIGMA_T_CGS,
    AM3_C_CGS,
    SecondarySpeciesDistribution,
    SecondaryTargetPhotonField,
    initialize_secondary_inverse_compton_kernel,
    initialize_secondary_synchrotron_kernel,
    muon_inverse_compton_rate,
    muon_synchrotron_rate,
    pion_inverse_compton_rate,
    pion_synchrotron_rate,
    solve_secondary_radiation_spectrum,
)


def test_secondary_radiation_smoke_and_regression() -> None:
    dln = 0.2
    hadron_energy_gev = np.exp(np.arange(8, dtype=float) * dln)
    photon_energy_gev = np.exp((np.arange(12, dtype=float) - 4.0) * dln)
    magnetic_field_g = 25.0

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

    synch_kernel = initialize_secondary_synchrotron_kernel(hadron_energy_gev, photon_energy_gev, magnetic_field_g)
    ic_kernel = initialize_secondary_inverse_compton_kernel(
        hadron_energy_gev,
        photon_energy_gev,
        target.ind_min_energy_pho_hadgrid,
    )
    pion_syn = pion_synchrotron_rate(hadron_energy_gev, species, photon_energy_gev, magnetic_field_g, kernel=synch_kernel)
    muon_syn = muon_synchrotron_rate(hadron_energy_gev, species, photon_energy_gev, magnetic_field_g, kernel=synch_kernel)
    pion_ic = pion_inverse_compton_rate(hadron_energy_gev, species, target, kernel=ic_kernel)
    muon_ic = muon_inverse_compton_rate(hadron_energy_gev, species, target, kernel=ic_kernel)
    total = solve_secondary_radiation_spectrum(
        hadron_energy_gev=hadron_energy_gev,
        species=species,
        target=target,
        magnetic_field_g=magnetic_field_g,
        synch_kernel=synch_kernel,
        ic_kernel=ic_kernel,
    )

    for arr in (pion_syn, muon_syn, pion_ic, muon_ic):
        assert arr.shape == photon_energy_gev.shape
        assert np.all(np.isfinite(arr))
        assert np.all(arr >= 0.0)

    assert np.allclose(total.pion_synch_rate_per_gev, pion_syn, rtol=1.0e-14, atol=0.0)
    assert np.allclose(total.muon_synch_rate_per_gev, muon_syn, rtol=1.0e-14, atol=0.0)
    assert np.allclose(total.pion_ic_rate_per_gev, pion_ic, rtol=1.0e-14, atol=0.0)
    assert np.allclose(total.muon_ic_rate_per_gev, muon_ic, rtol=1.0e-14, atol=0.0)

    pion_total = species.pion_plus_per_gev + species.pion_minus_per_gev
    muon_total = (
        species.muon_minus_left_per_gev
        + species.muon_minus_right_per_gev
        + species.muon_plus_left_per_gev
        + species.muon_plus_right_per_gev
    )
    expected_pion_syn = synch_kernel.dln_energy * (synch_kernel.kernel_pion @ pion_total)
    expected_muon_syn = synch_kernel.dln_energy * (synch_kernel.kernel_muon @ muon_total)
    expected_pion_ic = _manual_ic_channel(
        target.photons_on_had_grid_per_gev,
        pion_total,
        ic_kernel.delta_e_pi,
        ic_kernel.jmax_pi,
        ic_kernel.dln_energy,
        AM3_C_CGS * AM3_SIGMA_T_CGS / (AM3_MASS_PION_CHARGED_GEV / AM3_MASS_ELECTRON_GEV) ** 2,
    )
    expected_muon_ic = _manual_ic_channel(
        target.photons_on_had_grid_per_gev,
        muon_total,
        ic_kernel.delta_e_mu,
        ic_kernel.jmax_mu,
        ic_kernel.dln_energy,
        AM3_C_CGS * AM3_SIGMA_T_CGS / (AM3_MASS_MUON_GEV / AM3_MASS_ELECTRON_GEV) ** 2,
    )
    assert np.allclose(pion_syn, expected_pion_syn, rtol=1.0e-14, atol=0.0)
    assert np.allclose(muon_syn, expected_muon_syn, rtol=1.0e-14, atol=0.0)
    assert np.allclose(pion_ic, expected_pion_ic, rtol=1.0e-14, atol=0.0)
    assert np.allclose(muon_ic, expected_muon_ic, rtol=1.0e-14, atol=0.0)


def test_secondary_radiation_zero_inputs() -> None:
    dln = 0.25
    hadron_energy_gev = np.exp(np.arange(6, dtype=float) * dln)
    photon_energy_gev = np.exp((np.arange(10, dtype=float) - 3.0) * dln)
    species = SecondarySpeciesDistribution(
        pion_plus_per_gev=np.zeros_like(hadron_energy_gev),
        pion_minus_per_gev=np.zeros_like(hadron_energy_gev),
        muon_minus_left_per_gev=np.zeros_like(hadron_energy_gev),
        muon_minus_right_per_gev=np.zeros_like(hadron_energy_gev),
        muon_plus_left_per_gev=np.zeros_like(hadron_energy_gev),
        muon_plus_right_per_gev=np.zeros_like(hadron_energy_gev),
    )
    target = SecondaryTargetPhotonField(
        photon_energy_gev=photon_energy_gev,
        photons_on_had_grid_per_gev=np.zeros_like(photon_energy_gev),
        ind_min_energy_pho_hadgrid=2,
    )
    output = solve_secondary_radiation_spectrum(
        hadron_energy_gev=hadron_energy_gev,
        species=species,
        target=target,
        magnetic_field_g=10.0,
    )
    assert np.all(output.pion_synch_rate_per_gev == 0.0)
    assert np.all(output.muon_synch_rate_per_gev == 0.0)
    assert np.all(output.pion_ic_rate_per_gev == 0.0)
    assert np.all(output.muon_ic_rate_per_gev == 0.0)


def _manual_ic_channel(
    photons: np.ndarray,
    particles: np.ndarray,
    delta_e: np.ndarray,
    jmax: np.ndarray,
    dln: float,
    coeff: float,
) -> np.ndarray:
    out = np.zeros_like(photons)
    for j in range(photons.size):
        z = 0.0
        for i in range(particles.size):
            if j < int(delta_e[i]) or j > int(jmax[i]):
                continue
            z += photons[j - int(delta_e[i])] * particles[i]
        out[j] = z * dln * coeff
    return out


def main() -> None:
    test_secondary_radiation_smoke_and_regression()
    test_secondary_radiation_zero_inputs()


if __name__ == "__main__":
    main()
