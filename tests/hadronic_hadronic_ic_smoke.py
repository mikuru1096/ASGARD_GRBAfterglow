from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.hadronic_processes import (  # noqa: E402
    AM3_C_CGS,
    AM3_MASS_ELECTRON_GEV,
    AM3_MASS_MUON_GEV,
    AM3_MASS_PION_CHARGED_GEV,
    AM3_MASS_PROTON_GEV,
    AM3_SIGMA_T_CGS,
    initialize_hadronic_inverse_compton_kernel,
    solve_hadronic_inverse_compton,
)


def test_hadronic_ic_smoke_and_regression() -> None:
    dln = 0.25
    hadron_energy_gev = np.exp(np.arange(7, dtype=float) * dln)
    photon_energy_gev = np.exp((np.arange(11, dtype=float) - 4.0) * dln)
    ind_min = 4

    photons = np.linspace(1.0e-4, 3.0e-4, photon_energy_gev.size)
    protons = np.linspace(2.0e-10, 8.0e-10, hadron_energy_gev.size)
    pion_plus = np.linspace(1.0e-10, 2.0e-10, hadron_energy_gev.size)
    pion_minus = np.linspace(0.5e-10, 1.0e-10, hadron_energy_gev.size)
    muon_minus_left = np.linspace(0.8e-11, 2.2e-11, hadron_energy_gev.size)
    muon_minus_right = np.linspace(0.9e-11, 2.0e-11, hadron_energy_gev.size)
    muon_plus_left = np.linspace(1.2e-11, 1.8e-11, hadron_energy_gev.size)
    muon_plus_right = np.linspace(1.1e-11, 1.9e-11, hadron_energy_gev.size)

    kernel = initialize_hadronic_inverse_compton_kernel(hadron_energy_gev, photon_energy_gev, ind_min)
    output = solve_hadronic_inverse_compton(
        hadron_energy_gev=hadron_energy_gev,
        photon_energy_gev=photon_energy_gev,
        photons_on_had_grid_per_gev=photons,
        protons_per_gev=protons,
        pion_plus_per_gev=pion_plus,
        pion_minus_per_gev=pion_minus,
        muon_minus_left_per_gev=muon_minus_left,
        muon_minus_right_per_gev=muon_minus_right,
        muon_plus_left_per_gev=muon_plus_left,
        muon_plus_right_per_gev=muon_plus_right,
        ind_min_energy_pho_hadgrid=ind_min,
        kernel=kernel,
    )

    assert output.epsilon_p_ic.shape == photon_energy_gev.shape
    assert output.epsilon_pi_ic.shape == photon_energy_gev.shape
    assert output.epsilon_mu_ic.shape == photon_energy_gev.shape
    assert np.all(np.isfinite(output.epsilon_p_ic))
    assert np.all(np.isfinite(output.epsilon_pi_ic))
    assert np.all(np.isfinite(output.epsilon_mu_ic))
    assert np.all(output.epsilon_p_ic >= 0.0)
    assert np.all(output.epsilon_pi_ic >= 0.0)
    assert np.all(output.epsilon_mu_ic >= 0.0)

    coeff_p = AM3_C_CGS * AM3_SIGMA_T_CGS / (AM3_MASS_PROTON_GEV / AM3_MASS_ELECTRON_GEV) ** 2
    coeff_pi = AM3_C_CGS * AM3_SIGMA_T_CGS / (AM3_MASS_PION_CHARGED_GEV / AM3_MASS_ELECTRON_GEV) ** 2
    coeff_mu = AM3_C_CGS * AM3_SIGMA_T_CGS / (AM3_MASS_MUON_GEV / AM3_MASS_ELECTRON_GEV) ** 2

    expected_p = _manual_channel(photons, protons, kernel.delta_e_p, kernel.jmax_p, kernel.dln_energy, coeff_p)
    expected_pi = _manual_channel(
        photons,
        pion_plus + pion_minus,
        kernel.delta_e_pi,
        kernel.jmax_pi,
        kernel.dln_energy,
        coeff_pi,
    )
    expected_mu = _manual_channel(
        photons,
        muon_minus_left + muon_minus_right + muon_plus_left + muon_plus_right,
        kernel.delta_e_mu,
        kernel.jmax_mu,
        kernel.dln_energy,
        coeff_mu,
    )

    assert np.allclose(output.epsilon_p_ic, expected_p, rtol=1.0e-14, atol=0.0)
    assert np.allclose(output.epsilon_pi_ic, expected_pi, rtol=1.0e-14, atol=0.0)
    assert np.allclose(output.epsilon_mu_ic, expected_mu, rtol=1.0e-14, atol=0.0)


def test_hadronic_ic_zero_distributions() -> None:
    dln = 0.3
    hadron_energy_gev = np.exp(np.arange(5, dtype=float) * dln)
    photon_energy_gev = np.exp((np.arange(9, dtype=float) - 2.0) * dln)
    zeros_had = np.zeros_like(hadron_energy_gev)
    zeros_ph = np.zeros_like(photon_energy_gev)

    output = solve_hadronic_inverse_compton(
        hadron_energy_gev=hadron_energy_gev,
        photon_energy_gev=photon_energy_gev,
        photons_on_had_grid_per_gev=zeros_ph,
        protons_per_gev=zeros_had,
        pion_plus_per_gev=zeros_had,
        pion_minus_per_gev=zeros_had,
        muon_minus_left_per_gev=zeros_had,
        muon_minus_right_per_gev=zeros_had,
        muon_plus_left_per_gev=zeros_had,
        muon_plus_right_per_gev=zeros_had,
        ind_min_energy_pho_hadgrid=2,
    )

    assert np.all(output.epsilon_p_ic == 0.0)
    assert np.all(output.epsilon_pi_ic == 0.0)
    assert np.all(output.epsilon_mu_ic == 0.0)


def _manual_channel(
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
    test_hadronic_ic_smoke_and_regression()
    test_hadronic_ic_zero_distributions()


if __name__ == "__main__":
    main()
