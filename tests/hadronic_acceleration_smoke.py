from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.hadronic_acceleration import (
    ERG_PER_GEV,
    InjectionConfig,
    build_species_injection_operators,
    estimate_max_gamma,
    external_photon_cooling_timescale_s,
    proton_acceleration_timescale_s,
    species_properties,
)


def test_proton_acceleration_timescale() -> None:
    gamma = np.logspace(1.0, 7.0, 64)
    t_acc = proton_acceleration_timescale_s(gamma=gamma, b_field_g=30.0, eta_acc=3.0)
    ratio = t_acc / gamma
    assert np.all(np.isfinite(t_acc))
    assert np.all(np.diff(t_acc) > 0.0)
    assert np.allclose(ratio, ratio[0], rtol=1.0e-12, atol=0.0)


def test_species_aware_injection_operators() -> None:
    gamma = np.logspace(1.0, 9.0, 8192)
    l_inj = 3.5e40
    cfg_p = InjectionConfig(
        species="proton",
        luminosity_erg_s=l_inj,
        spectral_index=2.2,
        gamma_min=30.0,
        gamma_max=2.0e8,
        gamma_cut=8.0e7,
    )
    cfg_mu = InjectionConfig(
        species="muon_plus",
        luminosity_erg_s=l_inj,
        spectral_index=2.2,
        gamma_min=30.0,
        gamma_max=2.0e8,
        gamma_cut=8.0e7,
    )
    q_map = build_species_injection_operators(gamma, [cfg_p, cfg_mu])
    q_p = q_map["proton"]
    q_mu = q_map["muon_plus"]

    m_p_erg = species_properties("proton").mass_gev * ERG_PER_GEV
    m_mu_erg = species_properties("muon_plus").mass_gev * ERG_PER_GEV
    l_p = np.trapezoid(q_p * gamma * m_p_erg, gamma)
    l_mu = np.trapezoid(q_mu * gamma * m_mu_erg, gamma)
    assert np.isclose(l_p, l_inj, rtol=2.0e-3, atol=0.0)
    assert np.isclose(l_mu, l_inj, rtol=2.0e-3, atol=0.0)
    assert not np.allclose(q_p, q_mu)


def test_max_energy_estimator_with_external_cooling() -> None:
    def toy_external_rate(species: str, gamma: np.ndarray, context: dict[str, float]) -> np.ndarray:
        _ = species
        return context["k0"] * gamma**2

    gamma_sample = np.logspace(2.0, 8.0, 32)
    t_ext = external_photon_cooling_timescale_s(
        species="proton",
        gamma=gamma_sample,
        cooling_rate=toy_external_rate,
        context={"k0": 1.0e-8},
    )
    assert np.all(np.isfinite(t_ext))
    assert np.all(np.diff(t_ext) < 0.0)

    out = estimate_max_gamma(
        species="proton",
        b_field_g=120.0,
        radius_cm=6.0e15,
        gamma_bulk=200.0,
        eta_acc=5.0,
        cooling_rate=toy_external_rate,
        context={"k0": 1.0e-8},
    )
    assert out.gamma_ext is not None
    assert out.gamma_max > 1.0
    assert out.gamma_max <= out.gamma_dyn
    assert out.gamma_max <= out.gamma_syn


def main() -> None:
    test_proton_acceleration_timescale()
    test_species_aware_injection_operators()
    test_max_energy_estimator_with_external_cooling()


if __name__ == "__main__":
    main()
