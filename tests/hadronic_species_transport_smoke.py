from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.hadronic_species_transport import (
    ChargedMuonDistribution,
    ChargedPionDistribution,
    HadronicSpeciesSources,
    HadronicSpeciesState,
    NeutronDistribution,
    advance_species_transport_explicit,
    spherical_divergence_rate,
)


def main() -> None:
    gamma = np.logspace(1.0, 5.0, 96)

    neutron0 = 1.0e-10 * gamma ** (-2.1)
    pion0 = 4.0e-9 * np.exp(-((np.log10(gamma) - 2.8) / 0.55) ** 2)
    muon0 = 3.0e-9 * np.exp(-((np.log10(gamma) - 2.6) / 0.55) ** 2)

    state0 = HadronicSpeciesState(
        neutron=NeutronDistribution(gamma=gamma, density_per_gamma=neutron0),
        charged_pion=ChargedPionDistribution(
            gamma=gamma,
            plus_density_per_gamma=0.55 * pion0,
            minus_density_per_gamma=0.45 * pion0,
        ),
        charged_muon=ChargedMuonDistribution(
            gamma=gamma,
            minus_left_density_per_gamma=0.30 * muon0,
            minus_right_density_per_gamma=0.20 * muon0,
            plus_left_density_per_gamma=0.25 * muon0,
            plus_right_density_per_gamma=0.25 * muon0,
        ),
    )
    sources0 = HadronicSpeciesSources(
        neutron_per_gamma_s=np.zeros_like(gamma),
        charged_pion_plus_per_gamma_s=np.zeros_like(gamma),
        charged_pion_minus_per_gamma_s=np.zeros_like(gamma),
        charged_muon_minus_left_per_gamma_s=np.zeros_like(gamma),
        charged_muon_minus_right_per_gamma_s=np.zeros_like(gamma),
        charged_muon_plus_left_per_gamma_s=np.zeros_like(gamma),
        charged_muon_plus_right_per_gamma_s=np.zeros_like(gamma),
    )

    dt_decay_s = 5.0e-7
    div_rate = spherical_divergence_rate(radius_cm=1.0e15, expansion_speed_cm_s=5.0e8)
    state1 = advance_species_transport_explicit(
        state=state0,
        sources=sources0,
        dt_s=dt_decay_s,
        b_field_g=30.0,
        divergence_rate_s_inv=div_rate,
    )

    n0 = float(np.trapezoid(state0.neutron.density_per_gamma, gamma))
    pi0 = float(np.trapezoid(state0.charged_pion.plus_density_per_gamma + state0.charged_pion.minus_density_per_gamma, gamma))
    mu0 = float(
        np.trapezoid(
            state0.charged_muon.minus_left_density_per_gamma
            + state0.charged_muon.minus_right_density_per_gamma
            + state0.charged_muon.plus_left_density_per_gamma
            + state0.charged_muon.plus_right_density_per_gamma,
            gamma,
        )
    )
    n1 = float(np.trapezoid(state1.neutron.density_per_gamma, gamma))
    pi1 = float(np.trapezoid(state1.charged_pion.plus_density_per_gamma + state1.charged_pion.minus_density_per_gamma, gamma))
    mu1 = float(
        np.trapezoid(
            state1.charged_muon.minus_left_density_per_gamma
            + state1.charged_muon.minus_right_density_per_gamma
            + state1.charged_muon.plus_left_density_per_gamma
            + state1.charged_muon.plus_right_density_per_gamma,
            gamma,
        )
    )

    assert np.all(np.isfinite(state1.neutron.density_per_gamma))
    assert np.all(np.isfinite(state1.charged_pion.plus_density_per_gamma))
    assert np.all(np.isfinite(state1.charged_pion.minus_density_per_gamma))
    assert np.all(np.isfinite(state1.charged_muon.minus_left_density_per_gamma))
    assert np.all(np.isfinite(state1.charged_muon.minus_right_density_per_gamma))
    assert np.all(np.isfinite(state1.charged_muon.plus_left_density_per_gamma))
    assert np.all(np.isfinite(state1.charged_muon.plus_right_density_per_gamma))
    assert np.all(state1.neutron.density_per_gamma >= 0.0)
    assert np.all(state1.charged_pion.plus_density_per_gamma >= 0.0)
    assert np.all(state1.charged_pion.minus_density_per_gamma >= 0.0)
    assert np.all(state1.charged_muon.minus_left_density_per_gamma >= 0.0)
    assert np.all(state1.charged_muon.minus_right_density_per_gamma >= 0.0)
    assert np.all(state1.charged_muon.plus_left_density_per_gamma >= 0.0)
    assert np.all(state1.charged_muon.plus_right_density_per_gamma >= 0.0)

    assert n1 < n0
    assert pi1 < pi0
    assert mu1 < mu0
    assert (1.0 - pi1 / pi0) > (1.0 - n1 / n0)
    assert (1.0 - mu1 / mu0) > (1.0 - n1 / n0)

    inj = np.exp(-((np.log10(gamma) - 3.0) / 0.22) ** 2)
    sources1 = HadronicSpeciesSources(
        neutron_per_gamma_s=np.zeros_like(gamma),
        charged_pion_plus_per_gamma_s=0.55e-3 * inj,
        charged_pion_minus_per_gamma_s=0.45e-3 * inj,
        charged_muon_minus_left_per_gamma_s=2.4e-4 * inj,
        charged_muon_minus_right_per_gamma_s=1.6e-4 * inj,
        charged_muon_plus_left_per_gamma_s=2.0e-4 * inj,
        charged_muon_plus_right_per_gamma_s=2.0e-4 * inj,
    )
    state2 = advance_species_transport_explicit(
        state=state1,
        sources=sources1,
        dt_s=1.0e-8,
        b_field_g=10.0,
        divergence_rate_s_inv=div_rate,
    )
    pi2 = float(np.trapezoid(state2.charged_pion.plus_density_per_gamma + state2.charged_pion.minus_density_per_gamma, gamma))
    mu2 = float(
        np.trapezoid(
            state2.charged_muon.minus_left_density_per_gamma
            + state2.charged_muon.minus_right_density_per_gamma
            + state2.charged_muon.plus_left_density_per_gamma
            + state2.charged_muon.plus_right_density_per_gamma,
            gamma,
        )
    )

    assert pi2 > pi1
    assert mu2 > mu1
    print("hadronic_species_transport_smoke: ok")


if __name__ == "__main__":
    main()
