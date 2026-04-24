from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.hadronic_pp import PI0_MASS_GEV, PROTON_MASS_GEV, solve_pp_delta


def _pp_threshold_kinetic_energy_gev() -> float:
    return 2.0 * PI0_MASS_GEV + PI0_MASS_GEV * PI0_MASS_GEV / (2.0 * PROTON_MASS_GEV)


def test_pp_minimal_backend_smoke() -> None:
    proton_energy_gev = np.logspace(np.log10(PROTON_MASS_GEV + 1.0e-3), 6.0, 256)
    proton_density_per_gev = 1.0e-11 * proton_energy_gev ** (-2.2)

    neutrino_energy_gev = np.logspace(-3.0, 4.0, 220)
    pair_energy_gev = neutrino_energy_gev.copy()
    gamma_energy_gev = 2.0 * neutrino_energy_gev

    output = solve_pp_delta(
        proton_energy_gev=proton_energy_gev,
        proton_density_per_gev=proton_density_per_gev,
        target_proton_density_cm3=3.0e7,
        gamma_energy_gev=gamma_energy_gev,
        neutrino_energy_gev=neutrino_energy_gev,
        pair_energy_gev=pair_energy_gev,
    )

    assert output.gamma_rate_per_gev.shape == gamma_energy_gev.shape
    assert output.neutrino_rate_per_gev.shape == neutrino_energy_gev.shape
    assert output.pair_rate_per_gev.shape == pair_energy_gev.shape
    assert output.proton_loss_rate.shape == proton_energy_gev.shape

    assert np.all(np.isfinite(output.gamma_rate_per_gev))
    assert np.all(np.isfinite(output.neutrino_rate_per_gev))
    assert np.all(np.isfinite(output.pair_rate_per_gev))
    assert np.all(np.isfinite(output.proton_loss_rate))

    assert np.all(output.proton_loss_rate <= 0.0)
    assert float((-output.proton_loss_rate).max()) > 0.0

    threshold_parent = PROTON_MASS_GEV + _pp_threshold_kinetic_energy_gev()
    gamma_threshold = 0.5 * 0.17 * threshold_parent
    neutrino_threshold = 0.25 * 0.17 * threshold_parent
    pair_threshold = neutrino_threshold
    assert np.all(output.gamma_rate_per_gev[gamma_energy_gev < gamma_threshold] == 0.0)
    assert np.all(output.neutrino_rate_per_gev[neutrino_energy_gev < neutrino_threshold] == 0.0)
    assert np.all(output.pair_rate_per_gev[pair_energy_gev < pair_threshold] == 0.0)

    active_pair = output.pair_rate_per_gev > 0.0
    ratio_nu_pair = output.neutrino_rate_per_gev[active_pair] / output.pair_rate_per_gev[active_pair]
    assert np.allclose(ratio_nu_pair, 3.0, rtol=1.0e-12, atol=0.0)

    active_gamma = (output.gamma_rate_per_gev > 0.0) & (output.neutrino_rate_per_gev > 0.0)
    ratio_nu_gamma = output.neutrino_rate_per_gev[active_gamma] / output.gamma_rate_per_gev[active_gamma]
    assert np.allclose(ratio_nu_gamma, 6.0, rtol=1.0e-12, atol=0.0)


def main() -> None:
    test_pp_minimal_backend_smoke()


if __name__ == "__main__":
    main()
