from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.hadronic_pair_production import ELECTRON_MASS_GEV, solve_pair_production


def _aligned_grids() -> tuple[np.ndarray, np.ndarray]:
    dln = 0.1
    n_ph = 84
    n_e = 56
    ind_min = 28
    x_ph = np.exp((-ind_min + np.arange(n_ph, dtype=float)) * dln)
    gm_e = np.exp(np.arange(n_e, dtype=float) * dln)
    return x_ph * ELECTRON_MASS_GEV, gm_e * ELECTRON_MASS_GEV


def test_pair_production_smoke() -> None:
    photon_energy_gev, electron_energy_gev = _aligned_grids()
    x_ph = photon_energy_gev / ELECTRON_MASS_GEV
    photon_density_per_gev = 3.0e9 * np.power(x_ph, -2.1)

    out = solve_pair_production(
        photon_energy_gev=photon_energy_gev,
        photon_density_per_gev=photon_density_per_gev,
        electron_energy_gev=electron_energy_gev,
    )

    assert np.all(np.isfinite(out.photon_loss_rate))
    assert np.all(np.isfinite(out.pair_injection_rate_per_gev_per_species))
    assert np.all(np.isfinite(out.pair_injection_rate_per_gev_total))
    assert float(np.max(out.photon_loss_rate)) > 0.0
    assert float(np.max(out.pair_injection_rate_per_gev_total)) > 0.0


def test_pair_production_energy_closure() -> None:
    photon_energy_gev, electron_energy_gev = _aligned_grids()
    x_ph = photon_energy_gev / ELECTRON_MASS_GEV
    photon_density_per_gev = 4.0e9 * np.power(x_ph, -2.25)

    out = solve_pair_production(
        photon_energy_gev=photon_energy_gev,
        photon_density_per_gev=photon_density_per_gev,
        electron_energy_gev=electron_energy_gev,
    )

    assert out.absorbed_power_gev_per_cm3_s > 0.0
    assert out.injected_power_gev_per_cm3_s > 0.0
    assert np.isclose(
        out.injected_power_gev_per_cm3_s,
        out.absorbed_power_gev_per_cm3_s,
        rtol=5.0e-7,
        atol=0.0,
    )


def main() -> None:
    test_pair_production_smoke()
    test_pair_production_energy_closure()


if __name__ == "__main__":
    main()
