from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.hadronic_cascade import compute_iterative_pair_cascade
from asgard_core.hadronic_pair_production import ELECTRON_MASS_GEV


def _aligned_inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dln = 0.1
    n_ph = 84
    n_e = 56
    ind_min = 28
    x_ph = np.exp((-ind_min + np.arange(n_ph, dtype=float)) * dln)
    gm_e = np.exp(np.arange(n_e, dtype=float) * dln)
    photon_energy_gev = x_ph * ELECTRON_MASS_GEV
    electron_energy_gev = gm_e * ELECTRON_MASS_GEV
    photon_density_per_gev = 3.0e9 * np.power(x_ph, -2.1)
    return photon_energy_gev, photon_density_per_gev, electron_energy_gev


def test_pair_cascade_smoke() -> None:
    photon_energy_gev, photon_density_per_gev, electron_energy_gev = _aligned_inputs()

    out = compute_iterative_pair_cascade(
        photon_energy_gev=photon_energy_gev,
        photon_density_per_gev=photon_density_per_gev,
        electron_energy_gev=electron_energy_gev,
        radius_cm=1.0e16,
        gamma_bulk=100.0,
        b_field_g=1.0e10,
        max_iterations=2,
    )

    assert np.all(np.isfinite(out.pair_syn_luminosity_hz))
    assert np.all(out.pair_syn_luminosity_hz > 0.0)
    assert np.all(out.tau_pair_path >= 0.0)
    assert float(np.max(out.pair_syn_luminosity_hz)) > 0.0
    assert float(np.max(out.tau_pair_path)) > 0.0
    assert out.absorbed_power_gev_per_cm3_s > 0.0


def test_pair_cascade_zero_magnetic_field_boundary() -> None:
    photon_energy_gev, photon_density_per_gev, electron_energy_gev = _aligned_inputs()

    out = compute_iterative_pair_cascade(
        photon_energy_gev=photon_energy_gev,
        photon_density_per_gev=photon_density_per_gev,
        electron_energy_gev=electron_energy_gev,
        radius_cm=1.0e16,
        gamma_bulk=100.0,
        b_field_g=0.0,
        max_iterations=2,
    )

    assert np.all(out.pair_syn_luminosity_hz == 0.0)
    assert float(np.max(out.tau_pair_path)) > 0.0
    assert out.absorbed_power_gev_per_cm3_s > 0.0


def main() -> None:
    test_pair_cascade_smoke()
    test_pair_cascade_zero_magnetic_field_boundary()


if __name__ == "__main__":
    main()
