from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.hadronic_cascade import compute_iterative_pair_cascade, compute_time_dependent_pair_cascade_sequence
from asgard_core.hadronic_pair_production import ELECTRON_MASS_GEV
from src import constants


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


def test_time_dependent_pair_cascade_sequence() -> None:
    photon_energy_gev, photon_density_per_gev, electron_energy_gev = _aligned_inputs()
    primary = np.stack((photon_density_per_gev, 0.8 * photon_density_per_gev), axis=1)
    frequency_hz = photon_energy_gev * constants.para_gev2hz

    out = compute_time_dependent_pair_cascade_sequence(
        photon_energy_gev=photon_energy_gev,
        primary_photon_density_per_gev=primary,
        electron_energy_gev=electron_energy_gev,
        frequency_hz=frequency_hz,
        radius_cm=np.array([1.0e16, 1.2e16], dtype=float),
        gamma_bulk=np.array([100.0, 95.0], dtype=float),
        observer_time_s=np.array([100.0, 130.0], dtype=float),
        b_field_g=np.array([1.0e10, 8.0e9], dtype=float),
        num_threads=1,
        index_syn_integr=2,
        substeps_per_shell=2,
    )

    assert out.pair_density_per_gamma.shape == (electron_energy_gev.size, 2)
    assert np.all(np.isfinite(out.pair_syn_luminosity_hz))
    assert np.all(np.isfinite(out.pair_syn_seed_per_hz))
    assert np.all(out.tau_pair_path >= 0.0)
    assert float(np.max(out.pair_density_per_gamma[:, 1])) > 0.0
    assert float(np.max(out.pair_syn_luminosity_hz)) > 0.0


def main() -> None:
    test_pair_cascade_smoke()
    test_pair_cascade_zero_magnetic_field_boundary()
    test_time_dependent_pair_cascade_sequence()


if __name__ == "__main__":
    main()
