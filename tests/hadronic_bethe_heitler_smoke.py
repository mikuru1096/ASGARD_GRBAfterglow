from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.hadronic_bethe_heitler import solve_bethe_heitler
from src import Radiation as RadiationKernel


def test_direct_bh_kernel() -> None:
    proton_energy_gev = np.logspace(-1, 2, 9)
    proton_density_per_gev = 1.0e-10 * proton_energy_gev ** (-2.0)
    photon_energy_gev = np.logspace(-9, -3, 11)
    photon_density_per_gev = 1.0e8 * photon_energy_gev ** (-1.5)
    electron_energy_gev = np.logspace(-6, 0, 13)

    output = solve_bethe_heitler(
        proton_energy_gev,
        proton_density_per_gev,
        photon_energy_gev,
        photon_density_per_gev,
        electron_energy_gev,
    )
    assert np.all(np.isfinite(output.pair_rate_per_gev))
    assert np.all(np.isfinite(output.proton_loss_rate))
    assert float(output.pair_rate_per_gev.max()) > 0.0
    assert float((-output.proton_loss_rate).max()) > 0.0


def test_annihilation_tau_extra() -> None:
    nu_hz = np.logspace(10.0, 20.0, 8)
    radius_cm = np.logspace(15.0, 16.0, 5)
    gamma_bulk = np.full(5, 100.0)
    seed = np.full((8, 5), 1.0e-30)
    tau_extra = np.zeros((8, 5))
    tau_extra[3:, :] = 0.5

    absorption_no_extra = RadiationKernel.annihilation(
        gamma_bulk,
        radius_cm,
        nu_hz,
        seed,
        seed,
        np.zeros_like(tau_extra),
        1,
    )
    absorption_extra = RadiationKernel.annihilation(
        gamma_bulk,
        radius_cm,
        nu_hz,
        seed,
        seed,
        tau_extra,
        1,
    )
    assert np.all(absorption_extra <= absorption_no_extra)


def main() -> None:
    test_direct_bh_kernel()
    test_annihilation_tau_extra()


if __name__ == "__main__":
    main()
