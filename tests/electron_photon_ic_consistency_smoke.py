from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from src import Radiation, constants
from src.Electron import electron_forward_fullhide_1d as Electron1D


def _simpson_log_integral(values: np.ndarray, grid: np.ndarray) -> float:
    weights = np.ones(grid.size)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    return float((np.log(grid[1]) - np.log(grid[0])) / 3.0 * np.dot(weights, values))


def _ic_loss_power(gamma: np.ndarray, density: np.ndarray, frequency: np.ndarray, seed: np.ndarray) -> float:
    dot_gamma_over_gamma = Electron1D.fs_electron_ic_cooling_loss_shell(gamma, frequency, seed, 1)
    return _simpson_log_integral(
        density * dot_gamma_over_gamma * gamma**2 * constants.para_m_energy,
        gamma,
    )


def _ic_source_power(radius_cm: float, gamma: np.ndarray, density: np.ndarray, frequency: np.ndarray, seed: np.ndarray) -> float:
    luminosity, _ = Radiation.ssc_spec(
        np.array([radius_cm], dtype=float),
        gamma,
        density[:, None],
        frequency,
        seed[:, None],
        1,
    )
    return _simpson_log_integral(luminosity[:, 0], frequency)


def main() -> None:
    radius_cm = 3.0e16
    gamma = np.logspace(1.0, 5.0, 65)
    frequency = np.logspace(9.0, 23.0, 65)
    density = 1.0e42 * gamma ** -2.2
    seed = 1.0e-18 * (frequency / 1.0e14) ** -0.4

    loss = _ic_loss_power(gamma, density, frequency, seed)
    source = _ic_source_power(radius_cm, gamma, density, frequency, seed)
    loss_seed2 = _ic_loss_power(gamma, density, frequency, 2.0 * seed)
    source_seed2 = _ic_source_power(radius_cm, gamma, density, frequency, 2.0 * seed)
    loss_e2 = _ic_loss_power(gamma, 2.0 * density, frequency, seed)
    source_e2 = _ic_source_power(radius_cm, gamma, 2.0 * density, frequency, seed)

    assert loss > 0.0 and source > 0.0
    np.testing.assert_allclose(loss_seed2 / loss, 2.0, rtol=1.0e-12)
    np.testing.assert_allclose(source_seed2 / source, 2.0, rtol=1.0e-12)
    np.testing.assert_allclose(loss_e2 / loss, 2.0, rtol=1.0e-12)
    np.testing.assert_allclose(source_e2 / source, 2.0, rtol=1.0e-12)
    np.testing.assert_allclose(loss, source, rtol=5.0e-12)
    print("electron_photon_ic_consistency_smoke: ok")


if __name__ == "__main__":
    main()
