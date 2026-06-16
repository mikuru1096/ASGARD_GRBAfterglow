from __future__ import annotations

from _repo_path import ensure_repo_root_on_path
import numpy as np


ensure_repo_root_on_path()

from tests.public_api_builders import top_hat_model


def main() -> None:
    model = top_hat_model()

    scalar = model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total
    assert scalar.shape == (1,)
    assert np.all(np.isfinite(scalar))

    times = np.logspace(2.0, 6.0, 8)
    freqs = np.array([1.0e9, 1.0e14, 1.0e18])

    grid = model.flux_density_grid(times, freqs).total
    assert grid.shape == (freqs.size, times.size)
    assert np.all(np.isfinite(grid))

    pairwise = model.flux_density(times[: freqs.size], freqs).total
    assert pairwise.shape == (freqs.size,)
    assert np.all(np.isfinite(pairwise))

    exposure = model.flux_density_exposures(
        times[: freqs.size],
        np.full(freqs.size, 1.0e14),
        0.2 * times[: freqs.size],
    ).total
    assert exposure.shape == (freqs.size,)
    assert np.all(np.isfinite(exposure))


if __name__ == "__main__":
    main()
