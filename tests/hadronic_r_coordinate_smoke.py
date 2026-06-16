from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.asgard_runtime import _hadronic_shell_comoving_dt_from_radius
from src import constants


def test_first_shell_uses_radius_interval_not_absolute_radius() -> None:
    gamma_bulk = np.array([20.0, 20.0, 20.0], dtype=float)
    beta = np.sqrt(1.0 - 1.0 / gamma_bulk[0] ** 2)
    radius = np.array([1.0e16, 1.04e16, 1.10e16], dtype=float)
    expected = (radius[1] - radius[0]) / (beta * gamma_bulk[0] * constants.para_c)
    got = _hadronic_shell_comoving_dt_from_radius(radius, gamma_bulk, 0)
    np.testing.assert_allclose(got, expected, rtol=1.0e-14, atol=0.0)


def main() -> None:
    test_first_shell_uses_radius_interval_not_absolute_radius()
    print("hadronic_r_coordinate_smoke: ok")


if __name__ == "__main__":
    main()
