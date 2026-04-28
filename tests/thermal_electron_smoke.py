from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import run_fit
from asgard_core.asgard_config import FitConfig


def _config(*, thermal_electrons: bool, electron_solver: str = "fullhide_1d") -> FitConfig:
    return FitConfig(
        num_threads=1,
        electron_solver=electron_solver,
        thermal_electrons=thermal_electrons,
        index_y=0,
        num_gam_e=81,
        num_nu=61,
        num_r=48,
        num_theta=16,
        num_tobs=20,
        plot_lc=False,
        show_plots=False,
    )


def main() -> None:
    result_off = run_fit(_config(thermal_electrons=False))
    result_on = run_fit(_config(thermal_electrons=True))
    if not np.all(np.isfinite(result_off.bands_flux)):
        raise AssertionError("thermal off bands_flux contains non-finite values.")
    if not np.all(np.isfinite(result_on.bands_flux)):
        raise AssertionError("thermal on bands_flux contains non-finite values.")
    if np.max(np.abs(result_on.bands_flux - result_off.bands_flux)) <= 0.0:
        raise AssertionError("thermal electron switch did not change the forward spectrum.")

    try:
        run_fit(_config(thermal_electrons=True, electron_solver="weno5_1d"))
    except NotImplementedError:
        pass
    else:
        raise AssertionError("thermal electrons should reject unsupported electron solvers.")

    print("thermal_electron_smoke: ok")


if __name__ == "__main__":
    main()
