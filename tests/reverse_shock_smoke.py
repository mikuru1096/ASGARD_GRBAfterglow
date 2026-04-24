from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import run_fit
from asgard_core.asgard_config import FitConfig, ReverseShockConfig


def _config(index_y: int) -> FitConfig:
    return FitConfig(
        index_y=index_y,
        index_syn_integr=2,
        num_threads=1,
        num_gam_e=61,
        num_nu=41,
        num_r=40,
        num_theta=24,
        num_tobs=24,
        reverse=True,
        plot_lc=False,
        show_plots=False,
        reverse_shock=ReverseShockConfig(
            enabled=True,
            delta_t_s=10.0,
            epsilon_e=0.1,
            epsilon_b=1.0e-2,
            p=2.4,
            f_e=1.0,
        ),
    )


def _run(index_y: int) -> None:
    result = run_fit(_config(index_y))
    assert result.rs_nu_m is not None
    assert result.rs_nu_c is not None
    assert result.rs_nu_a is not None
    assert np.all(np.isfinite(result.bands_flux))
    assert np.all(np.isfinite(result.rs_nu_m))
    assert np.all(np.isfinite(result.rs_nu_c))
    assert np.all(np.isfinite(result.rs_nu_a))


def main() -> None:
    for index_y in (0, 3):
        _run(index_y)
    print("reverse-shock-smoke-ok")


if __name__ == "__main__":
    main()
