from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import run_fit
from asgard_core.asgard_config import FitConfig, ReverseShockConfig
from asgard_core.asgard_state import make_query_setup, solve_state_from_setup


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
    config = _config(index_y)
    result = run_fit(config)
    assert result.rs_nu_m is not None
    assert result.rs_nu_c is not None
    assert result.rs_nu_a is not None
    assert np.all(np.isfinite(result.bands_flux))
    assert np.all(np.isfinite(result.rs_nu_m))
    assert np.all(np.isfinite(result.rs_nu_c))
    assert np.all(np.isfinite(result.rs_nu_a))

    setup = make_query_setup(config, np.logspace(2.0, 5.0, 6), np.array([1.0e9, 1.0e14]))
    state = solve_state_from_setup(config, setup)
    rs = state.dynamics.reverse_shock
    assert rs is not None
    assert np.all(np.isfinite(rs.magnetic_field_g))
    assert np.all(np.isfinite(rs.internal_energy_erg))
    assert np.all(np.isfinite(rs.comoving_volume_cm3))
    assert np.all(np.isfinite(rs.gamma34))
    assert np.all(rs.internal_energy_erg > 0.0)
    assert np.all(rs.comoving_volume_cm3 > 0.0)
    assert np.all(rs.magnetic_field_g >= 0.0)
    assert np.all(rs.gamma34 >= 1.0)
    post = np.asarray(state.dynamics.radius, dtype=float) >= float(rs.r_cross)
    if np.count_nonzero(post) > 2:
        thermal_energy = rs.internal_energy_erg[post] / rs.swept_mass_g[post]
        assert np.all(np.diff(thermal_energy) <= 0.0)


def main() -> None:
    for index_y in (0, 3):
        _run(index_y)
    print("reverse-shock-smoke-ok")


if __name__ == "__main__":
    main()
