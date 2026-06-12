from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.asgard_config import FitConfig, ReverseShockConfig
from asgard_core.asgard_physics_utils import ambient_density
from asgard_core.asgard_state import make_query_setup, solve_state_from_setup


def _base_config() -> FitConfig:
    return FitConfig(
        index_y=0,
        index_syn_integr=2,
        num_threads=1,
        num_gam_e=81,
        num_nu=48,
        num_r=96,
        num_theta=16,
        num_tobs=64,
        t_obs_min_log10=2.0,
        t_obs_max_log10=6.0,
        eta_0=120.0,
        e_iso=1.0e53,
        d_ne=1.0,
        a_star=-1.0,
        electron_solver="fullhide_1d",
        reverse=True,
        include_forward_ssc=False,
        reverse_shock=ReverseShockConfig(
            enabled=True,
            delta_t_s=30.0,
            epsilon_e=0.12,
            epsilon_b=3.0e-3,
            p=2.35,
            f_e=1.0,
            include_ssc=False,
            include_cross_zone_ic=False,
        ),
    )


def _run_single_bump_equivalence() -> None:
    radius = np.logspace(15.0, 18.0, 64)
    scalar = _base_config()
    scalar.r_tr = 3.0e16
    scalar.f_jump = 6.0
    scalar.f_wide = 0.08
    arrays = _base_config()
    arrays.jump_r_cm = (scalar.r_tr,)
    arrays.jump_factor = (scalar.f_jump,)
    arrays.jump_width_log10 = (scalar.f_wide,)
    assert np.allclose(ambient_density(radius, scalar), ambient_density(radius, arrays), rtol=0.0, atol=0.0)


def _run_multi_bump_reverse() -> None:
    config = _base_config()
    config.jump_r_cm = (2.5e16, 7.0e16)
    config.jump_factor = (5.0, 8.0)
    config.jump_width_log10 = (0.10, 0.12)
    times = np.logspace(2.0, 6.0, 18)
    freqs = np.array([1.0e10, 1.0e14], dtype=float)
    setup = make_query_setup(config, times, freqs)
    state = solve_state_from_setup(config, setup)
    assert state.reverse_emission is not None
    assert state.reverse_emission.secondary_rs is not None
    secondary = state.reverse_emission.secondary_rs
    assert np.all(np.isfinite(state.reverse_emission.l_syn_spec))
    assert np.any(secondary.luminosity_syn > 0.0)
    for values in (
        secondary.gamma_contact,
        secondary.pressure_3,
        secondary.gamma_43,
        secondary.dissipated_energy_density,
        secondary.magnetic_field_g,
        secondary.nu_m,
        secondary.nu_c,
        secondary.nu_a,
    ):
        assert np.all(np.isfinite(values))
    active = secondary.magnetic_field_g > 0.0
    assert np.all(secondary.gamma_contact[active] > 1.0)
    assert np.all(secondary.gamma_43[active] >= 1.0)
    assert np.all(secondary.pressure_3[active] > 0.0)
    assert np.all(secondary.dissipated_energy_density[active] > 0.0)


def main() -> None:
    _run_single_bump_equivalence()
    _run_multi_bump_reverse()
    print("multi-density-reverse-shock-smoke-ok")


if __name__ == "__main__":
    main()
