from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_coupling import build_coupled_shock_geometry, build_cross_zone_seed_fields
from asgard_presets import build_reverse_demo_config
from asgard_setup import build_simulation_setup
from asgard_runtime import solve_dynamics, solve_electron, solve_reverse_shock_emission
from src import Radiation


def main() -> None:
    config = build_reverse_demo_config()
    config.num_gam_e = 101
    config.num_nu = 101
    config.num_r = 120
    config.num_theta = 180

    setup = build_simulation_setup(config)
    dynamics = solve_dynamics(setup.boundary, config)
    electron = solve_electron(setup.boundary, dynamics, setup.seed_frequency_hz, config)
    reverse_emission = solve_reverse_shock_emission(setup.boundary, dynamics, setup.seed_frequency_hz, config)

    geometry = build_coupled_shock_geometry(dynamics, config)
    seed_fs_to_rs, seed_rs_to_fs = build_cross_zone_seed_fields(electron.seed_syn, reverse_emission.seed_syn, geometry)

    l_cic_fs, _ = Radiation.ssc_spec(
        dynamics.radius,
        electron.gam_e,
        electron.d_n_gam_e,
        setup.seed_frequency_hz,
        seed_rs_to_fs,
        config.num_threads,
    )
    l_cic_rs, _ = Radiation.ssc_spec(
        dynamics.radius,
        dynamics.reverse_shock.gam_e,
        dynamics.reverse_shock.d_n_gam_e,
        setup.seed_frequency_hz,
        seed_fs_to_rs,
        config.num_threads,
    )

    assert np.all(np.isfinite(geometry.fs_width_cm))
    assert np.all(np.isfinite(geometry.rs_width_cm))
    assert np.all(geometry.fs_width_cm > 0.0)
    assert np.all(geometry.rs_width_cm >= 0.0)
    assert np.all(geometry.center_delay_s >= 0.0)
    assert np.all(geometry.proper_time_s[1:] >= geometry.proper_time_s[:-1])
    assert np.all(np.isfinite(l_cic_fs))
    assert np.all(np.isfinite(l_cic_rs))
    assert float(np.max(l_cic_fs)) > 0.0
    assert float(np.max(l_cic_rs)) > 0.0

    print("PASS: coupled IC check succeeded.")


if __name__ == "__main__":
    main()
