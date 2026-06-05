from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.hadronic_am3_solver import HUMMER_PROCESS_GROUP_LABELS, solve_hummer_2010_response_processes


def main() -> None:
    radius_cm = np.array([1.0e16], dtype=float)
    gam_p = np.logspace(4.5, 5.3, 48)
    d_n_gam_p = 3.0e28 * (gam_p / gam_p[0]) ** (-2.2)
    d_n_gam_p = d_n_gam_p[:, None]

    v_seed_hz = np.logspace(18.3, 18.9, 128)
    seed_target_hz = 1.0e-20 * (v_seed_hz / v_seed_hz[0]) ** (-1.4)
    seed_target_hz = seed_target_hz[:, None]

    output = solve_hummer_2010_response_processes(
        radius_cm,
        gam_p,
        d_n_gam_p,
        v_seed_hz,
        seed_target_hz,
        48,
        include_pg=True,
        include_neutrino=True,
    )

    assert output.l_had_pg_gamma.shape == (v_seed_hz.size, radius_cm.size)
    assert output.neutrino_frequency_hz.shape == (48,)
    assert output.neutrino_luminosity.shape == (48, radius_cm.size)
    assert output.am3_process_power.shape == (len(HUMMER_PROCESS_GROUP_LABELS), gam_p.size, radius_cm.size)

    assert np.all(np.isfinite(output.l_had_pg_gamma))
    assert np.all(np.isfinite(output.neutrino_luminosity))
    assert np.all(np.isfinite(output.am3_process_power))
    assert np.any(output.l_had_pg_gamma > 0.0)
    assert np.any(output.neutrino_luminosity > 0.0)
    assert np.any(output.am3_process_power[0] > 0.0)
    assert np.any(output.am3_process_power[1] > 0.0)
    assert np.any(output.am3_process_power[2] > 0.0)


if __name__ == "__main__":
    main()
