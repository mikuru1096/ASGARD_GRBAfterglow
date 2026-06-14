from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet


def _model() -> Model:
    return Model(
        jet=TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=80.0, duration=20.0),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-4, p=2.3, xi_N=0.2, ssc=False),
        rvs_rad=Radiation(eps_e=0.15, eps_B=3.0e-3, p=2.35, xi_N=1.0, ssc=False),
        setups=Setups(
            rvs_shock=True,
            reverse_delta_t_s=20.0,
            fwd_ssc=False,
            rvs_ssc=False,
            ssc_cooling=False,
            num_threads=1,
            num_r=48,
            num_theta=8,
            num_tobs=24,
            num_gam_e=41,
            num_nu=31,
            observer_time_min_s=1.0e2,
            observer_time_max_s=1.0e5,
            jump_r_cm=(3.0e16,),
            jump_factor=(6.0,),
            jump_width_log10=(0.10,),
            electron_solver="fullhide_1d",
        ),
    )


def main() -> None:
    times = np.logspace(2.0, 5.0, 8)
    freqs = np.array([1.0e10], dtype=float)
    model = _model()
    baseline = model.flux_density_grid(times, freqs)
    adaptive = model.flux_density_grid_adaptive(times, freqs)

    assert baseline.total.shape == (freqs.size, times.size)
    assert adaptive.total.shape[0] == freqs.size
    assert adaptive.total.shape[1] == adaptive.time_s.size
    assert adaptive.time_s.size > times.size
    assert np.all(np.diff(adaptive.time_s) > 0.0)
    assert np.isclose(adaptive.time_s[0], times[0])
    assert np.isclose(adaptive.time_s[-1], times[-1])
    assert np.all(np.isfinite(adaptive.total))
    for time_s in times:
        assert np.any(np.isclose(adaptive.time_s, time_s, rtol=1.0e-13, atol=0.0))
    assert np.any(adaptive.rev.sync > 0.0)
    details = model.details(times[0], times[-1])
    assert details.rev is not None
    assert details.rev.secondary_rs_event_active is not None
    assert np.any(details.rev.secondary_rs_event_active)
    assert details.rev.secondary_rs_start_radius is not None
    assert details.rev.secondary_rs_shock_end_radius is not None
    assert details.rev.secondary_rs_B is not None
    active = details.rev.secondary_rs_event_active
    assert np.all(details.rev.secondary_rs_shock_end_radius[active] >= details.rev.secondary_rs_start_radius[active])
    latest_end = float(np.max(details.rev.secondary_rs_shock_end_radius[active]))
    cooling_tail = (details.rev.radius > latest_end) & (details.rev.secondary_rs_B > 0.0)
    assert np.any(cooling_tail)
    tail_axis_times = details.rev.t_obs[cooling_tail]
    covered_tail = tail_axis_times[(tail_axis_times >= adaptive.time_s[0]) & (tail_axis_times <= adaptive.time_s[-1])]
    assert covered_tail.size > 0
    for time_s in covered_tail[: min(3, covered_tail.size)]:
        assert np.any(np.isclose(adaptive.time_s, time_s, rtol=1.0e-13, atol=0.0))
    print("adaptive-observer-time-grid-smoke-ok")


if __name__ == "__main__":
    main()
