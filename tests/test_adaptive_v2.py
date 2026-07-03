"""Test new flux-based adaptive theta integration."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_state import make_query_cfg, make_query_setup, solve_state_from_setup
from asgard_core.structured_jet_kernel import (
    _project_structured_chi_ring_state,
    _solve_project_structured_chi_ring_flux_adaptive,
)


class Jet:
    kind = "gaussian"
    theta_max = 0.1
    def energy_iso(self, phi, theta):
        return 1e53 * np.exp(-0.5 * (theta / 0.04) ** 2)
    def gamma0(self, phi, theta):
        return 100.0 * np.exp(-0.5 * (theta / 0.08) ** 2)


class Obs:
    theta_obs = 0.3; phi_obs = 0.0


class Rad:
    ssc = False


class Setups:
    structured_num_theta = 96; structured_num_phi = 8
    structured_parallel_mode = "outer"; num_threads = 1
    structured_outer_threads = None; structured_inner_threads = None
    geometry_kernel = "chi_eats_2d"; electron_solver = "fullhide_2d"
    rvs_shock = False; hadronic_enabled = False
    include_cross_zone_ic = False; pair_cascade_iterations = 1
    structured_adaptive_rtol = 0.01; structured_adaptive_max_depth = 4


class Model:
    jet = Jet(); observer = Obs(); fwd_rad = Rad(); setups = Setups()


def make_patch_config(m, **kw):
    c = RuntimeConfig()
    c.num_r = 50; c.num_gam_e = 96; c.num_nu = 64; c.num_threads = 1
    c.electron_solver = "fullhide_2d"; c.geometry_kernel = "chi_eats_2d"
    c.z = 0.4; c.opening_angle_jet = kw.get("opening_angle_jet", 0.1)
    c.theta_v = kw.get("theta_v", 0.0); c.e_iso = kw.get("e_iso", 1e53)
    c.eta_0 = kw.get("gamma0", 100.0); c.epsilon_e = 0.1
    c.epsilon_b = 1e-3; c.p = 2.5; c.f_e = 0.1; c.index_y = 2
    return c


def main():
    model = Model()
    theta_edges = np.linspace(0.0, 0.1, 97)
    solve_times = np.logspace(2, 8, 50)
    obs_times = np.array([1e5, 1e6])
    obs_freqs = np.array([3e9, 1e14])

    # Uniform reference
    print("=== Uniform reference (96 rings) ===")
    t0 = time.perf_counter()
    ref = np.zeros((2, 2))
    for i in range(96):
        tc = 0.5 * (theta_edges[i] + theta_edges[i + 1])
        g0 = model.jet.gamma0(0.0, tc)
        if g0 < 1.4:
            continue
        eiso = model.jet.energy_iso(0.0, tc)
        c = make_patch_config(model, opening_angle_jet=theta_edges[i + 1] - theta_edges[i],
                              e_iso=eiso, gamma0=g0)
        qc = make_query_cfg(c, solve_times)
        qc.num_r = max(int(qc.num_r), int(solve_times.size))
        stp = make_query_setup(qc, solve_times, obs_freqs)
        state = solve_state_from_setup(qc, stp, assemble_observer=False)
        ref += _project_structured_chi_ring_state(state, theta_edges[i], theta_edges[i + 1],
                                                   0.3, 8, obs_times, obs_freqs)
    ref_t = time.perf_counter() - t0
    ref_s = float(np.sum(np.abs(ref)))
    print("Ref flux=%.6e time=%.1fs" % (ref_s, ref_t))

    # Adaptive
    print("\n=== Adaptive (rtol varying) ===")
    active = np.ones(96, dtype=np.int32)
    e_iso_arr = np.array([model.jet.energy_iso(0.0, 0.5 * (theta_edges[i] + theta_edges[i + 1]))
                          for i in range(96)])
    gamma0_arr = np.array([model.jet.gamma0(0.0, 0.5 * (theta_edges[i] + theta_edges[i + 1]))
                           for i in range(96)])
    tc_arr = 0.5 * (theta_edges[:-1] + theta_edges[1:])

    for rtol in [0.1, 0.05, 0.02, 0.01, 0.005]:
        t0 = time.perf_counter()
        try:
            flux, _ = _solve_project_structured_chi_ring_flux_adaptive(
                model, make_patch_config, tc_arr, e_iso_arr, gamma0_arr, active,
                theta_edges, solve_times, obs_times, obs_freqs,
                adaptive_rtol=rtol, adaptive_max_depth=4,
            )
            dt = time.perf_counter() - t0
            fsum = float(np.sum(np.abs(flux)))
            err = abs(fsum - ref_s) / max(ref_s, 1e-100)
            spd = ref_t / max(dt, 1e-6)
            print("rtol=%.3f flux=%.6e err=%.1e time=%.1fs spd=%.1fx" % (rtol, fsum, err, dt, spd))
        except Exception as e:
            print("rtol=%.3f FAILED: %s" % (rtol, str(e)[:80]))


if __name__ == "__main__":
    main()
