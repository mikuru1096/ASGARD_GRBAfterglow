"""EATS convergence: vary num_r, measure per-time flux accuracy."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_state import make_query_cfg, make_query_setup, solve_state_from_setup
from asgard_core.structured_jet_kernel import _project_structured_chi_ring_state


def solve_one(nr, n_nu_seed, obs_t, obs_f, tlo, thi, tv, nphi):
    c = RuntimeConfig()
    c.num_r = nr; c.num_gam_e = 96; c.num_nu = n_nu_seed; c.num_threads = 1
    c.electron_solver = "fullhide_2d"; c.geometry_kernel = "chi_eats_2d"
    c.z = 0.4; c.opening_angle_jet = thi - tlo; c.theta_v = tv
    c.e_iso = 1e53; c.eta_0 = 100.0; c.epsilon_e = 0.1
    c.epsilon_b = 1e-3; c.p = 2.5; c.f_e = 0.1; c.index_y = 2
    st = np.logspace(2, 8, nr)
    qc = make_query_cfg(c, st); qc.num_r = max(int(qc.num_r), int(st.size))
    setup = make_query_setup(qc, st, obs_f)
    state = solve_state_from_setup(qc, setup, assemble_observer=False)
    return _project_structured_chi_ring_state(state, tlo, thi, tv, nphi, obs_t, obs_f)


def main():
    obs_times = np.array([1e4, 1e5, 1e6, 1e7])
    obs_freqs = np.array([3e9, 1e14])
    tlo, thi, tv, nphi = 0.01, 0.03, 0.3, 4

    # Reference: fine grid
    t0 = time.perf_counter()
    flux_ref = solve_one(256, 128, obs_times, obs_freqs, tlo, thi, tv, nphi)
    ref_time = time.perf_counter() - t0
    ref_sum = float(np.sum(np.abs(flux_ref)))
    print(f"Reference (nr=256): total_abs={ref_sum:.6e}  time={ref_time:.1f}s")
    for i, t in enumerate(obs_times):
        print(f"  t={t:.0e}: {flux_ref[0,i]:.4e} {flux_ref[1,i]:.4e}")

    # Vary num_r
    print(f"\n{'nr':>4s}  {'err':>8s}  {'time':>6s}  {'spd':>5s}  {'flux_per_t'}")
    best_nt, best_err = 256, 0
    for nr in [192, 128, 96, 64, 48, 40, 36, 32]:
        try:
            t0 = time.perf_counter()
            flux = solve_one(nr, 128, obs_times, obs_freqs, tlo, thi, tv, nphi)
            dt = time.perf_counter() - t0
            # Per-time relative error, take max
            per_t = []
            for i in range(len(obs_times)):
                r = abs(float(flux[0,i]) - float(flux_ref[0,i])) / max(float(flux_ref[0,i]), 1e-100)
                per_t.append(r)
            err = max(per_t)
            spd = ref_time / max(dt, 1e-6)
            ft = " ".join(f"{flux[0,i]:.4e}" for i in range(len(obs_times)))
            print(f"{nr:4d}  {err:8.1e}  {dt:5.1f}s  {spd:4.1f}x  {ft}")
            if err < 0.01 and nr < best_nt:
                best_nt, best_err = nr, err
        except Exception as e:
            print(f"{nr:4d}  FAILED: {str(e)[:60]}")

    print(f"\nBest converged: nr={best_nt} (err={best_err:.1e})")


if __name__ == "__main__":
    main()
