from __future__ import annotations

import numpy as np

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet


def _strong_overlap_model() -> Model:
    jump_r = (1.0e15, 1.25e15, 1.56e15)
    return Model(
        jet=TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=100.0, duration=10.0),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=0.3, eps_B=1.0e-5, p=2.3, xi_N=0.1, ssc=False),
        rvs_rad=Radiation(eps_e=0.3, eps_B=1.0e-5, p=2.4, xi_N=0.1, ssc=False),
        setups=Setups(
            rvs_shock=True,
            reverse_delta_t_s=10.0,
            reverse_sigma=0.1,
            fwd_ssc=False,
            rvs_ssc=False,
            ssc_cooling=False,
            num_threads=1,
            num_r=72,
            num_theta=8,
            num_tobs=16,
            num_gam_e=31,
            num_nu=21,
            initial_radius_cm=1.0e13,
            observer_time_min_s=1.0e-2,
            observer_time_max_s=1.0e6,
            electron_solver="fullhide_1d",
            jump_r_cm=jump_r,
            jump_factor=(80.0,) * len(jump_r),
            jump_width_log10=(0.18,) * len(jump_r),
        ),
    )


def main() -> None:
    details = _strong_overlap_model().details(1.0e-2, 1.0e6)
    assert details.rev is not None
    rev = details.rev
    branch_gamma = np.asarray(rev.secondary_rs_branch_gamma_contact, dtype=float)
    branch_gamma_43 = np.asarray(rev.secondary_rs_branch_gamma_43, dtype=float)
    assert branch_gamma.shape[0] == 3
    assert branch_gamma.shape == branch_gamma_43.shape
    assert branch_gamma.shape[1] >= 100
    samples = np.count_nonzero(branch_gamma > 1.0, axis=1)
    assert np.all(samples >= 20)
    assert rev.secondary_rs_start_tobs_axis is not None
    assert rev.secondary_rs_shock_end_tobs_axis is not None
    start_t = np.asarray(rev.secondary_rs_start_tobs_axis, dtype=float)
    end_t = np.asarray(rev.secondary_rs_shock_end_tobs_axis, dtype=float)
    assert start_t[1] < end_t[0]
    assert start_t[2] < end_t[1]
    fs_gamma = np.asarray(details.fwd.Gamma, dtype=float)
    ratio = branch_gamma / fs_gamma.reshape(1, -1)
    active = branch_gamma > 1.0
    assert np.max(ratio[active]) < 1.0
    assert np.all(branch_gamma_43[active] >= 1.0)
    print("secondary-reverse-branch-resolution-smoke-ok")


if __name__ == "__main__":
    main()
