from __future__ import annotations

import numpy as np

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet


def _model(jump_r: tuple[float, ...], jump_width: float) -> Model:
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
            jump_width_log10=(jump_width,) * len(jump_r),
        ),
    )


def _reacceleration_energy(model: Model) -> tuple[np.ndarray, np.ndarray]:
    details = model.details(1.0e-2, 1.0e6)
    assert details.rev is not None
    seed = np.asarray(details.rev.secondary_rs_branch_reacceleration_seed_energy_erg, dtype=float)
    reaccelerated = np.asarray(details.rev.secondary_rs_branch_reaccelerated_energy_erg, dtype=float)
    assert seed.shape == reaccelerated.shape
    assert np.all(np.isfinite(seed))
    assert np.all(np.isfinite(reaccelerated))
    return seed, reaccelerated


def _branch_compression(model: Model) -> np.ndarray:
    details = model.details(1.0e-2, 1.0e6)
    assert details.rev is not None
    compression = np.asarray(details.rev.secondary_rs_branch_compression, dtype=float)
    gamma43 = np.asarray(details.rev.secondary_rs_branch_gamma_43, dtype=float)
    assert np.all(np.isfinite(compression))
    assert np.all(np.isfinite(gamma43))
    return compression, gamma43


def main() -> None:
    separated_seed, separated_reaccelerated = _reacceleration_energy(
        _model((1.0e15, 4.0e15, 1.6e16), 0.08)
    )
    separated_seed_by_branch = np.sum(separated_seed, axis=1)
    separated_reaccelerated_by_branch = np.sum(separated_reaccelerated, axis=1)
    assert separated_seed_by_branch[0] == 0.0
    assert separated_reaccelerated_by_branch[0] == 0.0
    assert np.all(separated_seed_by_branch[1:] > 0.0)
    assert np.all(separated_reaccelerated_by_branch[1:] > separated_seed_by_branch[1:])

    overlap_model = _model((1.0e15, 1.25e15, 1.56e15), 0.18)
    overlap_seed, overlap_reaccelerated = _reacceleration_energy(overlap_model)
    overlap_compression, overlap_gamma43 = _branch_compression(overlap_model)
    seed_by_branch = np.sum(overlap_seed, axis=1)
    reaccelerated_by_branch = np.sum(overlap_reaccelerated, axis=1)
    active_reaccel = overlap_seed > 0.0
    assert seed_by_branch[0] == 0.0
    assert reaccelerated_by_branch[0] == 0.0
    assert np.all(seed_by_branch[1:] > 0.0)
    assert np.all(reaccelerated_by_branch[1:] > seed_by_branch[1:])
    assert np.all(overlap_compression[active_reaccel] > 0.0)
    assert np.all(overlap_gamma43[active_reaccel] > 1.0)
    print("secondary-reverse-reacceleration-smoke-ok")


if __name__ == "__main__":
    main()
