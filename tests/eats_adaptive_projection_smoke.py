from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.asgard_config import RuntimeConfig, SimulationSetup
from asgard_core.asgard_postprocess import observe_flux
from asgard_core.asgard_setup import build_simulation_setup
from src import Interpolation


def _projection_inputs(theta_v: float = 0.025):
    boundary = np.zeros(30, dtype=float)
    boundary[0] = 120.0
    boundary[3] = 1.0e15
    boundary[7] = 0.2
    boundary[8] = 0.08
    boundary[9] = theta_v
    radius = np.geomspace(1.0e16, 7.0e17, 32)
    gamma = np.geomspace(90.0, 8.0, radius.size)
    r_tobs = (1.0 + boundary[7]) * radius / (2.0 * gamma**2) / 2.99792458e10
    seed = np.geomspace(1.0e8, 1.0e20, 64)
    obs = np.array([1.0e10, 1.0e14, 1.0e18], dtype=float)
    theta_center = 0.5 * boundary[8]
    arrival = r_tobs + radius * (1.0 - np.cos(theta_center)) * (1.0 + boundary[7]) / 2.99792458e10
    tobs = np.geomspace(arrival[4], arrival[-5], 40)
    source = (seed[:, None] / 1.0e12) ** 0.25 * (radius[None, :] / radius[0]) ** -0.9
    return boundary, radius, gamma, r_tobs, seed, obs, tobs, source


def case_adaptive_depth_one_matches_uniform_double_theta():
    boundary, radius, gamma, r_tobs, seed, obs, tobs, source = _projection_inputs()
    reference = Interpolation.sed_interpolation(boundary, r_tobs, gamma, radius, source, seed, obs, tobs, 16, 8, 1)
    adaptive = Interpolation.sed_adaptive_theta(
        boundary,
        r_tobs,
        gamma,
        radius,
        source,
        seed,
        obs,
        tobs,
        0.0,
        1,
        8,
        8,
        1,
    )
    np.testing.assert_allclose(adaptive, reference, rtol=1.0e-12, atol=1.0e-12)
    return {"max_reference": float(np.max(reference)), "max_abs_diff": float(np.max(np.abs(adaptive - reference)))}


def case_adaptive_depth_zero_matches_legacy_baseline():
    boundary, radius, gamma, r_tobs, seed, obs, tobs, source = _projection_inputs()
    legacy = Interpolation.sed_interpolation(boundary, r_tobs, gamma, radius, source, seed, obs, tobs, 8, 8, 1)
    adaptive = Interpolation.sed_adaptive_theta(
        boundary,
        r_tobs,
        gamma,
        radius,
        source,
        seed,
        obs,
        tobs,
        1.0e-2,
        0,
        8,
        8,
        1,
    )
    np.testing.assert_allclose(adaptive, legacy, rtol=0.0, atol=0.0)
    return {"max_legacy": float(np.max(legacy))}


def case_postprocess_routes_adaptive_kernel():
    boundary, radius, gamma, r_tobs, seed, obs, tobs, source = _projection_inputs()
    setup = SimulationSetup(
        luminosity_distance_cm=1.0,
        boundary=boundary,
        seed_frequency_hz=seed,
        observer_time_s=tobs,
    )
    config = RuntimeConfig(
        geometry_kernel="sed_adaptive_theta",
        eats_num_theta=8,
        eats_num_phi=8,
        num_threads=1,
        projection_adaptive_rtol=1.0e-30,
        projection_adaptive_max_depth=1,
    )
    routed = observe_flux(setup, r_tobs, gamma, radius, source, obs, config)
    direct = Interpolation.sed_adaptive_theta(
        boundary,
        r_tobs,
        gamma,
        radius,
        source,
        seed,
        obs,
        tobs,
        1.0e-30,
        1,
        8,
        8,
        1,
    )
    np.testing.assert_allclose(routed, direct, rtol=0.0, atol=0.0)
    return {"routed_shape": list(routed.shape)}


def case_adaptive_projection_thread_invariant():
    boundary, radius, gamma, r_tobs, seed, obs, tobs, source = _projection_inputs()
    single = Interpolation.sed_adaptive_theta(
        boundary,
        r_tobs,
        gamma,
        radius,
        source,
        seed,
        obs,
        tobs,
        1.0e-4,
        3,
        10,
        8,
        1,
    )
    threaded = Interpolation.sed_adaptive_theta(
        boundary,
        r_tobs,
        gamma,
        radius,
        source,
        seed,
        obs,
        tobs,
        1.0e-4,
        3,
        10,
        8,
        4,
    )
    atol = 1.0e-12 * float(np.max(np.abs(single)))
    np.testing.assert_allclose(threaded, single, rtol=1.0e-12, atol=atol)
    return {"max_abs_diff": float(np.max(np.abs(threaded - single)))}


def case_setup_validates_adaptive_contract():
    build_simulation_setup(
        RuntimeConfig(
            geometry_kernel="sed_adaptive_theta",
            projection_adaptive_rtol=1.0e-2,
            projection_adaptive_max_depth=2,
        )
    )
    for kwargs, expected in (
        ({"projection_adaptive_rtol": 0.0}, "projection_adaptive_rtol"),
        ({"projection_adaptive_max_depth": -1}, "projection_adaptive_max_depth"),
    ):
        try:
            build_simulation_setup(RuntimeConfig(geometry_kernel="sed_adaptive_theta", **kwargs))
        except ValueError as exc:
            assert expected in str(exc)
        else:
            raise AssertionError(f"accepted invalid adaptive projection config: {kwargs}")
    return {"validated": True}


def main() -> None:
    cases = [
        ("[1/5] adaptive_theta:depth_zero_matches_legacy", case_adaptive_depth_zero_matches_legacy_baseline),
        ("[2/5] adaptive_theta:depth_one_matches_uniform_double_theta", case_adaptive_depth_one_matches_uniform_double_theta),
        ("[3/5] adaptive_theta:postprocess_routes_kernel", case_postprocess_routes_adaptive_kernel),
        ("[4/5] adaptive_theta:thread_invariant", case_adaptive_projection_thread_invariant),
        ("[5/5] adaptive_theta:setup_validation", case_setup_validates_adaptive_contract),
    ]
    for label, fn in cases:
        print(f"  {label} ...", flush=True)
        fn()


if __name__ == "__main__":
    main()
