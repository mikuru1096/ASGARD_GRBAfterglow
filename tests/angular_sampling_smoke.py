from __future__ import annotations

from pathlib import Path
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _make_model(
    *,
    backend: str = "python_patch",
    sampling: str = "uniform",
    theta_obs: float = 0.08,
    gamma0: float = 120.0,
    theta_c: float = 0.08,
    theta_max: float = 0.24,
):
    from ASGARD import GaussianJet, ISM, Model, Observer, Radiation, Setups

    return Model(
        GaussianJet(E_iso=1.0e52, Gamma0=gamma0, theta_c=theta_c, theta_max=theta_max),
        ISM(1.0),
        Observer(1.0e26, 0.1, theta_obs),
        Radiation(0.1, 1.0e-3, 2.3, ssc=False),
        setups=Setups(
            structured_backend=backend,
            patch_sampling=sampling,
            patch_theta=20,
            patch_phi=15,
            num_gam_e=21,
            num_nu=17,
            num_r=12,
            num_theta=8,
            num_phi=8,
            num_tobs=8,
            num_threads=1,
            electron_solver="fullhide_1d",
            patch_sampling_pilot_theta=8,
            patch_sampling_num_times=3,
            patch_sampling_beaming_factor=3.0,
            patch_sampling_beaming_resolution=8.0,
        ),
    )


def _assert_uniform_grid_matches_legacy() -> None:
    from asgard_core.angular_sampling import build_patch_grid

    model = _make_model(sampling="uniform")
    grid = build_patch_grid(model)
    theta_edges = np.linspace(0.0, model.jet.theta_max, model.setups.patch_theta + 1)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, model.setups.patch_phi + 1)
    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    domega = (np.cos(theta_edges[:-1]) - np.cos(theta_edges[1:]))[:, None] * np.diff(phi_edges)[None, :]
    np.testing.assert_allclose(grid.theta_edges, theta_edges, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(grid.phi_edges, phi_edges, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(grid.theta_centers, theta_centers, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(grid.phi_centers, phi_centers, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(grid.domega, domega, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(grid.patch_half_angle, np.sqrt(np.maximum(domega, 1.0e-12) / np.pi))


def _assert_adaptive_grid() -> None:
    from asgard_core.angular_sampling import build_patch_grid

    model = _make_model(sampling="dominant_region_ioka_v1")
    grid = build_patch_grid(model)
    if grid.domega.shape != (20, 15):
        raise AssertionError(f"unexpected adaptive domega shape: {grid.domega.shape}")
    expected_solid_angle = 2.0 * np.pi * (1.0 - np.cos(model.jet.theta_max))
    np.testing.assert_allclose(float(np.sum(grid.domega)), expected_solid_angle, rtol=1.0e-12)
    if not np.all(np.diff(grid.theta_edges) > 0.0):
        raise AssertionError("adaptive theta edges are not strictly increasing")
    np.testing.assert_allclose(grid.phi_centers[0], 0.0, rtol=0.0, atol=1.0e-14)
    np.testing.assert_allclose(grid.phi_centers[-1], np.pi, rtol=0.0, atol=1.0e-14)

    axis_model = _make_model(sampling="dominant_region_ioka_v1", theta_obs=0.0)
    axis_grid = build_patch_grid(axis_model)
    if axis_grid.domega.shape != (20, 15):
        raise AssertionError(f"unexpected on-axis adaptive domega shape: {axis_grid.domega.shape}")
    np.testing.assert_allclose(axis_grid.phi_centers[0], 0.0, rtol=0.0, atol=1.0e-14)
    np.testing.assert_allclose(axis_grid.phi_centers[-1], np.pi, rtol=0.0, atol=1.0e-14)
    np.testing.assert_allclose(float(np.sum(axis_grid.domega)), expected_solid_angle, rtol=1.0e-12)


def _assert_error_contracts() -> None:
    from asgard_core.angular_sampling import build_patch_grid

    try:
        build_patch_grid(_make_model(sampling="not_a_sampler"))
    except ValueError:
        pass
    else:
        raise AssertionError("unknown patch_sampling did not raise ValueError")

    try:
        _make_model(backend="fortran_1d", sampling="dominant_region_ioka_v1").flux_density_grid(
            np.array([1.0e4]),
            np.array([1.0e14]),
        )
    except NotImplementedError:
        pass
    else:
        raise AssertionError("fortran_1d accepted dominant_region_ioka_v1")


def _assert_generation_speed() -> None:
    from asgard_core.angular_sampling import build_patch_grid

    model = _make_model(sampling="dominant_region_ioka_v1")
    samples = []
    for _ in range(20):
        start = time.perf_counter()
        build_patch_grid(model)
        samples.append(time.perf_counter() - start)
    median = float(np.median(samples))
    if median >= 0.05:
        raise AssertionError(f"dominant_region_ioka_v1 grid generation median={median:.4f}s")
    print(f"dominant_region_ioka_v1 median generation time: {median:.4f}s")


def _assert_time_sampler_beaming_phi_resolution() -> None:
    from asgard_core.angular_sampling import _pilot_gamma_theta_time, _sampling_times, build_patch_grid

    model = _make_model(
        sampling="dominant_region_ioka_time_v1",
        theta_obs=0.04,
        gamma0=600.0,
        theta_c=0.04,
        theta_max=0.12,
    )
    model.setups.structured_outer_threads = 2
    times = np.geomspace(1.0, 1.0e3, 6)
    grid = build_patch_grid(model, times)
    if grid.theta_centers.size <= model.setups.patch_theta:
        raise AssertionError("time-dependent sampler did not increase theta resolution for Gamma0=600")
    if grid.phi_centers.size <= model.setups.patch_phi:
        raise AssertionError("time-dependent sampler did not increase phi resolution for Gamma0=600")
    theta_probe = np.linspace(0.0, model.jet.theta_max, max(4 * model.setups.patch_theta, 256))
    gamma_probe = _pilot_gamma_theta_time(model, theta_probe, _sampling_times(model, times))
    max_gamma_sin_theta = float(np.max(gamma_probe * np.sin(theta_probe)[:, None]))
    max_allowed_step = model.setups.patch_sampling_beaming_factor / (
        model.setups.patch_sampling_beaming_resolution * max_gamma_sin_theta
    )
    if np.pi / float(grid.phi_centers.size - 1) > max_allowed_step:
        raise AssertionError("time-dependent sampler phi spacing violates the beaming-resolution criterion")


def main() -> None:
    _assert_uniform_grid_matches_legacy()
    _assert_adaptive_grid()
    _assert_error_contracts()
    _assert_generation_speed()
    _assert_time_sampler_beaming_phi_resolution()
    print("angular_sampling_smoke passed")


if __name__ == "__main__":
    main()
