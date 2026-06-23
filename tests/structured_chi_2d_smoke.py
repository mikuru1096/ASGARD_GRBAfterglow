from __future__ import annotations

import time

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core import tabulated_angular_jet
from src import Interpolation
from tests.public_api_builders import numerics, radiation, solver_options, top_hat_model


def _structured_kernel_inputs():
    num_theta = 3
    num_phi = 4
    num_chi = 3
    num_r = 7
    seed = np.geomspace(1.0e8, 1.0e18, 12)
    obs = np.array([1.0e10, 1.0e14], dtype=float)
    boundary = np.zeros(30, dtype=float)
    boundary[7] = 0.08
    boundary[8] = 0.12
    boundary[9] = 0.19
    radius_base = np.geomspace(1.0e16, 5.0e17, num_r)
    theta_scale = 1.0 + 0.04 * np.arange(num_theta, dtype=float)
    radius = radius_base[:, None] * theta_scale[None, :]
    gamma = np.geomspace(70.0, 8.0, num_r)[:, None] * (1.0 - 0.06 * np.arange(num_theta)[None, :])
    r_tobs = (1.0 + boundary[7]) * radius / (4.0 * gamma**2) / 2.99792458e10
    chi_offset = np.linspace(0.0, 1.0, num_chi)
    chi_radius = radius[None, :, :] * (1.0 - 0.025 * chi_offset[:, None, None])
    chi_gamma = gamma[None, :, :] * (1.0 - 0.035 * chi_offset[:, None, None])
    chi_weight = np.array([0.55, 0.30, 0.15], dtype=float)[:, None, None] * np.ones((1, num_r, num_theta))
    source = (
        (seed[:, None, None, None] / 1.0e12) ** 0.15
        * (radius[None, None, :, :] / radius_base[0]) ** -0.7
        * (1.0 + 0.2 * chi_offset[None, :, None, None])
    )
    tau = 0.04 * (seed[:, None, None, None] / seed[0]) ** -0.1 * (1.0 + chi_offset[None, :, None, None])
    tau = tau * np.ones((1, 1, num_r, num_theta))
    theta_centers = (np.arange(num_theta, dtype=float) + 0.5) * boundary[8] / float(num_theta)
    phi_centers = (np.arange(num_phi, dtype=float) + 0.5) * (2.0 * np.pi / float(num_phi))
    arrival = r_tobs + (1.0 + boundary[7]) * (radius - radius * np.cos(theta_centers)[None, :]) / 2.99792458e10
    tobs = np.geomspace(float(np.min(arrival[1:-1, :])), float(np.max(arrival[1:-1, :])), 9)
    return boundary, r_tobs, radius, source, tau, chi_radius, chi_gamma, chi_weight, seed, obs, tobs, phi_centers


def case_structured_chi_batch_matches_direct_sum():
    boundary, r_tobs, radius, source, tau, chi_radius, chi_gamma, chi_weight, seed, obs, tobs, phi_centers = (
        _structured_kernel_inputs()
    )
    num_theta = r_tobs.shape[1]
    num_phi = phi_centers.size
    dtheta = boundary[8] / float(num_theta)
    dphi = 2.0 * np.pi / float(num_phi)

    start = time.perf_counter()
    batch = Interpolation.sed_interpolation_chi_structured_axisym(
        boundary,
        r_tobs,
        radius,
        source,
        tau,
        chi_radius,
        chi_gamma,
        chi_weight,
        seed,
        obs,
        tobs,
        num_phi,
        1,
    )
    batch_seconds = time.perf_counter() - start

    direct = np.zeros_like(batch)
    start = time.perf_counter()
    for i_theta in range(num_theta):
        theta_lo = dtheta * i_theta
        theta_hi = dtheta * (i_theta + 1)
        theta_center = dtheta * (i_theta + 0.5)
        domega = (np.cos(theta_lo) - np.cos(theta_hi)) * dphi
        for phi_center in phi_centers:
            mu = np.cos(boundary[9]) * np.cos(theta_center) + np.sin(boundary[9]) * np.sin(theta_center) * np.cos(phi_center)
            patch_boundary = np.array(boundary, copy=True)
            patch_boundary[9] = np.arccos(mu)
            direct += Interpolation.sed_interpolation_chi_surface_element(
                patch_boundary,
                r_tobs[:, i_theta],
                radius[:, i_theta],
                source[:, :, :, i_theta],
                tau[:, :, :, i_theta],
                chi_radius[:, :, i_theta],
                chi_gamma[:, :, i_theta],
                chi_weight[:, :, i_theta],
                seed,
                obs,
                tobs,
                domega,
            )
    direct_seconds = time.perf_counter() - start

    assert np.any(batch > 0.0)
    np.testing.assert_allclose(batch, direct, rtol=1.0e-12, atol=1.0e-12)
    return {
        "batch_calls": 1,
        "direct_calls": int(num_theta * num_phi),
        "batch_seconds": batch_seconds,
        "direct_seconds": direct_seconds,
    }


def case_structured_chi_electron_cached_ring_sum():
    boundary, r_tobs, radius, _source, _tau, chi_radius, chi_gamma, chi_weight, seed, obs, tobs, phi_centers = (
        _structured_kernel_inputs()
    )
    boundary = np.array(boundary, copy=True)
    boundary[12] = 1.0e27
    num_theta = r_tobs.shape[1]
    num_phi = phi_centers.size
    num_chi = chi_radius.shape[0]
    num_r = radius.shape[0]
    gam_e = np.geomspace(2.0, 2.0e4, 9)
    gam_shape = (gam_e[:, None, None, None] / gam_e[0]) ** -2.2
    chi_shape = 1.0 + 0.1 * np.arange(num_chi, dtype=float)[None, :, None, None]
    radius_shape = (radius[None, None, :, :] / radius[0, 0]) ** -0.5
    theta_shape = 1.0 + 0.03 * np.arange(num_theta, dtype=float)[None, None, None, :]
    dne_chi = 1.0e44 * gam_shape * chi_shape * radius_shape * theta_shape
    b_chi = 0.08 * (radius[None, :, :] / radius[0, 0]) ** -0.6
    b_chi = b_chi * (1.0 + 0.04 * np.arange(num_chi, dtype=float)[:, None, None])
    dtheta = boundary[8] / float(num_theta)

    batch = Interpolation.sed_interpolation_chi_structured_axisym_electron_cached(
        boundary,
        r_tobs,
        radius,
        dne_chi,
        b_chi,
        chi_radius,
        chi_gamma,
        chi_weight,
        gam_e,
        seed,
        obs,
        tobs,
        num_phi,
        1,
    )
    ring_sum = np.zeros_like(batch)
    for i_theta in range(num_theta):
        ring_sum += Interpolation.sed_interpolation_chi_structured_axisym_electron_cached_ring(
            boundary,
            r_tobs[:, i_theta],
            radius[:, i_theta],
            dne_chi[:, :, :, i_theta],
            b_chi[:, :, i_theta],
            chi_radius[:, :, i_theta],
            chi_gamma[:, :, i_theta],
            chi_weight[:, :, i_theta],
            gam_e,
            seed,
            obs,
            tobs,
            dtheta * i_theta,
            dtheta * (i_theta + 1),
            num_phi,
        )
    assert np.any(batch > 0.0)
    np.testing.assert_allclose(ring_sum, batch, rtol=1.0e-12, atol=1.0e-12)
    return {"max_flux": float(np.max(batch)), "num_r": int(num_r)}


def case_public_structured_fullhide_2d_chi_path():
    model = _public_structured_chi_model()
    times = np.array([1.0e3, 1.0e4], dtype=float)
    freq = np.array([1.0e9], dtype=float)
    flux = model.flux_density_grid(times, freq)
    assert flux.total.shape == (1, 2)
    assert np.all(np.isfinite(flux.total))
    assert np.max(flux.total) > 0.0
    details = model._last_details
    assert details is not None
    assert len(details.patches) == 8

    total_matrix = model._total_matrix(times, freq, projection_kind="lightcurve")
    assert total_matrix.shape == flux.total.shape
    assert np.all(np.isfinite(total_matrix))
    assert np.max(total_matrix) > 0.0
    return {"max_flux": float(np.max(flux.total)), "patch_count": len(details.patches)}


def case_public_structured_parallel_matches_serial():
    times = np.array([1.0e3, 1.0e4], dtype=float)
    freq = np.array([1.0e9], dtype=float)
    serial = _public_structured_chi_model(num_threads=1)
    parallel = _public_structured_chi_model(num_threads=2, structured_outer_threads=2, structured_inner_threads=1)
    serial_flux = serial.flux_density_grid(times, freq).total
    parallel_flux = parallel.flux_density_grid(times, freq).total
    np.testing.assert_allclose(parallel_flux, serial_flux, rtol=1.0e-12, atol=1.0e-12)
    assert parallel._last_details is not None
    assert len(parallel._last_details.patches) == 8
    return {"parallel_max_flux": float(np.max(parallel_flux)), "patch_count": len(parallel._last_details.patches)}



def _public_structured_chi_model(
    num_threads=1,
    structured_outer_threads=None,
    structured_inner_threads=None,
    patch_sampling="uniform",
):
    model = top_hat_model(
        jet=tabulated_angular_jet(
            theta_rad=np.array([0.0, 0.04, 0.10, 0.16], dtype=float),
            energy_iso_erg=np.array([1.0e51, 4.0e50, 8.0e49, 1.0e45], dtype=float),
            lorentz_factor=np.array([90.0, 18.0, 1.4, 1.2], dtype=float),
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        fwd_rad=radiation(include_ssc=False, epsilon_B=1.0e-3, magnetic_decay_alpha_t=0.4),
        rvs_rad=None,
        numerics=numerics(
            structured_num_theta=3,
            structured_num_phi=4,
            num_radius=8,
            eats_num_theta=4,
            eats_num_phi=2,
            num_observer_time=5,
            num_electron_gamma=8,
            num_photon_frequency=8,
            downstream_num_chi=3,
            num_threads=num_threads,
        ),
        solver_options=solver_options(
            electron_solver="fullhide_2d",
            geometry_projection="chi_eats_2d",
            structured_backend="fortran_1d",
            ssc_cooling_mode="nakar_y_thomson",
            synchrotron_integration="fixed_grid",
            structured_parallel_mode="outer",
            structured_outer_threads=structured_outer_threads,
            structured_inner_threads=structured_inner_threads,
            patch_sampling=patch_sampling,
        ),
    )
    model.observer.viewing_angle_rad = 0.2
    return model


def main() -> None:
    results = {
        "kernel_parity": case_structured_chi_batch_matches_direct_sum(),
        "electron_ring_parity": case_structured_chi_electron_cached_ring_sum(),
        "public_path": case_public_structured_fullhide_2d_chi_path(),
        "parallel_path": case_public_structured_parallel_matches_serial(),
    }
    print(
        "structured-chi-2d-smoke-ok "
        f"batch_calls={results['kernel_parity']['batch_calls']} "
        f"direct_calls={results['kernel_parity']['direct_calls']} "
        f"batch_seconds={results['kernel_parity']['batch_seconds']:.6g} "
        f"direct_seconds={results['kernel_parity']['direct_seconds']:.6g} "
        f"patch_count={results['public_path']['patch_count']}"
    )


if __name__ == "__main__":
    main()
