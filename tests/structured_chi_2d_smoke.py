from __future__ import annotations

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


def case_structured_chi_ring_precomputed_matches_single_ring_chi():
    boundary, r_tobs, radius, source, tau, chi_radius, chi_gamma, chi_weight, seed, obs, tobs, phi_centers = (
        _structured_kernel_inputs()
    )
    num_theta = r_tobs.shape[1]
    num_phi_half = phi_centers.size
    dtheta = boundary[8] / float(num_theta)
    ring_boundary = np.array(boundary, copy=True)
    ring_boundary[8] = dtheta

    precomputed = Interpolation.sed_interpolation_chi_structured_axisym_ring_precomputed(
        boundary,
        r_tobs[:, 0],
        radius[:, 0],
        source[:, :, :, 0],
        tau[:, :, :, 0],
        chi_radius[:, :, 0],
        chi_gamma[:, :, 0],
        chi_weight[:, :, 0],
        seed,
        obs,
        tobs,
        0.0,
        dtheta,
        2 * num_phi_half,
    )
    reference = Interpolation.sed_interpolation_chi(
        ring_boundary,
        r_tobs[:, 0],
        radius[:, 0],
        source[:, :, :, 0],
        tau[:, :, :, 0],
        chi_radius[:, :, 0],
        chi_gamma[:, :, 0],
        chi_weight[:, :, 0],
        seed,
        obs,
        tobs,
        1,
        num_phi_half,
        1,
    )
    assert np.any(precomputed > 0.0)
    np.testing.assert_allclose(precomputed, reference, rtol=1.0e-12, atol=1.0e-12)
    return {"max_flux": float(np.max(precomputed)), "num_phi_full": int(2 * num_phi_half)}


def case_chi_electron_cached_top_hat_path():
    boundary, r_tobs, radius, _source, _tau, chi_radius, chi_gamma, chi_weight, seed, obs, tobs, phi_centers = (
        _structured_kernel_inputs()
    )
    boundary = np.array(boundary, copy=True)
    boundary[12] = 1.0e27
    num_chi = chi_radius.shape[0]
    num_r = radius.shape[0]
    gam_e = np.geomspace(2.0, 2.0e4, 9)
    gam_shape = (gam_e[:, None, None] / gam_e[0]) ** -2.2
    chi_shape = 1.0 + 0.1 * np.arange(num_chi, dtype=float)[None, :, None]
    radius_shape = (radius[None, None, :, 0] / radius[0, 0]) ** -0.5
    dne_chi = 1.0e44 * gam_shape * chi_shape * radius_shape
    b_chi = 0.08 * (radius[None, :, 0] / radius[0, 0]) ** -0.6
    b_chi = b_chi * (1.0 + 0.04 * np.arange(num_chi, dtype=float)[:, None])

    flux = Interpolation.sed_interpolation_chi_electron_cached(
        boundary,
        r_tobs[:, 0],
        radius[:, 0],
        dne_chi,
        b_chi,
        chi_radius[:, :, 0],
        chi_gamma[:, :, 0],
        chi_weight[:, :, 0],
        gam_e,
        seed,
        obs,
        tobs,
        1,
        phi_centers.size,
        1,
    )
    assert np.any(flux > 0.0)
    assert np.all(np.isfinite(flux))
    return {"max_flux": float(np.max(flux)), "num_r": int(num_r)}


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
        "ring_precomputed_parity": case_structured_chi_ring_precomputed_matches_single_ring_chi(),
        "electron_cached_top_hat": case_chi_electron_cached_top_hat_path(),
        "public_path": case_public_structured_fullhide_2d_chi_path(),
        "parallel_path": case_public_structured_parallel_matches_serial(),
    }
    print(
        "structured-chi-2d-smoke-ok "
        f"ring_max_flux={results['ring_precomputed_parity']['max_flux']:.6g} "
        f"electron_cached_max_flux={results['electron_cached_top_hat']['max_flux']:.6g} "
        f"patch_count={results['public_path']['patch_count']}"
    )


if __name__ == "__main__":
    main()
