from __future__ import annotations

from types import SimpleNamespace

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core import Model
from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_state import _build_observer_setup_from_state, solve_state
from src import Interpolation
from tests.public_api_builders import numerics, solver_options, top_hat_model


NUM_GAM_E = 8
NUM_CHI = 4
NUM_NU = 8
NUM_R = 8
NUM_THETA = 8
NUM_TOBS = 4


def _build_model(solver: str) -> Model:
    return top_hat_model(
        numerics=numerics(
            num_electron_gamma=NUM_GAM_E,
            num_chi=NUM_CHI,
            num_photon_frequency=NUM_NU,
            num_radius=NUM_R,
            num_theta=NUM_THETA,
            num_observer_time=NUM_TOBS,
        ),
        solver_options=solver_options(electron_solver=solver),
    )


def _build_model_with_geometry(solver: str, geometry_kernel: str) -> Model:
    return top_hat_model(
        numerics=numerics(
            num_electron_gamma=NUM_GAM_E,
            num_chi=NUM_CHI,
            num_photon_frequency=NUM_NU,
            num_radius=NUM_R,
            num_theta=NUM_THETA,
            num_observer_time=NUM_TOBS,
        ),
        solver_options=solver_options(electron_solver=solver, geometry_projection=geometry_kernel),
    )


def case_basic_smoke():
    model = _build_model("fullhide_2d")
    flux = model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total
    assert flux.shape == (1,)
    assert np.all(np.isfinite(flux))
    return {"flux_shape": list(flux.shape), "flux": float(flux[0])}


def case_electron_grid():
    state = solve_state(
        RuntimeConfig(
            electron_solver="fullhide_2d",
            num_gam_e=NUM_GAM_E,
            num_chi=NUM_CHI,
            num_nu=NUM_NU,
            num_r=NUM_R,
            num_theta=NUM_THETA,
            num_tobs=NUM_TOBS,
        ),
        np.array([1.0e2, 1.1e2]),
    )
    electron = state.electron
    chi_grid = np.asarray(electron.chi_grid, dtype=float)
    d_n_gam_e_chi = np.asarray(electron.d_n_gam_e_chi, dtype=float)
    d_n_gam_e = np.asarray(electron.d_n_gam_e, dtype=float)
    l_syn_spec_chi = np.asarray(electron.l_syn_spec_chi, dtype=float)
    seed_syn_chi = np.asarray(electron.seed_syn_chi, dtype=float)
    tau_syn_chi = np.asarray(electron.tau_syn_chi, dtype=float)
    chi_radius_cm = np.asarray(electron.chi_radius_cm, dtype=float)
    chi_gamma_bulk = np.asarray(electron.chi_gamma_bulk, dtype=float)
    chi_dvolume_weight = np.asarray(electron.chi_dvolume_weight, dtype=float)

    assert chi_grid.ndim == 1
    assert chi_grid.shape == (NUM_CHI,)
    assert np.all(np.isfinite(chi_grid))

    assert d_n_gam_e_chi.ndim == 3
    assert d_n_gam_e_chi.shape[0] == NUM_GAM_E
    assert d_n_gam_e_chi.shape[1] == NUM_CHI
    assert d_n_gam_e_chi.shape[2] == d_n_gam_e.shape[1]
    assert np.all(np.isfinite(d_n_gam_e_chi))
    assert l_syn_spec_chi.shape == (NUM_NU, NUM_CHI, d_n_gam_e.shape[1])
    assert seed_syn_chi.shape == l_syn_spec_chi.shape
    assert tau_syn_chi.shape == l_syn_spec_chi.shape
    assert chi_radius_cm.shape == (NUM_CHI, d_n_gam_e.shape[1])
    assert chi_gamma_bulk.shape == chi_radius_cm.shape
    assert chi_dvolume_weight.shape == chi_radius_cm.shape
    assert np.all(np.isfinite(l_syn_spec_chi))
    assert np.all(np.isfinite(seed_syn_chi))
    assert np.all(np.isfinite(tau_syn_chi))
    assert np.all(np.isfinite(chi_radius_cm))
    assert np.all(np.isfinite(chi_gamma_bulk))
    assert np.all(np.isfinite(chi_dvolume_weight))
    assert np.all(chi_radius_cm > 0.0)
    q_active = 1.0 - (3.0 / 4.0) ** 4
    q_centers = (np.arange(NUM_CHI, dtype=float) + 0.5) * q_active / NUM_CHI
    expected_chi = (1.0 - q_centers) ** (-4.0 / 3.0)
    np.testing.assert_allclose(chi_grid, expected_chi, rtol=1.0e-14, atol=1.0e-14)
    np.testing.assert_allclose(chi_dvolume_weight, q_active / NUM_CHI, rtol=1.0e-14, atol=1.0e-14)
    np.testing.assert_allclose(np.sum(chi_dvolume_weight, axis=0), q_active, rtol=1.0e-14, atol=1.0e-14)
    radius = np.asarray(state.dynamics.radius, dtype=float)
    gamma_front = np.asarray(state.dynamics.r_gamma, dtype=float)
    assert np.all(chi_radius_cm <= radius[None, :])
    assert np.all(np.diff(chi_radius_cm, axis=0) < 0.0)
    assert np.all(chi_radius_cm[-1, :] > 0.5 * radius)
    assert np.all(chi_gamma_bulk <= gamma_front[None, :])
    assert np.all(chi_gamma_bulk >= 1.0)

    return {
        "chi_shape": list(chi_grid.shape),
        "d_n_gam_e_chi_shape": list(d_n_gam_e_chi.shape),
        "d_n_gam_e_shape": list(d_n_gam_e.shape),
        "l_syn_spec_chi_shape": list(l_syn_spec_chi.shape),
    }


def case_chi_eats_geometry_smoke():
    model = _build_model_with_geometry("fullhide_2d", "chi_eats_2d")
    flux = model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total
    assert flux.shape == (1,)
    assert np.all(np.isfinite(flux))
    return {"flux_shape": list(flux.shape), "flux": float(flux[0])}


def case_chi_eats_projection_kind_routes():
    model = _build_model_with_geometry("fullhide_2d", "chi_eats_2d")
    times = np.array([1.0e4, 3.0e4], dtype=float)
    lightcurve_freqs = np.array([1.0e10, 1.0e14], dtype=float)
    spectrum_freqs = np.array([1.0e10, 1.0e14, 1.0e18], dtype=float)
    lightcurve = model.flux_density_grid(times, lightcurve_freqs, projection_kind="lightcurve").total
    spectrum = model.spectrum(1.0e4, spectrum_freqs)
    spectrum_grid = model.flux_density_grid(np.array([1.0e4]), spectrum_freqs, projection_kind="sed").total[:, 0]
    assert lightcurve.shape == (lightcurve_freqs.size, times.size)
    assert spectrum.shape == spectrum_freqs.shape
    np.testing.assert_allclose(spectrum, spectrum_grid, rtol=0.0, atol=0.0)
    assert np.all(np.isfinite(lightcurve))
    assert np.all(np.isfinite(spectrum))
    try:
        model.flux_density_grid(times, lightcurve_freqs, projection_kind="bad")
    except ValueError as exc:
        assert "projection_kind" in str(exc)
    else:
        raise AssertionError("invalid projection_kind was accepted")
    return {"lightcurve_shape": list(lightcurve.shape), "spectrum_shape": list(spectrum.shape)}


def case_chi_eats_rejects_1d_solver():
    model = _build_model_with_geometry("fullhide_1d", "chi_eats_2d")
    try:
        model.flux_density(np.array([1.0e4]), np.array([1.0e14]))
    except ValueError as exc:
        assert "requires a 2d electron solver" in str(exc)
        return {"error": str(exc)}
    raise AssertionError("chi_eats_2d accepted a 1d electron solver")


def case_off_axis_phi_collapse_rejected():
    model = _build_model_with_geometry("fullhide_2d", "chi_eats_2d")
    model.observer.theta_obs = 0.03
    model.setups.num_phi = 1
    try:
        model.flux_density(np.array([1.0e4]), np.array([1.0e14]))
    except ValueError as exc:
        assert "off-axis EATS projection requires num_phi >= 2" in str(exc)
        return {"error": str(exc)}
    raise AssertionError("off-axis EATS projection accepted num_phi=1")


def case_project_flux_grid_syncs_observer_theta_boundary():
    boundary = np.zeros(30, dtype=float)
    boundary[8] = 0.1
    boundary[9] = 0.2
    state = SimpleNamespace(
        config=SimpleNamespace(opening_angle_jet=0.12, theta_v=0.5),
        setup=SimpleNamespace(
            luminosity_distance_cm=1.0e26,
            boundary=boundary,
            seed_frequency_hz=np.array([1.0e10, 1.0e14], dtype=float),
        ),
    )
    setup = _build_observer_setup_from_state(state, np.array([1.0e4], dtype=float))
    assert setup.boundary[8] == 0.12
    assert setup.boundary[9] == 0.5
    assert boundary[8] == 0.1
    assert boundary[9] == 0.2
    return {
        "observer_opening": float(setup.boundary[8]),
        "observer_theta": float(setup.boundary[9]),
        "state_opening": float(boundary[8]),
        "state_theta": float(boundary[9]),
    }


def case_model_cache_includes_observer_angle():
    model = top_hat_model(
        numerics=numerics(
            num_electron_gamma=NUM_GAM_E,
            num_chi=NUM_CHI,
            num_photon_frequency=NUM_NU,
            num_radius=NUM_R,
            num_theta=NUM_THETA,
            num_phi=12,
            num_observer_time=NUM_TOBS,
        ),
        solver_options=solver_options(electron_solver="fullhide_2d", geometry_projection="chi_eats_2d"),
    )
    times = np.geomspace(1.0e4, 1.0e7, 8)
    freq = np.array([1.0e14], dtype=float)
    model.observer.theta_obs = 0.2
    near = model.flux_density_grid(times, freq, projection_kind="lightcurve").total[0]
    model.observer.theta_obs = 0.5
    far = model.flux_density_grid(times, freq, projection_kind="lightcurve").total[0]
    mask = np.isfinite(near) & np.isfinite(far) & (near > 0.0) & (far > 0.0)
    assert np.any(mask)
    log_delta = np.max(np.abs(np.log10(far[mask] / near[mask])))
    assert log_delta > 0.5
    return {"max_abs_log10_delta": float(log_delta)}


def case_on_axis_phi_collapse_matches_explicit_phi():
    boundary = np.zeros(30, dtype=float)
    boundary[0] = 300.0
    boundary[3] = 1.0e15
    boundary[7] = 0.1
    boundary[8] = 0.1
    boundary[9] = 0.0
    radius = np.geomspace(1.0e16, 1.0e18, 16)
    gamma = np.geomspace(200.0, 20.0, radius.size)
    r_tobs = (1.0 + boundary[7]) * radius / (2.0 * gamma**2) / 2.99792458e10
    seed = np.geomspace(1.0e8, 1.0e20, 16)
    obs = np.array([1.0e14])
    tobs = np.geomspace(r_tobs[2], r_tobs[-2], 12)
    source = (seed[:, None] / 1.0e12) ** 0.3 * (radius[None, :] / radius[0]) ** -1.1
    collapsed = Interpolation.sed_interpolation(boundary, r_tobs, gamma, radius, source, seed, obs, tobs, 24, 1, 1)
    explicit = Interpolation.sed_interpolation(boundary, r_tobs, gamma, radius, source, seed, obs, tobs, 24, 24, 1)
    np.testing.assert_allclose(collapsed, explicit, rtol=1.0e-12, atol=1.0e-12)
    return {"collapsed_sum": float(np.sum(collapsed)), "explicit_sum": float(np.sum(explicit))}


def _thin_shell_projection_inputs():
    boundary = np.zeros(30, dtype=float)
    boundary[0] = 100.0
    boundary[3] = 1.0e15
    boundary[7] = 0.1
    boundary[8] = 0.08
    boundary[9] = 0.0
    radius = np.geomspace(1.0e16, 5.0e17, 40)
    gamma = np.geomspace(80.0, 12.0, radius.size)
    r_tobs = (1.0 + boundary[7]) * radius / (2.0 * gamma**2) / 2.99792458e10
    seed = np.geomspace(1.0e8, 1.0e20, 80)
    obs = np.array([1.0e10, 1.0e14, 1.0e18], dtype=float)
    theta_center = 0.5 * boundary[8]
    arrival = r_tobs + radius * (1.0 - np.cos(theta_center)) * (1.0 + boundary[7]) / 2.99792458e10
    tobs = np.geomspace(arrival[5], arrival[-6], 48)
    source = (seed[:, None] / 1.0e12) ** 0.2 * (radius[None, :] / radius[0]) ** -0.8
    return boundary, radius, gamma, r_tobs, seed, obs, tobs, source


def case_chi_projection_delta_layer_matches_thin_shell():
    boundary, radius, gamma, r_tobs, seed, obs, tobs, source = _thin_shell_projection_inputs()
    legacy = Interpolation.sed_interpolation(boundary, r_tobs, gamma, radius, source, seed, obs, tobs, 32, 24, 1)
    num_chi = 5
    source_chi = np.zeros((seed.size, num_chi, radius.size), dtype=float)
    source_chi[:, 0, :] = source
    tau_chi = np.zeros_like(source_chi)
    r_chi = np.tile(radius, (num_chi, 1))
    gamma_chi = np.tile(gamma, (num_chi, 1))
    chi_weight = np.ones((num_chi, radius.size), dtype=float)
    projected = Interpolation.sed_interpolation_chi(
        boundary,
        r_tobs,
        radius,
        source_chi,
        tau_chi,
        r_chi,
        gamma_chi,
        chi_weight,
        seed,
        obs,
        tobs,
        32,
        24,
        1,
    )
    assert np.any(legacy > 0.0)
    np.testing.assert_allclose(projected, legacy, rtol=1.0e-12, atol=1.0e-12)
    return {"max_flux": float(np.max(legacy)), "max_abs_diff": float(np.max(np.abs(projected - legacy)))}


def case_chi_projection_finite_width_converges_to_thin_shell():
    boundary, radius, gamma, r_tobs, seed, obs, tobs, source = _thin_shell_projection_inputs()
    legacy = Interpolation.sed_interpolation(boundary, r_tobs, gamma, radius, source, seed, obs, tobs, 32, 24, 1)
    num_chi = 5
    weights = np.array([0.52, 0.22, 0.13, 0.08, 0.05], dtype=float)[:, None]
    source_chi = np.repeat(source[:, None, :], num_chi, axis=1)
    tau_chi = np.zeros_like(source_chi)
    chi_weight = np.repeat(weights, radius.size, axis=1)
    offsets = np.linspace(0.0, 1.0, num_chi)[:, None]
    mask = legacy > np.max(legacy) * 1.0e-12
    errors = []
    for eps in (0.5, 0.25, 0.125, 0.0625, 0.0):
        r_chi = radius[None, :] * (1.0 - eps * 0.03 * offsets)
        gamma_chi = gamma[None, :] * (1.0 - eps * 0.08 * offsets)
        projected = Interpolation.sed_interpolation_chi(
            boundary,
            r_tobs,
            radius,
            source_chi,
            tau_chi,
            r_chi,
            gamma_chi,
            chi_weight,
            seed,
            obs,
            tobs,
            32,
            24,
            1,
        )
        errors.append(float(np.sqrt(np.mean((np.log(projected[mask]) - np.log(legacy[mask])) ** 2))))
    assert all(later < earlier for earlier, later in zip(errors[:-2], errors[1:-1]))
    assert errors[-1] < 1.0e-12
    return {"log_rms_errors": errors}


def case_chi_ssa_cell_split_invariance():
    boundary = np.zeros(30, dtype=float)
    boundary[0] = 100.0
    boundary[3] = 1.0e15
    boundary[7] = 0.0
    boundary[8] = 0.1
    boundary[9] = 0.0
    radius = np.linspace(1.0e16, 1.3e16, 4)
    r_tobs = np.linspace(1.0e2, 4.0e2, radius.size)
    theta_center = 0.5 * boundary[8]
    arrival = r_tobs + (radius - radius * np.cos(theta_center)) / 2.99792458e10
    seed = np.geomspace(1.0e6, 1.0e20, 64)
    obs = np.array([1.0e10])
    tobs = np.array([0.5 * (arrival[1] + arrival[2])])
    values = []
    for num_chi in (1, 2, 4, 8):
        source = np.ones((seed.size, num_chi, radius.size), dtype=float) / float(num_chi)
        tau = np.ones_like(source) * (2.0 / float(num_chi))
        r_chi = np.tile(radius, (num_chi, 1))
        gamma_chi = np.full((num_chi, radius.size), 30.0)
        chi_weight = np.ones((num_chi, radius.size), dtype=float)
        flux = Interpolation.sed_interpolation_chi(
            boundary,
            r_tobs,
            radius,
            source,
            tau,
            r_chi,
            gamma_chi,
            chi_weight,
            seed,
            obs,
            tobs,
            1,
            1,
            1,
        )
        values.append(float(flux[0, 0]))
    np.testing.assert_allclose(values, values[-1], rtol=1.0e-12, atol=1.0e-12)
    return {"flux_by_num_chi": values}


def case_chi_ssa_nonuniform_tau_matches_manual():
    boundary = np.ones(24, dtype=float)
    boundary[7] = 0.0
    boundary[9] = 0.0
    gamma = 12.0
    beta = np.sqrt(1.0 - gamma**-2)
    doppler_den = gamma * (1.0 - beta)
    radius = np.array([1.0e15, 1.2e15], dtype=float)
    r_tobs = np.array([1.0, 3.0], dtype=float)
    seed = np.array([1.0e9, 1.0e10, 1.0e11], dtype=float)
    observed = np.array([seed[1] / doppler_den], dtype=float)
    tobs = np.array([2.0], dtype=float)
    d_omega = 1.0e-4
    tau_profile = np.array([0.1, 1.3, 0.4, 2.2], dtype=float)
    num_chi = tau_profile.size
    source = np.ones((seed.size, num_chi, radius.size), dtype=float)
    tau = np.repeat(tau_profile[None, :, None], seed.size, axis=0)
    tau = np.repeat(tau, radius.size, axis=2)
    r_chi = np.tile(radius, (num_chi, 1))
    gamma_chi = np.full((num_chi, radius.size), gamma)
    chi_weight = np.full((num_chi, radius.size), 1.0 / num_chi)
    flux = Interpolation.sed_interpolation_chi_surface_element(
        boundary,
        r_tobs,
        radius,
        source,
        tau,
        r_chi,
        gamma_chi,
        chi_weight,
        seed,
        observed,
        tobs,
        d_omega,
    )
    front_tau = 0.0
    escaped = 0.0
    for tau_cell in tau_profile:
        escaped += (1.0 / num_chi) * np.exp(-front_tau) * (1.0 - np.exp(-tau_cell)) / tau_cell
        front_tau += tau_cell
    expected = escaped * d_omega / (4.0 * np.pi) * doppler_den**-3
    np.testing.assert_allclose(flux[0, 0], expected, rtol=1.0e-12, atol=0.0)
    return {"flux": float(flux[0, 0]), "expected": float(expected)}


def main() -> None:
    cases = [
        ("[1/13] fullhide_2d:basic_smoke", case_basic_smoke),
        ("[2/13] fullhide_2d:electron_grid", case_electron_grid),
        ("[3/13] fullhide_2d:chi_eats_geometry", case_chi_eats_geometry_smoke),
        ("[4/13] chi_eats_2d:projection_kind_routes", case_chi_eats_projection_kind_routes),
        ("[5/13] chi_eats_2d:rejects_1d_solver", case_chi_eats_rejects_1d_solver),
        ("[6/13] eats:rejects_off_axis_phi_collapse", case_off_axis_phi_collapse_rejected),
        ("[7/13] eats:syncs_observer_theta_boundary", case_project_flux_grid_syncs_observer_theta_boundary),
        ("[8/13] model_cache:observer_angle", case_model_cache_includes_observer_angle),
        ("[9/13] eats:on_axis_phi_collapse", case_on_axis_phi_collapse_matches_explicit_phi),
        ("[10/13] chi_eats_2d:delta_layer_thin_shell", case_chi_projection_delta_layer_matches_thin_shell),
        ("[11/13] chi_eats_2d:finite_width_converges", case_chi_projection_finite_width_converges_to_thin_shell),
        ("[12/13] chi_eats_2d:ssa_cell_split_invariance", case_chi_ssa_cell_split_invariance),
        ("[13/13] chi_eats_2d:ssa_nonuniform_tau", case_chi_ssa_nonuniform_tau_matches_manual),
    ]
    for label, fn in cases:
        print(f"  {label} ...", flush=True)
        fn()


if __name__ == "__main__":
    main()
