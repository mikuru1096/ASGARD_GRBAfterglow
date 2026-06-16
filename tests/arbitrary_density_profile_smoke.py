from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.api_model import (
    Hadronic,
    Model,
    Numerics,
    Observer,
    ObserverGrid,
    Radiation,
    ReverseShock,
    SolverOptions,
    TabulatedMedium,
    UniformMedium,
    top_hat_jet,
)
from asgard_core.api_model import _build_fit_config_for_patch
from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_physics_utils import ambient_density
from asgard_core.asgard_state import make_query_setup, solve_state_from_setup


def _profile_points() -> tuple[tuple[float, ...], tuple[float, ...]]:
    radius = (1.0e15, 3.0e15, 1.0e16, 3.0e16)
    density = (1.0e2, 1.0e1, 3.0e2, 2.0e0)
    return radius, density


def _run_fit_config_profile() -> None:
    radius, density = _profile_points()
    config = RuntimeConfig(
        num_threads=1,
        num_r=48,
        num_theta=12,
        num_tobs=32,
        num_gam_e=41,
        num_nu=32,
        t_obs_min_log10=1.0,
        t_obs_max_log10=5.0,
        eta_0=80.0,
        e_iso=1.0e52,
        d_ne=1.0,
        a_star=-1.0,
        include_forward_ssc=False,
        electron_solver="fullhide_1d",
        density_profile_radius_cm=radius,
        density_profile_n_cm3=density,
    )
    probe = np.array([1.0e15, np.sqrt(3.0e15 * 1.0e16), 3.0e16])
    expected = np.exp(np.interp(np.log(probe), np.log(radius), np.log(density)))
    assert np.allclose(ambient_density(probe, config), expected, rtol=1.0e-14, atol=0.0)

    times = np.logspace(1.0, 5.0, 12)
    freqs = np.array([1.0e10, 1.0e14], dtype=float)
    state = solve_state_from_setup(config, make_query_setup(config, times, freqs))
    assert np.all(np.isfinite(state.dynamics.radius))
    assert np.all(np.isfinite(state.dynamics.r_gamma))
    assert np.all(state.components.fwd_sync >= 0.0)


def _run_public_medium_profile() -> None:
    radius, density = _profile_points()
    numerics = Numerics(
        num_radius=48,
        num_theta=12,
        num_phi=1,
        num_observer_time=32,
        num_electron_gamma=41,
        num_photon_frequency=32,
        num_chi=None,
        num_threads=1,
        electron_adaptive_substeps=False,
        electron_substep_rtol=0.02,
        electron_substep_min=100,
        electron_substep_max=1000,
        initial_radius_cm=1.0e14,
    )
    fwd_rad = Radiation(
        epsilon_e=0.1,
        epsilon_B=1.0e-4,
        p=2.3,
        proton_energy_fraction=0.0,
        epsilon_b_floor=None,
        magnetic_decay_alpha_t=0.0,
        magnetic_decay_t0_s=1.0,
        accelerated_electron_fraction=0.1,
        thermal_electrons=False,
        include_ssc=False,
        include_kn_correction=False,
        proton_synch=True,
        include_pgamma=False,
        bethe_heitler=False,
        hadronic_inverse_compton=False,
        pp=False,
        neutrino=False,
        acceleration_efficiency=1.0,
        reverse_proton_energy_fraction=0.0,
        pgamma_scheme="disabled",
        pair_production=False,
    )
    solver_options = SolverOptions(
        electron_solver="fullhide_1d",
        dynamics_solver="forward_legacy",
        geometry_projection="sed_legacy",
        electron_photon_coupling="separated",
        ssc_cooling_mode="none",
        synchrotron_integration="fixed_grid",
        cooling_kernel="legacy",
        radiation_kernel="legacy",
        structured_backend="fortran_1d",
        patch_sampling="uniform",
        patch_projection="auto",
        patch_sampling_pilot_theta=0,
        patch_sampling_num_times=12,
        patch_sampling_beaming_factor=3.0,
        patch_sampling_beaming_resolution=8.0,
        structured_parallel_mode="outer",
        structured_outer_threads=None,
        structured_inner_threads=None,
        fullhide2d_transport_model="legacy",
        fullhide2d_stochastic_accel_norm=0.0,
        fullhide2d_escape_mode="closed",
    )
    model = Model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=80.0,
            opening_angle_rad=0.1,
            shell_duration_s=10.0,
            magnetar=None,
            spreading=False,
        ),
        medium=TabulatedMedium(radius, density, label="test_density_profile"),
        observer=Observer(z=0.1, viewing_angle_rad=0.0, viewing_azimuth_rad=0.0, luminosity_distance_cm=1.0e27),
        fwd_rad=fwd_rad,
        rvs_rad=None,
        numerics=numerics,
        observer_grid=ObserverGrid(time_min_s=10.0, time_max_s=1.0e5),
        solver_options=solver_options,
        reverse_shock=ReverseShock(enabled=False, shell_duration_s=10.0, upstream_sigma=0.0, include_cross_zone_ic=False, include_ssc=False),
        hadronic=Hadronic(enabled=False, solver="legacy_1d", num_proton_gamma=161, num_neutrino_frequency=121, pgamma_scheme="disabled", pair_cascade_iterations=1),
    )
    flux = model.flux_density_grid(np.logspace(1.0, 5.0, 10), np.array([1.0e10, 1.0e14]))
    assert np.all(np.isfinite(flux.total))
    assert np.any(flux.total > 0.0)

    setup_override = Model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=80.0,
            opening_angle_rad=0.1,
            shell_duration_s=10.0,
            magnetar=None,
            spreading=False,
        ),
        medium=UniformMedium(number_density_cm3=1.0),
        observer=Observer(z=0.1, viewing_angle_rad=0.0, viewing_azimuth_rad=0.0, luminosity_distance_cm=1.0e27),
        fwd_rad=fwd_rad,
        rvs_rad=None,
        numerics=numerics,
        observer_grid=ObserverGrid(time_min_s=10.0, time_max_s=1.0e5),
        solver_options=solver_options,
        reverse_shock=ReverseShock(enabled=False, shell_duration_s=10.0, upstream_sigma=0.0, include_cross_zone_ic=False, include_ssc=False),
        hadronic=Hadronic(enabled=False, solver="legacy_1d", num_proton_gamma=161, num_neutrino_frequency=121, pgamma_scheme="disabled", pair_cascade_iterations=1),
    )
    override_flux = setup_override.flux_density_grid(np.logspace(1.0, 5.0, 10), np.array([1.0e10, 1.0e14]))
    assert np.all(np.isfinite(override_flux.total))
    assert np.any(override_flux.total > 0.0)
    cfg_medium = _build_fit_config_for_patch(model, theta_v=0.0, opening_angle_jet=0.1, e_iso=1.0e52, gamma0=80.0)
    cfg_uniform = _build_fit_config_for_patch(setup_override, theta_v=0.0, opening_angle_jet=0.1, e_iso=1.0e52, gamma0=80.0)
    assert cfg_medium.density_profile_radius_cm == radius
    assert cfg_medium.density_profile_n_cm3 == density
    assert cfg_uniform.density_profile_radius_cm == ()
    assert cfg_uniform.density_profile_n_cm3 == ()


if __name__ == "__main__":
    _run_fit_config_profile()
    _run_public_medium_profile()
    print("arbitrary_density_profile_smoke: ok")
