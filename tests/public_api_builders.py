from __future__ import annotations

from asgard_core import (
    Hadronic,
    Model,
    Numerics,
    Observer,
    ObserverGrid,
    Radiation,
    ReverseShock,
    SolverOptions,
    UniformMedium,
    top_hat_jet,
)


def radiation(**updates) -> Radiation:
    return Radiation(**(dict(
        epsilon_e=0.1,
        epsilon_B=1.0e-3,
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
    ) | updates))


def numerics(**updates) -> Numerics:
    return Numerics(**(dict(
            structured_num_theta=12,
            structured_num_phi=24,
        num_radius=48,
        eats_num_theta=12,
        eats_num_phi=1,
        num_observer_time=24,
        num_electron_gamma=41,
        num_photon_frequency=31,
        downstream_num_chi=None,
        num_threads=1,
        electron_adaptive_substeps=False,
        electron_substep_rtol=0.02,
        electron_substep_min=100,
        electron_substep_max=1000,
        initial_radius_cm=1.0e14,
    ) | updates))


def observer_grid(**updates) -> ObserverGrid:
    return ObserverGrid(**(dict(time_min_s=1.0e2, time_max_s=1.0e6) | updates))


def solver_options(**updates) -> SolverOptions:
    return SolverOptions(**(dict(
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
        adaptive_grid="manual",
        projection_adaptive_rtol=2.0e-2,
        projection_adaptive_max_depth=4,
        fullhide2d_transport_model="legacy",
        fullhide2d_stochastic_accel_norm=0.0,
        fullhide2d_escape_mode="closed",
    ) | updates))


def reverse_shock(**updates) -> ReverseShock:
    return ReverseShock(**(dict(
        enabled=False,
        shell_duration_s=10.0,
        upstream_sigma=0.0,
        include_cross_zone_ic=False,
        include_ssc=False,
    ) | updates))


def hadronic(**updates) -> Hadronic:
    return Hadronic(**(dict(
        enabled=False,
        solver="legacy_1d",
        num_proton_gamma=161,
        num_neutrino_frequency=121,
        pgamma_scheme="disabled",
        pair_cascade_iterations=1,
    ) | updates))


def top_hat_model(**updates) -> Model:
    return Model(**(dict(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=300.0,
            opening_angle_rad=0.1,
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        medium=UniformMedium(number_density_cm3=1.0),
        observer=Observer(z=0.1, viewing_angle_rad=0.0, viewing_azimuth_rad=0.0, luminosity_distance_cm=1.0e26),
        fwd_rad=radiation(),
        rvs_rad=None,
        numerics=numerics(),
        observer_grid=observer_grid(),
        solver_options=solver_options(),
        reverse_shock=reverse_shock(),
        hadronic=hadronic(),
    ) | updates))
