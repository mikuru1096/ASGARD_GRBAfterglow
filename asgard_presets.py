from __future__ import annotations

from asgard_models import FitConfig, ReverseShockConfig, SpectrumOutputConfig, default_num_threads


def build_baseline_config(**overrides) -> FitConfig:
    config = FitConfig(
        num_threads=default_num_threads(),
        index_dyn=3,
        index_y=2,
        index_syn_integr=2,
        num_gam_e=201,
        num_nu=201,
        num_r=300,
        num_theta=300,
        num_phi=1,
        z=0.4,
        eta_0=1.0e2,
        epsilon_e=1.0e-1,
        epsilon_b=1.0e-3,
        p=2.5,
        opening_angle_jet=1.0e-1,
        theta_v=0.0,
        f_e=1.0e-1,
        e_iso=1.0e53,
        d_ne=1.0e-1,
        a_star=-1.0,
        r0=1.0e9,
        ebv=0.0,
        rv=2.93,
        lyman_ar=0.0,
        f_sys=-1.0,
        plot_lc=False,
    )
    for key, value in overrides.items():
        setattr(config, key, value)
    return config


def build_spectrum_demo_config(**overrides) -> FitConfig:
    return build_baseline_config(
        spectrum_output=SpectrumOutputConfig(enabled=True, num_nu_obs=200),
        **overrides,
    )


def build_reverse_demo_config(**overrides) -> FitConfig:
    return build_baseline_config(
        reverse_shock=ReverseShockConfig(
            enabled=True,
            delta_t_s=10.0,
            epsilon_e=1.0e-1,
            epsilon_b=1.0e-2,
            p=2.4,
            f_e=1.0,
        ),
        **overrides,
    )
