from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet, observe
from asgard_core.asgard_plot import plot_light_curve, plot_spectrum
from asgard_core.asgard_presets import build_spectrum_demo_config


def main() -> None:
    config = build_spectrum_demo_config()
    model = Model(
        jet=TophatJet(E_iso=config.e_iso, Gamma0=config.eta_0, theta_j=config.opening_angle_jet),
        medium=ISM(n_ism=config.d_ne),
        observer=Observer(z=config.z, theta_obs=config.theta_v),
        fwd_rad=Radiation(eps_e=config.epsilon_e, eps_B=config.epsilon_b, p=config.p, xi_N=config.f_e, ssc=True),
        setups=Setups(
            num_r=config.num_r,
            num_theta=config.num_theta,
            num_phi=config.num_phi,
            num_gam_e=config.num_gam_e,
            num_nu=config.num_nu,
            num_tobs=config.num_tobs,
            observer_time_min_s=10**config.t_obs_min_log10,
            observer_time_max_s=10**config.t_obs_max_log10,
            ssc_cooling=True,
        ),
    )
    result = observe(model, config=config, spectrum_output=config.spectrum_output)
    plot_light_curve(result, outfile="Radiation_Lightcurves.pdf", show=False)
    plot_spectrum(result, times_s=[1e3, 1e4, 1e5], quantity="nufnu", outfile="Radiation_Spectra.pdf", show=False)


if __name__ == "__main__":
    main()
