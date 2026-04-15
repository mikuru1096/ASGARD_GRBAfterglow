from VegasAfterglow import (
    ISM, Model, Observer, Radiation, Setups, TophatJet, observe,
)
from asgard_models import SpectrumOutputConfig
from asgard_plot import plot_light_curve, plot_spectrum
from asgard_presets import build_spectrum_demo_config


def build_demo_model():
    c = build_spectrum_demo_config()
    jet = TophatJet(E_iso=c.e_iso, Gamma0=c.eta_0, theta_j=c.opening_angle_jet)
    medium = ISM(n_ism=c.d_ne)
    fwd_rad = Radiation(eps_e=c.epsilon_e, eps_B=c.epsilon_b, p=c.p, xi_N=c.f_e, ssc=True)
    observer = Observer(z=c.z, theta_obs=c.theta_v)
    setups = Setups(
        num_r=c.num_r, num_theta=c.num_theta, num_phi=c.num_phi,
        num_gam_e=c.num_gam_e, num_nu=c.num_nu, num_tobs=c.num_tobs,
        observer_time_min_s=10**c.t_obs_min_log10,
        observer_time_max_s=10**c.t_obs_max_log10,
        ssc_cooling=True,
    )
    return Model(medium=medium, jet=jet, observer=observer, fwd_rad=fwd_rad, setups=setups), c


def main() -> None:
    model, config = build_demo_model()
    result = observe(model, config=config, spectrum_output=config.spectrum_output)
    plot_light_curve(result, outfile="Radiation_Lightcurves.pdf", show=False)
    plot_spectrum(result, times_s=[1e3, 1e4, 1e5], quantity="nufnu", outfile="Radiation_Spectra.pdf", show=False)


if __name__ == "__main__":
    main()
