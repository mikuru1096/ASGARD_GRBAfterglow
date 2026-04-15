from VegasAfterglow import (
    ISM, Model, Observer, Radiation, Setups, TophatJet, observe,
)
from asgard_plot import plot_characteristic_frequencies, plot_light_curve
from asgard_presets import build_reverse_demo_config


def build_demo_model():
    c = build_reverse_demo_config()
    jet = TophatJet(E_iso=c.e_iso, Gamma0=c.eta_0, theta_j=c.opening_angle_jet)
    medium = ISM(n_ism=c.d_ne)
    fwd_rad = Radiation(eps_e=c.epsilon_e, eps_B=c.epsilon_b, p=c.p, xi_N=c.f_e, ssc=True)
    rvs_rad = Radiation(eps_e=c.reverse_shock.epsilon_e, eps_B=c.reverse_shock.epsilon_b, p=c.reverse_shock.p, xi_N=c.reverse_shock.f_e)
    observer = Observer(z=c.z, theta_obs=c.theta_v)
    setups = Setups(
        num_r=c.num_r, num_theta=c.num_theta, num_phi=c.num_phi,
        num_gam_e=c.num_gam_e, num_nu=c.num_nu, num_tobs=c.num_tobs,
        observer_time_min_s=10**c.t_obs_min_log10,
        observer_time_max_s=10**c.t_obs_max_log10,
        ssc_cooling=True,
        rvs_shock=True,
    )
    return Model(medium=medium, jet=jet, observer=observer, fwd_rad=fwd_rad, rvs_rad=rvs_rad, setups=setups), c


def main() -> None:
    model, config = build_demo_model()
    result = observe(model, config=config)
    plot_light_curve(result, outfile="Radiation_Lightcurves_RS.pdf", show=False)
    plot_characteristic_frequencies(result, include_reverse=True, outfile="Characteristic_Frequencies_RS.pdf", show=False)


if __name__ == "__main__":
    main()
