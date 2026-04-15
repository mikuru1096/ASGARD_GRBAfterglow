from __future__ import annotations

from pathlib import Path
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_models import ReverseShockConfig, SpectrumOutputConfig
from asgard_presets import build_baseline_config
from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet, observe


def build_contract_model_and_config():
    c = build_baseline_config(
        num_gam_e=161,
        num_nu=161,
        num_r=240,
        num_theta=240,
        reverse_shock=ReverseShockConfig(
            enabled=True,
            delta_t_s=10.0,
            epsilon_e=1.0e-1,
            epsilon_b=1.0e-2,
            p=2.4,
            f_e=1.0,
        ),
        spectrum_output=SpectrumOutputConfig(enabled=True, num_nu_obs=96),
    )
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
        ssc_cooling=True, rvs_shock=True,
    )
    model = Model(medium=medium, jet=jet, observer=observer, fwd_rad=fwd_rad, rvs_rad=rvs_rad, setups=setups)
    return model, c


def main() -> None:
    model, config = build_contract_model_and_config()
    result = observe(model, config=config, spectrum_output=config.spectrum_output)

    assert result.characteristic_time_s.shape == result.nu_m.shape
    assert result.characteristic_time_s.shape == result.nu_c.shape
    assert result.characteristic_time_s.shape == result.nu_a.shape
    assert result.rs_nu_m is not None
    assert result.rs_nu_c is not None
    assert result.rs_nu_a is not None
    assert result.characteristic_time_s.shape == result.rs_nu_m.shape
    assert result.characteristic_time_s.shape == result.rs_nu_c.shape
    assert result.characteristic_time_s.shape == result.rs_nu_a.shape
    assert result.spectrum_freq_hz is not None
    assert result.spectrum_fnu is not None
    assert result.spectrum_fnu.shape == (result.spectrum_freq_hz.shape[0], result.t_obs_s.shape[0])

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        from asgard_plot import plot_characteristic_frequencies, plot_spectrum
        freq_plot = tmpdir_path / "characteristic.pdf"
        spec_plot = tmpdir_path / "spectrum.pdf"

        plot_characteristic_frequencies(result, include_reverse=True, outfile=str(freq_plot), show=False)
        plot_spectrum(result, times_s=[1.0e3, 1.0e4], quantity="nufnu", outfile=str(spec_plot), show=False)

        assert freq_plot.exists() and freq_plot.stat().st_size > 0
        assert spec_plot.exists() and spec_plot.stat().st_size > 0

    print("PASS: API contract check succeeded.")


if __name__ == "__main__":
    main()
