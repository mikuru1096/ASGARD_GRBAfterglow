from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet
from src.Electron.electron_radiation import electron_radiation_kernel
from src import constants
import src.Hadronic.FS_hadronic_1d as hadronic_fortran_module


def _base_setups() -> Setups:
    return Setups(
        electron_solver="fullhide_1d",
        num_threads=1,
        num_gam_e=24,
        num_gam_p=24,
        num_nu=32,
        num_nu_nu=16,
        num_r=28,
        num_theta=8,
        num_tobs=28,
        patch_theta=3,
        patch_phi=12,
    )


def _model(theta_obs: float, *, reverse: bool = False, hadronic: bool = False, rs_hadronic: bool = False) -> Model:
    return Model(
        TophatJet(0.1, 1.0e52, 300.0),
        ISM(1.0),
        Observer(1.0e26, 0.1, theta_obs),
        Radiation(
            0.1,
            1.0e-3,
            2.3,
            epsilon_p=0.3 if hadronic else 0.0,
            proton_synch=True,
            reverse_epsilon_p=0.3 if rs_hadronic else 0.0,
        ),
        rvs_rad=Radiation(0.1, 1.0e-3, 2.4) if reverse else None,
        setups=Setups(
            **{
                **_base_setups().__dict__,
                "rvs_shock": reverse,
                "reverse_delta_t_s": 10.0,
                "hadronic_enabled": hadronic,
                "hadronic_solver": "legacy_1d",
            }
        ),
    )


def test_on_axis_shock_random_cancels() -> None:
    times = np.array([1.0e3, 1.0e4, 1.0e5], dtype=float)
    freqs = np.array([1.0e10, 1.0e14], dtype=float)
    pol = _model(0.0).polarization(times, freqs, magnetic_geometry="shock_random")
    assert pol.I_sync.shape == (freqs.size, times.size)
    assert np.all(np.isfinite(pol.I_sync))
    assert np.all(pol.linear_polarization < 1.0e-10)


def test_off_axis_is_finite_and_smooth() -> None:
    times = np.logspace(3.0, 6.0, 6)
    freqs = np.array([1.0e14], dtype=float)
    pol = _model(0.12).polarization(times, freqs, magnetic_geometry="shock_random")
    values = pol.linear_polarization[0]
    assert np.all(np.isfinite(values))
    assert np.all(values >= 0.0)
    assert np.all(values <= 1.0)
    assert np.max(np.abs(np.diff(values))) < 0.6


def test_toroidal_differs_from_shock_random() -> None:
    times = np.array([1.0e4, 1.0e5], dtype=float)
    freqs = np.array([1.0e14], dtype=float)
    model = _model(0.12)
    shock = model.polarization(times, freqs, magnetic_geometry="shock_random")
    toroidal = model.polarization(times, freqs, magnetic_geometry="toroidal")
    assert np.all(np.isfinite(toroidal.polarization_angle_rad))
    delta = np.max(np.abs(toroidal.Q - shock.Q) + np.abs(toroidal.U - shock.U))
    assert float(delta) > 0.0


def test_reverse_and_hadronic_synchrotron_components() -> None:
    times = np.array([1.0e4], dtype=float)
    freqs = np.array([1.0e14, 1.0e18], dtype=float)
    rs_pol = _model(0.12, reverse=True).polarization(times, freqs)
    had_pol = _model(0.12, hadronic=True).polarization(times, freqs)
    rs_had_pol = _model(0.12, reverse=True, rs_hadronic=True).polarization(times, freqs)
    assert np.any(rs_pol.components["rev_sync"]["I"] > 0.0)
    assert np.any(had_pol.components["hadronic_sync"]["I"] > 0.0)
    assert np.any(rs_had_pol.components["rev_hadronic_sync"]["I"] > 0.0)


def test_fortran_polarization_kernel_conserves_intensity() -> None:
    gamma = np.logspace(1.0, 4.0, 32)
    dnde = gamma ** -2.3
    freq = np.logspace(8.0, 18.0, 24)
    p_perp, p_parallel, pi_nu = electron_radiation_kernel.get_syn_polarization_selected(
        1,
        1.0e16,
        1.0,
        1,
        gamma,
        dnde,
        freq,
        2.3,
    )
    p_syn, _ = electron_radiation_kernel.get_syn_selected(1, 1.0e16, 1.0, 1, gamma, dnde, freq)
    assert np.allclose(p_perp + p_parallel, p_syn, rtol=1.0e-13, atol=0.0)
    expected = (2.3 + 1.0) / (2.3 + 7.0 / 3.0)
    assert np.allclose(pi_nu, expected, rtol=1.0e-15, atol=0.0)


def test_fortran_polarization_kernel_tracks_curved_spectrum() -> None:
    gamma = np.logspace(1.0, 6.0, 96)
    dnde = gamma ** -1.8 * np.exp(-gamma / 3.0e4)
    freq = np.logspace(8.0, 20.0, 48)
    _, _, pi_nu = electron_radiation_kernel.get_syn_polarization_selected(
        1,
        1.0e16,
        1.0,
        1,
        gamma,
        dnde,
        freq,
        1.8,
    )
    assert np.all(np.isfinite(pi_nu))
    assert np.max(pi_nu) - np.min(pi_nu) > 1.0e-3


def test_hadronic_polarization_kernel_tracks_curved_spectrum() -> None:
    energy = np.logspace(0.0, 7.0, 96) * constants.para_m_p_gev
    density = energy ** -1.8 * np.exp(-energy / 3.0e4)
    freq = np.logspace(10.0, 28.0, 48)
    pi_nu = hadronic_fortran_module.fs_hadronic_syn_polarization_shell(
        energy,
        density,
        freq,
        constants.para_m_p_gev,
        30.0,
        1.8,
    )
    assert np.all(np.isfinite(pi_nu))
    assert np.max(pi_nu) - np.min(pi_nu) > 1.0e-3


def main() -> None:
    test_on_axis_shock_random_cancels()
    test_off_axis_is_finite_and_smooth()
    test_toroidal_differs_from_shock_random()
    test_reverse_and_hadronic_synchrotron_components()
    test_fortran_polarization_kernel_conserves_intensity()
    test_fortran_polarization_kernel_tracks_curved_spectrum()
    test_hadronic_polarization_kernel_tracks_curved_spectrum()
    print("polarization_smoke: ok")


if __name__ == "__main__":
    main()
