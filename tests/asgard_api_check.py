from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from VegasAfterglow import GaussianJet, ISM, Medium, Model, Observer, Radiation, Setups, TophatJet, Wind


def main() -> None:
    times_s = np.logspace(2.0, 5.0, 16)
    frequencies_hz = np.array([4.63e14, 9.0e9], dtype=float)
    observer = Observer(z=0.4, theta_obs=0.0)
    direct_setups = Setups(num_threads=8, num_r=80, num_theta=64, num_tobs=24, patch_theta=1, patch_phi=1)
    structured_setups = Setups(num_threads=8, num_r=80, num_theta=64, num_tobs=24, patch_theta=2, patch_phi=3)

    tophat = Model(
        medium=ISM(0.1),
        jet=TophatJet(E_iso=1.0e53, lf=1.0e2, theta_j=1.0e-1),
        observer=observer,
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.5, xi_N=1.0e-1, ssc=True, kn=False),
        setups=direct_setups,
    )
    grid = tophat.flux_density_grid(times_s, frequencies_hz)
    assert grid.shape == (frequencies_hz.shape[0], times_s.shape[0])
    spec = tophat.spectrum(1.0e3, frequencies_hz)
    assert spec.shape == (frequencies_hz.shape[0],)
    details = tophat.details()
    assert details.fwd.t_obs.shape == details.fwd.nu_m.shape

    gaussian = Model(
        medium=ISM(0.1),
        jet=GaussianJet(E_iso=1.0e53, lf=1.0e2, theta_c=5.0e-2, theta_max=2.0e-1),
        observer=Observer(z=0.4, theta_obs=8.0e-2),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.5, xi_N=1.0e-1, ssc=False, kn=False),
        setups=structured_setups,
    )
    gaussian_grid = gaussian.flux_density_grid(times_s, frequencies_hz)
    assert np.all(np.isfinite(gaussian_grid))
    assert gaussian_grid.shape == (frequencies_hz.shape[0], times_s.shape[0])

    custom_ism = Medium(
        rho=lambda _phi, _theta, radius_cm: np.full(np.asarray(radius_cm, dtype=float).shape, 0.1, dtype=float),
        label="custom_ism",
    )
    gaussian_custom_ism = Model(
        medium=custom_ism,
        jet=GaussianJet(E_iso=1.0e53, lf=1.0e2, theta_c=5.0e-2, theta_max=2.0e-1),
        observer=Observer(z=0.4, theta_obs=8.0e-2),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.5, xi_N=1.0e-1, ssc=False, kn=False),
        setups=structured_setups,
    )
    custom_grid = gaussian_custom_ism.flux_density_grid(times_s, frequencies_hz).total
    rel_custom_ism = np.max(np.abs(custom_grid - gaussian_grid.total) / np.maximum(gaussian_grid.total, 1.0e-99))
    assert rel_custom_ism <= 3.0e-2

    custom_wind_floor = Medium(
        rho=lambda _phi, _theta, radius_cm: np.maximum(0.1, 1.0 * 3.0e35 / np.asarray(radius_cm, dtype=float) ** 2),
        label="custom_wind_floor",
    )
    wind_floor_model = Model(
        medium=Wind(A_star=1.0, n_ism=0.1),
        jet=GaussianJet(E_iso=1.0e53, lf=1.0e2, theta_c=5.0e-2, theta_max=2.0e-1),
        observer=Observer(z=0.4, theta_obs=8.0e-2),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.5, xi_N=1.0e-1, ssc=False, kn=False),
        setups=structured_setups,
    )
    custom_wind_floor_model = Model(
        medium=custom_wind_floor,
        jet=GaussianJet(E_iso=1.0e53, lf=1.0e2, theta_c=5.0e-2, theta_max=2.0e-1),
        observer=Observer(z=0.4, theta_obs=8.0e-2),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.5, xi_N=1.0e-1, ssc=False, kn=False),
        setups=structured_setups,
    )
    wind_floor_grid = wind_floor_model.flux_density_grid(times_s, frequencies_hz).total
    custom_wind_floor_grid = custom_wind_floor_model.flux_density_grid(times_s, frequencies_hz).total
    rel_custom_wind_floor = np.max(np.abs(custom_wind_floor_grid - wind_floor_grid) / np.maximum(wind_floor_grid, 1.0e-99))
    assert rel_custom_wind_floor <= 3.0e-2

    custom_wind_inner = Medium(
        rho=lambda _phi, _theta, radius_cm: np.minimum(1.0e3, 1.0 * 3.0e35 / np.asarray(radius_cm, dtype=float) ** 2),
        label="custom_wind_inner",
    )
    wind_inner_model = Model(
        medium=Wind(A_star=1.0, n_ism=0.0, n0=1.0e3),
        jet=GaussianJet(E_iso=1.0e53, lf=1.0e2, theta_c=5.0e-2, theta_max=2.0e-1),
        observer=Observer(z=0.4, theta_obs=8.0e-2),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.5, xi_N=1.0e-1, ssc=False, kn=False),
        setups=structured_setups,
    )
    custom_wind_inner_model = Model(
        medium=custom_wind_inner,
        jet=GaussianJet(E_iso=1.0e53, lf=1.0e2, theta_c=5.0e-2, theta_max=2.0e-1),
        observer=Observer(z=0.4, theta_obs=8.0e-2),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.5, xi_N=1.0e-1, ssc=False, kn=False),
        setups=structured_setups,
    )
    wind_inner_grid = wind_inner_model.flux_density_grid(times_s, frequencies_hz).total
    custom_wind_inner_grid = custom_wind_inner_model.flux_density_grid(times_s, frequencies_hz).total
    rel_custom_wind_inner = np.max(np.abs(custom_wind_inner_grid - wind_inner_grid) / np.maximum(wind_inner_grid, 1.0e-99))
    assert rel_custom_wind_inner <= 3.0e-2

    print("PASS: VegasAfterglow API check succeeded.")


if __name__ == "__main__":
    main()
