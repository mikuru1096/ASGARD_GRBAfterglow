from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import Ejecta, ISM, Medium, Model, Observer, Radiation, Setups, TophatJet, Wind
from asgard_core.asgard_config import FitConfig
from asgard_core.asgard_runtime import solve_hadronic


def _small_setups(**updates) -> Setups:
    values = {
        "num_threads": 1,
        "num_gam_e": 24,
        "num_nu": 24,
        "num_r": 24,
        "num_theta": 6,
        "num_tobs": 12,
        "observer_time_min_s": 1.0e2,
        "observer_time_max_s": 1.0e4,
        "electron_solver": "fullhide_1d",
    }
    values.update(updates)
    return Setups(**values)


def _base_model(*, medium=None, jet=None, radiation=None, setups=None) -> Model:
    return Model(
        jet=jet if jet is not None else TophatJet(theta_c=0.1, E_iso=1.0e51, Gamma0=80.0),
        medium=medium if medium is not None else ISM(1.0),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.08),
        fwd_rad=radiation if radiation is not None else Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.3),
        setups=setups if setups is not None else _small_setups(),
    )


def _expect_not_implemented(label: str, action) -> None:
    try:
        action()
    except NotImplementedError:
        return
    raise AssertionError(f"{label} should raise NotImplementedError.")


def _test_custom_medium_rejected() -> None:
    custom = Medium(lambda phi, theta, radius: np.ones_like(np.asarray(radius, dtype=float)), kind="custom")
    model = _base_model(medium=custom)
    _expect_not_implemented(
        "custom Medium kernel dispatch",
        lambda: model.flux_density_grid(np.array([1.0e3]), np.array([1.0e10])),
    )


def _test_jet_spreading_rejected() -> None:
    model = _base_model(jet=TophatJet(theta_c=0.1, E_iso=1.0e51, Gamma0=80.0, spreading=True))
    _expect_not_implemented(
        "jet spreading dynamics",
        lambda: model.flux_density_grid(np.array([1.0e3]), np.array([1.0e10])),
    )


def _test_wind_k_rejected() -> None:
    _expect_not_implemented("wind k != 2", lambda: Wind(Astar=0.1, k=1.5))


def _test_thermal_solver_rejected() -> None:
    model = _base_model(
        radiation=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.3, thermal_electrons=True),
        setups=_small_setups(electron_solver="weno5_1d"),
    )
    _expect_not_implemented(
        "thermal electrons outside fullhide_1d",
        lambda: model.flux_density_grid(np.array([1.0e3]), np.array([1.0e10])),
    )


def _test_nonaxisymmetric_toroidal_rejected() -> None:
    jet = Ejecta(
        E_iso=lambda phi, theta: 1.0e51 * (1.0 + 0.1 * np.cos(phi)),
        Gamma0=lambda phi, theta: 80.0,
        theta_max=0.2,
    )
    model = _base_model(jet=jet)
    _expect_not_implemented(
        "toroidal polarization on non-axisymmetric ejecta",
        lambda: model.polarization(np.array([1.0e3]), np.array([1.0e10]), magnetic_geometry="toroidal"),
    )


def _test_chi_hadronic_rejected() -> None:
    config = FitConfig(electron_solver="fullhide_2d")
    config.hadronic.enabled = True
    config.hadronic.epsilon_p = 0.3
    config.hadronic.include_proton_synch = True
    config.hadronic.num_gam_p = 16
    config.hadronic.num_nu_nu = 8
    _expect_not_implemented(
        "2D/chi-resolved hadronic transport",
        lambda: solve_hadronic(None, None, None, np.array([1.0e10]), np.zeros((1, 1)), config),
    )


def main() -> None:
    _test_custom_medium_rejected()
    _test_jet_spreading_rejected()
    _test_wind_k_rejected()
    _test_thermal_solver_rejected()
    _test_nonaxisymmetric_toroidal_rejected()
    _test_chi_hadronic_rejected()
    print("public-backend-limits-smoke-ok")


if __name__ == "__main__":
    main()
