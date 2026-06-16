from __future__ import annotations

import sys

from _repo_path import ensure_repo_root_on_path
import numpy as np


ensure_repo_root_on_path()

from asgard_core import (
    Numerics,
    Observer,
    ObserverGrid,
    WindMedium,
    top_hat_jet,
)
from tests.public_api_builders import (
    hadronic as _base_hadronic,
    radiation as _radiation,
    reverse_shock as _reverse_shock,
    solver_options as _base_solver_options,
    top_hat_model as _base_top_hat_model,
)

MODE = "high" if "--high" in sys.argv else "quick"
GRID = {
    "quick": {"lc": 24, "spec": 32, "pair": 24, "expo": 12, "gam": 81, "nu": 49, "r": 80, "theta": 80, "tobs": 48},
    "high": {"lc": 100, "spec": 100, "pair": 200, "expo": 200, "gam": 201, "nu": 201, "r": 300, "theta": 300, "tobs": 200},
}[MODE]
SOLVERS = ("fullhide_1d", "fullhide_2d", "charint_1d", "charint_2d")


def _numerics(solver: str, *, num_phi: int = 1) -> Numerics:
    return Numerics(
        num_electron_gamma=GRID["gam"],
        num_photon_frequency=GRID["nu"],
        num_radius=GRID["r"],
        num_theta=GRID["theta"],
        num_phi=num_phi,
        num_observer_time=GRID["tobs"],
        num_chi=8 if solver.endswith("_2d") else None,
        num_threads=1,
        electron_adaptive_substeps=False,
        electron_substep_rtol=0.02,
        electron_substep_min=100,
        electron_substep_max=1000,
        initial_radius_cm=1.0e14,
    )


def _solver_options(solver: str):
    return _base_solver_options(electron_solver=solver, ssc_cooling_mode="nakar_y_thomson")


def _hadronic(enabled: bool):
    return _base_hadronic(
        enabled=enabled,
        num_proton_gamma=41 if enabled else 161,
        num_neutrino_frequency=31 if enabled else 121,
    )


def _build_readme_model(solver: str):
    return _base_top_hat_model(
        numerics=_numerics(solver),
        observer_grid=ObserverGrid(time_min_s=1.0e2, time_max_s=1.0e8),
        solver_options=_solver_options(solver),
        hadronic=_hadronic(False),
    )


def case_quickstart(solver: str):
    model = _build_readme_model(solver)
    value = model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total
    assert value.shape == (1,)
    assert np.all(np.isfinite(value))


def case_lightcurve_grid(solver: str):
    model = _build_readme_model(solver)
    times = np.logspace(2.0, 8.0, GRID["lc"])
    bands = np.array([1.0e9, 1.0e14, 1.0e18])
    grid = model.flux_density_grid(times, bands).total
    assert grid.shape == (3, GRID["lc"])
    assert np.all(np.isfinite(grid))


def case_spectrum_grid(solver: str):
    model = _build_readme_model(solver)
    times = np.array([1.0e3, 1.0e5, 1.0e7])
    freqs = np.logspace(9.0, 22.0, GRID["spec"])
    grid = model.flux_density_grid(times, freqs).total
    assert grid.shape == (GRID["spec"], 3)
    assert np.all(np.isfinite(grid))


def case_pairs(solver: str):
    model = _build_readme_model(solver)
    times = np.logspace(2.0, 8.0, GRID["pair"])
    freqs = np.logspace(9.0, 18.0, GRID["pair"])
    values = model.flux_density(times, freqs).total
    assert values.shape == (GRID["pair"],)
    assert np.all(np.isfinite(values))


def case_exposures(solver: str):
    model = _build_readme_model(solver)
    times = np.logspace(2.0, 6.0, GRID["expo"])
    freqs = np.full(times.shape, 1.0e14)
    exposures = 0.2 * times
    values = model.flux_density_exposures(times, freqs, exposures).total
    assert values.shape == (GRID["expo"],)
    assert np.all(np.isfinite(values))


def case_new_public_configs():
    model = _base_top_hat_model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=120.0,
            opening_angle_rad=0.12,
            shell_duration_s=50.0,
            magnetar=None,
            spreading=False,
        ),
        medium=WindMedium(a_star=0.1, density_floor_cm3=1.0e-3, density_cap_cm3=10.0),
        observer=Observer(z=0.2, viewing_angle_rad=0.02, viewing_azimuth_rad=0.0, luminosity_distance_cm=None),
        fwd_rad=_radiation(include_ssc=True, proton_energy_fraction=0.01),
        rvs_rad=_radiation(epsilon_B=1.0e-2, p=2.4, accelerated_electron_fraction=1.0),
        numerics=_numerics("fullhide_1d", num_phi=8),
        observer_grid=ObserverGrid(time_min_s=1.0e2, time_max_s=1.0e6),
        solver_options=_solver_options("fullhide_1d"),
        reverse_shock=_reverse_shock(enabled=True),
        hadronic=_hadronic(True),
    )
    value = model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total
    assert value.shape == (1,)
    assert np.all(np.isfinite(value))


def main() -> None:
    cases = [
        ("quickstart", case_quickstart),
        ("lightcurve_grid", case_lightcurve_grid),
        ("spectrum_grid", case_spectrum_grid),
        ("pair_points", case_pairs),
        ("exposures", case_exposures),
    ]
    total = len(SOLVERS) * len(cases)
    done = 0
    for solver in SOLVERS:
        for name, fn in cases:
            done += 1
            label = f"[{done}/{total}] {solver}:{name}"
            print(f"  {label} ...", flush=True)
            fn(solver)
    print("  [extra] new_public_configs ...", flush=True)
    case_new_public_configs()


if __name__ == "__main__":
    main()
