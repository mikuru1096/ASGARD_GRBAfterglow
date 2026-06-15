from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import (
    Hadronic,
    Model,
    Numerics,
    Observer,
    ObserverGrid,
    Radiation,
    ReverseShock,
    SolverOptions,
    UniformMedium,
    WindMedium,
    top_hat_jet,
)
from tests._bench_common import run_case

MODE = "high" if "--high" in sys.argv else "quick"
GRID = {
    "quick": {"lc": 24, "spec": 32, "pair": 24, "expo": 12, "gam": 81, "nu": 49, "r": 80, "theta": 80, "tobs": 48},
    "high": {"lc": 100, "spec": 100, "pair": 200, "expo": 200, "gam": 201, "nu": 201, "r": 300, "theta": 300, "tobs": 200},
}[MODE]
SOLVERS = ("fullhide_1d", "fullhide_2d", "charint_1d", "charint_2d")


def _radiation(**updates) -> Radiation:
    values = dict(
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
    )
    values.update(updates)
    return Radiation(**values)


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


def _solver_options(solver: str) -> SolverOptions:
    return SolverOptions(
        electron_solver=solver,
        dynamics_solver="forward_legacy",
        geometry_projection="sed_legacy",
        electron_photon_coupling="separated",
        ssc_cooling_mode="nakar_y_thomson",
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
        fullhide2d_transport_model="legacy",
        fullhide2d_stochastic_accel_norm=0.0,
        fullhide2d_escape_mode="closed",
    )


def _reverse_shock(enabled: bool = False) -> ReverseShock:
    return ReverseShock(
        enabled=enabled,
        shell_duration_s=10.0,
        upstream_sigma=0.0,
        include_cross_zone_ic=False,
        include_ssc=False,
    )


def _hadronic(enabled: bool = False) -> Hadronic:
    return Hadronic(
        enabled=enabled,
        solver="legacy_1d",
        num_proton_gamma=41 if enabled else 161,
        num_neutrino_frequency=31 if enabled else 121,
        pgamma_scheme="disabled",
        pair_cascade_iterations=1,
    )


def _build_readme_model(solver: str) -> Model:
    return Model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=300.0,
            opening_angle_rad=0.1,
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        medium=UniformMedium(number_density_cm3=1.0),
        observer=Observer(luminosity_distance_cm=1.0e26, z=0.1, viewing_angle_rad=0.0, viewing_azimuth_rad=0.0),
        fwd_rad=_radiation(),
        numerics=_numerics(solver),
        observer_grid=ObserverGrid(time_min_s=1.0e2, time_max_s=1.0e8),
        solver_options=_solver_options(solver),
        reverse_shock=_reverse_shock(False),
        hadronic=_hadronic(False),
    )


def case_quickstart(solver: str):
    model = _build_readme_model(solver)
    value = model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total
    assert value.shape == (1,)
    assert np.all(np.isfinite(value))
    return {"solver": solver, "shape": list(value.shape), "value": float(value[0])}


def case_lightcurve_grid(solver: str):
    model = _build_readme_model(solver)
    times = np.logspace(2.0, 8.0, GRID["lc"])
    bands = np.array([1.0e9, 1.0e14, 1.0e18])
    grid = model.flux_density_grid(times, bands).total
    assert grid.shape == (3, GRID["lc"])
    assert np.all(np.isfinite(grid))
    return {"solver": solver, "shape": list(grid.shape), "peak": float(np.max(grid))}


def case_spectrum_grid(solver: str):
    model = _build_readme_model(solver)
    times = np.array([1.0e3, 1.0e5, 1.0e7])
    freqs = np.logspace(9.0, 22.0, GRID["spec"])
    grid = model.flux_density_grid(times, freqs).total
    assert grid.shape == (GRID["spec"], 3)
    assert np.all(np.isfinite(grid))
    return {"solver": solver, "shape": list(grid.shape), "peak": float(np.max(grid))}


def case_pairs(solver: str):
    model = _build_readme_model(solver)
    times = np.logspace(2.0, 8.0, GRID["pair"])
    freqs = np.logspace(9.0, 18.0, GRID["pair"])
    values = model.flux_density(times, freqs).total
    assert values.shape == (GRID["pair"],)
    assert np.all(np.isfinite(values))
    return {"solver": solver, "shape": list(values.shape), "median": float(np.median(values))}


def case_exposures(solver: str):
    model = _build_readme_model(solver)
    times = np.logspace(2.0, 6.0, GRID["expo"])
    freqs = np.full(times.shape, 1.0e14)
    exposures = 0.2 * times
    values = model.flux_density_exposures(times, freqs, exposures).total
    assert values.shape == (GRID["expo"],)
    assert np.all(np.isfinite(values))
    return {"solver": solver, "shape": list(values.shape), "median": float(np.median(values))}


def case_new_public_configs():
    model = Model(
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
        reverse_shock=ReverseShock(
            enabled=True,
            shell_duration_s=10.0,
            upstream_sigma=0.0,
            include_cross_zone_ic=False,
            include_ssc=False,
        ),
        hadronic=_hadronic(True),
    )
    value = model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total
    assert value.shape == (1,)
    assert np.all(np.isfinite(value))
    return {"shape": list(value.shape), "value": float(value[0])}


def _finite_relmax(lhs: np.ndarray, rhs: np.ndarray) -> float:
    lhs = np.asarray(lhs, dtype=float)
    rhs = np.asarray(rhs, dtype=float)
    mask = np.isfinite(lhs) & np.isfinite(rhs)
    if not np.any(mask):
        return float("nan")
    scale = np.maximum(np.abs(lhs[mask]), 1.0e-300)
    return float(np.max(np.abs(lhs[mask] - rhs[mask]) / scale))


def main() -> None:
    cases = [
        ("quickstart", case_quickstart),
        ("lightcurve_grid", case_lightcurve_grid),
        ("spectrum_grid", case_spectrum_grid),
        ("pair_points", case_pairs),
        ("exposures", case_exposures),
    ]
    results = []
    total = len(SOLVERS) * len(cases)
    done = 0
    for solver in SOLVERS:
        for name, fn in cases:
            done += 1
            label = f"[{done}/{total}] {solver}:{name}"
            results.append(run_case(label, lambda fn=fn, solver=solver: fn(solver)))
    results.append(run_case("[extra] new_public_configs", case_new_public_configs))

    failed = [item for item in results if item["status"] == "FAIL"]
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
