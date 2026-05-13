from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet
from tests._bench_common import run_case


MODE = "high" if "--high" in sys.argv else "quick"
GRID = {
    "quick": {"times": 8, "gam": 32, "nu": 36, "r": 32, "tobs": 32, "patch_theta": 3, "patch_phi": 12},
    "high": {"times": 32, "gam": 81, "nu": 81, "r": 96, "tobs": 96, "patch_theta": 5, "patch_phi": 24},
}[MODE]


def _setups(*, reverse: bool = False, hadronic: bool = False) -> Setups:
    return Setups(
        electron_solver="fullhide_1d",
        num_threads=1,
        num_gam_e=GRID["gam"],
        num_gam_p=GRID["gam"],
        num_nu=GRID["nu"],
        num_nu_nu=16,
        num_r=GRID["r"],
        num_tobs=GRID["tobs"],
        num_theta=8,
        patch_theta=GRID["patch_theta"],
        patch_phi=GRID["patch_phi"],
        rvs_shock=reverse,
        reverse_delta_t_s=10.0,
        hadronic_enabled=hadronic,
        hadronic_solver="legacy_1d",
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
            epsilon_p=0.25 if hadronic else 0.0,
            proton_synch=True,
            reverse_epsilon_p=0.25 if rs_hadronic else 0.0,
        ),
        rvs_rad=Radiation(0.1, 1.0e-3, 2.4) if reverse else None,
        setups=_setups(reverse=reverse, hadronic=hadronic),
    )


def _unwrap_pa(pa_rad: np.ndarray) -> np.ndarray:
    return 0.5 * np.unwrap(2.0 * np.asarray(pa_rad, dtype=float))


def _assert_physical_series(pi: np.ndarray, pa_rad: np.ndarray) -> dict[str, float]:
    values = np.asarray(pi, dtype=float)
    pa = _unwrap_pa(pa_rad)
    assert np.all(np.isfinite(values))
    assert np.all(np.isfinite(pa))
    assert np.all(values >= 0.0)
    assert np.all(values <= 1.0)
    max_pi_step = float(np.max(np.abs(np.diff(values)))) if values.size > 1 else 0.0
    max_pa_step = float(np.max(np.abs(np.diff(pa)))) if pa.size > 1 else 0.0
    assert max_pi_step < 0.5
    assert max_pa_step < np.pi / 2.0
    return {"max_pi_step": max_pi_step, "max_pa_step_rad": max_pa_step}


def case_on_axis_cancellation() -> dict[str, float]:
    times = np.logspace(3.0, 6.0, GRID["times"])
    pol = _model(0.0).polarization(times, np.array([1.0e14]), magnetic_geometry="shock_random")
    values = pol.linear_polarization[0]
    assert np.all(np.isfinite(pol.I_sync))
    assert float(np.max(values)) < 1.0e-10
    return {"max_linear_polarization": float(np.max(values))}


def case_off_axis_geometries() -> dict[str, float]:
    times = np.logspace(3.0, 6.0, GRID["times"])
    freqs = np.array([1.0e14])
    model = _model(0.12)
    shock = model.polarization(times, freqs, magnetic_geometry="shock_random")
    toroidal = model.polarization(times, freqs, magnetic_geometry="toroidal")
    shock_metrics = _assert_physical_series(shock.linear_polarization[0], shock.polarization_angle_rad[0])
    toroidal_metrics = _assert_physical_series(toroidal.linear_polarization[0], toroidal.polarization_angle_rad[0])
    stokes_delta = float(np.max(np.abs(toroidal.Q - shock.Q) + np.abs(toroidal.U - shock.U)))
    assert stokes_delta > 0.0
    return {
        "shock_max_pi": float(np.max(shock.linear_polarization[0])),
        "toroidal_max_pi": float(np.max(toroidal.linear_polarization[0])),
        "stokes_delta": stokes_delta,
        **{f"shock_{key}": value for key, value in shock_metrics.items()},
        **{f"toroidal_{key}": value for key, value in toroidal_metrics.items()},
    }


def case_component_paths() -> dict[str, float]:
    times = np.logspace(3.8, 5.2, max(3, GRID["times"] // 2))
    freqs = np.array([1.0e14, 1.0e18])
    pol = _model(0.12, reverse=True, hadronic=True, rs_hadronic=True).polarization(times, freqs)
    assert np.all(np.isfinite(pol.I_sync))
    assert np.any(pol.components["fwd_sync"]["I"] > 0.0)
    assert np.any(pol.components["rev_sync"]["I"] > 0.0)
    assert np.any(pol.components["hadronic_sync"]["I"] > 0.0)
    assert np.any(pol.components["rev_hadronic_sync"]["I"] > 0.0)
    metrics = _assert_physical_series(pol.linear_polarization[0], pol.polarization_angle_rad[0])
    return {
        "max_total_pi": float(np.max(pol.linear_polarization)),
        "max_i_sync": float(np.max(pol.I_sync)),
        **metrics,
    }


def main() -> None:
    cases = [
        run_case("polarization:on_axis_cancellation", case_on_axis_cancellation),
        run_case("polarization:off_axis_geometries", case_off_axis_geometries),
        run_case("polarization:component_paths", case_component_paths),
    ]
    failed = [item for item in cases if item["status"] == "FAIL"]
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
