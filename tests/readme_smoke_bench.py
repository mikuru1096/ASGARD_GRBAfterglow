from __future__ import annotations

from pathlib import Path
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet

MODE = "high" if "--high" in sys.argv else "quick"
GRID = {
    "quick": {"lc": 24, "spec": 32, "pair": 24, "expo": 12, "gam": 81, "nu": 49, "r": 80, "theta": 80, "tobs": 48},
    "high": {"lc": 100, "spec": 100, "pair": 200, "expo": 200, "gam": 201, "nu": 201, "r": 300, "theta": 300, "tobs": 200},
}[MODE]
SOLVERS = ("fullhide", "slc1", "charint")


def _run_case(name, fn):
    print(f"  {name} ...", flush=True)
    t0 = time.perf_counter()
    try:
        payload = fn()
        item = {"name": name, "status": "PASS", "seconds": time.perf_counter() - t0, "payload": payload}
    except NotImplementedError as exc:
        item = {"name": name, "status": "UNSUPPORTED", "seconds": time.perf_counter() - t0, "payload": str(exc)}
    except Exception as exc:
        item = {"name": name, "status": "FAIL", "seconds": time.perf_counter() - t0, "payload": f"{type(exc).__name__}: {exc}"}
    print(f"  {name}: {item['status']} {item['seconds']:.2f}s", flush=True)
    return item


def _build_readme_model(solver: str) -> Model:
    return Model(
        TophatJet(0.1, 1.0e52, 300.0),
        ISM(1.0),
        Observer(1.0e26, 0.1, 0.0),
        Radiation(0.1, 1.0e-3, 2.3),
        setups=Setups(
            electron_solver=solver,
            num_gam_e=GRID["gam"],
            num_nu=GRID["nu"],
            num_r=GRID["r"],
            num_theta=GRID["theta"],
            num_tobs=GRID["tobs"],
            electron_adaptive_substeps=False,
        ),
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
            results.append(_run_case(label, lambda fn=fn, solver=solver: fn(solver)))

    failed = [item for item in results if item["status"] == "FAIL"]
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
