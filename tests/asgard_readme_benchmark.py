from __future__ import annotations

from pathlib import Path
import json
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, TophatJet


def _run_case(name, fn):
    t0 = time.perf_counter()
    try:
        payload = fn()
        return {"name": name, "status": "PASS", "seconds": time.perf_counter() - t0, "payload": payload}
    except NotImplementedError as exc:
        return {"name": name, "status": "UNSUPPORTED", "seconds": time.perf_counter() - t0, "payload": str(exc)}
    except Exception as exc:
        return {"name": name, "status": "FAIL", "seconds": time.perf_counter() - t0, "payload": f"{type(exc).__name__}: {exc}"}


def _build_readme_model() -> Model:
    return Model(
        TophatJet(0.1, 1.0e52, 300.0),
        ISM(1.0),
        Observer(1.0e26, 0.1, 0.0),
        Radiation(0.1, 1.0e-3, 2.3),
    )


def case_quickstart():
    model = _build_readme_model()
    value = model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total
    assert value.shape == (1,)
    assert np.all(np.isfinite(value))
    return {"shape": list(value.shape), "value": float(value[0])}


def case_lightcurve_grid():
    model = _build_readme_model()
    times = np.logspace(2.0, 8.0, 100)
    bands = np.array([1.0e9, 1.0e14, 1.0e18])
    grid = model.flux_density_grid(times, bands).total
    assert grid.shape == (3, 100)
    assert np.all(np.isfinite(grid))
    return {"shape": list(grid.shape), "peak": float(np.max(grid))}


def case_spectrum_grid():
    model = _build_readme_model()
    times = np.array([1.0e3, 1.0e5, 1.0e7])
    freqs = np.logspace(9.0, 22.0, 100)
    grid = model.flux_density_grid(times, freqs).total
    assert grid.shape == (100, 3)
    assert np.all(np.isfinite(grid))
    return {"shape": list(grid.shape), "peak": float(np.max(grid))}


def case_pairs():
    model = _build_readme_model()
    times = np.logspace(2.0, 8.0, 128)
    freqs = np.logspace(9.0, 18.0, 128)
    values = model.flux_density(times, freqs).total
    assert values.shape == (128,)
    assert np.all(np.isfinite(values))
    return {"shape": list(values.shape), "median": float(np.median(values))}


def case_exposures():
    model = _build_readme_model()
    times = np.logspace(2.0, 6.0, 32)
    freqs = np.full(times.shape, 1.0e14)
    exposures = 0.2 * times
    values = model.flux_density_exposures(times, freqs, exposures).total
    assert values.shape == (32,)
    assert np.all(np.isfinite(values))
    return {"shape": list(values.shape), "median": float(np.median(values))}


def main() -> None:
    cases = [
        ("quickstart", case_quickstart),
        ("lightcurve_grid", case_lightcurve_grid),
        ("spectrum_grid", case_spectrum_grid),
        ("pair_points", case_pairs),
        ("exposures", case_exposures),
    ]
    results = [_run_case(name, fn) for name, fn in cases]
    print(json.dumps(results, indent=2))

    failed = [item for item in results if item["status"] == "FAIL"]
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
