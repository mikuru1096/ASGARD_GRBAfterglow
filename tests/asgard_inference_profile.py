from __future__ import annotations

import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_paths import ASGARD_DOC_DIR
from ASGARD import Fitter, ISM, Model, ObsData, Observer, ParamDef, Radiation, Scale, Setups, TophatJet
from asgard_inference import compile_inference_problem, evaluate_compiled_loglike


OUTPUT_JSON = ASGARD_DOC_DIR / "inference_fastpath_profile.json"
NUM_EVALS = 4


def _build_case(ssc: bool) -> tuple[Fitter, list[dict[str, float]]]:
    model = Model(
        medium=ISM(0.3),
        jet=TophatJet(E_iso=2.0e52, Gamma0=120.0, theta_c=0.1),
        observer=Observer(z=0.3, theta_obs=0.01, lumi_dist=2.0e27),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.35, xi_N=0.12, ssc=ssc),
        setups=Setups(num_threads=2, num_r=120, num_theta=120, num_tobs=48, num_gam_e=81, num_nu=81, fwd_ssc=ssc),
    )
    times = np.logspace(2.0, 4.0, 8)
    freqs = np.array([1.0e9, 3.0e9, 1.0e10, 3.0e10, 1.0e14, 3.0e14, 1.0e17, 1.0e18], dtype=float)
    flux = model.flux_density(times, freqs).total
    data = ObsData()
    data.add_flux_density(times, freqs, flux, 0.05 * np.maximum(flux, 1.0e-99))
    fitter = Fitter(
        model=model,
        data=data,
        params=[
            ParamDef("log10_Eiso", "jet.E_iso", 52.0, 53.5, Scale.LOG10),
            ParamDef("p", "fwd_rad.p", 2.1, 2.6, Scale.LINEAR),
        ],
    )
    values = [
        {"log10_Eiso": 52.2, "p": 2.25},
        {"log10_Eiso": 52.5, "p": 2.30},
        {"log10_Eiso": 52.8, "p": 2.35},
        {"log10_Eiso": 53.0, "p": 2.40},
    ]
    return fitter, values


def _legacy_loglike(fitter: Fitter, values: dict[str, float]) -> float:
    model = fitter.build_model(values)
    loglike = 0.0
    for dataset in fitter.data.flux_density:
        pred = model.flux_density(dataset["times_s"], dataset["frequencies_hz"]).total
        resid = (pred - dataset["flux"]) / dataset["flux_err"]
        loglike -= 0.5 * float(np.sum(resid * resid))
    return loglike


def _profile_case(ssc: bool) -> dict[str, float]:
    fitter, values = _build_case(ssc)
    problem = compile_inference_problem(fitter.data, fitter.model, params=fitter.params)

    t0 = perf_counter()
    legacy = [_legacy_loglike(fitter, value) for value in values[:NUM_EVALS]]
    legacy_seconds = perf_counter() - t0

    timings: dict[str, float] = {}
    t1 = perf_counter()
    fast = [evaluate_compiled_loglike(problem, value, timings=timings) for value in values[:NUM_EVALS]]
    fast_seconds = perf_counter() - t1

    if not np.allclose(legacy, fast, rtol=5.0e-7, atol=5.0e-7):
        raise AssertionError("Compiled inference profile case diverged from legacy results.")

    fortran_labels = [key for key in timings if key.startswith("Dynamics.") or key.startswith("Electron.") or key.startswith("Radiation.")]
    fortran_seconds = float(sum(timings[key] for key in fortran_labels))
    python_seconds = max(fast_seconds - fortran_seconds, 0.0)
    ratio = python_seconds / fortran_seconds if fortran_seconds > 0.0 else np.inf
    return {
        "legacy_seconds": legacy_seconds,
        "compiled_seconds": fast_seconds,
        "fortran_seconds": fortran_seconds,
        "python_seconds": python_seconds,
        "python_over_fortran": ratio,
    }


def main() -> None:
    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "sync_only": _profile_case(ssc=False),
        "with_ssc": _profile_case(ssc=True),
    }
    with OUTPUT_JSON.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
    print(OUTPUT_JSON)
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
