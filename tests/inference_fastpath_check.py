from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import Fitter, ISM, Model, ObsData, Observer, ParamDef, Radiation, Scale, Setups, TophatJet
from ASGARD.api import _compute_total_matrix
from asgard_inference import compile_inference_problem, evaluate_compiled_loglike


def _legacy_loglike(fitter: Fitter, values: dict[str, float]) -> float:
    model = fitter.build_model(values)
    loglike = 0.0
    for dataset in fitter.data.flux_density:
        pred = model.flux_density(dataset["times_s"], dataset["frequencies_hz"]).total
        resid = (pred - dataset["flux"]) / dataset["flux_err"]
        loglike -= 0.5 * float(np.sum(resid * resid))
    for dataset in fitter.data.spectrum:
        pred = model.spectrum(dataset["time_s"], dataset["frequencies_hz"])
        resid = (pred - dataset["flux"]) / dataset["flux_err"]
        loglike -= 0.5 * float(np.sum(resid * resid))
    for dataset in fitter.data.flux:
        pred = model.flux(
            dataset["time_s"],
            dataset["nu_min_hz"],
            dataset["nu_max_hz"],
            num_points=dataset["num_points"],
        )
        resid = (pred - dataset["flux"]) / dataset["flux_err"]
        loglike -= 0.5 * float(resid * resid)
    return loglike


def main() -> None:
    setups = Setups(num_threads=2, num_r=120, num_theta=120, num_tobs=48, num_gam_e=81, num_nu=81, fwd_ssc=True)
    truth_model = Model(
        medium=ISM(0.2),
        jet=TophatJet(E_iso=3.0e52, Gamma0=150.0, theta_c=0.12),
        observer=Observer(z=0.2, theta_obs=0.02, lumi_dist=3.0e27),
        fwd_rad=Radiation(eps_e=0.08, eps_B=5.0e-4, p=2.4, xi_N=0.15, ssc=True),
        setups=setups,
    )

    pair_times = np.logspace(2.0, 4.0, 6)
    pair_freqs = np.array([1.0e9, 3.0e9, 1.0e14, 1.0e15, 1.0e17, 1.0e18], dtype=float)
    pair_flux = truth_model.flux_density(pair_times, pair_freqs).total
    spectrum_freqs = np.logspace(9.0, 18.0, 12)
    spectrum_flux = truth_model.spectrum(2.0e4, spectrum_freqs)
    band_flux = float(truth_model.flux(5.0e4, 1.0e14, 5.0e14, num_points=24))

    data = ObsData()
    data.add_flux_density(pair_times, pair_freqs, pair_flux, 0.05 * np.maximum(pair_flux, 1.0e-99))
    data.add_spectrum(2.0e4, spectrum_freqs, spectrum_flux, 0.05 * np.maximum(spectrum_flux, 1.0e-99))
    data.add_flux(1.0e14, 5.0e14, 5.0e4, band_flux, max(0.05 * abs(band_flux), 1.0e-99), num_points=24)

    fitter = Fitter(
        model=truth_model,
        data=data,
        params=[
            ParamDef("log10_Eiso", "jet.E_iso", 52.0, 53.0, Scale.LOG10),
            ParamDef("p", "fwd_rad.p", 2.2, 2.6, Scale.LINEAR),
        ],
    )
    values = {"log10_Eiso": np.log10(3.0e52), "p": 2.4}

    compiled_problem = compile_inference_problem(data, truth_model, params=fitter.params)
    compiled_loglike = evaluate_compiled_loglike(compiled_problem, values)
    legacy_loglike = _legacy_loglike(fitter, values)
    if not np.isclose(compiled_loglike, legacy_loglike, rtol=5.0e-7, atol=5.0e-7):
        raise AssertionError(f"compiled_loglike={compiled_loglike}, legacy_loglike={legacy_loglike}")

    compare_times = np.logspace(2.0, 5.0, 10)
    compare_freqs = np.array([1.0e9, 1.0e14, 1.0e18], dtype=float)
    full_total = truth_model.flux_density_grid(compare_times, compare_freqs).total
    fast_total = _compute_total_matrix(truth_model, compare_times, compare_freqs)
    if not np.allclose(full_total, fast_total, rtol=5.0e-7, atol=5.0e-12):
        raise AssertionError("total_only observer path deviates from the legacy total flux grid.")

    print("PASS: inference fast-path matches legacy loglike and total flux grids.")


if __name__ == "__main__":
    main()
