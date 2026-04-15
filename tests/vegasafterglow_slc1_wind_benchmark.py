from __future__ import annotations

from pathlib import Path
import json
import sys
import time

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from VegasAfterglow import Fitter, Model, ObsData, Observer, ParamDef, Radiation, Scale, Setups, TophatJet, Wind, units


OUTPUT_DIR = ROOT / "output" / "vegasafterglow_doc"
OUTPUT_JSON = OUTPUT_DIR / "slc1_wind_benchmark_compare.json"
OUTPUT_PNG = OUTPUT_DIR / "slc1_wind_benchmark_compare.png"
OUTPUT_PDF = OUTPUT_DIR / "slc1_wind_benchmark_compare.pdf"


def _build_wind_model(solver: str) -> Model:
    return Model(
        TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0),
        Wind(A_star=0.1, n0=1.0),
        Observer(1.0e26, 0.1, 0.0),
        Radiation(0.1, 1.0e-3, 2.3),
        setups=Setups(electron_solver=solver),
    )


def _run_case(name, fn):
    t0 = time.perf_counter()
    try:
        payload = fn()
        return {"name": name, "status": "PASS", "seconds": time.perf_counter() - t0, "payload": payload}
    except NotImplementedError as exc:
        return {"name": name, "status": "UNSUPPORTED", "seconds": time.perf_counter() - t0, "payload": str(exc)}
    except Exception as exc:
        return {"name": name, "status": "FAIL", "seconds": time.perf_counter() - t0, "payload": f"{type(exc).__name__}: {exc}"}


def _cases_for_solver(solver: str) -> list[dict]:
    def case_quickstart():
        model = _build_wind_model(solver)
        value = model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total
        assert value.shape == (1,)
        assert np.all(np.isfinite(value))
        return {"shape": list(value.shape), "value": float(value[0])}

    def case_lightcurve_grid():
        model = _build_wind_model(solver)
        times = np.logspace(2.0, 8.0, 100)
        bands = np.array([1.0e9, 1.0e14, 1.0e18])
        grid = model.flux_density_grid(times, bands).total
        assert grid.shape == (3, 100)
        assert np.all(np.isfinite(grid))
        return {"shape": list(grid.shape), "peak": float(np.max(grid))}

    def case_spectrum_grid():
        model = _build_wind_model(solver)
        times = np.array([1.0e3, 1.0e5, 1.0e7])
        freqs = np.logspace(9.0, 22.0, 100)
        grid = model.flux_density_grid(times, freqs).total
        assert grid.shape == (100, 3)
        assert np.all(np.isfinite(grid))
        return {"shape": list(grid.shape), "peak": float(np.max(grid))}

    def case_pairs():
        model = _build_wind_model(solver)
        times = np.logspace(2.0, 8.0, 128)
        freqs = np.logspace(9.0, 18.0, 128)
        values = model.flux_density(times, freqs).total
        assert values.shape == (128,)
        assert np.all(np.isfinite(values))
        return {"shape": list(values.shape), "median": float(np.median(values))}

    def case_exposures():
        model = _build_wind_model(solver)
        times = np.logspace(2.0, 6.0, 32)
        freqs = np.full(times.shape, 1.0e14)
        exposures = 0.2 * times
        values = model.flux_density_exposures(times, freqs, exposures).total
        assert values.shape == (32,)
        assert np.all(np.isfinite(values))
        return {"shape": list(values.shape), "median": float(np.median(values))}

    def case_details():
        model = _build_wind_model(solver)
        details = model.details(t_min=1.0e2, t_max=1.0e6)
        assert details.fwd.t_obs.ndim == 1
        assert np.all(np.isfinite(details.fwd.t_obs))
        assert np.all(np.isfinite(details.fwd.nu_m))
        return {"n_t": int(details.fwd.t_obs.shape[0]), "n_patches": len(details.patches)}

    def case_fitter_cfg():
        truth_model = _build_wind_model(solver)
        times = np.logspace(2.0, 4.0, 8)
        freqs = np.full(times.shape, 1.0e14)
        flux = truth_model.flux_density(times, freqs).total
        data = ObsData()
        data.add_flux_density(t=times, nu=freqs, f_nu=flux, err=0.05 * np.maximum(flux, 1.0e-99))
        fitter = Fitter(
            data=data,
            z=0.1,
            lumi_dist=1.0e26,
            jet="tophat",
            medium="wind",
            E_iso=1.0e52,
            Gamma0=300.0,
            theta_c=0.1,
            A_star=0.1,
            n0=1.0,
            eps_e=0.1,
            eps_B=1.0e-3,
            p=2.3,
            electron_solver=solver,
        )
        result = fitter.fit(
            [
                ParamDef("E_iso", 52.0, 52.0, Scale.log),
                ParamDef("p", 2.3, 2.3, Scale.fixed),
            ],
            nsteps=12,
            nburn=4,
        )
        assert np.isfinite(result.best_loglike)
        return {"best_loglike": float(result.best_loglike), "best_params": result.best_params}

    def case_sky_image():
        model = _build_wind_model(solver)
        t_obs = np.array([1.0e6])
        nu_obs = 1.0e9
        img = model.sky_image(t_obs, nu_obs=nu_obs, fov=500.0 * units.uas, npixel=64)
        flux_from_image = img.image.sum(axis=(1, 2)) * img.pixel_solid_angle
        flux_direct = model.flux_density_grid(t_obs, np.array([nu_obs])).total[0, :]
        ratio = flux_from_image / flux_direct
        assert img.image.shape == (1, 64, 64)
        assert np.all(np.isfinite(img.image))
        assert np.isfinite(ratio[0])
        return {"shape": list(img.image.shape), "ratio": float(ratio[0])}

    cases = [
        ("quickstart", case_quickstart),
        ("lightcurve_grid", case_lightcurve_grid),
        ("spectrum_grid", case_spectrum_grid),
        ("pair_points", case_pairs),
        ("exposures", case_exposures),
        ("details", case_details),
        ("fitter_cfg", case_fitter_cfg),
        ("sky_image", case_sky_image),
    ]
    return [_run_case(name, fn) for name, fn in cases]


def _plot(results: dict[str, list[dict]]) -> None:
    names = [item["name"] for item in results["fullhide"]]
    fullhide_seconds = np.array([item["seconds"] for item in results["fullhide"]], dtype=float)
    slc1_seconds = np.array([item["seconds"] for item in results["slc1"]], dtype=float)
    x = np.arange(len(names), dtype=float)
    width = 0.38

    fig, ax = plt.subplots(figsize=(12, 5.5), constrained_layout=True)
    ax.bar(x - width / 2.0, fullhide_seconds, width=width, label="fullhide")
    ax.bar(x + width / 2.0, slc1_seconds, width=width, label="slc1")
    ax.set_xticks(x, names, rotation=20, ha="right")
    ax.set_ylabel("seconds")
    ax.set_title("ASGARD wind benchmark: fullhide vs slc1")
    ax.grid(axis="y", alpha=0.3)
    ax.legend()

    for xi, s0, s1 in zip(x, fullhide_seconds, slc1_seconds):
        if s1 > 0.0:
            ax.text(xi + width / 2.0, s1, f"x{s0/s1:.2f}", ha="center", va="bottom", fontsize=8)

    fig.savefig(OUTPUT_PNG, dpi=200)
    fig.savefig(OUTPUT_PDF)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {
        "fullhide": _cases_for_solver("fullhide"),
        "slc1": _cases_for_solver("slc1"),
    }

    with OUTPUT_JSON.open("w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False)

    _plot(results)

    print(OUTPUT_JSON)
    print(OUTPUT_PNG)
    print(OUTPUT_PDF)
    print(json.dumps(results, indent=2, ensure_ascii=False))

    failed = [
        (solver, item["name"], item["payload"])
        for solver, solver_results in results.items()
        for item in solver_results
        if item["status"] == "FAIL"
    ]
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
