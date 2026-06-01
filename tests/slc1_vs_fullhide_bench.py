from __future__ import annotations

from pathlib import Path
import sys
import time

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.asgard_paths import ASGARD_DOC_DIR
from ASGARD import Fitter, ISM, Model, Observer, ParamDef, Radiation, Scale, Setups, TophatJet, Wind, units, make_empty_obs, make_flux_density_entry

MODE = "high" if "--high" in sys.argv else "quick"
GRID = {
    "quick": {"lc": 24, "spec": 32, "pair": 24, "expo": 12, "fit_steps": 4, "fit_burn": 1, "pix": 24, "gam": 81, "nu": 49, "r": 80, "theta": 80, "tobs": 48},
    "high": {"lc": 100, "spec": 100, "pair": 200, "expo": 200, "fit_steps": 12, "fit_burn": 4, "pix": 64, "gam": 201, "nu": 201, "r": 300, "theta": 300, "tobs": 200},
}[MODE]


OUTPUT_DIR = ASGARD_DOC_DIR
OUTPUT_PNG = OUTPUT_DIR / "fullhide_1d_charint_1d_benchmark_compare.png"
OUTPUT_PDF = OUTPUT_DIR / "fullhide_1d_charint_1d_benchmark_compare.pdf"


def _build_solver_model(solver: str) -> Model:
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


def _cases_for_solver(solver: str) -> list[dict]:
    def case_quickstart():
        model = _build_solver_model(solver)
        value = model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total
        assert value.shape == (1,)
        assert np.all(np.isfinite(value))
        return {"shape": list(value.shape), "value": float(value[0])}

    def case_lightcurve_grid():
        model = _build_solver_model(solver)
        times = np.logspace(2.0, 8.0, GRID["lc"])
        bands = np.array([1.0e9, 1.0e14, 1.0e18])
        grid = model.flux_density_grid(times, bands).total
        assert grid.shape == (3, GRID["lc"])
        assert np.all(np.isfinite(grid))
        return {"shape": list(grid.shape), "peak": float(np.max(grid))}

    def case_spectrum_grid():
        model = _build_solver_model(solver)
        times = np.array([1.0e3, 1.0e5, 1.0e7])
        freqs = np.logspace(9.0, 22.0, GRID["spec"])
        grid = model.flux_density_grid(times, freqs).total
        assert grid.shape == (GRID["spec"], 3)
        assert np.all(np.isfinite(grid))
        return {"shape": list(grid.shape), "peak": float(np.max(grid))}

    def case_pairs():
        model = _build_solver_model(solver)
        times = np.logspace(2.0, 8.0, GRID["pair"])
        freqs = np.logspace(9.0, 18.0, GRID["pair"])
        values = model.flux_density(times, freqs).total
        assert values.shape == (GRID["pair"],)
        assert np.all(np.isfinite(values))
        return {"shape": list(values.shape), "median": float(np.median(values))}

    def case_exposures():
        model = _build_solver_model(solver)
        times = np.logspace(2.0, 6.0, GRID["expo"])
        freqs = np.full(times.shape, 1.0e14)
        exposures = 0.2 * times
        values = model.flux_density_exposures(times, freqs, exposures).total
        assert values.shape == (GRID["expo"],)
        assert np.all(np.isfinite(values))
        return {"shape": list(values.shape), "median": float(np.median(values))}

    def case_details():
        model = _build_solver_model(solver)
        details = model.details(t_min=1.0e2, t_max=1.0e6)
        assert details.fwd.t_obs.ndim == 1
        assert np.all(np.isfinite(details.fwd.t_obs))
        assert np.all(np.isfinite(details.fwd.nu_m))
        return {"n_t": int(details.fwd.t_obs.shape[0]), "n_patches": len(details.patches)}

    def case_fitter_cfg():
        truth_model = Model(
            medium=Wind(A_star=0.1, n0=1.0),
            jet=TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0),
            observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
            fwd_rad=Radiation(0.1, 1.0e-3, 2.3),
            setups=Setups(electron_solver=solver),
        )
        times = np.logspace(2.0, 4.0, 8)
        freqs = np.full(times.shape, 1.0e14)
        flux = truth_model.flux_density(times, freqs).total
        data = make_empty_obs()
        data["flux_density"].append(make_flux_density_entry(times_s=times, frequencies_hz=freqs, flux=flux, flux_err=0.05 * np.maximum(flux, 1.0e-99)))
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
                ParamDef("E_iso", 52.0, 52.0, Scale.LOG),
                ParamDef("p", 2.3, 2.3, Scale.FIXED),
            ],
            nsteps=GRID["fit_steps"],
            nburn=GRID["fit_burn"],
        )
        assert np.isfinite(result.best_loglike)
        return {"best_loglike": float(result.best_loglike), "best_params": result.best_params}

    def case_sky_image():
        model = _build_solver_model(solver)
        t_obs = np.array([1.0e6])
        nu_obs = 1.0e9
        img = model.sky_image(t_obs, nu_obs=nu_obs, fov=500.0 * units.uas, npixel=GRID["pix"])
        flux_from_image = img.image.sum(axis=(1, 2)) * img.pixel_solid_angle
        flux_direct = model.flux_density_grid(t_obs, np.array([nu_obs])).total[0, :]
        ratio = flux_from_image / flux_direct
        assert img.image.shape == (1, GRID["pix"], GRID["pix"])
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
    results = []
    total = len(cases)
    for idx, (name, fn) in enumerate(cases, start=1):
        results.append(_run_case(f"[{idx}/{total}] {solver}:{name}", lambda fn=fn, solver=solver: fn()))
    return results


def _plot_compare(results: dict[str, list[dict]]) -> None:
    names = [item["name"] for item in results["fullhide_1d"]]
    fullhide_seconds = np.array([item["seconds"] for item in results["fullhide_1d"]], dtype=float)
    charint_seconds = np.array([item["seconds"] for item in results["charint_1d"]], dtype=float)
    x = np.arange(len(names), dtype=float)
    width = 0.38

    fig, ax = plt.subplots(figsize=(12, 5.5), constrained_layout=True)
    ax.bar(x - width / 2.0, fullhide_seconds, width=width, label="fullhide_1d")
    ax.bar(x + width / 2.0, charint_seconds, width=width, label="charint_1d")
    ax.set_xticks(x, names, rotation=20, ha="right")
    ax.set_ylabel("seconds")
    ax.set_title("ASGARD benchmark: fullhide_1d vs charint_1d")
    ax.grid(axis="y", alpha=0.3)
    ax.legend()

    for xi, s0, s1 in zip(x, fullhide_seconds, charint_seconds):
        if s1 > 0.0:
            ax.text(xi + width / 2.0, s1, f"x{s0/s1:.2f}", ha="center", va="bottom", fontsize=8)

    fig.savefig(OUTPUT_PNG, dpi=200)
    fig.savefig(OUTPUT_PDF)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {"fullhide_1d": _cases_for_solver("fullhide_1d"), "charint_1d": _cases_for_solver("charint_1d")}

    _plot_compare(results)

    print(OUTPUT_PNG)
    print(OUTPUT_PDF)

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
