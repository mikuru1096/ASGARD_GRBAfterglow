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

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet, Wind, units


OUTPUT_DIR = ROOT / "output" / "benchmark_exp_tail"
BENCHMARK_JSON = OUTPUT_DIR / "benchmark_compare.json"
BENCHMARK_PNG = OUTPUT_DIR / "benchmark_compare.png"
BENCHMARK_PDF = OUTPUT_DIR / "benchmark_compare.pdf"
MULTIBAND_NPZ = OUTPUT_DIR / "multiband_afterglow_data.npz"
MULTIBAND_PNG = OUTPUT_DIR / "multiband_afterglow_compare.png"
MULTIBAND_PDF = OUTPUT_DIR / "multiband_afterglow_compare.pdf"

SOLVERS = ("fullhide", "slc1", "slc1_mmg2")
MEDIA = ("ism", "wind")
NUM_GAM_BY_SOLVER = {"fullhide": 121, "slc1": 32, "slc1_mmg2": 32}
IC_EPSILON_E = 0.2
IC_EPSILON_B = 1.0e-5
IC_P = 2.3
IC_E_ISO = 3.0e52
IC_GAMMA0 = 250.0
IC_N_ISM = 10.0
IC_A_STAR = 1.0
IC_WIND_N0 = 1.0e-6
BANDS = {
    "Radio": 3.0e9,
    "Optical": 4.6e14,
    "X-ray": 2.417989e17,
    "GeV": 2.417989e23,
    "TeV": 2.417989e26,
}


def _build_model(solver: str, medium_name: str) -> Model:
    jet = TophatJet(theta_c=0.1, E_iso=IC_E_ISO, Gamma0=IC_GAMMA0)
    if medium_name == "ism":
        medium = ISM(n_ism=IC_N_ISM)
    elif medium_name == "wind":
        medium = Wind(A_star=IC_A_STAR, n0=IC_WIND_N0)
    else:
        raise ValueError(medium_name)
    return Model(
        jet=jet,
        medium=medium,
        observer=Observer(1.0e26, 0.1, 0.0),
        fwd_rad=Radiation(IC_EPSILON_E, IC_EPSILON_B, IC_P, ssc=True, kn=True),
        setups=Setups(electron_solver=solver, num_gam_e=NUM_GAM_BY_SOLVER[solver], num_nu=121, num_r=160),
    )


def _tight_log_ylim(ax, arrays: list[np.ndarray], lower_q: float = 0.02, upper_q: float = 0.98) -> None:
    positive = []
    for arr in arrays:
        values = np.asarray(arr, dtype=float)
        values = values[np.isfinite(values) & (values > 0.0)]
        if values.size:
            positive.append(values)
    if not positive:
        return
    merged = np.concatenate(positive)
    y_lo = np.quantile(merged, lower_q)
    y_hi = np.quantile(merged, upper_q)
    if y_lo <= 0.0 or y_hi <= y_lo:
        return
    ax.set_ylim(10.0 ** (np.log10(y_lo) - 0.35), 10.0 ** (np.log10(y_hi) + 0.35))


def _save_figure(fig, png_path: Path, pdf_path: Path, dpi: int) -> None:
    fig.savefig(png_path, dpi=dpi)
    try:
        fig.savefig(pdf_path)
    except PermissionError:
        print(f"skip locked pdf: {pdf_path}", flush=True)


def _run_case(name, fn):
    t0 = time.perf_counter()
    try:
        payload = fn()
        return {"name": name, "status": "PASS", "seconds": time.perf_counter() - t0, "payload": payload}
    except NotImplementedError as exc:
        return {"name": name, "status": "UNSUPPORTED", "seconds": time.perf_counter() - t0, "payload": str(exc)}
    except Exception as exc:
        return {"name": name, "status": "FAIL", "seconds": time.perf_counter() - t0, "payload": f"{type(exc).__name__}: {exc}"}


def _benchmark_cases(solver: str, medium_name: str) -> list[dict]:
    print(f"benchmark {medium_name} {solver}", flush=True)

    def case_quickstart():
        model = _build_model(solver, medium_name)
        value = model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total
        assert value.shape == (1,)
        assert np.all(np.isfinite(value))
        return {"shape": list(value.shape), "value": float(value[0])}

    def case_lightcurve_grid():
        model = _build_model(solver, medium_name)
        times = np.logspace(2.0, 8.0, 60)
        bands = np.array([1.0e9, 1.0e14, 1.0e18])
        grid = model.flux_density_grid(times, bands).total
        assert grid.shape == (3, 60)
        assert np.all(np.isfinite(grid))
        return {"shape": list(grid.shape), "peak": float(np.max(grid))}

    def case_spectrum_grid():
        model = _build_model(solver, medium_name)
        times = np.array([1.0e3, 1.0e5, 1.0e7])
        freqs = np.logspace(9.0, 26.5, 72)
        grid = model.flux_density_grid(times, freqs).total
        assert grid.shape == (72, 3)
        assert np.all(np.isfinite(grid))
        return {"shape": list(grid.shape), "peak": float(np.max(grid))}

    def case_pairs():
        model = _build_model(solver, medium_name)
        times = np.logspace(2.0, 8.0, 64)
        freqs = np.logspace(9.0, 25.0, 64)
        values = model.flux_density(times, freqs).total
        assert values.shape == (64,)
        assert np.all(np.isfinite(values))
        return {"shape": list(values.shape), "median": float(np.median(values))}

    def case_exposures():
        model = _build_model(solver, medium_name)
        times = np.logspace(2.0, 6.0, 24)
        freqs = np.full(times.shape, 1.0e14)
        exposures = 0.2 * times
        values = model.flux_density_exposures(times, freqs, exposures).total
        assert values.shape == (24,)
        assert np.all(np.isfinite(values))
        return {"shape": list(values.shape), "median": float(np.median(values))}

    def case_details():
        model = _build_model(solver, medium_name)
        details = model.details(t_min=1.0e2, t_max=1.0e6)
        assert details.fwd.t_obs.ndim == 1
        assert np.all(np.isfinite(details.fwd.nu_m))
        return {"n_t": int(details.fwd.t_obs.shape[0]), "n_patches": len(details.patches)}

    def case_sky_image():
        model = _build_model(solver, medium_name)
        t_obs = np.array([1.0e6])
        nu_obs = 1.0e9
        img = model.sky_image(t_obs, nu_obs=nu_obs, fov=500.0 * units.uas, npixel=64)
        flux_from_image = img.image.sum(axis=(1, 2)) * img.pixel_solid_angle
        flux_direct = model.flux_density_grid(t_obs, np.array([nu_obs])).total[0, :]
        ratio = flux_from_image / flux_direct
        assert img.image.shape == (1, 64, 64)
        assert np.all(np.isfinite(img.image))
        return {"shape": list(img.image.shape), "ratio": float(ratio[0])}

    cases = [
        ("quickstart", case_quickstart),
        ("lightcurve_grid", case_lightcurve_grid),
        ("spectrum_grid", case_spectrum_grid),
        ("pair_points", case_pairs),
        ("exposures", case_exposures),
        ("details", case_details),
        ("sky_image", case_sky_image),
    ]
    return [_run_case(name, fn) for name, fn in cases]


def _plot_benchmark(results: dict[str, dict[str, list[dict]]]) -> None:
    names = [item["name"] for item in results["ism"]["fullhide"]]
    x = np.arange(len(names), dtype=float)
    width = 0.35
    colors = {"fullhide": "#1f77b4", "slc1": "#2ca02c", "slc1_mmg2": "#d62728"}

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True, sharey=True)
    for ax, medium_name in zip(axes, MEDIA):
        offsets = np.linspace(-width, width, len(SOLVERS))
        for offset, solver in zip(offsets, SOLVERS):
            seconds = np.array([item["seconds"] for item in results[medium_name][solver]], dtype=float)
            ax.bar(x + offset, seconds, width=0.28, label=solver, color=colors[solver])
        ax.set_xticks(x, names, rotation=20, ha="right")
        ax.set_title(medium_name.upper())
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("seconds")
    axes[0].legend()
    fig.suptitle("Benchmark with exponential-cutoff injection tail")
    _save_figure(fig, BENCHMARK_PNG, BENCHMARK_PDF, dpi=200)
    plt.close(fig)


def _multiband_curves() -> dict[str, dict[str, dict[str, np.ndarray]]]:
    times = np.logspace(2.0, 8.0, 96)
    results: dict[str, dict[str, dict[str, np.ndarray]]] = {}
    for medium_name in MEDIA:
        results[medium_name] = {}
        for solver in SOLVERS:
            print(f"multiband {medium_name} {solver}", flush=True)
            model = _build_model(solver, medium_name)
            freqs = np.array(list(BANDS.values()), dtype=float)
            flux = model.flux_density_grid(times, freqs).total
            results[medium_name][solver] = {"times": times, "flux": flux}
    return results


def _plot_multiband(curves: dict[str, dict[str, dict[str, np.ndarray]]]) -> None:
    colors = {"fullhide": "#1f77b4", "slc1": "#2ca02c", "slc1_mmg2": "#d62728"}
    fig, axes = plt.subplots(len(BANDS), len(MEDIA), figsize=(11, 13), constrained_layout=True, sharex=True)
    for col, medium_name in enumerate(MEDIA):
        for row, (band_name, freq) in enumerate(BANDS.items()):
            ax = axes[row, col]
            series = []
            for solver in SOLVERS:
                times = curves[medium_name][solver]["times"]
                flux = curves[medium_name][solver]["flux"][row, :]
                series.append(flux)
                ax.loglog(times, flux, label=solver, color=colors[solver])
            if row == 0:
                ax.set_title(f"{medium_name.upper()} forward-shock + SSC")
            if col == 0:
                ax.set_ylabel(f"{band_name}\nFlux")
            if row == len(BANDS) - 1:
                ax.set_xlabel("Time (s)")
            ax.text(0.03, 0.90, f"{band_name}: {freq:.2e} Hz", transform=ax.transAxes, fontsize=8)
            ax.grid(True, which="both", alpha=0.2)
            _tight_log_ylim(ax, series)
    axes[0, 0].legend(loc="lower left", fontsize=8)
    fig.suptitle("Radio-to-TeV multi-band afterglow with SSC")
    _save_figure(fig, MULTIBAND_PNG, MULTIBAND_PDF, dpi=220)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    benchmark_results: dict[str, dict[str, list[dict]]] = {}
    for medium_name in MEDIA:
        benchmark_results[medium_name] = {}
        for solver in SOLVERS:
            benchmark_results[medium_name][solver] = _benchmark_cases(solver, medium_name)

    with BENCHMARK_JSON.open("w", encoding="utf-8") as fh:
        json.dump(benchmark_results, fh, indent=2, ensure_ascii=False)
    _plot_benchmark(benchmark_results)

    curves = _multiband_curves()
    np.savez(
        MULTIBAND_NPZ,
        bands=np.array(list(BANDS.values()), dtype=float),
        band_names=np.array(list(BANDS.keys())),
        media=np.array(MEDIA),
        solvers=np.array(SOLVERS),
        ism_times=curves["ism"]["fullhide"]["times"],
        wind_times=curves["wind"]["fullhide"]["times"],
        ism_fullhide=curves["ism"]["fullhide"]["flux"],
        ism_slc1=curves["ism"]["slc1"]["flux"],
        ism_slc1_mmg2=curves["ism"]["slc1_mmg2"]["flux"],
        wind_fullhide=curves["wind"]["fullhide"]["flux"],
        wind_slc1=curves["wind"]["slc1"]["flux"],
        wind_slc1_mmg2=curves["wind"]["slc1_mmg2"]["flux"],
    )
    _plot_multiband(curves)

    failed = [
        (medium_name, solver, item["name"], item["payload"])
        for medium_name, medium_results in benchmark_results.items()
        for solver, solver_results in medium_results.items()
        for item in solver_results
        if item["status"] == "FAIL"
    ]
    print(BENCHMARK_JSON)
    print(BENCHMARK_PNG)
    print(MULTIBAND_PNG)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
