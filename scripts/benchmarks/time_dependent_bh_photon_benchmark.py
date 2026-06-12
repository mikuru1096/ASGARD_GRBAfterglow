from __future__ import annotations

import argparse
import csv
import json
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
import sys

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet, Wind
from ASGARD.api_observe import _build_fit_config_for_patch, _solve_patch_state
from asgard_core.asgard_paths import ASGARD_DOC_DIR


OUTPUT_DIR = ASGARD_DOC_DIR / "time_dependent_bh_photon_benchmark"
ASGARD_COLOR = "#1f77b4"
JOINT_COLOR = "#d62728"
BAND_COLORS = ("#1f77b4", "#2ca02c", "#d62728")
GRID_STYLE = {"alpha": 0.25, "linestyle": ":", "linewidth": 0.5}
BANDS_HZ = np.array([1.0e9, 1.0e14, 1.0e18], dtype=float)
BAND_LABELS = ["1 GHz", "1e14 Hz", "1e18 Hz"]
EPOCHS_S = np.array([1.0e3, 1.0e5, 1.0e7], dtype=float)

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.dpi": 180,
        "savefig.dpi": 220,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linestyle": ":",
        "grid.linewidth": 0.5,
    }
)


@dataclass(frozen=True)
class BenchmarkScenario:
    name: str
    medium: str
    n_ism: float
    e_iso: float
    gamma0: float
    eps_e: float
    eps_b: float
    epsilon_p: float
    a_star: float = 0.0
    wind_n_ism: float = 1.0e-15
    wind_n0: float | None = None
    theta_j: float = 0.1
    theta_obs: float = 0.0
    luminosity_distance_cm: float = 1.0e26
    z: float = 0.1
    p: float = 2.3


SCENARIOS = (
    BenchmarkScenario(
        name="weak_feedback",
        medium="ism",
        n_ism=1.0,
        e_iso=1.0e52,
        gamma0=120.0,
        eps_e=0.1,
        eps_b=1.0e-3,
        epsilon_p=0.03,
    ),
    BenchmarkScenario(
        name="bh_active",
        medium="ism",
        n_ism=10.0,
        e_iso=3.0e53,
        gamma0=180.0,
        eps_e=0.03,
        eps_b=1.0e-2,
        epsilon_p=0.5,
    ),
    BenchmarkScenario(
        name="strong_wind_bh",
        medium="wind",
        n_ism=1.0e-15,
        e_iso=1.0e54,
        gamma0=300.0,
        eps_e=0.1,
        eps_b=3.0e-2,
        epsilon_p=1.0,
        a_star=0.03,
        wind_n_ism=1.0e-15,
        luminosity_distance_cm=1.0e26,
    ),
)


@dataclass(frozen=True)
class GridMode:
    num_gam_e: int
    num_gam_p: int
    num_nu: int
    num_nu_nu: int
    num_r: int
    num_theta: int
    num_tobs: int
    lc_points: int
    sed_points: int


GRID_MODES = {
    "quick": GridMode(20, 24, 24, 16, 56, 10, 18, 36, 56),
    "formal": GridMode(48, 64, 64, 32, 80, 40, 80, 120, 180),
}


@dataclass
class BenchmarkRun:
    scenario: BenchmarkScenario
    coupling: str
    model: Model
    lightcurve_times_s: np.ndarray
    lightcurve_flux: np.ndarray
    sed_frequency_hz: np.ndarray
    sed_flux: np.ndarray
    state: object
    elapsed_s: float


def _save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=220, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {path}")
    print(f"[save] {path.with_suffix('.pdf')}")


def _git_text(args: list[str]) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True, encoding="utf-8").strip()


def _model(scenario: BenchmarkScenario, coupling: str, grid: GridMode) -> Model:
    medium = (
        ISM(scenario.n_ism)
        if scenario.medium == "ism"
        else Wind(A_star=scenario.a_star, n_ism=scenario.wind_n_ism, n0=scenario.wind_n0)
    )
    return Model(
        TophatJet(scenario.theta_j, scenario.e_iso, scenario.gamma0),
        medium,
        Observer(scenario.luminosity_distance_cm, scenario.z, scenario.theta_obs),
        Radiation(
            scenario.eps_e,
            scenario.eps_b,
            scenario.p,
            epsilon_p=scenario.epsilon_p,
            ssc=True,
            kn=True,
            proton_synch=True,
            bethe_heitler=True,
            pgamma_scheme="hummer_2010_response",
        ),
        setups=Setups(
            electron_solver="fullhide_1d",
            electron_photon_coupling=coupling,
            num_gam_e=grid.num_gam_e,
            num_gam_p=grid.num_gam_p,
            num_nu=grid.num_nu,
            num_nu_nu=grid.num_nu_nu,
            num_r=grid.num_r,
            num_theta=grid.num_theta,
            num_tobs=grid.num_tobs,
            hadronic_enabled=True,
            hadronic_solver="am3_1d",
            pgamma_scheme="hummer_2010_response",
        ),
    )


def _state_for_model(model: Model, times_s: np.ndarray, frequencies_hz: np.ndarray):
    config = _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )
    return _solve_patch_state(model, config, times_s, frequencies_hz)


def _run_case(scenario: BenchmarkScenario, coupling: str, grid: GridMode) -> BenchmarkRun:
    model = _model(scenario, coupling, grid)
    times = np.logspace(2.0, 8.0, grid.lc_points)
    sed_frequency = np.logspace(8.0, 24.0, grid.sed_points)
    t0 = time.perf_counter()
    lightcurve = model.component_fluxes(times, BANDS_HZ).total
    sed = model.flux_density_grid(EPOCHS_S, sed_frequency, projection_kind="sed").total
    state = _state_for_model(model, times, BANDS_HZ)
    elapsed = time.perf_counter() - t0
    _assert_case_physics(scenario, coupling, lightcurve, sed, state)
    return BenchmarkRun(scenario, coupling, model, times, lightcurve, sed_frequency, sed, state, elapsed)


def _assert_case_physics(
    scenario: BenchmarkScenario,
    coupling: str,
    lightcurve: np.ndarray,
    sed: np.ndarray,
    state,
) -> None:
    if not np.all(np.isfinite(lightcurve)) or not np.all(np.isfinite(sed)):
        raise RuntimeError(f"{coupling} benchmark produced non-finite flux.")
    if float(np.max(lightcurve)) <= 0.0 or float(np.max(sed)) <= 0.0:
        raise RuntimeError(f"{coupling} benchmark produced no positive flux.")
    if state.hadronic is None or state.hadronic.tau_bh is None or state.hadronic.d_n_gam_e_bh is None:
        raise RuntimeError(f"{coupling} benchmark requires BH diagnostics.")
    if not np.all(np.isfinite(state.hadronic.tau_bh)):
        raise RuntimeError(f"{coupling} benchmark produced non-finite tau_bh.")
    if np.any(np.asarray(state.hadronic.tau_bh, dtype=float) < 0.0):
        raise RuntimeError(f"{coupling} benchmark produced negative tau_bh.")
    if coupling == "joint" and state.adapter_reports["electron"].grid_semantics != "log-gamma-1d-joint-cooling":
        raise RuntimeError("joint benchmark did not use the joint cooling electron adapter.")
    _assert_smooth_shell_profile(
        f"{coupling} electron number",
        np.trapezoid(state.electron.d_n_gam_e, state.electron.gam_e, axis=0),
        limit=5.0e3,
        startup_skip=3 if scenario.medium == "wind" else 0,
    )


def _assert_smooth_shell_profile(name: str, profile: np.ndarray, *, limit: float, startup_skip: int = 0) -> None:
    positive = np.asarray(profile, dtype=float)
    positive = positive[np.isfinite(positive) & (positive > 0.0)]
    positive = positive[startup_skip:]
    if positive.size < 4:
        raise RuntimeError(f"{name} has too few positive shell samples.")
    adjacent = np.maximum(positive[1:] / positive[:-1], positive[:-1] / positive[1:])
    if float(np.max(adjacent)) > limit:
        raise RuntimeError(f"{name} adjacent shell ratio is too large: {float(np.max(adjacent)):.3g}")


def _run_benchmark(grid: GridMode) -> dict[tuple[str, str], BenchmarkRun]:
    runs: dict[tuple[str, str], BenchmarkRun] = {}
    for scenario in SCENARIOS:
        for coupling in ("separated", "joint"):
            print(f"[run] {scenario.name} {coupling}")
            runs[(scenario.name, coupling)] = _run_case(scenario, coupling, grid)
    return runs


def _ratio(joint: np.ndarray, separated: np.ndarray) -> np.ndarray:
    return np.asarray(joint, dtype=float) / np.maximum(np.asarray(separated, dtype=float), np.finfo(float).tiny)


def _plot_lightcurves(runs: dict[tuple[str, str], BenchmarkRun]) -> None:
    fig, axes = plt.subplots(len(SCENARIOS), 2, figsize=(10.8, 3.2 * len(SCENARIOS)), sharex=True)
    for i, scenario in enumerate(SCENARIOS):
        sep = runs[(scenario.name, "separated")]
        joint = runs[(scenario.name, "joint")]
        ax_lc, ax_ratio = axes[i]
        for j, (label, color) in enumerate(zip(BAND_LABELS, BAND_COLORS)):
            ax_lc.loglog(sep.lightcurve_times_s, sep.lightcurve_flux[j], color=color, alpha=0.55, linestyle="--")
            ax_lc.loglog(joint.lightcurve_times_s, joint.lightcurve_flux[j], color=color, alpha=0.95, label=label)
            ax_ratio.semilogx(
                joint.lightcurve_times_s,
                _ratio(joint.lightcurve_flux[j], sep.lightcurve_flux[j]),
                color=color,
                label=label,
            )
        ax_lc.set_title(scenario.name)
        ax_lc.set_ylabel(r"$F_\nu$ [mJy]")
        ax_lc.grid(**GRID_STYLE)
        ax_ratio.axhline(1.0, color="0.35", linewidth=0.9)
        ax_ratio.set_ylabel("joint / separated")
        ax_ratio.grid(**GRID_STYLE)
    axes[-1, 0].set_xlabel(r"$t_{\rm obs}$ [s]")
    axes[-1, 1].set_xlabel(r"$t_{\rm obs}$ [s]")
    band_handles = [Line2D([0], [0], color=color, label=label) for label, color in zip(BAND_LABELS, BAND_COLORS)]
    coupling_handles = [
        Line2D([0], [0], color="0.2", linestyle="--", label="separated"),
        Line2D([0], [0], color="0.2", linestyle="-", label="joint"),
    ]
    axes[0, 0].legend(handles=band_handles + coupling_handles, ncol=2)
    axes[0, 1].legend()
    _save(fig, OUTPUT_DIR / "time_dependent_bh_photon_lightcurves.png")


def _plot_seds(runs: dict[tuple[str, str], BenchmarkRun]) -> None:
    fig, axes = plt.subplots(1, len(SCENARIOS), figsize=(5.4 * len(SCENARIOS), 4.6), sharey=True)
    for ax, scenario in zip(axes, SCENARIOS):
        sep = runs[(scenario.name, "separated")]
        joint = runs[(scenario.name, "joint")]
        for i, (epoch, color) in enumerate(zip(EPOCHS_S, BAND_COLORS)):
            ax.loglog(sep.sed_frequency_hz, sep.sed_flux[:, i], color=color, alpha=0.55, linestyle="--")
            ax.loglog(joint.sed_frequency_hz, joint.sed_flux[:, i], color=color, alpha=0.95, label=f"{epoch:.0e} s")
        ax.set_title(scenario.name)
        ax.set_xlabel(r"$\nu$ [Hz]")
        ax.grid(**GRID_STYLE)
    axes[0].set_ylabel(r"$F_\nu$ [mJy]")
    axes[0].legend(title="epochs")
    _save(fig, OUTPUT_DIR / "time_dependent_bh_photon_seds.png")


def _shell_integrals(run: BenchmarkRun) -> dict[str, np.ndarray]:
    state = run.state
    tau_bh = np.asarray(state.hadronic.tau_bh, dtype=float)
    return {
        "radius_cm": np.asarray(state.dynamics.radius, dtype=float),
        "electron_number": np.trapezoid(state.electron.d_n_gam_e, state.electron.gam_e, axis=0),
        "bh_pair_number": np.trapezoid(state.hadronic.d_n_gam_e_bh, state.electron.gam_e, axis=0),
        "target_seed_integral": np.trapezoid(
            state.photon_field.hadronic_target_seed,
            state.photon_field.seed_frequency_hz,
            axis=0,
        ),
        "tau_bh_max": np.max(tau_bh, axis=0),
        "tau_bh_logmean": np.trapezoid(tau_bh, np.log(state.photon_field.seed_frequency_hz), axis=0)
        / np.log(state.photon_field.seed_frequency_hz[-1] / state.photon_field.seed_frequency_hz[0]),
    }


def _plot_shell_diagnostics(runs: dict[tuple[str, str], BenchmarkRun]) -> None:
    fig, axes = plt.subplots(len(SCENARIOS), 4, figsize=(14.0, 3.0 * len(SCENARIOS)), sharex="row")
    fields = [
        ("electron_number", r"$N_e$"),
        ("bh_pair_number", r"$N_{e,\rm BH}$"),
        ("target_seed_integral", r"$\int n_\gamma d\nu$"),
        ("tau_bh_max", r"$\max_\nu \tau_{\rm BH}$"),
    ]
    for i, scenario in enumerate(SCENARIOS):
        for coupling, color, linestyle in (("separated", ASGARD_COLOR, "--"), ("joint", JOINT_COLOR, "-")):
            values = _shell_integrals(runs[(scenario.name, coupling)])
            plot_slice = slice(8, None) if scenario.medium == "wind" else slice(None)
            radius = values["radius_cm"][plot_slice]
            for ax, (key, label) in zip(axes[i], fields):
                ax.loglog(radius, values[key][plot_slice], color=color, linestyle=linestyle, label=coupling)
                ax.set_ylabel(label)
                ax.grid(**GRID_STYLE)
        axes[i, 0].set_title(scenario.name)
    for ax in axes[-1]:
        ax.set_xlabel(r"$R$ [cm]")
    axes[0, -1].legend()
    _save(fig, OUTPUT_DIR / "time_dependent_bh_photon_shell_diagnostics.png")


def _summary_rows(runs: dict[tuple[str, str], BenchmarkRun]) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for scenario in SCENARIOS:
        sep = runs[(scenario.name, "separated")]
        joint = runs[(scenario.name, "joint")]
        ratio = _ratio(joint.lightcurve_flux, sep.lightcurve_flux)
        rows.append(
            {
                "scenario": scenario.name,
                "elapsed_separated_s": sep.elapsed_s,
                "elapsed_joint_s": joint.elapsed_s,
                "max_tau_bh_joint": float(np.max(joint.state.hadronic.tau_bh)),
                "median_lightcurve_ratio_joint_over_separated": float(np.median(ratio)),
                "max_lightcurve_ratio_joint_over_separated": float(np.max(ratio)),
            }
        )
    return rows


def _plot_summary(runs: dict[tuple[str, str], BenchmarkRun]) -> None:
    rows = _summary_rows(runs)
    scenarios = [str(row["scenario"]) for row in rows]
    scenario_labels = [scenario.replace("_", "\n") for scenario in scenarios]
    x = np.arange(len(rows), dtype=float)
    width = 0.34
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8))
    axes[0].bar(x - 0.5 * width, [float(row["elapsed_separated_s"]) for row in rows], width, color=ASGARD_COLOR, label="separated")
    axes[0].bar(x + 0.5 * width, [float(row["elapsed_joint_s"]) for row in rows], width, color=JOINT_COLOR, label="joint")
    axes[0].set_ylabel("elapsed [s]")
    axes[0].set_xticks(x, scenario_labels)
    axes[0].legend()
    axes[1].semilogy(x, [float(row["max_tau_bh_joint"]) for row in rows], "o-", color=JOINT_COLOR)
    axes[1].set_ylabel(r"$\max \tau_{\rm BH}$")
    axes[1].set_xticks(x, scenario_labels)
    axes[2].plot(x, [float(row["median_lightcurve_ratio_joint_over_separated"]) for row in rows], "o-", label="median", color="#9467bd")
    axes[2].plot(x, [float(row["max_lightcurve_ratio_joint_over_separated"]) for row in rows], "s--", label="max", color="#8c564b")
    axes[2].axhline(1.0, color="0.35", linewidth=0.9)
    axes[2].set_ylabel("LC ratio joint / separated")
    axes[2].set_xticks(x, scenario_labels)
    axes[2].legend()
    for ax in axes:
        ax.grid(**GRID_STYLE)
    _save(fig, OUTPUT_DIR / "time_dependent_bh_photon_coupling_summary.png")


def _write_lightcurve_csv(runs: dict[tuple[str, str], BenchmarkRun]) -> None:
    path = OUTPUT_DIR / "time_dependent_bh_photon_lightcurves.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["scenario", "coupling", "frequency_hz", "time_s", "flux_mjy"])
        for run in runs.values():
            for i_freq, frequency in enumerate(BANDS_HZ):
                for time_s, flux in zip(run.lightcurve_times_s, run.lightcurve_flux[i_freq]):
                    writer.writerow([run.scenario.name, run.coupling, frequency, time_s, flux])
    print(f"[save] {path}")


def _write_sed_csv(runs: dict[tuple[str, str], BenchmarkRun]) -> None:
    path = OUTPUT_DIR / "time_dependent_bh_photon_seds.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["scenario", "coupling", "epoch_s", "frequency_hz", "flux_mjy"])
        for run in runs.values():
            for i_epoch, epoch in enumerate(EPOCHS_S):
                for frequency, flux in zip(run.sed_frequency_hz, run.sed_flux[:, i_epoch]):
                    writer.writerow([run.scenario.name, run.coupling, epoch, frequency, flux])
    print(f"[save] {path}")


def _write_shell_csv(runs: dict[tuple[str, str], BenchmarkRun]) -> None:
    path = OUTPUT_DIR / "time_dependent_bh_photon_shell_diagnostics.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "scenario",
                "coupling",
                "radius_cm",
                "electron_number",
                "bh_pair_number",
                "target_seed_integral",
                "tau_bh_max",
                "tau_bh_logmean",
            ]
        )
        for run in runs.values():
            values = _shell_integrals(run)
            for i in range(values["radius_cm"].size):
                writer.writerow([run.scenario.name, run.coupling, *[values[key][i] for key in values]])
    print(f"[save] {path}")


def _write_summary(runs: dict[tuple[str, str], BenchmarkRun], mode: str, command: str) -> None:
    path = OUTPUT_DIR / "time_dependent_bh_photon_summary.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "scenario",
                "elapsed_separated_s",
                "elapsed_joint_s",
                "max_tau_bh_joint",
                "median_lightcurve_ratio_joint_over_separated",
                "max_lightcurve_ratio_joint_over_separated",
            ]
        )
        for row in _summary_rows(runs):
            writer.writerow([row[key] for key in row])
    metadata = {
        "mode": mode,
        "command": command,
        "git_head": _git_text(["rev-parse", "--short", "HEAD"]),
        "git_status_short_branch": _git_text(["status", "--short", "--branch"]),
        "git_diff_stat": _git_text(["diff", "--stat"]),
        "scenarios": [scenario.__dict__ for scenario in SCENARIOS],
        "shell_plot_startup_skip_by_medium": {"ism": 0, "wind": 8},
        "bands_hz": BANDS_HZ.tolist(),
        "epochs_s": EPOCHS_S.tolist(),
        "outputs": sorted(path.name for path in OUTPUT_DIR.glob("time_dependent_bh_photon_*")),
    }
    meta_path = OUTPUT_DIR / "time_dependent_bh_photon_metadata.json"
    meta_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"[save] {path}")
    print(f"[save] {meta_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark separated vs joint time-dependent BH photon coupling.")
    parser.add_argument("--mode", choices=sorted(GRID_MODES), default="quick")
    args = parser.parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    grid = GRID_MODES[args.mode]
    command = " ".join(sys.argv)
    runs = _run_benchmark(grid)
    _plot_lightcurves(runs)
    _plot_seds(runs)
    _plot_shell_diagnostics(runs)
    _plot_summary(runs)
    _write_lightcurve_csv(runs)
    _write_sed_csv(runs)
    _write_shell_csv(runs)
    _write_summary(runs, args.mode, command)
    print(f"time_dependent_bh_photon_benchmark: ok ({OUTPUT_DIR})")


if __name__ == "__main__":
    main()
