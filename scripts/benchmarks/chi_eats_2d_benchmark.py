from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core import Model, Observer, UniformMedium, top_hat_jet
from asgard_core.api_model import _solve_patch_state
from asgard_core.asgard_paths import DOC_ROOT
from asgard_core.asgard_setup import build_setup
from asgard_core.asgard_state import query_cfg, project_flux, solve_setup
from scripts.benchmarks.benchmark_common import DATA_ROOT, FIGURE_ROOT, environment, joint_significant, plot_style, save_figure, write_json
from tests.public_api_builders import hadronic, numerics, observer_grid, radiation, reverse_shock, solver_options


OUTPUT_DIR = DOC_ROOT / "chi_eats_2d_benchmark"
LIGHTCURVE_TIMES = np.logspace(2.0, 9.0, 128)
SED_TIMES = np.array([1.0e3, 1.0e5, 1.0e7], dtype=float)
LIGHTCURVE_BANDS = np.array([1.0e9, 1.0e14, 1.0e18], dtype=float)
SED_FREQUENCIES = np.logspace(6.0, 28.0, 240)
CHI_GRID_SCAN_TIMES = np.logspace(2.0, 9.0, 64)
THETA_J = 0.1
CHI_GRID_SCAN_THETA_RATIO = 0.5
THETA_RATIO_CASES = np.array([0.0, 0.5, 1.0, 1.5, 3.0, 5.0], dtype=float)
E_ISO = 1.0e53
GAMMA0 = 300.0
EPSILON_E = 0.1
EPSILON_B = 1.0e-3
EPSILON_B_FLOOR = EPSILON_B
MAGNETIC_DECAY_ALPHA_T = 0.0
MAGNETIC_DECAY_T0_S = 1.0
ELECTRON_P = 2.3
LUMINOSITY_DISTANCE_CM = 1.0e26
REDSHIFT = 0.1
ISM_N = 1.0
WIND_ASTAR = 0.1
WIND_N_ISM = 1.0e-15

QUICK_GRID = dict(num_gam_e=31, num_chi=24, num_nu=41, num_r=300, num_theta=300, num_phi=50, num_tobs=128)
FORMAL_GRID = dict(num_gam_e=61, num_chi=24, num_nu=81, num_r=300, num_theta=300, num_phi=50, num_tobs=128)
THETA_SCAN_QUICK_GRID = QUICK_GRID
THETA_SCAN_FORMAL_GRID = FORMAL_GRID
CHI_GRID_SCAN_VALUES = np.array([32, 64, 128, 256, 512], dtype=int)
CHI_GRID_SCAN_GRID = dict(num_gam_e=16, num_chi=32, num_nu=21, num_r=150, num_theta=150, num_phi=25, num_tobs=64)

PLOT_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "pdf.fonttype": 42,
    "font.size": 7.5,
    "axes.labelsize": 8.0,
    "axes.titlesize": 8.6,
    "legend.fontsize": 6.7,
    "legend.title_fontsize": 7.0,
    "xtick.labelsize": 7.0,
    "ytick.labelsize": 7.0,
    "axes.linewidth": 0.75,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
    "figure.dpi": 140,
}
GEOMETRY_STYLES = {
    "sed_legacy": ("#111111", "-"),
    "chi_eats_2d": ("#0072B2", "--"),
}
BAND_COLORS = {1.0e9: "#0072B2", 1.0e14: "#D55E00", 1.0e18: "#009E73"}
TIME_COLORS = {1.0e3: "#0072B2", 1.0e5: "#D55E00", 1.0e7: "#009E73"}
SED_MODEL_STYLES = {
    "1d": ("#111111", "-"),
    "2d_legacy": ("#D55E00", ":"),
    "2d_chi": ("#0072B2", "--"),
}
METHOD_STYLES = {
    "legacy": ("sed legacy", "#111111", "-"),
    "chi": ("chi EATS", "#0072B2", "--"),
    "one_d": ("1D", "#111111", "-"),
    "two_d_chi": ("2D chi EATS", "#0072B2", "--"),
}


@dataclass(frozen=True)
class CaseSpec:
    solver: str
    geometry_kernel: str

    @property
    def name(self) -> str:
        return f"{self.solver}_{self.geometry_kernel}"


@dataclass(frozen=True)
class MediumSpec:
    name: str
    label: str
    n_ism: float
    a_star: float


ISM_SPEC = MediumSpec(name="ism", label=fr"ISM, $n={ISM_N:g}$ cm$^{{-3}}$", n_ism=ISM_N, a_star=-1.0)
WIND_SPEC = MediumSpec(
    name="wind",
    label=fr"Wind, $A_\ast={WIND_ASTAR:g}$, $n_\mathrm{{ISM}}={WIND_N_ISM:.0e}$ cm$^{{-3}}$",
    n_ism=WIND_N_ISM,
    a_star=WIND_ASTAR,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Chi-resolved 2D EATS benchmark plots with shell-level SSC included.")
    parser.add_argument("--mode", choices=("quick", "formal"), default="quick")
    parser.add_argument("--solver", choices=("fullhide_2d",), default="fullhide_2d")
    parser.add_argument("--medium", choices=("ism", "wind", "both"), default="ism")
    parser.add_argument("--only", choices=("all", "theta-scan", "chi-grid-scan", "convergence"), default="all")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path)
    return parser.parse_args()


def _git_head() -> str:
    result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True, encoding="utf-8")
    return result.stdout.strip()


def _git_status() -> str:
    result = subprocess.run(["git", "status", "--short", "--branch"], cwd=ROOT, check=True, capture_output=True, text=True, encoding="utf-8")
    return result.stdout.strip()


def _git_diff_stat() -> str:
    result = subprocess.run(["git", "diff", "--stat"], cwd=ROOT, check=True, capture_output=True, text=True, encoding="utf-8")
    return result.stdout.strip()


def _one_d_solver(two_d_solver: str) -> str:
    if two_d_solver == "fullhide_2d":
        return "fullhide_1d"
    raise ValueError(f"unsupported 2d solver for comparison: {two_d_solver}")


def _selected_media(choice: str) -> tuple[MediumSpec, ...]:
    if choice == "ism":
        return (ISM_SPEC,)
    if choice == "wind":
        return (WIND_SPEC,)
    return (ISM_SPEC, WIND_SPEC)


def _medium_model(spec: MediumSpec):
    if spec.name == "ism":
        return ISM(n_ism=spec.n_ism)
    return Wind(A_star=spec.a_star, n_ism=spec.n_ism)


def _projection_grid(grid: dict[str, int], theta_v: float) -> dict[str, int]:
    out = dict(grid)
    if theta_v == 0.0:
        out["num_phi"] = 1
    return out


def _build_model(
    case: CaseSpec,
    grid: dict[str, int],
    *,
    medium: MediumSpec = ISM_SPEC,
    theta_j: float = THETA_J,
    theta_v: float = 0.0,
) -> Model:
    grid = _projection_grid(grid, theta_v)
    return Model(
        jet=TophatJet(theta_c=theta_j, E_iso=E_ISO, Gamma0=GAMMA0),
        medium=_medium_model(medium),
        observer=Observer(lumi_dist=LUMINOSITY_DISTANCE_CM, z=REDSHIFT, theta_obs=theta_v),
        fwd_rad=Radiation(
            eps_e=EPSILON_E,
            eps_B=EPSILON_B,
            p=ELECTRON_P,
            epsilon_b_floor=EPSILON_B_FLOOR,
            magnetic_decay_alpha_t=MAGNETIC_DECAY_ALPHA_T,
            magnetic_decay_t0_s=MAGNETIC_DECAY_T0_S,
            ssc=True,
            kn=True,
        ),
        setups=Setups(
            electron_solver=case.solver,
            geometry_kernel=case.geometry_kernel,
            num_gam_e=grid["num_gam_e"],
            num_chi=grid["num_chi"],
            num_nu=grid["num_nu"],
            num_r=grid["num_r"],
            num_theta=grid["num_theta"],
            num_phi=grid.get("num_phi", 1),
            num_tobs=grid["num_tobs"],
            electron_adaptive_substeps=False,
        ),
    )


def _build_fit_config(
    case: CaseSpec,
    grid: dict[str, int],
    *,
    medium: MediumSpec = ISM_SPEC,
    theta_j: float = THETA_J,
    theta_v: float = 0.0,
) -> FitConfig:
    grid = _projection_grid(grid, theta_v)
    return FitConfig(
        electron_solver=case.solver,
        geometry_kernel=case.geometry_kernel,
        include_forward_ssc=True,
        num_gam_e=grid["num_gam_e"],
        num_chi=grid["num_chi"],
        num_nu=grid["num_nu"],
        num_r=grid["num_r"],
        num_theta=grid["num_theta"],
        num_phi=grid.get("num_phi", 1),
        num_tobs=grid["num_tobs"],
        t_obs_min_log10=float(np.log10(LIGHTCURVE_TIMES[0])),
        t_obs_max_log10=float(np.log10(LIGHTCURVE_TIMES[-1])),
        eta_0=GAMMA0,
        e_iso=E_ISO,
        d_ne=medium.n_ism,
        a_star=medium.a_star,
        z=REDSHIFT,
        opening_angle_jet=theta_j,
        theta_v=theta_v,
        epsilon_e=EPSILON_E,
        epsilon_b=EPSILON_B,
        epsilon_b_floor=EPSILON_B_FLOOR,
        magnetic_decay_alpha_t=MAGNETIC_DECAY_ALPHA_T,
        magnetic_decay_t0_s=MAGNETIC_DECAY_T0_S,
        p=ELECTRON_P,
        electron_adaptive_substeps=False,
        luminosity_distance_cm_override=LUMINOSITY_DISTANCE_CM,
    )


def _collect_case(case: CaseSpec, grid: dict[str, int], medium: MediumSpec) -> dict[str, Any]:
    print(f"  collect {medium.name} {case.name}", flush=True)
    model = _build_model(case, grid, medium=medium)
    lightcurve_fnu = np.asarray(
        model.flux_density_grid(LIGHTCURVE_TIMES, LIGHTCURVE_BANDS, projection_kind="lightcurve").total,
        dtype=float,
    )
    sed_fnu = np.asarray(model.flux_density_grid(SED_TIMES, SED_FREQUENCIES, projection_kind="sed").total, dtype=float)
    detail_config = query_cfg(_build_fit_config(case, grid, medium=medium), LIGHTCURVE_TIMES)
    detail_setup = build_setup(detail_config, SED_FREQUENCIES, observer_time_s=LIGHTCURVE_TIMES)
    state = solve_setup(detail_config, detail_setup, requested_frequencies_hz=SED_FREQUENCIES)
    out: dict[str, Any] = {
        "lightcurve_nufnu": LIGHTCURVE_BANDS[:, None] * lightcurve_fnu,
        "sed_nufnu": SED_FREQUENCIES[:, None] * sed_fnu,
        "detail_time": np.asarray(state.components.fwd.characteristic_time_s, dtype=float),
        "nu_m": np.asarray(state.electron.nu_m, dtype=float),
        "nu_c": np.asarray(state.electron.nu_c, dtype=float),
        "nu_a": np.asarray(state.electron.nu_a, dtype=float),
    }
    if case.geometry_kernel == "chi_eats_2d":
        out["chi_radius_cm"] = np.asarray(state.electron.chi_radius_cm, dtype=float)
        out["chi_gamma_bulk"] = np.asarray(state.electron.chi_gamma_bulk, dtype=float)
        out["chi_dvolume_weight"] = np.asarray(state.electron.chi_dvolume_weight, dtype=float)
        out["l_syn_spec_chi"] = np.asarray(state.electron.l_syn_spec_chi, dtype=float)
        out["tau_syn_chi"] = np.asarray(state.electron.tau_syn_chi, dtype=float)
        out["seed_frequency_hz"] = np.asarray(state.setup.seed_frequency_hz, dtype=float)
    return out


def _collect_sed_case(case: CaseSpec, grid: dict[str, int], medium: MediumSpec) -> dict[str, np.ndarray]:
    print(f"  collect SED {medium.name} {case.name}", flush=True)
    model = _build_model(case, grid, medium=medium)
    sed_fnu = np.asarray(model.flux_density_grid(SED_TIMES, SED_FREQUENCIES, projection_kind="sed").total, dtype=float)
    return {"sed_nufnu": SED_FREQUENCIES[:, None] * sed_fnu}


def _collect_lightcurve_case(
    case: CaseSpec,
    grid: dict[str, int],
    *,
    medium: MediumSpec,
    theta_j: float,
    theta_v: float,
) -> dict[str, np.ndarray]:
    print(f"  collect LC {medium.name} {case.name} theta_v/theta_j={theta_v / theta_j:.2f}", flush=True)
    model = _build_model(case, grid, medium=medium, theta_j=theta_j, theta_v=theta_v)
    lightcurve_fnu = np.asarray(
        model.flux_density_grid(LIGHTCURVE_TIMES, LIGHTCURVE_BANDS, projection_kind="lightcurve").total,
        dtype=float,
    )
    return {"lightcurve_nufnu": LIGHTCURVE_BANDS[:, None] * lightcurve_fnu}


def _collect_theta_scan_base_state(
    case: CaseSpec,
    grid: dict[str, int],
    *,
    medium: MediumSpec,
    theta_j: float,
) -> dict[str, Any]:
    print(f"  solve theta-scan base {medium.name} {case.name}", flush=True)
    model = _build_model(case, grid, medium=medium, theta_j=theta_j, theta_v=0.0)
    base_config = _build_fit_config(case, grid, medium=medium, theta_j=theta_j, theta_v=0.0)
    state = _solve_patch_state(model, base_config, LIGHTCURVE_TIMES, LIGHTCURVE_BANDS)
    return {"state": state, "case": case, "grid": dict(grid)}


def _collect_lightcurve_from_state(
    base: dict[str, Any],
    *,
    medium: MediumSpec,
    theta_j: float,
    theta_v: float,
) -> dict[str, np.ndarray]:
    state = base["state"]
    theta_config = query_cfg(
        _build_fit_config(base["case"], base["grid"], medium=medium, theta_j=theta_j, theta_v=theta_v),
        LIGHTCURVE_TIMES,
    )
    theta_setup = build_setup(theta_config, LIGHTCURVE_BANDS, observer_time_s=LIGHTCURVE_TIMES)
    theta_state = replace(state, config=theta_config, setup=theta_setup)
    projected = project_flux(theta_state, LIGHTCURVE_TIMES, LIGHTCURVE_BANDS, mode="total_only")
    return {"lightcurve_nufnu": LIGHTCURVE_BANDS[:, None] * np.asarray(projected.components["total"], dtype=float)}


def _safe_ratio(numer: np.ndarray, denom: np.ndarray) -> np.ndarray:
    numer = np.asarray(numer, dtype=float)
    denom = np.asarray(denom, dtype=float)
    return np.divide(numer, denom, out=np.full_like(numer, np.nan), where=(numer > 0.0) & (denom > 0.0))


def _finite_positive(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr) & (arr > 0.0)]


def _positive_limits(arrays: list[np.ndarray], pad: float = 0.25) -> tuple[float, float] | None:
    values = [_finite_positive(arr) for arr in arrays]
    values = [arr for arr in values if arr.size]
    if not values:
        return None
    merged = np.concatenate(values)
    return 10.0 ** (np.log10(float(np.min(merged))) - pad), 10.0 ** (np.log10(float(np.max(merged))) + pad)


def _ratio_limits(arrays: list[np.ndarray], pad: float = 0.10) -> tuple[float, float] | None:
    values = [_finite_positive(arr) for arr in arrays]
    values = [arr for arr in values if arr.size]
    if not values:
        return None
    merged = np.concatenate(values)
    span = max(abs(float(np.log10(np.min(merged)))), abs(float(np.log10(np.max(merged))))) + pad
    return 10.0 ** (-span), 10.0**span


def _display_limits_quantile(arrays: list[np.ndarray], low: float = 0.03, high: float = 0.97, pad: float = 0.28) -> tuple[float, float] | None:
    values = [_finite_positive(arr) for arr in arrays]
    values = [arr for arr in values if arr.size]
    if not values:
        return None
    merged = np.concatenate(values)
    lo = float(np.quantile(merged, low))
    hi = float(np.quantile(merged, high))
    if lo <= 0.0 or hi <= lo:
        return _positive_limits(arrays, pad=pad)
    return 10.0 ** (np.log10(lo) - pad), 10.0 ** (np.log10(hi) + pad)


def _bounded_display_limits(
    arrays: list[np.ndarray],
    low: float = 0.05,
    high: float = 0.985,
    pad: float = 0.25,
    max_decades: float = 8.0,
) -> tuple[float, float] | None:
    limits = _display_limits_quantile(arrays, low=low, high=high, pad=pad)
    if limits is None:
        return None
    ymin, ymax = limits
    return max(ymin, ymax / 10.0**max_decades), ymax


def _ratio_display_limits(arrays: list[np.ndarray], pad: float = 0.15) -> tuple[float, float] | None:
    values = [_finite_positive(arr) for arr in arrays]
    values = [arr for arr in values if arr.size]
    if not values:
        return None
    merged = np.concatenate(values)
    lo = float(np.quantile(merged, 0.08))
    hi = float(np.quantile(merged, 0.92))
    if lo <= 0.0 or hi <= lo:
        return _ratio_limits(arrays, pad=pad)
    return 10.0 ** (np.log10(lo) - pad), 10.0 ** (np.log10(hi) + pad)


def _log10_ratio(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return np.log10(arr, out=np.full_like(arr, np.nan), where=(np.isfinite(arr) & (arr > 0.0)))


def _style_dex_ratio_axis(ax, xlabel: str, ylabel: str) -> None:
    ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.axhline(0.0, color="#222222", lw=0.7, ls=":")
    ax.set_ylim(-3.0, 3.0)
    ax.set_yticks(np.arange(-3.0, 3.1, 1.0))
    ax.tick_params(which="both", direction="out", length=2.8, width=0.65, top=False, right=False)
    ax.tick_params(which="minor", length=1.8, width=0.45)
    ax.grid(which="major", alpha=0.16, lw=0.45)
    ax.grid(which="minor", axis="x", alpha=0.05, lw=0.3)


def _setup_log_axes(ax, xlabel: str | None = None, ylabel: str | None = None) -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.tick_params(which="both", direction="out", length=2.8, width=0.65, top=False, right=False)
    ax.tick_params(which="minor", length=1.8, width=0.45)
    ax.grid(which="major", alpha=0.16, lw=0.45)
    ax.grid(which="minor", alpha=0.05, lw=0.3)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)


def _label_frequency(value: float) -> str:
    exp = int(np.floor(np.log10(value)))
    base = value / 10.0**exp
    return fr"${base:.0f}\times10^{{{exp}}}$ Hz"


def _save(fig, output_dir: Path, stem: str) -> list[Path]:
    paths = [output_dir / f"{stem}.png", output_dir / f"{stem}.pdf"]
    fig.savefig(paths[0], dpi=360, bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight")
    plt.close(fig)
    for path in paths:
        print(path)
    return paths


def _band_handles() -> list[Line2D]:
    return [
        Line2D([0], [0], color=BAND_COLORS[float(band)], lw=1.8, label=_label_frequency(float(band)))
        for band in LIGHTCURVE_BANDS
    ]


def _time_handles() -> list[Line2D]:
    return [Line2D([0], [0], color=TIME_COLORS[float(t_obs)], lw=1.8, label=fr"$t={t_obs:.0e}$ s") for t_obs in SED_TIMES]


def _method_handles(methods: tuple[str, ...]) -> list[Line2D]:
    handles = []
    for method in methods:
        label, color, linestyle = METHOD_STYLES[method]
        handles.append(Line2D([0], [0], color=color, ls=linestyle, lw=1.8, label=label))
    return handles


def _style_ratio_axis(ax, xlabel: str, ylabel: str) -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.axhline(1.0, color="#222222", lw=0.7, ls=":")
    ax.tick_params(which="both", direction="out", length=2.8, width=0.65, top=False, right=False)
    ax.tick_params(which="minor", length=1.8, width=0.45)
    ax.grid(which="major", alpha=0.16, lw=0.45)
    ax.grid(which="minor", alpha=0.05, lw=0.3)


def _plot_lightcurve_solver(solver: str, data: dict[str, dict[str, Any]], output_dir: Path, medium: MediumSpec) -> list[Path]:
    legacy = data[f"{solver}_sed_legacy"]
    chi = data[f"{solver}_chi_eats_2d"]
    fig, axes = plt.subplots(2, 1, figsize=(6.9, 5.2), constrained_layout=True, sharex=True, height_ratios=(1.45, 1.0))
    flux_values: list[np.ndarray] = []
    ratio_values: list[np.ndarray] = []
    for band in LIGHTCURVE_BANDS:
        idx = int(np.where(LIGHTCURVE_BANDS == band)[0][0])
        color = BAND_COLORS[float(band)]
        legacy_curve = legacy["lightcurve_nufnu"][idx]
        chi_curve = chi["lightcurve_nufnu"][idx]
        ratio = _safe_ratio(chi_curve, legacy_curve)
        flux_values.extend([legacy_curve, chi_curve])
        ratio_values.append(ratio)
        axes[0].loglog(LIGHTCURVE_TIMES, legacy_curve, color=color, ls="-", lw=1.45)
        axes[0].loglog(LIGHTCURVE_TIMES, chi_curve, color=color, ls="--", lw=1.45)
        mask = np.isfinite(ratio) & (ratio > 0.0)
        if np.any(mask):
            axes[1].semilogx(LIGHTCURVE_TIMES[mask], ratio[mask], color=color, lw=1.5)
    _setup_log_axes(axes[0], ylabel=r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
    _style_ratio_axis(axes[1], "Observer time [s]", "chi EATS / legacy")
    flux_lim = _positive_limits(flux_values)
    ratio_lim = _ratio_limits(ratio_values)
    if flux_lim:
        axes[0].set_ylim(*flux_lim)
    if ratio_lim:
        axes[1].set_ylim(*ratio_lim)
    leg1 = axes[0].legend(handles=_band_handles(), title="Band", loc="lower left", ncol=1)
    axes[0].add_artist(leg1)
    axes[0].legend(handles=_method_handles(("legacy", "chi")), title="Projection", loc="upper right")
    fig.suptitle(f"{solver}: total light-curve comparison with chi-resolved FS synchrotron, {medium.label}")
    return _save(fig, output_dir, f"{solver}_{medium.name}_lightcurve_chi_vs_legacy")


def _plot_sed_solver(solver: str, data: dict[str, dict[str, Any]], output_dir: Path, medium: MediumSpec) -> list[Path]:
    legacy = data[f"{solver}_sed_legacy"]
    chi = data[f"{solver}_chi_eats_2d"]
    fig, axes = plt.subplots(2, 1, figsize=(6.9, 5.2), constrained_layout=True, sharex=True, height_ratios=(1.45, 1.0))
    flux_values: list[np.ndarray] = []
    ratio_values: list[np.ndarray] = []
    for t_obs in SED_TIMES:
        idx = int(np.where(SED_TIMES == t_obs)[0][0])
        color = TIME_COLORS[float(t_obs)]
        legacy_sed = legacy["sed_nufnu"][:, idx]
        chi_sed = chi["sed_nufnu"][:, idx]
        ratio = _safe_ratio(chi_sed, legacy_sed)
        flux_values.extend([legacy_sed, chi_sed])
        ratio_values.append(ratio)
        axes[0].loglog(SED_FREQUENCIES, legacy_sed, color=color, ls="-", lw=1.45)
        axes[0].loglog(SED_FREQUENCIES, chi_sed, color=color, ls="--", lw=1.45)
        mask = np.isfinite(ratio) & (ratio > 0.0)
        if np.any(mask):
            axes[1].semilogx(SED_FREQUENCIES[mask], ratio[mask], color=color, lw=1.6, label=fr"$t={t_obs:.0e}$ s")
    _setup_log_axes(axes[0], ylabel=r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
    _style_ratio_axis(axes[1], "Observed frequency [Hz]", "chi EATS / legacy")
    flux_lim = _positive_limits(flux_values)
    ratio_lim = _ratio_limits(ratio_values)
    if flux_lim:
        axes[0].set_ylim(*flux_lim)
    if ratio_lim:
        axes[1].set_ylim(*ratio_lim)
    leg1 = axes[0].legend(handles=_time_handles(), title="Epoch", loc="lower left")
    axes[0].add_artist(leg1)
    axes[0].legend(handles=_method_handles(("legacy", "chi")), title="Projection", loc="upper right")
    fig.suptitle(f"{solver}: total SED comparison with chi-resolved FS synchrotron, {medium.label}")
    return _save(fig, output_dir, f"{solver}_{medium.name}_sed_chi_vs_legacy")


def _plot_1d_2d_sed_solver(
    solver: str,
    sed_data: dict[str, dict[str, np.ndarray]],
    output_dir: Path,
    medium: MediumSpec,
) -> list[Path]:
    one_d_solver = _one_d_solver(solver)
    one_d = sed_data[f"{one_d_solver}_sed_legacy"]["sed_nufnu"]
    two_d_legacy = sed_data[f"{solver}_sed_legacy"]["sed_nufnu"]
    two_d_chi = sed_data[f"{solver}_chi_eats_2d"]["sed_nufnu"]
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.4), constrained_layout=True, sharex=True, height_ratios=(1.45, 1.0))
    flux_values: list[np.ndarray] = []
    ratio_values: list[np.ndarray] = []
    model_curves = (
        ("1d", one_d, f"{one_d_solver}"),
        ("2d_legacy", two_d_legacy, f"{solver} legacy"),
        ("2d_chi", two_d_chi, f"{solver} chi EATS"),
    )
    for i_time, t_obs in enumerate(SED_TIMES):
        color = TIME_COLORS[float(t_obs)]
        for model_key, spectra, label in model_curves:
            _, linestyle = SED_MODEL_STYLES[model_key]
            y = spectra[:, i_time]
            flux_values.append(y)
            axes[0].loglog(SED_FREQUENCIES, y, color=color, ls=linestyle, lw=1.35)
        ratio_legacy = _safe_ratio(two_d_legacy[:, i_time], one_d[:, i_time])
        ratio_chi = _safe_ratio(two_d_chi[:, i_time], one_d[:, i_time])
        ratio_values.extend([ratio_legacy, ratio_chi])
        mask_legacy = np.isfinite(ratio_legacy) & (ratio_legacy > 0.0)
        mask_chi = np.isfinite(ratio_chi) & (ratio_chi > 0.0)
        if np.any(mask_legacy):
            axes[1].semilogx(SED_FREQUENCIES[mask_legacy], ratio_legacy[mask_legacy], color=color, ls=":", lw=1.5)
        if np.any(mask_chi):
            axes[1].semilogx(SED_FREQUENCIES[mask_chi], ratio_chi[mask_chi], color=color, ls="--", lw=1.6)
    _setup_log_axes(axes[0], ylabel=r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
    _style_ratio_axis(axes[1], "Observed frequency [Hz]", "2D / 1D")
    flux_lim = _positive_limits(flux_values)
    ratio_lim = _ratio_limits(ratio_values)
    if flux_lim:
        axes[0].set_ylim(*flux_lim)
    if ratio_lim:
        axes[1].set_ylim(*ratio_lim)
    leg1 = axes[0].legend(handles=_time_handles(), title="Epoch", loc="lower left")
    axes[0].add_artist(leg1)
    axes[0].legend(
        handles=[
            Line2D([0], [0], color="black", ls="-", lw=1.5, label=one_d_solver),
            Line2D([0], [0], color="black", ls=":", lw=1.5, label=f"{solver} legacy"),
            Line2D([0], [0], color="black", ls="--", lw=1.5, label=f"{solver} chi EATS"),
        ],
        title="Model",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0,
    )
    axes[1].legend(
        handles=[
            Line2D([0], [0], color="black", ls=":", lw=1.5, label="2D legacy / 1D"),
            Line2D([0], [0], color="black", ls="--", lw=1.5, label="2D chi EATS / 1D"),
        ],
        loc="upper left",
    )
    fig.suptitle(f"{solver}: 2D/1D total SED comparison, {medium.label}")
    return _save(fig, output_dir, f"{solver}_{medium.name}_sed_2d_vs_1d")


def _plot_tophat_theta_scan_solver(
    solver: str,
    theta_data: dict[str, dict[float, dict[str, np.ndarray]]],
    output_dir: Path,
    medium: MediumSpec,
) -> list[Path]:
    one_d_solver = _one_d_solver(solver)
    fig, axes = plt.subplots(
        2,
        THETA_RATIO_CASES.size,
        figsize=(15.5, 5.7),
        constrained_layout=True,
        sharex="col",
        height_ratios=(1.4, 1.0),
    )
    for col, theta_ratio in enumerate(THETA_RATIO_CASES):
        ax_flux = axes[0, col]
        ax_ratio = axes[1, col]
        one_d = theta_data[one_d_solver][float(theta_ratio)]["lightcurve_nufnu"]
        two_d = theta_data[solver][float(theta_ratio)]["lightcurve_nufnu"]
        panel_flux_values: list[np.ndarray] = []
        for band in LIGHTCURVE_BANDS:
            idx = int(np.where(LIGHTCURVE_BANDS == band)[0][0])
            color = BAND_COLORS[float(band)]
            one_curve = one_d[idx]
            two_curve = two_d[idx]
            ratio = _safe_ratio(two_curve, one_curve)
            dex_ratio = _log10_ratio(ratio)
            panel_flux_values.extend([one_curve, two_curve])
            ax_flux.loglog(LIGHTCURVE_TIMES, one_curve, color=color, ls="-", lw=1.25)
            ax_flux.loglog(LIGHTCURVE_TIMES, two_curve, color=color, ls="--", lw=1.25)
            mask = np.isfinite(dex_ratio)
            if np.any(mask):
                ax_ratio.semilogx(LIGHTCURVE_TIMES[mask], dex_ratio[mask], color=color, lw=1.25)
        ax_flux.set_title(fr"$\theta_v/\theta_j={theta_ratio:g}$")
        _setup_log_axes(ax_flux, ylabel=r"$\nu F_\nu$" if col == 0 else None)
        _style_dex_ratio_axis(ax_ratio, "Observer time [s]", r"$\log_{10}$(2D chi / 1D)" if col == 0 else "")
        if col == 0:
            ax_ratio.set_ylabel(r"$\log_{10}$(2D chi / 1D)")
        flux_lim = _bounded_display_limits(panel_flux_values)
        if flux_lim:
            ax_flux.set_ylim(*flux_lim)
    fig.legend(
        handles=_band_handles() + _method_handles(("one_d", "two_d_chi")),
        loc="lower center",
        bbox_to_anchor=(0.5, -0.015),
        ncol=5,
        columnspacing=1.5,
        handlelength=2.2,
    )
    fig.suptitle(f"{one_d_solver} vs {solver} chi-EATS total top-hat light curves, {medium.label}")
    return _save(fig, output_dir, f"{solver}_{medium.name}_tophat_theta_scan_2d_vs_1d")


def _collect_chi_grid_scan_case(solver: str, num_chi: int) -> dict[str, Any]:
    grid = dict(CHI_GRID_SCAN_GRID)
    grid["num_chi"] = int(num_chi)
    theta_v = CHI_GRID_SCAN_THETA_RATIO * THETA_J
    case = CaseSpec(solver=solver, geometry_kernel="chi_eats_2d")
    print(f"  collect chi-grid scan ISM {case.name} num_chi={num_chi}", flush=True)
    model = _build_model(case, grid, medium=ISM_SPEC, theta_j=THETA_J, theta_v=theta_v)
    tic = time.perf_counter()
    lightcurve_fnu = np.asarray(
        model.flux_density_grid(CHI_GRID_SCAN_TIMES, LIGHTCURVE_BANDS, projection_kind="lightcurve").total,
        dtype=float,
    )
    runtime_s = time.perf_counter() - tic
    return {"lightcurve_nufnu": LIGHTCURVE_BANDS[:, None] * lightcurve_fnu, "runtime_s": runtime_s, "grid": grid}


def _warm_chi_grid_scan_runtime(solver: str) -> None:
    warm_chi = int(CHI_GRID_SCAN_VALUES[0])
    grid = dict(CHI_GRID_SCAN_GRID)
    grid["num_chi"] = warm_chi
    theta_v = CHI_GRID_SCAN_THETA_RATIO * THETA_J
    case = CaseSpec(solver=solver, geometry_kernel="chi_eats_2d")
    print(f"  warm chi-grid scan ISM {case.name} num_chi={warm_chi}", flush=True)
    model = _build_model(case, grid, medium=ISM_SPEC, theta_j=THETA_J, theta_v=theta_v)
    model.flux_density_grid(CHI_GRID_SCAN_TIMES, LIGHTCURVE_BANDS, projection_kind="lightcurve")


def _plot_chi_grid_scan_solver(solver: str, scan_data: dict[int, dict[str, Any]], output_dir: Path) -> list[Path]:
    reference_chi = int(CHI_GRID_SCAN_VALUES[-1])
    reference = scan_data[reference_chi]["lightcurve_nufnu"]
    colors = mpl.colormaps["viridis"](np.linspace(0.12, 0.82, CHI_GRID_SCAN_VALUES.size))
    color_by_chi = {int(num_chi): colors[i] for i, num_chi in enumerate(CHI_GRID_SCAN_VALUES)}
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.0), constrained_layout=True)
    axes = axes.reshape(2, 2)

    runtime_values = np.array([scan_data[int(num_chi)]["runtime_s"] for num_chi in CHI_GRID_SCAN_VALUES], dtype=float)
    axes[0, 0].loglog(CHI_GRID_SCAN_VALUES, runtime_values, color="#222222", marker="o", ms=4.2, lw=1.35)
    axes[0, 0].set_xlabel(r"$N_\chi$")
    axes[0, 0].set_ylabel("Wall time [s]")
    axes[0, 0].set_title("Runtime")
    axes[0, 0].grid(which="major", alpha=0.18, lw=0.5)
    axes[0, 0].grid(which="minor", alpha=0.06, lw=0.35)

    ratio_axes = [axes[0, 1], axes[1, 0], axes[1, 1]]
    for i_band, band in enumerate(LIGHTCURVE_BANDS):
        ax = ratio_axes[i_band]
        panel_dex: list[np.ndarray] = []
        for num_chi in CHI_GRID_SCAN_VALUES[:-1]:
            data = scan_data[int(num_chi)]["lightcurve_nufnu"][i_band]
            ratio = _safe_ratio(data, reference[i_band])
            dex_ratio = _log10_ratio(ratio)
            mask = np.isfinite(dex_ratio)
            panel_dex.append(dex_ratio[mask])
            ax.semilogx(
                CHI_GRID_SCAN_TIMES[mask],
                dex_ratio[mask],
                color=color_by_chi[int(num_chi)],
                marker="o",
                markevery=8,
                ms=2.4,
                lw=1.15,
                label=fr"$N_\chi={int(num_chi)}$",
            )
        ax.axhline(0.0, color="#222222", lw=0.75, ls=":")
        ax.set_title(_label_frequency(float(band)))
        ax.set_xlabel("Observer time [s]")
        ax.set_ylabel(r"$\log_{10}(F_\chi/F_{512})$")
        finite_panel = np.concatenate([values for values in panel_dex if values.size])
        ymax = max(0.08, float(np.max(finite_panel)), 0.0)
        ymin = min(-0.08, float(np.min(finite_panel)), 0.0)
        margin = 0.08 * max(ymax - ymin, 1.0)
        ax.set_ylim(ymin - margin, ymax + margin)
        ax.grid(which="major", alpha=0.18, lw=0.5)
        ax.grid(which="minor", axis="x", alpha=0.06, lw=0.35)
    fig.legend(
        handles=[
            Line2D([0], [0], color=color_by_chi[int(num_chi)], marker="o", ms=3.0, lw=1.2, label=fr"$N_\chi={int(num_chi)}$")
            for num_chi in CHI_GRID_SCAN_VALUES[:-1]
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.015),
        ncol=4,
        columnspacing=1.1,
        handlelength=1.8,
    )
    fig.suptitle(
        rf"{solver} ISM chi-grid convergence, $\theta_v/\theta_j={CHI_GRID_SCAN_THETA_RATIO:g}$, "
        rf"half grid except $N_\chi$"
    )
    return _save(fig, output_dir, f"{solver}_ism_chi_grid_scan")


def _chi_grid_scan_rows(solver: str, scan_data: dict[int, dict[str, Any]]) -> list[dict[str, Any]]:
    reference_chi = int(CHI_GRID_SCAN_VALUES[-1])
    reference = scan_data[reference_chi]["lightcurve_nufnu"]
    rows: list[dict[str, Any]] = []
    flux_rows: list[dict[str, Any]] = []
    for num_chi in CHI_GRID_SCAN_VALUES:
        data = scan_data[int(num_chi)]["lightcurve_nufnu"]
        for i_band, band in enumerate(LIGHTCURVE_BANDS):
            ratio = _safe_ratio(data[i_band], reference[i_band])
            finite = ratio[np.isfinite(ratio) & (ratio > 0.0)]
            rows.append(
                {
                    "medium": "ism",
                    "solver": solver,
                    "num_chi": int(num_chi),
                    "band_hz": float(band),
                    "runtime_s": float(scan_data[int(num_chi)]["runtime_s"]),
                    "ratio_to_chi512_median": float(np.median(finite)),
                    "ratio_to_chi512_min": float(np.min(finite)),
                    "ratio_to_chi512_max": float(np.max(finite)),
                    "smoothness": _smoothness_metric(data[i_band]),
                }
            )
    return rows


def _chi_grid_scan_lightcurve_rows(solver: str, scan_data: dict[int, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for num_chi in CHI_GRID_SCAN_VALUES:
        data = scan_data[int(num_chi)]["lightcurve_nufnu"]
        runtime_s = float(scan_data[int(num_chi)]["runtime_s"])
        for i_band, band in enumerate(LIGHTCURVE_BANDS):
            for i_time, t_obs in enumerate(CHI_GRID_SCAN_TIMES):
                rows.append(
                    {
                        "medium": "ism",
                        "solver": solver,
                        "num_chi": int(num_chi),
                        "band_hz": float(band),
                        "observer_time_s": float(t_obs),
                        "nu_fnu": float(data[i_band, i_time]),
                        "runtime_s": runtime_s,
                    }
                )
    return rows


def _write_chi_grid_scan_summary(
    output_dir: Path,
    rows: list[dict[str, Any]],
    lightcurve_rows: list[dict[str, Any]],
    metadata: dict[str, Any],
) -> tuple[Path, Path, Path]:
    csv_path = output_dir / "chi_grid_scan_ism_summary.csv"
    raw_csv_path = output_dir / "chi_grid_scan_ism_lightcurves.csv"
    json_path = output_dir / "chi_grid_scan_ism_metadata.json"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    with raw_csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(lightcurve_rows[0].keys()))
        writer.writeheader()
        writer.writerows(lightcurve_rows)
    json_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
    print(csv_path)
    print(raw_csv_path)
    print(json_path)
    return csv_path, raw_csv_path, json_path


def _plot_chi_diagnostics(solver: str, data: dict[str, dict[str, Any]], output_dir: Path, medium: MediumSpec) -> list[Path]:
    chi = data[f"{solver}_chi_eats_2d"]
    t = chi["detail_time"]
    radius = chi["chi_radius_cm"]
    gamma = chi["chi_gamma_bulk"]
    weight = chi["chi_dvolume_weight"]
    tau = chi["tau_syn_chi"]
    emissivity = chi["l_syn_spec_chi"]
    seed_frequency = chi["seed_frequency_hz"]
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.6), constrained_layout=True)
    cmap = mpl.colormaps["viridis"]
    norm = mpl.colors.Normalize(vmin=1.0, vmax=float(radius.shape[0]))
    for idx in range(radius.shape[0]):
        color = cmap(norm(float(idx + 1)))
        active = weight[idx] > 0.0
        if np.any(active):
            axes[0, 0].loglog(t[active], radius[idx, active], color=color, lw=0.85)
            axes[0, 1].loglog(t[active], gamma[idx, active], color=color, lw=0.85)
            axes[1, 0].loglog(t[active], weight[idx, active], color=color, lw=0.85)
    tau_shell = np.sum(tau, axis=1)
    emiss_shell = np.sum(emissivity * weight[None, :, :], axis=1)
    for t_obs in SED_TIMES:
        i_shell = int(np.argmin(np.abs(t - t_obs)))
        axes[1, 1].loglog(seed_frequency, tau_shell[:, i_shell], color=TIME_COLORS[float(t_obs)], ls="-", lw=1.5)
        axes[1, 1].loglog(seed_frequency, emiss_shell[:, i_shell], color=TIME_COLORS[float(t_obs)], ls="--", lw=1.5)
    _setup_log_axes(axes[0, 0], ylabel=r"$r(\chi)$ [cm]")
    _setup_log_axes(axes[0, 1], ylabel=r"$\Gamma(\chi)$")
    _setup_log_axes(axes[1, 0], "Observer time [s]", r"$\chi\,\Delta\eta\,\ln 10$")
    _setup_log_axes(axes[1, 1], "Frequency [Hz]", "cell sums [arb.]")
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=[axes[0, 0], axes[0, 1], axes[1, 0]], fraction=0.035, pad=0.02)
    cbar.set_label("chi bin")
    cbar.ax.tick_params(length=2.5, width=0.6, labelsize=6.5)
    axes[1, 1].legend(
        handles=[
            Line2D([0], [0], color="black", ls="-", lw=1.3, label=r"$\sum_\chi\tau_\nu$"),
            Line2D([0], [0], color="black", ls="--", lw=1.3, label=r"$\sum_\chi P_\nu w_\chi$"),
        ],
        loc="upper right",
    )
    fig.suptitle(f"{solver}: chi-resolved geometry and radiation diagnostics, {medium.label}")
    return _save(fig, output_dir, f"{solver}_{medium.name}_chi_diagnostics")


def _smoothness_metric(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        arr = arr.reshape(-1)
    mask = np.isfinite(arr) & (arr > 0.0)
    if np.count_nonzero(mask) < 3:
        return float("nan")
    log_values = np.log(arr[mask])
    curvature = np.diff(log_values, n=2)
    return float(np.nanmax(np.abs(curvature)))


def _summary_rows(solver: str, data: dict[str, dict[str, Any]], medium: MediumSpec) -> list[dict[str, Any]]:
    legacy = data[f"{solver}_sed_legacy"]
    chi = data[f"{solver}_chi_eats_2d"]
    rows: list[dict[str, Any]] = []
    for band in LIGHTCURVE_BANDS:
        idx = int(np.where(LIGHTCURVE_BANDS == band)[0][0])
        ratio = _safe_ratio(chi["lightcurve_nufnu"][idx], legacy["lightcurve_nufnu"][idx])
        finite = ratio[np.isfinite(ratio) & (ratio > 0.0)]
        rows.append(
            {
                "medium": medium.name,
                "solver": solver,
                "quantity": "lightcurve",
                "coordinate": float(band),
                "ratio_median": float(np.median(finite)),
                "ratio_min": float(np.min(finite)),
                "ratio_max": float(np.max(finite)),
                "smoothness_legacy": _smoothness_metric(legacy["lightcurve_nufnu"][idx]),
                "smoothness_chi_eats_2d": _smoothness_metric(chi["lightcurve_nufnu"][idx]),
            }
        )
    for t_obs in SED_TIMES:
        idx = int(np.where(SED_TIMES == t_obs)[0][0])
        ratio = _safe_ratio(chi["sed_nufnu"][:, idx], legacy["sed_nufnu"][:, idx])
        finite = ratio[np.isfinite(ratio) & (ratio > 0.0)]
        rows.append(
            {
                "medium": medium.name,
                "solver": solver,
                "quantity": "sed",
                "coordinate": float(t_obs),
                "ratio_median": float(np.median(finite)),
                "ratio_min": float(np.min(finite)),
                "ratio_max": float(np.max(finite)),
                "smoothness_legacy": _smoothness_metric(legacy["sed_nufnu"][:, idx]),
                "smoothness_chi_eats_2d": _smoothness_metric(chi["sed_nufnu"][:, idx]),
            }
        )
    return rows


def _sed_2d_1d_summary_rows(
    solver: str,
    sed_data: dict[str, dict[str, np.ndarray]],
    medium: MediumSpec,
) -> list[dict[str, Any]]:
    one_d_solver = _one_d_solver(solver)
    one_d = sed_data[f"{one_d_solver}_sed_legacy"]["sed_nufnu"]
    two_d_legacy = sed_data[f"{solver}_sed_legacy"]["sed_nufnu"]
    two_d_chi = sed_data[f"{solver}_chi_eats_2d"]["sed_nufnu"]
    rows: list[dict[str, Any]] = []
    for t_obs in SED_TIMES:
        idx = int(np.where(SED_TIMES == t_obs)[0][0])
        for label, spectra in (("sed_2d_legacy_over_1d", two_d_legacy), ("sed_2d_chi_over_1d", two_d_chi)):
            ratio = _safe_ratio(spectra[:, idx], one_d[:, idx])
            finite = ratio[np.isfinite(ratio) & (ratio > 0.0)]
            rows.append(
                {
                    "medium": medium.name,
                    "solver": solver,
                    "quantity": label,
                    "coordinate": float(t_obs),
                    "ratio_median": float(np.median(finite)),
                    "ratio_min": float(np.min(finite)),
                    "ratio_max": float(np.max(finite)),
                    "smoothness_legacy": _smoothness_metric(one_d[:, idx]),
                    "smoothness_chi_eats_2d": _smoothness_metric(spectra[:, idx]),
                }
            )
    return rows


def _theta_scan_summary_rows(
    solver: str,
    theta_data: dict[str, dict[float, dict[str, np.ndarray]]],
    medium: MediumSpec,
) -> list[dict[str, Any]]:
    one_d_solver = _one_d_solver(solver)
    rows: list[dict[str, Any]] = []
    for theta_ratio in THETA_RATIO_CASES:
        one_d = theta_data[one_d_solver][float(theta_ratio)]["lightcurve_nufnu"]
        two_d = theta_data[solver][float(theta_ratio)]["lightcurve_nufnu"]
        for band in LIGHTCURVE_BANDS:
            idx = int(np.where(LIGHTCURVE_BANDS == band)[0][0])
            ratio = _safe_ratio(two_d[idx], one_d[idx])
            finite = ratio[np.isfinite(ratio) & (ratio > 0.0)]
            rows.append(
                {
                    "medium": medium.name,
                    "solver": solver,
                    "quantity": f"theta_scan_band_{band:.0e}",
                    "coordinate": float(theta_ratio),
                    "ratio_median": float(np.median(finite)),
                    "ratio_min": float(np.min(finite)),
                    "ratio_max": float(np.max(finite)),
                    "smoothness_legacy": _smoothness_metric(one_d[idx]),
                    "smoothness_chi_eats_2d": _smoothness_metric(two_d[idx]),
                }
            )
    return rows


def _write_summary(output_dir: Path, rows: list[dict[str, Any]], metadata: dict[str, Any]) -> tuple[Path, Path]:
    csv_path = output_dir / "chi_eats_2d_summary.csv"
    json_path = output_dir / "chi_eats_2d_metadata.json"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")
    print(csv_path)
    print(json_path)
    return csv_path, json_path


def _verify_outputs(paths: list[Path]) -> None:
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"benchmark output was not written: {path}")


def _frontier_values(mode: str) -> dict[str, tuple[int, ...]]:
    if mode == "formal":
        return {"num_chi": (16, 32, 64, 128), "num_theta": (40, 80, 160, 320), "num_phi": (8, 16, 32, 64)}
    return {"num_chi": (8, 16), "num_theta": (24, 48), "num_phi": (4, 8)}


def _frontier_grid(mode: str) -> dict[str, int]:
    if mode == "formal":
        return dict(num_gam_e=41, num_chi=32, num_nu=51, num_r=240, num_theta=80, num_phi=16, num_tobs=80)
    return dict(num_gam_e=21, num_chi=16, num_nu=25, num_r=80, num_theta=48, num_phi=8, num_tobs=32)


def _frontier_flux(grid: dict[str, int], theta_ratio: float, solver: str = "fullhide_2d") -> tuple[np.ndarray, float]:
    two_d = solver.endswith("_2d")
    model = Model(
        jet=top_hat_jet(energy_iso_erg=E_ISO, initial_lorentz_factor=GAMMA0, opening_angle_rad=THETA_J, shell_duration_s=None, magnetar=None, spreading=False),
        medium=UniformMedium(number_density_cm3=ISM_N),
        observer=Observer(z=REDSHIFT, viewing_angle_rad=theta_ratio * THETA_J, viewing_azimuth_rad=0.0, luminosity_distance_cm=LUMINOSITY_DISTANCE_CM),
        fwd_rad=radiation(epsilon_e=EPSILON_E, epsilon_B=EPSILON_B, p=ELECTRON_P, include_ssc=True),
        rvs_rad=None,
        numerics=numerics(num_radius=grid["num_r"], eats_num_theta=grid["num_theta"], eats_num_phi=grid["num_phi"] if theta_ratio else 1, num_observer_time=grid["num_tobs"], num_electron_gamma=grid["num_gam_e"], num_photon_frequency=grid["num_nu"], downstream_num_chi=grid["num_chi"] if two_d else None, num_threads=1, electron_adaptive_substeps=False, initial_radius_cm=1.0e14),
        observer_grid=observer_grid(time_min_s=float(CHI_GRID_SCAN_TIMES[0]), time_max_s=float(CHI_GRID_SCAN_TIMES[-1])),
        solver_options=solver_options(electron_solver=solver, geometry_projection="chi_eats_2d" if two_d else "sed_legacy", ssc_cooling_mode="numeric_ic_kn"),
        reverse_shock=reverse_shock(),
        hadronic=hadronic(),
    )
    start = time.perf_counter()
    flux = model.flux_density_grid(CHI_GRID_SCAN_TIMES, LIGHTCURVE_BANDS, projection_kind="lightcurve").total
    return LIGHTCURVE_BANDS[:, None] * np.asarray(flux), time.perf_counter() - start


def _frontier_error(values: np.ndarray, reference: np.ndarray) -> tuple[float, float, int]:
    mask = joint_significant(values, reference, 1.0e-10)
    error = np.abs(np.log10(values[mask] / reference[mask]))
    return float(np.median(error)), float(np.percentile(error, 95.0)), int(error.size)


def _run_frontier(mode: str, output_dir: Path | None, figure_dir: Path | None) -> None:
    base = _frontier_grid(mode)
    scans = _frontier_values(mode)
    data_dir = output_dir or DATA_ROOT / "chi_eats_limits"
    figure_dir = figure_dir or FIGURE_ROOT / "chi_eats_limits"
    rows: list[dict[str, Any]] = []
    flux_rows: list[dict[str, Any]] = []
    raw: dict[str, dict[int, tuple[np.ndarray, float]]] = {}
    for key, values in scans.items():
        raw[key] = {}
        for value in values:
            grid = dict(base)
            grid[key] = value
            raw[key][value] = _frontier_flux(grid, 0.5)
        reference = raw[key][values[-1]][0]
        uncertainty = _frontier_error(raw[key][values[-2]][0], reference)[1]
        for value in values:
            median, p95, count = _frontier_error(raw[key][value][0], reference)
            rows.append(dict(scan=key, value=value, theta_ratio=0.5, band="all", median_abs_log10=median, p95_abs_log10=p95, valid_count=count, runtime_s=raw[key][value][1], reference_uncertainty_p95=uncertainty))
            for band, curve in zip(LIGHTCURVE_BANDS, raw[key][value][0]):
                flux_rows.extend(dict(scan=key, value=value, theta_ratio=0.5, band_hz=float(band), observer_time_s=float(tobs), nu_fnu=float(flux)) for tobs, flux in zip(CHI_GRID_SCAN_TIMES, curve))
    reference_grid = dict(base)
    two_d, runtime = _frontier_flux(reference_grid, 0.0)
    one_grid = dict(reference_grid)
    one_grid["num_phi"] = 1
    one_d, _ = _frontier_flux(one_grid, 0.0, "fullhide_1d")
    median, p95, count = _frontier_error(two_d, one_d)
    rows.append(dict(scan="one_d_limit", value=base["num_chi"], theta_ratio=0.0, band="all", median_abs_log10=median, p95_abs_log10=p95, valid_count=count, runtime_s=runtime, reference_uncertainty_p95=float("nan")))
    reference_grid["num_theta"] = scans["num_theta"][-2]
    reference_grid["num_phi"] = scans["num_phi"][-2]
    for theta_ratio in (0.0, 0.5, 1.5):
        flux, runtime = _frontier_flux(reference_grid, theta_ratio)
        coarse_grid = dict(reference_grid)
        coarse_grid["num_theta"] = scans["num_theta"][-3]
        coarse_grid["num_phi"] = scans["num_phi"][-3]
        coarse, _ = _frontier_flux(coarse_grid, theta_ratio)
        radio = np.s_[0:1, :]
        median, p95, count = _frontier_error(coarse[radio], flux[radio])
        rows.append(dict(scan="view_angle", value=base["num_theta"], theta_ratio=theta_ratio, band="radio_ssa", median_abs_log10=median, p95_abs_log10=p95, valid_count=count, runtime_s=runtime, reference_uncertainty_p95=float("nan")))
    data_dir.mkdir(parents=True, exist_ok=True)
    csv_path = data_dir / "chi_eats_convergence.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (data_dir / "chi_eats_lightcurves.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(flux_rows[0]))
        writer.writeheader()
        writer.writerows(flux_rows)
    write_json(data_dir / "metadata.json", environment(mode, threads=1, grid=base, repeats=1) | {"total_flux_mask": "both arrays exceed 1e-10 of their own peak", "reference": "largest grid value; uncertainty is the p95 difference between the two largest values", "ssa_band_hz": float(LIGHTCURVE_BANDS[0])})
    plt.rcParams.update(plot_style())
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.0), constrained_layout=True)
    styles = {"num_chi": ("o", "-"), "num_theta": ("s", "--"), "num_phi": ("^", ":")}
    for key, (marker, linestyle) in styles.items():
        subset = [row for row in rows if row["scan"] == key]
        axes[0].loglog([row["value"] for row in subset], [row["p95_abs_log10"] for row in subset], marker=marker, linestyle=linestyle, label=key.replace("num_", r"$N_\mathrm{") + "}$")
        axes[1].loglog([row["runtime_s"] for row in subset], [row["p95_abs_log10"] for row in subset], marker=marker, linestyle=linestyle)
    axes[0].set(xlabel="Grid cells", ylabel=r"p95 $|\log_{10}(F/F_\mathrm{ref})|$")
    axes[1].set(xlabel="Wall time [s]", ylabel=r"p95 $|\log_{10}(F/F_\mathrm{ref})|$")
    for axis in axes:
        axis.grid(alpha=0.18)
    axes[0].legend()
    save_figure(fig, figure_dir / "chi_eats_convergence")
    plt.close(fig)


def main() -> None:
    args = _parse_args()
    if args.only == "convergence":
        _run_frontier(args.mode, args.output_dir if args.output_dir != OUTPUT_DIR else None, args.figure_dir)
        return
    grid = QUICK_GRID if args.mode == "quick" else FORMAL_GRID
    theta_grid = THETA_SCAN_QUICK_GRID if args.mode == "quick" else THETA_SCAN_FORMAL_GRID
    solvers = (args.solver,)
    media = _selected_media(args.medium)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(PLOT_STYLE)
    written: list[Path] = []
    rows: list[dict[str, Any]] = []
    if args.only == "chi-grid-scan":
        _warm_chi_grid_scan_runtime(args.solver)
        scan_data = {int(num_chi): _collect_chi_grid_scan_case(args.solver, int(num_chi)) for num_chi in CHI_GRID_SCAN_VALUES}
        written.extend(_plot_chi_grid_scan_solver(args.solver, scan_data, output_dir))
        scan_rows = _chi_grid_scan_rows(args.solver, scan_data)
        scan_metadata = {
            "git_head": _git_head(),
            "git_status": _git_status(),
            "git_diff_stat": _git_diff_stat(),
            "mode": args.mode,
            "only": args.only,
            "solver": args.solver,
            "medium": "ism",
            "geometry_kernel": "chi_eats_2d",
            "reference_num_chi": int(CHI_GRID_SCAN_VALUES[-1]),
            "runtime_warmup_num_chi": int(CHI_GRID_SCAN_VALUES[0]),
            "num_chi_values": [int(value) for value in CHI_GRID_SCAN_VALUES],
            "grid_except_num_chi": {key: int(value) for key, value in CHI_GRID_SCAN_GRID.items() if key != "num_chi"},
            "theta_v_over_theta_j": CHI_GRID_SCAN_THETA_RATIO,
            "physical_parameters": {
                "E_iso": E_ISO,
                "Gamma0": GAMMA0,
                "epsilon_e": EPSILON_E,
                "epsilon_B": EPSILON_B,
                "epsilon_B_floor": EPSILON_B_FLOOR,
                "magnetic_decay_alpha_t": MAGNETIC_DECAY_ALPHA_T,
                "magnetic_decay_t0_s": MAGNETIC_DECAY_T0_S,
                "p": ELECTRON_P,
                "theta_j": THETA_J,
                "z": REDSHIFT,
                "luminosity_distance_cm": LUMINOSITY_DISTANCE_CM,
                "ism_n": ISM_N,
            },
            "lightcurve_time_s": [float(CHI_GRID_SCAN_TIMES[0]), float(CHI_GRID_SCAN_TIMES[-1]), int(CHI_GRID_SCAN_TIMES.size)],
            "bands_hz": [float(value) for value in LIGHTCURVE_BANDS],
            "command": "rtk bash -lc 'source ~/.wsl_env && cd "
            + f'"{ROOT}" && uv run python scripts/benchmarks/chi_eats_2d_benchmark.py --mode {args.mode} --solver {args.solver} '
            + "--medium ism --only chi-grid-scan'",
            "acceptance": (
                "Fixed ISM fullhide_2d chi_eats_2d convergence scan. All non-chi grid dimensions are reduced "
                "relative to the quick benchmark: num_gam_e=16, num_nu=21, num_r=150, num_theta=150, "
                "num_phi=25, num_tobs=64. The plotted ratios use num_chi=512 as the reference and preserve "
                "raw ratio ranges in the CSV. Runtime measurements exclude one unplotted num_chi=32 warm-up run."
            ),
        }
        raw_scan_rows = _chi_grid_scan_lightcurve_rows(args.solver, scan_data)
        written.extend(_write_chi_grid_scan_summary(output_dir, scan_rows, raw_scan_rows, scan_metadata))
        _verify_outputs(written)
        return
    for medium in media:
        for solver in solvers:
            data: dict[str, dict[str, Any]] = {}
            sed_comparison_data: dict[str, dict[str, np.ndarray]] = {}
            theta_scan_data: dict[str, dict[float, dict[str, np.ndarray]]] = {}
            one_d_solver = _one_d_solver(solver)
            one_d_case = CaseSpec(solver=one_d_solver, geometry_kernel="sed_legacy")
            if args.only == "all":
                for geometry in ("sed_legacy", "chi_eats_2d"):
                    case = CaseSpec(solver=solver, geometry_kernel=geometry)
                    data[case.name] = _collect_case(case, grid, medium)
                    sed_comparison_data[case.name] = {"sed_nufnu": data[case.name]["sed_nufnu"]}
                sed_comparison_data[one_d_case.name] = _collect_sed_case(one_d_case, grid, medium)
                written.extend(_plot_lightcurve_solver(solver, data, output_dir, medium))
                written.extend(_plot_sed_solver(solver, data, output_dir, medium))
                written.extend(_plot_chi_diagnostics(solver, data, output_dir, medium))
                written.extend(_plot_1d_2d_sed_solver(solver, sed_comparison_data, output_dir, medium))
                rows.extend(_summary_rows(solver, data, medium))
                rows.extend(_sed_2d_1d_summary_rows(solver, sed_comparison_data, medium))

            theta_scan_data.setdefault(one_d_solver, {})
            theta_scan_data.setdefault(solver, {})
            two_d_theta_base = _collect_theta_scan_base_state(
                CaseSpec(solver=solver, geometry_kernel="chi_eats_2d"),
                theta_grid,
                medium=medium,
                theta_j=THETA_J,
            )
            for theta_ratio in THETA_RATIO_CASES:
                theta_v = float(theta_ratio) * THETA_J
                theta_scan_data[one_d_solver][float(theta_ratio)] = _collect_lightcurve_case(
                    one_d_case,
                    theta_grid,
                    medium=medium,
                    theta_j=THETA_J,
                    theta_v=theta_v,
                )
                print(f"  project LC {medium.name} {solver}_chi_eats_2d theta_v/theta_j={theta_v / THETA_J:.2f}", flush=True)
                theta_scan_data[solver][float(theta_ratio)] = _collect_lightcurve_from_state(
                    two_d_theta_base,
                    medium=medium,
                    theta_j=THETA_J,
                    theta_v=theta_v,
                )
            written.extend(_plot_tophat_theta_scan_solver(solver, theta_scan_data, output_dir, medium))
            rows.extend(_theta_scan_summary_rows(solver, theta_scan_data, medium))
    metadata = {
        "git_head": _git_head(),
        "git_status": _git_status(),
        "git_diff_stat": _git_diff_stat(),
        "mode": args.mode,
        "only": args.only,
        "solvers": list(solvers),
        "media": [medium.name for medium in media],
        "physical_parameters": {
            "E_iso": E_ISO,
            "Gamma0": GAMMA0,
            "epsilon_e": EPSILON_E,
            "epsilon_B": EPSILON_B,
            "epsilon_B_floor": EPSILON_B_FLOOR,
            "magnetic_decay_alpha_t": MAGNETIC_DECAY_ALPHA_T,
            "magnetic_decay_t0_s": MAGNETIC_DECAY_T0_S,
            "p": ELECTRON_P,
            "theta_j": THETA_J,
            "z": REDSHIFT,
            "luminosity_distance_cm": LUMINOSITY_DISTANCE_CM,
            "ism_n": ISM_N,
            "wind_Astar": WIND_ASTAR,
            "wind_n_ism": WIND_N_ISM,
        },
        "grid": grid,
        "theta_scan_grid": theta_grid,
        "on_axis_num_phi": 1,
        "lightcurve_time_s": [float(LIGHTCURVE_TIMES[0]), float(LIGHTCURVE_TIMES[-1]), int(LIGHTCURVE_TIMES.size)],
        "sed_frequency_hz": [float(SED_FREQUENCIES[0]), float(SED_FREQUENCIES[-1]), int(SED_FREQUENCIES.size)],
        "theta_v_over_theta_j": [float(value) for value in THETA_RATIO_CASES],
        "command": "rtk bash -lc 'source ~/.wsl_env && cd "
        + f'"{ROOT}" && uv run python scripts/benchmarks/chi_eats_2d_benchmark.py --mode {args.mode} --solver {args.solver} '
        + f"--medium {args.medium} --only {args.only}'",
        "acceptance": (
            "Total plotted flux includes forward SSC through the existing shell-level path. The chi_eats_2d change applies "
            "only to FS synchrotron+SSA; it may differ from sed_legacy near SSA and angular/finite-thickness sensitive "
            "epochs, but plotted total fluxes, break tracks, chi radius, local Gamma, and chi-weighted emission must remain "
            "finite and smooth. Light curves cover t_obs=1e2-1e9 s, SED covers 1e6-1e28 Hz, and the 2d/1d top-hat scan "
            "compares theta_v/theta_j = 0, 0.5, 1, 1.5, 3, 5. The 2D shell magnetic field is fixed to be constant "
            "inside each shell by setting epsilon_B_floor=epsilon_B and magnetic_decay_alpha_t=0."
        ),
    }
    summary_paths = _write_summary(output_dir, rows, metadata)
    written.extend(summary_paths)
    _verify_outputs(written)


if __name__ == "__main__":
    main()
