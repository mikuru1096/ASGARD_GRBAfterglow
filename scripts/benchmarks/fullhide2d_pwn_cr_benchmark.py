from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet, Wind
from asgard_core.asgard_paths import ASGARD_DOC_DIR

from VegasAfterglow import ISM as VegasISM
from VegasAfterglow import Model as VegasModel
from VegasAfterglow import Observer as VegasObserver
from VegasAfterglow import Radiation as VegasRadiation
from VegasAfterglow import TophatJet as VegasTophatJet
from VegasAfterglow import Wind as VegasWind


OUTPUT_DIR = ASGARD_DOC_DIR / "fullhide2d_pwn_cr_benchmark"

TIMES = np.array([1.0e3, 1.0e5, 1.0e7], dtype=float)
BANDS = np.array([1.0e9, 1.0e14, 1.0e18, 1.0e22], dtype=float)
LIGHTCURVE_TIMES = np.logspace(0.0, 8.0, 112)
FREQS = np.logspace(7.0, 30.0, 220)
ELECTRON_TARGET_GAMMA = np.logspace(np.log10(3.0), 8.0, 240)

EPSILON_E = 0.2
EPSILON_B = 1.0e-5
P_INDEX = 2.3
E_ISO = 3.0e52
GAMMA0 = 250.0
THETA_J = 0.1
LUMI_DIST = 1.0e26
REDSHIFT = 0.1
THETA_OBS = 0.0
N_ISM = 10.0
A_STAR = 1.0
WIND_N_ISM = 1.0e-15

FORMAL_GRID = dict(num_gam_e=121, num_chi=16, num_nu=121, num_r=160, num_theta=160, num_tobs=200)
QUICK_GRID = dict(num_gam_e=41, num_chi=6, num_nu=41, num_r=32, num_theta=32, num_tobs=32)

PLOT_STYLE = {
    "font.size": 10.0,
    "axes.labelsize": 10.5,
    "axes.titlesize": 11.0,
    "legend.fontsize": 8.0,
    "legend.title_fontsize": 8.0,
    "xtick.labelsize": 9.0,
    "ytick.labelsize": 9.0,
    "axes.linewidth": 0.95,
}
TIME_COLORS = {1.0e3: "#0072B2", 1.0e5: "#D55E00", 1.0e7: "#009E73"}
BAND_COLORS = {1.0e9: "#0072B2", 1.0e14: "#D55E00", 1.0e18: "#009E73", 1.0e22: "#CC79A7"}
SOLVER_STYLES = {
    "fullhide_1d": ("#111111", "-"),
    "fullhide_2d_legacy": ("#0072B2", ":"),
    "fullhide_2d_pwn_cr_v1": ("#D55E00", "--"),
    "charint_2d": ("#009E73", "-."),
    "VegasAfterglow": ("#666666", (0, (5, 2))),
}
SCAN_DECAY_ALPHA_T = -0.5
SCAN_DECAY_T0_S = 1.0
SCAN_EPSILON_B_FLOOR = 1.0e-8
SCAN_STOCHASTIC_ACCEL_NORM = 1.0e-1
SCAN_CASE_STYLES = {
    "pwn_closed": ("#0072B2", "-"),
    "pwn_free": ("#D55E00", "-"),
    "pwn_decay": ("#009E73", "--"),
    "pwn_free_decay": ("#CC79A7", "--"),
    "pwn_accel": ("#56B4E9", "-."),
    "pwn_free_accel": ("#E69F00", "-."),
    "pwn_decay_accel": ("#000000", ":"),
    "pwn_free_decay_accel": ("#999999", ":"),
}
OVERLAY_COLORS = {
    "VegasAfterglow": "#111111",
    "fullhide_2d_legacy": "#0072B2",
    "fullhide_2d_pwn_cr_v1": "#D55E00",
}
OVERLAY_LINESTYLES = {
    "VegasAfterglow": "-",
    "fullhide_2d_legacy": "--",
    "fullhide_2d_pwn_cr_v1": ":",
}
OVERLAY_MARKERS = {
    "VegasAfterglow": "o",
    "fullhide_2d_legacy": "s",
    "fullhide_2d_pwn_cr_v1": "^",
}
BAND_LINESTYLES = {1.0e9: "-", 1.0e18: "--"}
TIME_LINESTYLES = {1.0e3: "-", 1.0e7: "--"}


@dataclass(frozen=True)
class CaseSpec:
    name: str
    solver: str
    transport_model: str = "legacy"
    num_chi: int | None = None
    escape_mode: str = "closed"
    stochastic_accel_norm: float = 0.0
    magnetic_decay_alpha_t: float = 0.0
    magnetic_decay_t0_s: float = 1.0
    epsilon_b_floor: float | None = None
    label: str | None = None


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="fullhide_2d PWN/CR physical benchmark plots.")
    parser.add_argument("--mode", choices=("quick", "formal"), default="formal")
    parser.add_argument(
        "--targets",
        choices=("all", "legacy", "solver", "vegas", "scan"),
        default="all",
        help="plot group to generate; all collects every curve once and writes every figure",
    )
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    return parser.parse_args()


def _git_head() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return result.stdout.strip()


def _build_asgard_model(case: CaseSpec, medium_name: str, grid: dict[str, int]) -> Model:
    jet = TophatJet(theta_c=THETA_J, E_iso=E_ISO, Gamma0=GAMMA0)
    medium = ISM(n_ism=N_ISM) if medium_name == "ism" else Wind(A_star=A_STAR, n_ism=WIND_N_ISM)
    setups = Setups(
        electron_solver=case.solver,
        num_gam_e=grid["num_gam_e"],
        num_nu=grid["num_nu"],
        num_r=grid["num_r"],
        num_theta=grid["num_theta"],
        num_tobs=grid["num_tobs"],
        num_chi=case.num_chi,
        fullhide2d_transport_model=case.transport_model,
        fullhide2d_stochastic_accel_norm=case.stochastic_accel_norm,
        fullhide2d_escape_mode=case.escape_mode,
        electron_adaptive_substeps=False,
    )
    return Model(
        jet=jet,
        medium=medium,
        observer=Observer(LUMI_DIST, REDSHIFT, THETA_OBS),
        fwd_rad=Radiation(
            EPSILON_E,
            EPSILON_B,
            P_INDEX,
            epsilon_b_floor=case.epsilon_b_floor,
            magnetic_decay_alpha_t=case.magnetic_decay_alpha_t,
            magnetic_decay_t0_s=case.magnetic_decay_t0_s,
            ssc=True,
            kn=True,
        ),
        setups=setups,
    )


def _build_vegas_model(medium_name: str) -> VegasModel:
    jet = VegasTophatJet(theta_c=THETA_J, E_iso=E_ISO, Gamma0=GAMMA0)
    if medium_name == "ism":
        medium = VegasISM(n_ism=N_ISM)
    else:
        medium = VegasWind(A_star=A_STAR, n_ism=WIND_N_ISM)
    return VegasModel(
        jet=jet,
        medium=medium,
        observer=VegasObserver(lumi_dist=LUMI_DIST, z=REDSHIFT, theta_obs=THETA_OBS),
        fwd_rad=VegasRadiation(eps_e=EPSILON_E, eps_B=EPSILON_B, p=P_INDEX, ssc=True, kn=True),
        resolutions=(0.12, 0.25, 8),
    )


def _interp_positive_loglog(x_new: np.ndarray, x_old: np.ndarray, y_old: np.ndarray) -> np.ndarray:
    x_old = np.asarray(x_old, dtype=float)
    y_old = np.asarray(y_old, dtype=float)
    mask = np.isfinite(x_old) & np.isfinite(y_old) & (x_old > 0.0) & (y_old > 0.0)
    out = np.full_like(x_new, np.nan, dtype=float)
    if np.count_nonzero(mask) < 2:
        return out
    lx_new = np.log10(np.asarray(x_new, dtype=float))
    ly_new = np.interp(lx_new, np.log10(x_old[mask]), np.log10(y_old[mask]), left=np.nan, right=np.nan)
    finite = np.isfinite(ly_new)
    out[finite] = 10.0**ly_new[finite]
    return out


def _collect_asgard_case(model: Model) -> dict[str, np.ndarray]:
    print(f"  collect ASGARD {model.medium.kind}/{model.setups.electron_solver}/{model.setups.fullhide2d_transport_model}", flush=True)
    lc_fnu = np.asarray(model.flux_density_grid(LIGHTCURVE_TIMES, BANDS).total, dtype=float)
    sed_fnu = np.asarray(model.flux_density_grid(TIMES, FREQS).total, dtype=float)
    details = model.details(float(LIGHTCURVE_TIMES[0]), float(LIGHTCURVE_TIMES[-1]))
    t_detail = np.asarray(details.fwd.t_obs, dtype=float)
    gamma = np.asarray(details.fwd.gamma_e, dtype=float)
    dnde = np.asarray(details.fwd.dN_dgamma_e, dtype=float)
    electron = np.full((ELECTRON_TARGET_GAMMA.size, TIMES.size), np.nan, dtype=float)
    for i_time, t_obs in enumerate(TIMES):
        i_shell = int(np.argmin(np.abs(t_detail - t_obs)))
        electron[:, i_time] = _interp_positive_loglog(ELECTRON_TARGET_GAMMA, gamma, dnde[:, i_shell])
    return {
        "lightcurve_nufnu": BANDS[:, None] * lc_fnu,
        "sed_nufnu": FREQS[:, None] * sed_fnu,
        "electron_gamma": ELECTRON_TARGET_GAMMA.copy(),
        "electron_dnde": electron,
        "detail_time": t_detail,
        "nu_m": np.asarray(details.fwd.nu_m, dtype=float),
        "nu_c": np.asarray(details.fwd.nu_c, dtype=float),
        "nu_a": np.asarray(details.fwd.nu_a, dtype=float),
    }


def _collect_vegas_case(model: VegasModel, medium_name: str) -> dict[str, np.ndarray]:
    print(f"  collect VegasAfterglow {medium_name}", flush=True)
    lc_fnu = np.asarray(model.flux_density_grid(LIGHTCURVE_TIMES, BANDS).total, dtype=float)
    sed_fnu = np.asarray(model.flux_density_grid(TIMES, FREQS).total, dtype=float)
    return {
        "lightcurve_nufnu": BANDS[:, None] * lc_fnu,
        "sed_nufnu": FREQS[:, None] * sed_fnu,
    }


def _positive_limits(arrays: list[np.ndarray], pad: float = 0.35) -> tuple[float, float] | None:
    positive = []
    for arr in arrays:
        values = np.asarray(arr, dtype=float)
        values = values[np.isfinite(values) & (values > 0.0)]
        if values.size:
            panel_max = float(np.max(values))
            positive.append(values[values > panel_max * 1.0e-12])
    positive = [values for values in positive if values.size]
    if not positive:
        return None
    merged = np.concatenate(positive)
    lo = 10.0 ** (np.log10(float(np.min(merged))) - pad)
    hi = 10.0 ** (np.log10(float(np.max(merged))) + pad)
    return lo, hi


def _ratio_limits(arrays: list[np.ndarray], pad: float = 0.12) -> tuple[float, float] | None:
    positive = []
    for arr in arrays:
        values = np.asarray(arr, dtype=float)
        values = values[np.isfinite(values) & (values > 0.0)]
        if values.size:
            positive.append(values)
    if not positive:
        return None
    merged = np.concatenate(positive)
    span = max(abs(float(np.log10(np.min(merged)))), abs(float(np.log10(np.max(merged))))) + pad
    return 10.0 ** (-span), 10.0**span


def _electron_display_limits(arrays: list[np.ndarray], pad: float = 0.25, dynamic_dex: float = 8.0) -> tuple[float, float] | None:
    positive = []
    for arr in arrays:
        values = np.asarray(arr, dtype=float)
        values = values[np.isfinite(values) & (values > 0.0)]
        if values.size:
            positive.append(values)
    if not positive:
        return None
    merged = np.concatenate(positive)
    hi = float(np.quantile(merged, 0.995))
    if hi <= 0.0:
        return None
    lo_quantile = float(np.quantile(merged, 0.08))
    lo = max(lo_quantile, hi * 10.0 ** (-dynamic_dex))
    return 10.0 ** (np.log10(lo) - pad), 10.0 ** (np.log10(hi) + pad)


def _setup_log_axes(ax, xlabel: str | None = None, ylabel: str | None = None) -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.xaxis.set_major_locator(LogLocator(base=10.0))
    ax.yaxis.set_major_locator(LogLocator(base=10.0))
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.grid(which="major", alpha=0.18, lw=0.55)
    ax.grid(which="minor", alpha=0.08, lw=0.4)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)


def _safe_ratio(numer: np.ndarray, denom: np.ndarray) -> np.ndarray:
    numer = np.asarray(numer, dtype=float)
    denom = np.asarray(denom, dtype=float)
    return np.divide(numer, denom, out=np.full_like(numer, np.nan), where=(numer > 0.0) & (denom > 0.0))


def _label_frequency(value: float) -> str:
    exp = int(np.floor(np.log10(value)))
    base = value / 10.0**exp
    return fr"${base:.0f}\times10^{{{exp}}}$ Hz"


def _save(fig, output_dir: Path, stem: str) -> list[Path]:
    paths = [output_dir / f"{stem}.png", output_dir / f"{stem}.pdf"]
    fig.savefig(paths[0], dpi=220, bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight")
    plt.close(fig)
    for path in paths:
        print(path)
    return paths


def _plot_legacy_vs_pwncr(data: dict[str, dict[str, dict[str, np.ndarray]]], output_dir: Path) -> list[Path]:
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 8.2), constrained_layout=True)
    styles = {"fullhide_2d_legacy": "-", "fullhide_2d_pwn_cr_v1": "--"}
    labels = {"fullhide_2d_legacy": "legacy", "fullhide_2d_pwn_cr_v1": "pwn_cr_v1"}
    for col, medium_name in enumerate(("ism", "wind")):
        sed_ax = axes[0, col]
        ele_ax = axes[1, col]
        sed_values: list[np.ndarray] = []
        ele_values: list[np.ndarray] = []
        for i_time, t_obs in enumerate(TIMES):
            for case_name in ("fullhide_2d_legacy", "fullhide_2d_pwn_cr_v1"):
                color = TIME_COLORS[float(t_obs)]
                case = data[medium_name][case_name]
                sed = case["sed_nufnu"][:, i_time]
                electron = case["electron_dnde"][:, i_time]
                sed_values.append(sed)
                ele_values.append(electron)
                sed_ax.loglog(FREQS, sed, color=color, ls=styles[case_name], lw=1.8)
                ele_ax.loglog(case["electron_gamma"], electron, color=color, ls=styles[case_name], lw=1.8)
        sed_ax.set_title(f"{medium_name.upper()} SED")
        ele_ax.set_title(f"{medium_name.upper()} electron spectrum")
        _setup_log_axes(sed_ax, "Observed frequency [Hz]", r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
        _setup_log_axes(ele_ax, r"Electron Lorentz factor $\gamma_e$", r"$dN_e/d\gamma_e$")
        sed_lim = _positive_limits(sed_values)
        ele_lim = _electron_display_limits(ele_values)
        if sed_lim:
            sed_ax.set_ylim(*sed_lim)
        if ele_lim:
            ele_ax.set_ylim(*ele_lim)
    time_handles = [Line2D([0], [0], color=TIME_COLORS[float(t)], lw=2.0, label=fr"$t={t:.0e}$ s") for t in TIMES]
    model_handles = [Line2D([0], [0], color="black", lw=2.0, ls=styles[name], label=labels[name]) for name in styles]
    axes[0, 0].legend(handles=time_handles, title="Epoch", frameon=False, loc="lower left")
    axes[0, 1].legend(handles=model_handles, title="Transport", frameon=False, loc="lower left")
    fig.suptitle("fullhide_2d legacy vs PWN/CR transport")
    return _save(fig, output_dir, "legacy_vs_pwncr_sed_electron")


def _plot_lightcurve_breaks(data: dict[str, dict[str, dict[str, np.ndarray]]], output_dir: Path) -> list[Path]:
    fig, axes = plt.subplots(3, 2, figsize=(12.8, 10.0), constrained_layout=True, sharex="col")
    for col, medium_name in enumerate(("ism", "wind")):
        legacy = data[medium_name]["fullhide_2d_legacy"]
        pwn = data[medium_name]["fullhide_2d_pwn_cr_v1"]
        ax_lc = axes[0, col]
        ax_break = axes[1, col]
        ax_ratio = axes[2, col]
        for band in BANDS:
            idx = int(np.where(BANDS == band)[0][0])
            color = BAND_COLORS[float(band)]
            label = _label_frequency(float(band))
            ax_lc.loglog(LIGHTCURVE_TIMES, legacy["lightcurve_nufnu"][idx], color=color, ls="-", lw=1.7, label=f"legacy {label}")
            ax_lc.loglog(LIGHTCURVE_TIMES, pwn["lightcurve_nufnu"][idx], color=color, ls="--", lw=1.5, label=f"pwn_cr {label}")
            ratio = _safe_ratio(pwn["lightcurve_nufnu"][idx], legacy["lightcurve_nufnu"][idx])
            mask = np.isfinite(ratio) & (ratio > 0.0)
            if np.any(mask):
                ax_ratio.semilogx(LIGHTCURVE_TIMES[mask], ratio[mask], color=color, lw=1.4, label=label)
        for key, color, label in (("nu_m", "#0072B2", r"$\nu_m$"), ("nu_c", "#D55E00", r"$\nu_c$"), ("nu_a", "#009E73", r"$\nu_a$")):
            ax_break.loglog(legacy["detail_time"], legacy[key], color=color, ls="-", lw=1.7, label=f"legacy {label}")
            ax_break.loglog(pwn["detail_time"], pwn[key], color=color, ls="--", lw=1.5, label=f"pwn_cr {label}")
        ax_lc.set_title(f"{medium_name.upper()} light curves")
        ax_break.set_title(f"{medium_name.upper()} spectral breaks")
        _setup_log_axes(ax_lc, ylabel=r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
        _setup_log_axes(ax_break, ylabel="Frequency [Hz]")
        ax_ratio.axhline(1.0, color="black", lw=0.9, ls=":")
        ratio_values = [
            _safe_ratio(pwn["lightcurve_nufnu"][int(np.where(BANDS == band)[0][0])], legacy["lightcurve_nufnu"][int(np.where(BANDS == band)[0][0])])
            for band in BANDS
        ]
        ratio_lim = _positive_limits(ratio_values, pad=0.18)
        ax_ratio.set_yscale("log")
        if ratio_lim:
            ax_ratio.set_ylim(*ratio_lim)
        ax_ratio.set_xscale("log")
        ax_ratio.set_xlabel("Observer time [s]")
        ax_ratio.set_ylabel("pwn_cr / legacy")
        ax_ratio.grid(which="major", alpha=0.2)
        ax_lc.legend(frameon=False, fontsize=6.0, ncol=2)
        ax_break.legend(frameon=False, fontsize=7.0, ncol=2)
        ax_ratio.legend(frameon=False, fontsize=7.0, ncol=2)
    fig.suptitle("fullhide_2d transport light-curve and break-frequency diagnostics")
    return _save(fig, output_dir, "legacy_vs_pwncr_lightcurve_breaks")


def _plot_solver_compare(data: dict[str, dict[str, dict[str, np.ndarray]]], output_dir: Path) -> list[Path]:
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.2), constrained_layout=True)
    case_names = ("fullhide_1d", "fullhide_2d_legacy", "fullhide_2d_pwn_cr_v1", "charint_2d")
    for col, medium_name in enumerate(("ism", "wind")):
        sed_ax = axes[0, col]
        ele_ax = axes[1, col]
        sed_values: list[np.ndarray] = []
        ele_values: list[np.ndarray] = []
        for case_name in case_names:
            color, style = SOLVER_STYLES[case_name]
            sed = data[medium_name][case_name]["sed_nufnu"][:, 1]
            electron = data[medium_name][case_name]["electron_dnde"][:, 1]
            if case_name != "charint_2d":
                sed_values.append(sed)
                ele_values.append(electron)
            sed_ax.loglog(FREQS, sed, color=color, ls=style, lw=1.8)
            ele_ax.loglog(
                data[medium_name][case_name]["electron_gamma"],
                electron,
                color=color,
                ls=style,
                lw=1.8,
            )
        sed_ax.set_title(f"{medium_name.upper()} SED at 1e5 s")
        ele_ax.set_title(f"{medium_name.upper()} electron spectrum at 1e5 s")
        _setup_log_axes(sed_ax, "Observed frequency [Hz]", r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
        _setup_log_axes(ele_ax, r"Electron Lorentz factor $\gamma_e$", r"$dN_e/d\gamma_e$")
        sed_lim = _positive_limits(sed_values)
        ele_lim = _electron_display_limits(ele_values)
        if sed_lim:
            sed_ax.set_ylim(*sed_lim)
        if ele_lim:
            ele_ax.set_ylim(*ele_lim)
    handles = [Line2D([0], [0], color=SOLVER_STYLES[name][0], ls=SOLVER_STYLES[name][1], lw=2.0, label=name) for name in case_names]
    axes[0, 0].legend(handles=handles, frameon=False, loc="lower left")
    fig.suptitle("Solver comparison with PWN/CR transport overlay")
    return _save(fig, output_dir, "solver_transport_compare")


def _plot_vegas_overlay(
    data: dict[str, dict[str, dict[str, np.ndarray]]],
    vegas: dict[str, dict[str, np.ndarray]],
    output_dir: Path,
) -> list[Path]:
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.2), constrained_layout=True)
    line_defs = (
        ("VegasAfterglow", "VegasAfterglow", "Vegas"),
        ("fullhide_2d_legacy", "fullhide_2d_legacy", "legacy"),
        ("fullhide_2d_pwn_cr_v1", "fullhide_2d_pwn_cr_v1", "pwn_cr"),
    )
    for col, medium_name in enumerate(("ism", "wind")):
        lc_ax = axes[0, col]
        sed_ax = axes[1, col]
        for band in (1.0e9, 1.0e18):
            idx = int(np.where(BANDS == band)[0][0])
            for source_name, style_name, label_name in line_defs:
                source = vegas[medium_name] if source_name == "VegasAfterglow" else data[medium_name][source_name]
                lc_ax.loglog(
                    LIGHTCURVE_TIMES,
                    source["lightcurve_nufnu"][idx],
                    color=OVERLAY_COLORS[style_name],
                    ls=BAND_LINESTYLES[float(band)],
                    marker=OVERLAY_MARKERS[style_name],
                    markevery=14,
                    ms=3.2,
                    lw=1.55,
                    label=f"{label_name} {_label_frequency(float(band))}",
                )
        for i_time, t_obs in enumerate((1.0e3, 1.0e7)):
            for source_name, style_name, label_name in line_defs:
                source = vegas[medium_name] if source_name == "VegasAfterglow" else data[medium_name][source_name]
                sed_ax.loglog(
                    FREQS,
                    source["sed_nufnu"][:, i_time * 2],
                    color=OVERLAY_COLORS[style_name],
                    ls=TIME_LINESTYLES[float(t_obs)],
                    marker=OVERLAY_MARKERS[style_name],
                    markevery=32,
                    ms=3.2,
                    lw=1.55,
                    label=f"{label_name} t={t_obs:.0e}s",
                )
        lc_ax.set_title(f"{medium_name.upper()} light curve overlay")
        sed_ax.set_title(f"{medium_name.upper()} SED overlay")
        _setup_log_axes(lc_ax, "Observer time [s]", r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
        _setup_log_axes(sed_ax, "Observed frequency [Hz]", r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
        lc_ax.legend(frameon=False, fontsize=6.5, ncol=2)
        sed_ax.legend(frameon=False, fontsize=6.5, ncol=2)
    fig.suptitle("VegasAfterglow overlay for fullhide_2d transport")
    return _save(fig, output_dir, "vegas_overlay_transport")


def _scan_label(case: CaseSpec) -> str:
    if case.label:
        return case.label
    parts = [case.escape_mode]
    if case.magnetic_decay_alpha_t < 0.0:
        parts.append(fr"$\alpha_B={case.magnetic_decay_alpha_t:g}$")
    if case.stochastic_accel_norm > 0.0:
        parts.append(fr"$D_\gamma={case.stochastic_accel_norm:.0e}$")
    return ", ".join(parts)


def _plot_scan_lightcurves_absolute_group(
    data: dict[str, dict[str, dict[str, np.ndarray]]],
    cases: tuple[CaseSpec, ...],
    output_dir: Path,
    case_names: tuple[str, ...],
    stem: str,
    title: str,
) -> list[Path]:
    fig, axes = plt.subplots(len(BANDS), 2, figsize=(13.2, 11.8), constrained_layout=True, sharex="col")
    case_lookup = {case.name: case for case in cases}
    for col, medium_name in enumerate(("ism", "wind")):
        for row, band in enumerate(BANDS):
            ax = axes[row, col]
            idx = int(np.where(BANDS == band)[0][0])
            panel_values: list[np.ndarray] = []
            for case_name in case_names:
                if case_name == "fullhide_2d_legacy":
                    color, style, marker = "#111111", "-", None
                    label = "legacy"
                else:
                    color, style = SCAN_CASE_STYLES[case_name]
                    marker = "o"
                    label = _scan_label(case_lookup[case_name])
                values = data[medium_name][case_name]["lightcurve_nufnu"][idx]
                panel_values.append(values)
                ax.loglog(
                    LIGHTCURVE_TIMES,
                    values,
                    color=color,
                    ls=style,
                    lw=1.65,
                    marker=marker,
                    markevery=18,
                    ms=2.5,
                    label=label,
                )
            ax.set_title(f"{medium_name.upper()} {_label_frequency(float(band))}")
            _setup_log_axes(ax, ylabel=r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
            ylim = _positive_limits(panel_values, pad=0.22)
            if ylim:
                ax.set_ylim(*ylim)
            if row == len(BANDS) - 1:
                ax.set_xlabel("Observer time [s]")
    handles = []
    for case_name in case_names:
        if case_name == "fullhide_2d_legacy":
            handles.append(Line2D([0], [0], color="#111111", ls="-", lw=2.0, label="legacy"))
        else:
            handles.append(
                Line2D(
                    [0],
                    [0],
                    color=SCAN_CASE_STYLES[case_name][0],
                    ls=SCAN_CASE_STYLES[case_name][1],
                    marker="o",
                    ms=3.0,
                    lw=2.0,
                    label=_scan_label(case_lookup[case_name]),
                )
            )
    axes[0, 0].legend(handles=handles, frameon=False, fontsize=7.0, ncol=2, loc="best")
    fig.suptitle(title)
    return _save(fig, output_dir, stem)


def _plot_scan_lightcurves_absolute(
    data: dict[str, dict[str, dict[str, np.ndarray]]],
    cases: tuple[CaseSpec, ...],
    output_dir: Path,
) -> list[Path]:
    paths: list[Path] = []
    paths.extend(
        _plot_scan_lightcurves_absolute_group(
            data,
            cases,
            output_dir,
            ("fullhide_2d_legacy", "pwn_closed", "pwn_free", "pwn_accel", "pwn_free_accel"),
            "parameter_scan_lightcurves_absolute_transport_accel",
            "fullhide_2d PWN/CR parameter scan: absolute light curves, transport and acceleration",
        )
    )
    paths.extend(
        _plot_scan_lightcurves_absolute_group(
            data,
            cases,
            output_dir,
            (
                "fullhide_2d_legacy",
                "pwn_decay",
                "pwn_free_decay",
                "pwn_decay_accel",
                "pwn_free_decay_accel",
            ),
            "parameter_scan_lightcurves_absolute_decay",
            "fullhide_2d PWN/CR parameter scan: absolute light curves, microturbulence decay",
        )
    )
    return paths


def _plot_scan_lightcurve_ratios(
    data: dict[str, dict[str, dict[str, np.ndarray]]],
    cases: tuple[CaseSpec, ...],
    output_dir: Path,
) -> list[Path]:
    fig, axes = plt.subplots(len(BANDS), 2, figsize=(13.2, 11.8), constrained_layout=True, sharex="col")
    scan_cases = tuple(case for case in cases if case.name != "fullhide_2d_legacy")
    for col, medium_name in enumerate(("ism", "wind")):
        legacy = data[medium_name]["fullhide_2d_legacy"]["lightcurve_nufnu"]
        for row, band in enumerate(BANDS):
            ax = axes[row, col]
            idx = int(np.where(BANDS == band)[0][0])
            ratio_values: list[np.ndarray] = []
            for case in scan_cases:
                color, style = SCAN_CASE_STYLES[case.name]
                ratio = _safe_ratio(data[medium_name][case.name]["lightcurve_nufnu"][idx], legacy[idx])
                ratio_values.append(ratio)
                mask = np.isfinite(ratio) & (ratio > 0.0)
                if np.any(mask):
                    ax.semilogx(
                        LIGHTCURVE_TIMES[mask],
                        ratio[mask],
                        color=color,
                        ls=style,
                        lw=1.55,
                        marker="o",
                        markevery=18,
                        ms=2.4,
                    )
            ax.axhline(1.0, color="black", lw=0.85, ls=":")
            ax.set_yscale("log")
            ax.set_title(f"{medium_name.upper()} {_label_frequency(float(band))}")
            ax.set_ylabel("variant / legacy")
            ax.grid(which="major", alpha=0.2)
            ax.grid(which="minor", alpha=0.08)
            ylim = _ratio_limits(ratio_values)
            if ylim:
                ax.set_ylim(*ylim)
            if row == len(BANDS) - 1:
                ax.set_xlabel("Observer time [s]")
    handles = [
        Line2D([0], [0], color=SCAN_CASE_STYLES[case.name][0], ls=SCAN_CASE_STYLES[case.name][1], lw=2.0, label=_scan_label(case))
        for case in scan_cases
    ]
    axes[0, 0].legend(handles=handles, frameon=False, fontsize=7.0, ncol=2, loc="best")
    fig.suptitle("fullhide_2d PWN/CR parameter scan: light-curve ratios")
    return _save(fig, output_dir, "parameter_scan_lightcurve_ratios")


def _plot_scan_sed_electron_ratios(
    data: dict[str, dict[str, dict[str, np.ndarray]]],
    cases: tuple[CaseSpec, ...],
    output_dir: Path,
) -> list[Path]:
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.2), constrained_layout=True)
    scan_cases = tuple(case for case in cases if case.name != "fullhide_2d_legacy")
    time_index = int(np.where(TIMES == 1.0e5)[0][0])
    for col, medium_name in enumerate(("ism", "wind")):
        sed_ax = axes[0, col]
        ele_ax = axes[1, col]
        legacy_sed = data[medium_name]["fullhide_2d_legacy"]["sed_nufnu"][:, time_index]
        legacy_ele = data[medium_name]["fullhide_2d_legacy"]["electron_dnde"][:, time_index]
        sed_ratios: list[np.ndarray] = []
        ele_ratios: list[np.ndarray] = []
        for case in scan_cases:
            color, style = SCAN_CASE_STYLES[case.name]
            sed_ratio = _safe_ratio(data[medium_name][case.name]["sed_nufnu"][:, time_index], legacy_sed)
            ele_ratio = _safe_ratio(data[medium_name][case.name]["electron_dnde"][:, time_index], legacy_ele)
            sed_ratios.append(sed_ratio)
            ele_ratios.append(ele_ratio)
            sed_mask = np.isfinite(sed_ratio) & (sed_ratio > 0.0)
            ele_mask = np.isfinite(ele_ratio) & (ele_ratio > 0.0)
            sed_ax.loglog(FREQS[sed_mask], sed_ratio[sed_mask], color=color, ls=style, lw=1.6)
            ele_ax.loglog(
                data[medium_name][case.name]["electron_gamma"][ele_mask],
                ele_ratio[ele_mask],
                color=color,
                ls=style,
                lw=1.6,
            )
        for ax in (sed_ax, ele_ax):
            ax.axhline(1.0, color="black", lw=0.85, ls=":")
            ax.grid(which="major", alpha=0.2)
            ax.grid(which="minor", alpha=0.08)
            ax.tick_params(which="both", direction="in", top=True, right=True)
        sed_ax.set_title(f"{medium_name.upper()} SED ratio at 1e5 s")
        ele_ax.set_title(f"{medium_name.upper()} electron spectrum ratio at 1e5 s")
        sed_ax.set_xlabel("Observed frequency [Hz]")
        sed_ax.set_ylabel("variant / legacy")
        ele_ax.set_xlabel(r"Electron Lorentz factor $\gamma_e$")
        ele_ax.set_ylabel("variant / legacy")
        sed_ylim = _ratio_limits(sed_ratios)
        ele_ylim = _ratio_limits(ele_ratios)
        if sed_ylim:
            sed_ax.set_ylim(*sed_ylim)
        if ele_ylim:
            ele_ax.set_ylim(*ele_ylim)
    handles = [
        Line2D([0], [0], color=SCAN_CASE_STYLES[case.name][0], ls=SCAN_CASE_STYLES[case.name][1], lw=2.0, label=_scan_label(case))
        for case in scan_cases
    ]
    axes[0, 0].legend(handles=handles, frameon=False, fontsize=7.0, ncol=2, loc="best")
    fig.suptitle("fullhide_2d PWN/CR parameter scan: SED and electron ratios")
    return _save(fig, output_dir, "parameter_scan_sed_electron_ratios")


def _plot_scan_break_ratios(
    data: dict[str, dict[str, dict[str, np.ndarray]]],
    cases: tuple[CaseSpec, ...],
    output_dir: Path,
) -> list[Path]:
    break_defs = (("nu_m", r"$\nu_m$"), ("nu_c", r"$\nu_c$"), ("nu_a", r"$\nu_a$"))
    fig, axes = plt.subplots(3, 2, figsize=(13.0, 9.0), constrained_layout=True, sharex="col")
    scan_cases = tuple(case for case in cases if case.name != "fullhide_2d_legacy")
    for col, medium_name in enumerate(("ism", "wind")):
        legacy = data[medium_name]["fullhide_2d_legacy"]
        for row, (key, label) in enumerate(break_defs):
            ax = axes[row, col]
            ratio_values: list[np.ndarray] = []
            for case in scan_cases:
                color, style = SCAN_CASE_STYLES[case.name]
                payload = data[medium_name][case.name]
                ratio = _safe_ratio(payload[key], legacy[key])
                ratio_values.append(ratio)
                mask = np.isfinite(ratio) & (ratio > 0.0)
                if np.any(mask):
                    ax.semilogx(payload["detail_time"][mask], ratio[mask], color=color, ls=style, lw=1.55)
            ax.axhline(1.0, color="black", lw=0.85, ls=":")
            ax.set_yscale("log")
            ax.set_title(f"{medium_name.upper()} {label} ratio")
            ax.set_ylabel("variant / legacy")
            ax.grid(which="major", alpha=0.2)
            ax.grid(which="minor", alpha=0.08)
            ylim = _ratio_limits(ratio_values)
            if ylim:
                ax.set_ylim(*ylim)
            if row == 2:
                ax.set_xlabel("Observer time [s]")
    handles = [
        Line2D([0], [0], color=SCAN_CASE_STYLES[case.name][0], ls=SCAN_CASE_STYLES[case.name][1], lw=2.0, label=_scan_label(case))
        for case in scan_cases
    ]
    axes[0, 0].legend(handles=handles, frameon=False, fontsize=7.0, ncol=2, loc="best")
    fig.suptitle("fullhide_2d PWN/CR parameter scan: break-frequency ratios")
    return _save(fig, output_dir, "parameter_scan_break_ratios")


def _array_is_finite(values: np.ndarray) -> bool:
    arr = np.asarray(values, dtype=float)
    return bool(np.all(np.isfinite(arr)))


def _write_outputs(
    output_dir: Path,
    mode: str,
    data: dict[str, dict[str, dict[str, np.ndarray]]],
    vegas: dict[str, dict[str, np.ndarray]] | None,
    paths: list[Path],
) -> None:
    arrays: dict[str, np.ndarray] = {}
    summary: dict[str, object] = {
        "git_head": _git_head(),
        "mode": mode,
        "parameters": {
            "times_s": TIMES.tolist(),
            "bands_hz": BANDS.tolist(),
            "lightcurve_times_s": LIGHTCURVE_TIMES.tolist(),
            "frequency_minmax_hz": [float(FREQS[0]), float(FREQS[-1])],
            "electron_gamma_minmax": [float(ELECTRON_TARGET_GAMMA[0]), float(ELECTRON_TARGET_GAMMA[-1])],
            "physics": {
                "epsilon_e": EPSILON_E,
                "epsilon_b": EPSILON_B,
                "p": P_INDEX,
                "E_iso": E_ISO,
                "Gamma0": GAMMA0,
                "n_ism": N_ISM,
                "A_star": A_STAR,
                "wind_n_ism_floor": WIND_N_ISM,
            },
        },
        "figures": [str(path) for path in paths],
        "finite": {},
        "ratios": {},
    }
    finite_summary: dict[str, bool] = {}
    for medium_name, cases in data.items():
        for case_name, payload in cases.items():
            for key, values in payload.items():
                array_key = f"{medium_name}__{case_name}__{key}"
                arrays[array_key] = np.asarray(values)
                finite_summary[array_key] = _array_is_finite(values)
        pwn = cases["fullhide_2d_pwn_cr_v1"]["lightcurve_nufnu"]
        legacy = cases["fullhide_2d_legacy"]["lightcurve_nufnu"]
        ratio = _safe_ratio(pwn, legacy)
        arrays[f"{medium_name}__pwncr_over_legacy__lightcurve_nufnu"] = ratio
        finite_summary[f"{medium_name}__pwncr_over_legacy__lightcurve_nufnu"] = _array_is_finite(ratio[np.isfinite(ratio)])
        positive_ratio = ratio[np.isfinite(ratio) & (ratio > 0.0)]
        summary["ratios"][medium_name] = {
            "pwncr_over_legacy_lightcurve_min": float(np.min(positive_ratio)),
            "pwncr_over_legacy_lightcurve_max": float(np.max(positive_ratio)),
        }
    if vegas is not None:
        for medium_name, payload in vegas.items():
            for key, values in payload.items():
                array_key = f"{medium_name}__VegasAfterglow__{key}"
                arrays[array_key] = np.asarray(values)
                finite_summary[array_key] = _array_is_finite(values)
    summary["finite"] = finite_summary
    summary_path = output_dir / "benchmark_summary.json"
    arrays_path = output_dir / "benchmark_arrays.npz"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    np.savez_compressed(arrays_path, **arrays)
    print(summary_path)
    print(arrays_path)


def _ratio_stats(numer: np.ndarray, denom: np.ndarray) -> dict[str, float]:
    ratio = _safe_ratio(numer, denom)
    values = ratio[np.isfinite(ratio) & (ratio > 0.0)]
    return {
        "min": float(np.min(values)),
        "p05": float(np.percentile(values, 5.0)),
        "median": float(np.median(values)),
        "p95": float(np.percentile(values, 95.0)),
        "max": float(np.max(values)),
        "max_abs_log10": float(np.max(np.abs(np.log10(values)))),
    }


def _write_scan_outputs(
    output_dir: Path,
    mode: str,
    cases: tuple[CaseSpec, ...],
    data: dict[str, dict[str, dict[str, np.ndarray]]],
    paths: list[Path],
) -> None:
    arrays: dict[str, np.ndarray] = {}
    finite_summary: dict[str, bool] = {}
    summary_cases: dict[str, dict[str, Any]] = {}
    for case in cases:
        summary_cases[case.name] = {
            "label": _scan_label(case),
            "solver": case.solver,
            "transport_model": case.transport_model,
            "escape_mode": case.escape_mode,
            "stochastic_accel_norm": case.stochastic_accel_norm,
            "magnetic_decay_alpha_t": case.magnetic_decay_alpha_t,
            "magnetic_decay_t0_s": case.magnetic_decay_t0_s,
            "epsilon_b_floor": case.epsilon_b_floor,
        }
    summary: dict[str, Any] = {
        "git_head": _git_head(),
        "mode": mode,
        "scan_parameters": {
            "decay_alpha_t": SCAN_DECAY_ALPHA_T,
            "decay_t0_s": SCAN_DECAY_T0_S,
            "epsilon_b_floor": SCAN_EPSILON_B_FLOOR,
            "stochastic_accel_norm": SCAN_STOCHASTIC_ACCEL_NORM,
            "escape_modes": ["closed", "free_outer"],
            "wind_n_ism_floor": WIND_N_ISM,
            "lightcurve_start_s": float(LIGHTCURVE_TIMES[0]),
            "times_s": TIMES.tolist(),
            "bands_hz": BANDS.tolist(),
        },
        "cases": summary_cases,
        "figures": [str(path) for path in paths],
        "finite": finite_summary,
        "ratios": {},
    }
    for medium_name, cases_payload in data.items():
        summary["ratios"][medium_name] = {}
        legacy = cases_payload["fullhide_2d_legacy"]
        for case_name, payload in cases_payload.items():
            for key, values in payload.items():
                array_key = f"{medium_name}__{case_name}__{key}"
                arrays[array_key] = np.asarray(values)
                finite_summary[array_key] = _array_is_finite(values)
            if case_name == "fullhide_2d_legacy":
                continue
            summary["ratios"][medium_name][case_name] = {
                "lightcurve_nufnu": _ratio_stats(payload["lightcurve_nufnu"], legacy["lightcurve_nufnu"]),
                "sed_nufnu": _ratio_stats(payload["sed_nufnu"], legacy["sed_nufnu"]),
                "electron_dnde": _ratio_stats(payload["electron_dnde"], legacy["electron_dnde"]),
                "nu_m": _ratio_stats(payload["nu_m"], legacy["nu_m"]),
                "nu_c": _ratio_stats(payload["nu_c"], legacy["nu_c"]),
                "nu_a": _ratio_stats(payload["nu_a"], legacy["nu_a"]),
            }
            arrays[f"{medium_name}__{case_name}__over_legacy__lightcurve_nufnu"] = _safe_ratio(
                payload["lightcurve_nufnu"],
                legacy["lightcurve_nufnu"],
            )
            arrays[f"{medium_name}__{case_name}__over_legacy__sed_nufnu"] = _safe_ratio(
                payload["sed_nufnu"],
                legacy["sed_nufnu"],
            )
            arrays[f"{medium_name}__{case_name}__over_legacy__electron_dnde"] = _safe_ratio(
                payload["electron_dnde"],
                legacy["electron_dnde"],
            )
    summary_path = output_dir / "parameter_scan_summary.json"
    arrays_path = output_dir / "parameter_scan_arrays.npz"
    summary["finite"] = finite_summary
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    np.savez_compressed(arrays_path, **arrays)
    print(summary_path)
    print(arrays_path)


def _required_cases(grid: dict[str, int]) -> tuple[CaseSpec, ...]:
    return (
        CaseSpec("fullhide_1d", "fullhide_1d", num_chi=None),
        CaseSpec("fullhide_2d_legacy", "fullhide_2d", "legacy", grid["num_chi"]),
        CaseSpec("fullhide_2d_pwn_cr_v1", "fullhide_2d", "pwn_cr_v1", grid["num_chi"]),
        CaseSpec("charint_2d", "charint_2d", "legacy", grid["num_chi"]),
    )


def _scan_cases(grid: dict[str, int]) -> tuple[CaseSpec, ...]:
    n_chi = grid["num_chi"]
    return (
        CaseSpec("fullhide_2d_legacy", "fullhide_2d", "legacy", n_chi, label="legacy"),
        CaseSpec("pwn_closed", "fullhide_2d", "pwn_cr_v1", n_chi, label="closed"),
        CaseSpec("pwn_free", "fullhide_2d", "pwn_cr_v1", n_chi, escape_mode="free_outer", label="free_outer"),
        CaseSpec(
            "pwn_decay",
            "fullhide_2d",
            "pwn_cr_v1",
            n_chi,
            magnetic_decay_alpha_t=SCAN_DECAY_ALPHA_T,
            magnetic_decay_t0_s=SCAN_DECAY_T0_S,
            epsilon_b_floor=SCAN_EPSILON_B_FLOOR,
            label=fr"closed, $\alpha_B={SCAN_DECAY_ALPHA_T:g}$",
        ),
        CaseSpec(
            "pwn_free_decay",
            "fullhide_2d",
            "pwn_cr_v1",
            n_chi,
            escape_mode="free_outer",
            magnetic_decay_alpha_t=SCAN_DECAY_ALPHA_T,
            magnetic_decay_t0_s=SCAN_DECAY_T0_S,
            epsilon_b_floor=SCAN_EPSILON_B_FLOOR,
            label=fr"free, $\alpha_B={SCAN_DECAY_ALPHA_T:g}$",
        ),
        CaseSpec(
            "pwn_accel",
            "fullhide_2d",
            "pwn_cr_v1",
            n_chi,
            stochastic_accel_norm=SCAN_STOCHASTIC_ACCEL_NORM,
            label=fr"closed, $D_\gamma={SCAN_STOCHASTIC_ACCEL_NORM:.0e}$",
        ),
        CaseSpec(
            "pwn_free_accel",
            "fullhide_2d",
            "pwn_cr_v1",
            n_chi,
            escape_mode="free_outer",
            stochastic_accel_norm=SCAN_STOCHASTIC_ACCEL_NORM,
            label=fr"free, $D_\gamma={SCAN_STOCHASTIC_ACCEL_NORM:.0e}$",
        ),
        CaseSpec(
            "pwn_decay_accel",
            "fullhide_2d",
            "pwn_cr_v1",
            n_chi,
            stochastic_accel_norm=SCAN_STOCHASTIC_ACCEL_NORM,
            magnetic_decay_alpha_t=SCAN_DECAY_ALPHA_T,
            magnetic_decay_t0_s=SCAN_DECAY_T0_S,
            epsilon_b_floor=SCAN_EPSILON_B_FLOOR,
            label=fr"closed, decay+$D_\gamma$",
        ),
        CaseSpec(
            "pwn_free_decay_accel",
            "fullhide_2d",
            "pwn_cr_v1",
            n_chi,
            escape_mode="free_outer",
            stochastic_accel_norm=SCAN_STOCHASTIC_ACCEL_NORM,
            magnetic_decay_alpha_t=SCAN_DECAY_ALPHA_T,
            magnetic_decay_t0_s=SCAN_DECAY_T0_S,
            epsilon_b_floor=SCAN_EPSILON_B_FLOOR,
            label=fr"free, decay+$D_\gamma$",
        ),
    )


def main() -> None:
    args = _parse_args()
    plt.rcParams.update(PLOT_STYLE)
    output_dir = args.output_dir
    if args.targets == "scan" and output_dir == OUTPUT_DIR:
        output_dir = OUTPUT_DIR / "parameter_scan"
    output_dir.mkdir(parents=True, exist_ok=True)
    grid = FORMAL_GRID if args.mode == "formal" else QUICK_GRID
    media = ("ism", "wind")
    if args.targets == "scan":
        scan_cases = _scan_cases(grid)
        scan_data: dict[str, dict[str, dict[str, np.ndarray]]] = {medium_name: {} for medium_name in media}
        for medium_name in media:
            for case in scan_cases:
                scan_data[medium_name][case.name] = _collect_asgard_case(_build_asgard_model(case, medium_name, grid))
        scan_paths: list[Path] = []
        scan_paths.extend(_plot_scan_lightcurves_absolute(scan_data, scan_cases, output_dir))
        scan_paths.extend(_plot_scan_lightcurve_ratios(scan_data, scan_cases, output_dir))
        scan_paths.extend(_plot_scan_sed_electron_ratios(scan_data, scan_cases, output_dir))
        scan_paths.extend(_plot_scan_break_ratios(scan_data, scan_cases, output_dir))
        _write_scan_outputs(output_dir, args.mode, scan_cases, scan_data, scan_paths)
        return
    data: dict[str, dict[str, dict[str, np.ndarray]]] = {medium_name: {} for medium_name in media}
    for medium_name in media:
        for case in _required_cases(grid):
            data[medium_name][case.name] = _collect_asgard_case(_build_asgard_model(case, medium_name, grid))
    vegas_data = None
    if args.targets in ("all", "vegas"):
        vegas_data = {medium_name: _collect_vegas_case(_build_vegas_model(medium_name), medium_name) for medium_name in media}

    paths: list[Path] = []
    if args.targets in ("all", "legacy"):
        paths.extend(_plot_legacy_vs_pwncr(data, output_dir))
        paths.extend(_plot_lightcurve_breaks(data, output_dir))
    if args.targets in ("all", "solver"):
        paths.extend(_plot_solver_compare(data, output_dir))
    if args.targets in ("all", "vegas"):
        if vegas_data is None:
            raise RuntimeError("vegas target requires VegasAfterglow data.")
        paths.extend(_plot_vegas_overlay(data, vegas_data, output_dir))
    _write_outputs(output_dir, args.mode, data, vegas_data, paths)


if __name__ == "__main__":
    main()
