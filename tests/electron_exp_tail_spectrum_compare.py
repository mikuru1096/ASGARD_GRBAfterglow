from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from VegasAfterglow import ISM, Model, Observer, Radiation, Setups, TophatJet, Wind
from VegasAfterglow.api import _build_fit_config_for_patch, _resolve_patch_state


OUTPUT_DIR = ROOT / "output" / "benchmark_exp_tail"
OUTPUT_NPZ = OUTPUT_DIR / "spectrum_compare_data.npz"
OUTPUT_PNG = OUTPUT_DIR / "spectrum_compare.png"
OUTPUT_PDF = OUTPUT_DIR / "spectrum_compare.pdf"

SOLVERS = ("fullhide", "slc1", "slc1_mmg2")
MEDIA = ("ism", "wind")
TIMES = np.array([1.0e3, 1.0e5, 1.0e7], dtype=float)
FREQS = np.logspace(7.0, 30.0, 220)
ELECTRON_TARGET_GAMMA = np.logspace(np.log10(3.0), 8.0, 240)
NUM_GAM_BY_SOLVER = {"fullhide": 121, "slc1": 32, "slc1_mmg2": 32}
IC_EPSILON_E = 0.2
IC_EPSILON_B = 1.0e-5
IC_P = 2.3
IC_E_ISO = 3.0e52
IC_GAMMA0 = 250.0
IC_N_ISM = 10.0
IC_A_STAR = 1.0
IC_WIND_N0 = 1.0e-6
PLOT_STYLE = {
    "font.size": 10.5,
    "axes.labelsize": 11.0,
    "axes.titlesize": 11.5,
    "legend.fontsize": 8.5,
    "legend.title_fontsize": 8.5,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "axes.linewidth": 0.95,
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


def _tight_log_ylim(ax, arrays: list[np.ndarray], lower_q: float = 0.03, upper_q: float = 0.97) -> None:
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
    y_max = np.max(merged)
    if y_lo <= 0.0 or y_hi <= y_lo:
        return
    ax.set_ylim(10.0 ** (np.log10(y_lo) - 0.35), max(10.0 ** (np.log10(y_hi) + 0.35), 1.0e2 * y_max))


def _apply_academic_log_axes(ax, xlim: tuple[float, float], ylabel: str | None = None) -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.xaxis.set_major_locator(LogLocator(base=10.0))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_major_locator(LogLocator(base=10.0))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.tick_params(which="both", direction="in", top=True, right=True, length=5.5, width=0.9, pad=4.0)
    ax.tick_params(which="minor", length=3.0, width=0.75)
    ax.grid(which="major", alpha=0.18, lw=0.55)
    ax.grid(which="minor", alpha=0.08, lw=0.4)
    if ylabel is not None:
        ax.set_ylabel(ylabel)


def _interp_positive_loglog(x_new: np.ndarray, x_old: np.ndarray, y_old: np.ndarray) -> np.ndarray:
    x_old = np.asarray(x_old, dtype=float)
    y_old = np.asarray(y_old, dtype=float)
    mask = np.isfinite(x_old) & np.isfinite(y_old) & (x_old > 0.0) & (y_old > 0.0)
    if np.count_nonzero(mask) < 2:
        return np.zeros_like(x_new, dtype=float)
    lx_old = np.log10(x_old[mask])
    ly_old = np.log10(y_old[mask])
    lx_new = np.log10(np.asarray(x_new, dtype=float))
    ly_new = np.interp(lx_new, lx_old, ly_old, left=np.nan, right=np.nan)
    out = np.zeros_like(x_new, dtype=float)
    finite = np.isfinite(ly_new)
    out[finite] = 10.0 ** ly_new[finite]
    return out


def _collect_electron(model: Model) -> tuple[np.ndarray, np.ndarray]:
    config = _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )
    state = _resolve_patch_state(model, config, TIMES, FREQS)
    electron_grid = np.zeros((ELECTRON_TARGET_GAMMA.shape[0], TIMES.shape[0]), dtype=float)
    for time_index, t_obs in enumerate(TIMES):
        shell_index = int(np.argmin(np.abs(state.dynamics.r_tobs - t_obs)))
        electron_grid[:, time_index] = _interp_positive_loglog(
            ELECTRON_TARGET_GAMMA,
            np.asarray(state.electron.gam_e, dtype=float),
            np.asarray(state.electron.d_n_gam_e[:, shell_index], dtype=float),
        )
    return ELECTRON_TARGET_GAMMA, electron_grid


def _collect() -> dict[str, dict[str, np.ndarray]]:
    data: dict[str, dict[str, np.ndarray]] = {}
    for medium_name in MEDIA:
        data[medium_name] = {}
        for solver in SOLVERS:
            model = _build_model(solver, medium_name)
            sed_grid = model.flux_density_grid(TIMES, FREQS).total
            gam_grid, electron_grid = _collect_electron(model)
            data[medium_name][solver] = {
                "sed": FREQS[:, None] * sed_grid,
                "gamma": gam_grid,
                "electron": electron_grid,
            }
    return data


def _plot(data: dict[str, dict[str, np.ndarray | dict[str, np.ndarray]]]) -> None:
    plt.rcParams.update(PLOT_STYLE)
    colors = {1.0e3: "#1f77b4", 1.0e5: "#ff7f0e", 1.0e7: "#2ca02c"}
    linestyles = {"fullhide": "-", "slc1": "--", "slc1_mmg2": ":"}

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12.5, 8.2),
        constrained_layout=True,
    )
    time_handles = [
        Line2D([0], [0], color=colors[float(t_obs)], lw=2.0, ls="-", label=fr"$t_{{\rm obs}}={t_obs:.0e}\,\rm s$")
        for t_obs in TIMES
    ]
    solver_handles = [
        Line2D([0], [0], color="black", lw=2.0, ls=linestyles[solver], label=solver)
        for solver in SOLVERS
    ]
    for column, medium_name in enumerate(MEDIA):
        sed_ax = axes[0, column]
        electron_ax = axes[1, column]
        sed_series = []
        electron_series = []
        for time_index, t_obs in enumerate(TIMES):
            for solver in SOLVERS:
                sed_values = data[medium_name][solver]["sed"][:, time_index]
                electron_values = data[medium_name][solver]["electron"][:, time_index]
                gam_values = data[medium_name][solver]["gamma"]
                sed_series.append(sed_values)
                electron_series.append(electron_values)
                sed_ax.loglog(
                    FREQS,
                    sed_values,
                    color=colors[float(t_obs)],
                    ls=linestyles[solver],
                    lw=1.8,
                    label=f"{solver}, t={t_obs:.0e}s",
                )
                electron_ax.loglog(
                    gam_values,
                    electron_values,
                    color=colors[float(t_obs)],
                    ls=linestyles[solver],
                    lw=1.9,
                )
        sed_ax.set_title(f"{medium_name.upper()} Forward Shock")
        electron_ax.set_xlabel(r"Electron Lorentz factor $\gamma_e$")
        _tight_log_ylim(sed_ax, sed_series)
        _tight_log_ylim(electron_ax, electron_series)
        _apply_academic_log_axes(sed_ax, (FREQS[0], FREQS[-1]))
        _apply_academic_log_axes(electron_ax, (ELECTRON_TARGET_GAMMA[0], ELECTRON_TARGET_GAMMA[-1]))
        sed_ax.set_xlabel(r"Observed frequency $\nu$ (Hz)")
    panel_labels = ("(a)", "(b)", "(c)", "(d)")
    for label, ax in zip(panel_labels, axes.ravel(order="C")):
        ax.text(
            0.04,
            0.95,
            label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=11.0,
        )
    axes[0, 0].set_ylabel(r"$\nu F_\nu\ \mathrm{(erg\ cm^{-2}\ s^{-1})}$")
    axes[1, 0].set_ylabel(r"$dN_e/d\gamma_e$")
    axes[0, 0].legend(handles=time_handles, loc="lower left", fontsize=8, frameon=False, title="Time")
    axes[0, 1].legend(handles=solver_handles, loc="lower left", fontsize=8, frameon=False, title="Solver")
    fig.suptitle("Forward-shock SED and electron-spectrum comparison", fontsize=13)
    fig.savefig(OUTPUT_PNG, dpi=220)
    fig.savefig(OUTPUT_PDF)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data = _collect()
    np.savez(
        OUTPUT_NPZ,
        times=TIMES,
        freqs=FREQS,
        electron_gamma=ELECTRON_TARGET_GAMMA,
        ism_fullhide=data["ism"]["fullhide"]["sed"],
        ism_slc1=data["ism"]["slc1"]["sed"],
        ism_slc1_mmg2=data["ism"]["slc1_mmg2"]["sed"],
        wind_fullhide=data["wind"]["fullhide"]["sed"],
        wind_slc1=data["wind"]["slc1"]["sed"],
        wind_slc1_mmg2=data["wind"]["slc1_mmg2"]["sed"],
        ism_electron_fullhide=data["ism"]["fullhide"]["electron"],
        ism_electron_slc1=data["ism"]["slc1"]["electron"],
        ism_electron_slc1_mmg2=data["ism"]["slc1_mmg2"]["electron"],
        wind_electron_fullhide=data["wind"]["fullhide"]["electron"],
        wind_electron_slc1=data["wind"]["slc1"]["electron"],
        wind_electron_slc1_mmg2=data["wind"]["slc1_mmg2"]["electron"],
    )
    _plot(data)
    print(OUTPUT_NPZ)
    print(OUTPUT_PNG)
    print(OUTPUT_PDF)


if __name__ == "__main__":
    main()
