from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from VegasAfterglow import ISM, Model, Observer, Radiation, Setups, TophatJet, Wind


OUTPUT_DIR = ROOT / "output" / "benchmark_exp_tail"

BANDS = {
    "Radio": 3.0e9,
    "Optical": 4.6e14,
    "X-ray": 2.417989e17,
    "GeV": 2.417989e23,
    "TeV": 2.417989e26,
}
SOLVERS = ("fullhide", "slc1", "slc1_mmg2")
LINESTYLES = {"fullhide": "-", "slc1": "--", "slc1_mmg2": ":"}
NUM_GAM_BY_SOLVER = {"fullhide": 121, "slc1": 32, "slc1_mmg2": 32}
IC_EPSILON_E = 0.2
IC_EPSILON_B = 1.0e-5
IC_P = 2.3
IC_E_ISO = 3.0e52
IC_GAMMA0 = 250.0
IC_N_ISM = 10.0
IC_A_STAR = 1.0
IC_WIND_N0 = 1.0e-6


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


def _plot_multiband_lightcurves() -> Path:
    times = np.logspace(0.0, 8.0, 220)
    freqs = np.array(list(BANDS.values()), dtype=float)
    colors = plt.cm.viridis(np.linspace(0.0, 1.0, len(freqs)))

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=220, constrained_layout=True, sharey=True)
    for ax, medium_name in zip(axes, ("ism", "wind")):
        series = []
        for solver in SOLVERS:
            model = _build_model(solver, medium_name)
            result = model.flux_density_grid(times, freqs).total
            for i, (label, nu) in enumerate(BANDS.items()):
                series.append(result[i, :])
                ax.loglog(times, result[i, :], color=colors[i], ls=LINESTYLES[solver], lw=1.8,
                          label=f"{label} ({solver})")
        ax.set_title(f"{medium_name.upper()} IC-dominated Multi-band Light Curves")
        ax.set_xlabel("Time (s)")
        ax.grid(True, which="both", alpha=0.2)
        _tight_log_ylim(ax, series)
    axes[0].set_ylabel("Flux Density (erg/cm$^2$/s/Hz)")
    axes[0].legend(fontsize=7, ncol=2)
    out = OUTPUT_DIR / "docstyle_multiband_lightcurves.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def _plot_wind_lightcurves() -> Path:
    times = np.logspace(0.0, 8.0, 240)
    freqs = np.array(list(BANDS.values()), dtype=float)
    model = _build_model("fullhide", "wind")
    result = model.flux_density_grid(times, freqs).total

    plt.figure(figsize=(4.8, 3.6), dpi=220)
    series = []
    for i, (label, nu) in enumerate(BANDS.items()):
        series.append(result[i, :])
        plt.loglog(times, result[i, :], label=f"{label} ({nu:.1e} Hz)")
    _tight_log_ylim(plt.gca(), series)
    plt.xlabel("Time (s)")
    plt.ylabel("Flux Density (erg/cm$^2$/s/Hz)")
    plt.title("Wind IC-dominated Light Curves")
    plt.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    out = OUTPUT_DIR / "wind_quick_lc.png"
    plt.savefig(out, dpi=300)
    plt.close()
    return out


def _plot_wind_spectra() -> Path:
    freqs = np.logspace(7.0, 30.0, 240)
    epochs = np.array([1.0e3, 1.0e5, 1.0e7], dtype=float)
    model = _build_model("fullhide", "wind")
    result = model.flux_density_grid(epochs, freqs).total
    sed = freqs[:, None] * result
    colors = plt.cm.viridis(np.linspace(0.0, 1.0, len(epochs)))

    plt.figure(figsize=(4.8, 3.6), dpi=220)
    series = []
    for i, t_obs in enumerate(epochs):
        series.append(sed[:, i])
        plt.loglog(freqs, sed[:, i], color=colors[i], label=f"{t_obs:.0e} s")
    for label, nu in BANDS.items():
        plt.axvline(nu, ls="--", lw=0.9, alpha=0.45)
    _tight_log_ylim(plt.gca(), series)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("SED $\\nu F_\\nu$ (erg/cm$^2$/s)")
    plt.title("Wind IC-dominated SED")
    plt.legend(fontsize=8)
    plt.tight_layout()
    out = OUTPUT_DIR / "wind_quick_spec.png"
    plt.savefig(out, dpi=300)
    plt.close()
    return out


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outputs = [
        _plot_multiband_lightcurves(),
        _plot_wind_lightcurves(),
        _plot_wind_spectra(),
    ]
    for out in outputs:
        print(out)


if __name__ == "__main__":
    main()
