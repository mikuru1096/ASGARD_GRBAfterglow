from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from ASGARD import Fitter, ISM, Model, Observer, Param, Radiation, Scale, Setups, TophatJet


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "doc" / "assets" / "tutorials"


def tutorial_model() -> Model:
    return Model(
        jet=TophatJet(E_iso=1.0e52, Gamma0=300.0, theta_j=0.1),
        medium=ISM(n_ism=1.0),
        observer=Observer(z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.3, xi_N=0.1, ssc=True),
        setups=Setups(
            electron_solver="fullhide_1d",
            num_r=72,
            num_tobs=72,
            num_gam_e=81,
            num_nu=81,
            observer_time_min_s=1.0e2,
            observer_time_max_s=1.0e7,
            ssc_cooling=True,
        ),
    )


def save_figure(fig: plt.Figure, name: str) -> None:
    path = OUT / name
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(path.relative_to(ROOT))


def positive_for_log(values: np.ndarray) -> np.ndarray:
    return np.where(values > 0.0, values, np.nan)


def plot_quick_light_curves(model: Model) -> None:
    times = np.logspace(2, 7, 80)
    freqs = np.array([1.0e9, 1.0e14, 1.0e18])
    result = model.flux_density_grid(times, freqs)

    fig, ax = plt.subplots(figsize=(5.4, 3.7))
    for i, nu in enumerate(freqs):
        ax.loglog(times, positive_for_log(result.total[i]), label=rf"${nu:.0e}\,\mathrm{{Hz}}$")
    ax.set_xlabel(r"$t_{\rm obs}\ {\rm (s)}$")
    ax.set_ylabel(r"$F_\nu\ {\rm (erg\ cm^{-2}\ s^{-1}\ Hz^{-1})}$")
    ax.set_title("ASGARD multi-band light curves")
    ax.legend(frameon=False)
    ax.grid(True, which="both", alpha=0.22)
    save_figure(fig, "quick_light_curves.png")


def plot_quick_spectra(model: Model) -> None:
    times = np.array([1.0e3, 1.0e5, 1.0e7])
    freqs = np.logspace(8, 22, 110)
    result = model.flux_density_grid(times, freqs, projection_kind="sed")

    fig, ax = plt.subplots(figsize=(5.4, 3.7))
    for i, t_obs in enumerate(times):
        ax.loglog(freqs, positive_for_log(result.total[:, i]), label=rf"$t={t_obs:.0e}\,\mathrm{{s}}$")
    for nu in (1.0e9, 1.0e14, 1.0e18):
        ax.axvline(nu, color="0.55", lw=0.8, ls="--", alpha=0.6)
    ax.set_xlabel(r"$\nu\ {\rm (Hz)}$")
    ax.set_ylabel(r"$F_\nu\ {\rm (erg\ cm^{-2}\ s^{-1}\ Hz^{-1})}$")
    ax.set_title("ASGARD broadband spectra")
    ax.legend(frameon=False)
    ax.grid(True, which="both", alpha=0.22)
    save_figure(fig, "quick_spectra.png")


def plot_component_breakdown(model: Model) -> None:
    times = np.logspace(2, 7, 80)
    freq = np.array([1.0e14])
    result = model.flux_density_grid(times, freq)

    fig, ax = plt.subplots(figsize=(5.4, 3.7))
    ax.loglog(times, positive_for_log(result.total[0]), label="total", lw=2.1, color="black")
    ax.loglog(times, positive_for_log(result.fwd.sync[0]), label="FS synchrotron")
    ax.loglog(times, positive_for_log(result.fwd.ssc[0]), label="FS SSC")
    ax.set_xlabel(r"$t_{\rm obs}\ {\rm (s)}$")
    ax.set_ylabel(r"$F_\nu(10^{14}{\rm Hz})$")
    ax.set_title("Radiation component breakdown")
    ax.legend(frameon=False)
    ax.grid(True, which="both", alpha=0.22)
    save_figure(fig, "component_breakdown.png")


def plot_internal_quantities(model: Model) -> None:
    details = model.details(1.0e2, 1.0e7)
    track = details.fwd

    fig, axes = plt.subplots(2, 1, figsize=(5.8, 5.8), sharex=True)
    axes[0].loglog(track.radius, positive_for_log(track.Gamma), label=r"$\Gamma$")
    axes[0].loglog(track.radius, positive_for_log(track.B_comv), label=r"$B'$")
    axes[0].set_ylabel(r"$\Gamma,\ B'{\rm (G)}$")
    axes[0].legend(frameon=False)
    axes[0].grid(True, which="both", alpha=0.22)

    axes[1].loglog(track.radius, positive_for_log(track.nu_m), label=r"$\nu_m$")
    axes[1].loglog(track.radius, positive_for_log(track.nu_c), label=r"$\nu_c$")
    axes[1].loglog(track.radius, positive_for_log(track.nu_a), label=r"$\nu_a$")
    axes[1].set_xlabel(r"$R\ {\rm (cm)}$")
    axes[1].set_ylabel(r"$\nu\ {\rm (Hz)}$")
    axes[1].legend(frameon=False)
    axes[1].grid(True, which="both", alpha=0.22)
    fig.suptitle("Forward-shock internal evolution")
    save_figure(fig, "internal_quantities.png")


def plot_synthetic_fit(model: Model) -> None:
    times = np.logspace(3, 6, 12)
    freq = np.full_like(times, 1.0e14)
    truth = model.flux_density(times, freq).total
    phase = np.linspace(0.0, 2.0 * np.pi, times.size)
    observed = truth * (1.0 + 0.12 * np.sin(phase))
    sigma = 0.15 * truth

    fitter = Fitter(model=model)
    fitter.add_flux_density(times, freq, observed, sigma)
    fitter.params = [
        Param("logE", "jet.E_iso", 50.0, 54.0, Scale.LOG10),
        Param("logn", "medium.n_ism", -4.0, 2.0, Scale.LOG10),
        Param("logepsB", "fwd_rad.eps_B", -6.0, -1.0, Scale.LOG10),
        Param("p", "fwd_rad.p", 2.0, 2.8, Scale.LINEAR),
    ]
    loglike = fitter.loglike({"logE": 52.0, "logn": 0.0, "logepsB": -3.0, "p": 2.3})

    dense_times = np.logspace(3, 6, 90)
    dense_freq = np.full_like(dense_times, 1.0e14)
    dense_flux = model.flux_density(dense_times, dense_freq).total

    fig, ax = plt.subplots(figsize=(5.4, 3.7))
    ax.errorbar(times, observed, yerr=sigma, fmt="o", ms=4, capsize=2, label="synthetic data")
    ax.loglog(dense_times, positive_for_log(dense_flux), label=fr"baseline, $\ln\mathcal{{L}}={loglike:.1f}$")
    ax.set_xlabel(r"$t_{\rm obs}\ {\rm (s)}$")
    ax.set_ylabel(r"$F_\nu(10^{14}{\rm Hz})$")
    ax.set_title("Synthetic fitting check")
    ax.legend(frameon=False)
    ax.grid(True, which="both", alpha=0.22)
    save_figure(fig, "synthetic_fit.png")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    model = tutorial_model()
    plot_quick_light_curves(model)
    plot_quick_spectra(model)
    plot_component_breakdown(model)
    plot_internal_quantities(model)
    plot_synthetic_fit(model)


if __name__ == "__main__":
    main()
