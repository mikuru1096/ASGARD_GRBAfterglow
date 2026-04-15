from __future__ import annotations

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matplotlib.colors import LogNorm

from VegasAfterglow import ISM, Model, Observer, Radiation, TophatJet, units
from asgard_presets import build_baseline_config
from asgard_runtime import solve_dynamics
from asgard_setup import build_simulation_setup
from src import constants
from tests.vegasafterglow_comprehensive_validation import REGIME_TEST_CONFIGS, _make_model, fit_powerlaw


OUTPUT_DIR = ROOT / "output" / "vegasafterglow_doc"


def _label_from_value(value: float, unit: str) -> str:
    exp = int(np.floor(np.log10(value)))
    base = value / 10**exp
    return fr"${base:.1f} \times 10^{{{exp}}}$ {unit}"


def _base_model() -> Model:
    medium = ISM(n_ism=1.0)
    jet = TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0)
    obs = Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0)
    rad = Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.3)
    return Model(jet=jet, medium=medium, observer=obs, fwd_rad=rad)


def plot_lightcurves() -> Path:
    model = _base_model()
    times = np.logspace(2.0, 8.0, 200)
    bands = np.array([1.0e9, 1.0e14, 1.0e17])
    results = model.flux_density_grid(times, bands)

    plt.figure(figsize=(4.8, 3.6), dpi=200)
    for i, nu in enumerate(bands):
        plt.loglog(times, results.total[i, :], label=_label_from_value(nu, "Hz"))
    plt.annotate("jet break", xy=(3.0e4, 1.0e-26), xytext=(3.0e3, 5.0e-28), arrowprops=dict(arrowstyle="->"))
    plt.annotate(r"$\nu_m=\nu_a$", xy=(6.0e5, 3.0e-25), xytext=(7.5e4, 5.0e-24), arrowprops=dict(arrowstyle="->"))
    plt.annotate(r"$\nu=\nu_a$", xy=(1.5e6, 4.0e-25), xytext=(7.5e5, 5.0e-24), arrowprops=dict(arrowstyle="->"))
    plt.xlabel("Time (s)")
    plt.ylabel("Flux Density (erg/cm²/s/Hz)")
    plt.legend()
    plt.title("Light Curves")
    plt.tight_layout()
    out = OUTPUT_DIR / "quick-lc.png"
    plt.savefig(out, dpi=300)
    plt.close()
    return out


def plot_spectra() -> Path:
    model = _base_model()
    bands = np.array([1.0e9, 1.0e14, 1.0e17])
    frequencies = np.logspace(5.0, 22.0, 200)
    epochs = np.array([1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8])
    results = model.flux_density_grid(epochs, frequencies)

    plt.figure(figsize=(4.8, 3.6), dpi=200)
    colors = plt.cm.viridis(np.linspace(0.0, 1.0, len(epochs)))
    for i, time_s in enumerate(epochs):
        plt.loglog(frequencies, results.total[:, i], color=colors[i], label=_label_from_value(time_s, "s"))
    for i, band in enumerate(bands):
        plt.axvline(band, ls="--", color=f"C{i}")
    plt.xlabel("frequency (Hz)")
    plt.ylabel("flux density (erg/cm²/s/Hz)")
    plt.legend(ncol=2)
    plt.title("Synchrotron Spectra")
    plt.tight_layout()
    out = OUTPUT_DIR / "quick-spec.png"
    plt.savefig(out, dpi=300)
    plt.close()
    return out


def plot_reverse_shock() -> Path:
    medium = ISM(n_ism=1.0)
    jet = TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0, duration=100.0)
    obs = Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0)
    fwd_rad = Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.3)
    rvs_rad = Radiation(eps_e=1.0e-2, eps_B=1.0e-4, p=2.4)
    model = Model(jet=jet, medium=medium, observer=obs, fwd_rad=fwd_rad, rvs_rad=rvs_rad, resolutions=(0.15, 0.5, 10))

    times = np.logspace(2.0, 8.0, 200)
    bands = np.array([1.0e9, 1.0e14, 1.0e17])
    results = model.flux_density_grid(times, bands)

    plt.figure(figsize=(4.8, 3.6), dpi=200)
    for i, nu in enumerate(bands):
        label = _label_from_value(nu, "Hz")
        plt.loglog(times, results.fwd.sync[i, :], label=f"{label} (fwd)")
        plt.loglog(times, results.rvs.sync[i, :], ls="--", label=f"{label} (rvs)")
    plt.xlabel("Time (s)")
    plt.ylabel("Flux Density (erg/cm²/s/Hz)")
    plt.legend(ncol=2, fontsize=8)
    plt.title("Reverse Shock Synchrotron")
    plt.tight_layout()
    out = OUTPUT_DIR / "reverse-shock-sync.png"
    plt.savefig(out, dpi=300)
    plt.close()
    return out


def plot_ssc() -> Path:
    medium = ISM(n_ism=1.0)
    jet = TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0)
    obs = Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0)
    rad = Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.3, ssc=True, kn=True)
    model = Model(jet=jet, medium=medium, observer=obs, fwd_rad=rad)

    times = np.logspace(2.0, 8.0, 200)
    bands = np.array([1.0e9, 1.0e14, 1.0e17])
    results = model.flux_density_grid(times, bands)

    plt.figure(figsize=(4.8, 3.6), dpi=200)
    for i, nu in enumerate(bands):
        label = _label_from_value(nu, "Hz")
        plt.loglog(times, results.fwd.sync[i, :], label=f"{label} (sync)")
        plt.loglog(times, results.fwd.ssc[i, :], ls="--", label=f"{label} (SSC)")
    plt.xlabel("Time (s)")
    plt.ylabel("Flux Density (erg/cm²/s/Hz)")
    plt.legend(ncol=2, fontsize=8)
    plt.title("Self-Synchrotron Compton")
    plt.tight_layout()
    out = OUTPUT_DIR / "ssc-components.png"
    plt.savefig(out, dpi=300)
    plt.close()
    return out


def plot_pairs_and_exposures() -> Path:
    model = _base_model()
    times = np.logspace(2.0, 8.0, 128)
    freqs = np.logspace(9.0, 17.0, 128)
    exposures = 0.2 * times
    pair_flux = model.flux_density(times, freqs)
    expo_flux = model.flux_density_exposures(times, freqs, exposures, num_subsamples=8)

    fig, axes = plt.subplots(2, 1, figsize=(4.8, 5.6), dpi=200, sharex=True)
    axes[0].loglog(times, pair_flux.total, label="paired points")
    axes[0].loglog(times, expo_flux.total, "--", label="exposure averaged")
    axes[0].set_ylabel("Flux Density (erg/cm²/s/Hz)")
    axes[0].legend()
    axes[0].set_title("Pair vs Exposure-Averaged Flux")

    rel = np.abs(expo_flux.total - pair_flux.total) / np.maximum(pair_flux.total, 1.0e-99)
    axes[1].semilogx(times, rel)
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Relative Difference")
    axes[1].grid(True, which="both", alpha=0.25)
    plt.tight_layout()
    out = OUTPUT_DIR / "pair-exposure-comparison.png"
    plt.savefig(out, dpi=300)
    plt.close()
    return out


def plot_supported_internal_evolution() -> Path:
    config = build_baseline_config(
        z=0.1,
        d_ne=1.0,
        opening_angle_jet=0.3,
        theta_v=0.0,
        e_iso=1.0e52,
        eta_0=100.0,
        epsilon_e=1.0e-1,
        epsilon_b=1.0e-3,
        p=2.3,
        num_tobs=300,
        num_r=300,
    )
    setup = build_simulation_setup(config)
    dynamics = solve_dynamics(setup.boundary, config)

    gamma = dynamics.r_gamma
    radius = dynamics.radius
    beta = np.sqrt(np.maximum(1.0 - gamma**-2, 1.0e-30))
    dr = np.diff(radius, prepend=radius[0])
    t_src = np.cumsum(dr / (beta * constants.para_c))
    t_comv = np.cumsum(dr / (beta * gamma * constants.para_c))
    t_obs = dynamics.r_tobs

    b_comv = 0.39 * np.sqrt(config.epsilon_b * config.d_ne * gamma * np.maximum(gamma - 1.0, 0.0))
    n_p = dynamics.swept_mass_g / constants.para_m_p
    n_e = config.f_e * n_p
    quantities = [
        (gamma, r"$\Gamma$"),
        (b_comv, r"$B^\prime$ [G]"),
        (n_p, r"$N_p$"),
        (radius, r"$r$ [cm]"),
        (n_e, r"$N_e$"),
    ]
    frames = [
        (t_src, r"$t_{\rm src}$ [s]", "source frame"),
        (t_comv, r"$t^\prime$ [s]", "comoving frame"),
        (t_obs, r"$t_{\rm obs}$ [s]", "observer frame"),
    ]
    fig = plt.figure(figsize=(12.6, 12.0), dpi=180)
    for i, (time_axis, xlabel, title) in enumerate(frames):
        for j, (values, ylabel) in enumerate(quantities):
            ax = plt.subplot(len(quantities), len(frames), j * len(frames) + i + 1)
            if j == 0:
                ax.set_title(title)
            ax.loglog(time_axis, values, color=f"C{i}", lw=1.8)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.grid(True, which="both", alpha=0.2)
    plt.tight_layout()
    out = OUTPUT_DIR / "shock-evolution-supported.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    return out


def plot_sky_image_single() -> Path:
    model = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=200.0),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.3),
        resolutions=(0.08, 0.15, 10),
    )
    img = model.sky_image([1.0e6], nu_obs=1.0e9, fov=500.0 * units.uas, npixel=128)
    fig, ax = plt.subplots(figsize=(4.5, 3.8), dpi=200)
    extent = img.extent / units.uas
    im = ax.imshow(img.image[0].T, origin="lower", extent=extent, cmap="inferno", norm=LogNorm())
    ax.set_xlabel(r"$\Delta x$ ($\mu$as)")
    ax.set_ylabel(r"$\Delta y$ ($\mu$as)")
    ax.set_title(r"$t_{\rm obs}=10^6$ s, $\nu=1$ GHz")
    fig.colorbar(im, label=r"Surface brightness (erg/cm$^2$/s/Hz/sr)")
    plt.tight_layout()
    out = OUTPUT_DIR / "sky_image_single.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    return out


def plot_sky_image_offaxis() -> Path:
    model = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=200.0),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.4),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.3),
        resolutions=(0.08, 0.15, 10),
    )
    times = np.array([1.0e5, 1.0e6, 1.0e7])
    imgs = model.sky_image(times, nu_obs=1.0e9, fov=5000.0 * units.uas, npixel=128)
    extent = imgs.extent / units.uas
    vmin = imgs.image[imgs.image > 0.0].min()
    vmax = imgs.image.max()
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8), dpi=200, constrained_layout=True)
    for i, (ax, time_s) in enumerate(zip(axes, times)):
        im = ax.imshow(
            imgs.image[i].T,
            origin="lower",
            extent=extent,
            cmap="inferno",
            norm=LogNorm(vmin=vmin, vmax=vmax),
        )
        ax.set_title(fr"$t_{{\rm obs}}=10^{{{int(np.log10(time_s))}}}$ s")
        ax.set_xlabel(r"$\Delta x$ ($\mu$as)")
        if i == 0:
            ax.set_ylabel(r"$\Delta y$ ($\mu$as)")
    fig.colorbar(im, ax=axes, label=r"erg/cm$^2$/s/Hz/sr", shrink=0.85)
    out = OUTPUT_DIR / "sky_image_offaxis.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_sky_image_flux_comparison() -> Path:
    model = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=200.0),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.3),
        resolutions=(0.08, 0.15, 10),
    )
    t_obs = np.logspace(3.0, 8.0, 30)
    nu_obs = 1.0e9
    img = model.sky_image(t_obs, nu_obs=nu_obs, fov=2000.0 * units.uas, npixel=128)
    flux_from_image = img.image.sum(axis=(1, 2)) * img.pixel_solid_angle
    flux_direct = model.flux_density_grid(t_obs, np.array([nu_obs])).total[0, :]

    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        figsize=(4.8, 5.0),
        dpi=200,
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )
    ax1.loglog(t_obs, flux_direct, "k-", label="flux_density_grid")
    ax1.loglog(t_obs, flux_from_image, "o", ms=4, color="C1", label="sky_image (integrated)")
    ax1.set_ylabel(r"Flux density (erg/cm$^2$/s/Hz)")
    ax1.legend()
    ratio = flux_from_image / flux_direct
    ax2.semilogx(t_obs, ratio, "o-", ms=4, color="C1")
    ax2.axhline(1.0, color="k", ls="--", lw=0.8)
    ax2.set_ylabel("image / direct")
    ax2.set_xlabel("Observer time (s)")
    ax2.set_ylim(0.90, 1.10)
    plt.tight_layout()
    out = OUTPUT_DIR / "sky_image_flux_comparison.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    return out


def plot_radio_slope_diagnostics() -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.8), dpi=200, constrained_layout=True)
    for ax, regime in zip(axes, ("I", "II")):
        cfg = REGIME_TEST_CONFIGS[regime]
        model = _make_model(
            cfg["medium"],
            resolutions=(0.3, 2.0, 15),
            eps_e=cfg.get("eps_e", 0.01),
            eps_B=cfg.get("eps_B", 0.01),
            xi_e=cfg.get("xi_e", 1.0),
            n_ism=cfg.get("n_ism", 1.0),
            A_star=cfg.get("A_star", 0.1),
        )
        t_obs = cfg["t"]
        details = model.details(t_obs * 0.9, t_obs * 1.1)
        idx = int(np.argmin(np.abs(np.asarray(details.fwd.t_obs) - t_obs)))
        nu_a = float(np.asarray(details.fwd.nu_a)[idx])
        nu_m = float(np.asarray(details.fwd.nu_m)[idx])

        if regime == "I":
            freq = np.logspace(np.log10(nu_a / 500.0), np.log10(nu_a / 20.0), 24)
            expected = 2.0
            label = r"$\nu < \nu_a$"
        else:
            freq = np.logspace(np.log10(nu_m * 6.0), np.log10(nu_a / 10.0), 24)
            expected = 2.5
            label = r"$\nu_m < \nu < \nu_a$"

        flux = model.flux_density(np.full_like(freq, t_obs), freq).total
        slope = fit_powerlaw(freq, flux)
        anchor = np.median(flux[(flux > 0.0) & np.isfinite(flux)])
        fit_line = anchor * (freq / np.median(freq)) ** expected

        ax.loglog(freq, flux, "o-", color="C0", label=f"measured, slope={slope:.3f}")
        ax.loglog(freq, fit_line, "--", color="C1", label=f"expected={expected:.1f}")
        ax.axvline(nu_a, color="k", ls="--", lw=0.8, label=r"$\nu_a$")
        if regime == "II":
            ax.axvline(nu_m, color="gray", ls=":", lw=0.8, label=r"$\nu_m$")
        ax.set_title(f"Regime {regime}: {label}")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel(r"Flux density (erg/cm$^2$/s/Hz)")
        ax.legend(fontsize=8)
        ax.grid(True, which="both", alpha=0.2)

    out = OUTPUT_DIR / "radio_slope_diagnostics.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outputs = [
        plot_lightcurves(),
        plot_spectra(),
        plot_reverse_shock(),
        plot_ssc(),
        plot_pairs_and_exposures(),
        plot_supported_internal_evolution(),
        plot_sky_image_single(),
        plot_sky_image_offaxis(),
        plot_sky_image_flux_comparison(),
        plot_radio_slope_diagnostics(),
    ]
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
