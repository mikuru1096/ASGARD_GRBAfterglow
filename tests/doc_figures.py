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

from asgard_paths import ASGARD_DOC_DIR
from matplotlib.colors import LogNorm

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet, units

MODE = "high" if "--high" in sys.argv else "quick"
GRID = {
    "quick": {"lc": 48, "spec": 64, "pair": 32, "expo": 32, "sky_pix": 48, "sky_times": 8, "sky_flux_times": 12, "gam": 81, "nu": 49, "r": 80, "theta": 80, "tobs": 48},
    "high": {"lc": 100, "spec": 100, "pair": 200, "expo": 200, "sky_pix": 128, "sky_times": 3, "sky_flux_times": 30, "gam": 201, "nu": 201, "r": 300, "theta": 300, "tobs": 200},
}[MODE]


OUTPUT_DIR = ASGARD_DOC_DIR


def _label_from_value(value: float, unit: str) -> str:
    exp = int(np.floor(np.log10(value)))
    base = value / 10**exp
    return fr"${base:.1f} \times 10^{{{exp}}}$ {unit}"


def _base_model() -> Model:
    medium = ISM(n_ism=1.0)
    jet = TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0)
    obs = Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0)
    rad = Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.3)
    return Model(
        jet=jet,
        medium=medium,
        observer=obs,
        fwd_rad=rad,
        setups=Setups(
            num_gam_e=GRID["gam"],
            num_nu=GRID["nu"],
            num_r=GRID["r"],
            num_theta=GRID["theta"],
            num_tobs=GRID["tobs"],
            electron_adaptive_substeps=False,
        ),
    )


def _res() -> tuple[float, float, int]:
    return (0.08, 0.15, 10) if MODE == "high" else (0.12, 0.25, 8)


def plot_lightcurves() -> Path:
    model = _base_model()
    times = np.logspace(2.0, 8.0, GRID["lc"])
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
    frequencies = np.logspace(5.0, 22.0, GRID["spec"])
    epochs = np.array([1.0e2, 1.0e4, 1.0e6, 1.0e8]) if MODE == "quick" else np.array([1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8])
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
    model = Model(
        jet=jet,
        medium=medium,
        observer=obs,
        fwd_rad=fwd_rad,
        rvs_rad=rvs_rad,
        setups=Setups(
            num_gam_e=GRID["gam"],
            num_nu=GRID["spec"],
            num_r=GRID["r"],
            num_theta=GRID["theta"],
            num_tobs=GRID["tobs"],
            electron_adaptive_substeps=False,
        ),
        resolutions=_res(),
    )

    times = np.logspace(2.0, 8.0, GRID["lc"])
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
    model = Model(
        jet=jet,
        medium=medium,
        observer=obs,
        fwd_rad=rad,
        setups=Setups(
            num_gam_e=GRID["gam"],
            num_nu=GRID["spec"],
            num_r=GRID["r"],
            num_theta=GRID["theta"],
            num_tobs=GRID["tobs"],
            electron_adaptive_substeps=False,
        ),
        resolutions=_res(),
    )

    times = np.logspace(2.0, 8.0, GRID["lc"])
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
    times = np.logspace(2.0, 8.0, GRID["pair"])
    freqs = np.logspace(9.0, 17.0, GRID["pair"])
    exposures = 0.2 * times
    pair_flux = model.flux_density(times, freqs)
    expo_flux = model.flux_density_exposures(times, freqs, exposures, num_subsamples=4 if MODE == "quick" else 8)

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


def plot_sky_image_single() -> Path:
    model = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=200.0),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.3),
        setups=Setups(
            num_gam_e=GRID["gam"],
            num_nu=GRID["spec"],
            num_r=GRID["r"],
            num_theta=GRID["theta"],
            num_tobs=GRID["sky_times"],
            electron_adaptive_substeps=False,
        ),
        resolutions=_res(),
    )
    img = model.sky_image([1.0e6], nu_obs=1.0e9, fov=500.0 * units.uas, npixel=GRID["sky_pix"])
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
        setups=Setups(
            num_gam_e=GRID["gam"],
            num_nu=GRID["spec"],
            num_r=GRID["r"],
            num_theta=GRID["theta"],
            num_tobs=GRID["sky_times"],
            electron_adaptive_substeps=False,
        ),
        resolutions=_res(),
    )
    times = np.array([1.0e5, 1.0e6, 1.0e7]) if MODE == "high" else np.array([1.0e5, 1.0e6, 1.0e7])
    imgs = model.sky_image(times, nu_obs=1.0e9, fov=5000.0 * units.uas, npixel=GRID["sky_pix"])
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
        setups=Setups(
            num_gam_e=GRID["gam"],
            num_nu=GRID["spec"],
            num_r=GRID["r"],
            num_theta=GRID["theta"],
            num_tobs=GRID["sky_flux_times"],
            electron_adaptive_substeps=False,
        ),
        resolutions=_res(),
    )
    t_obs = np.logspace(3.0, 8.0, GRID["sky_flux_times"])
    nu_obs = 1.0e9
    img = model.sky_image(t_obs, nu_obs=nu_obs, fov=2000.0 * units.uas, npixel=GRID["sky_pix"])
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


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plots = [
        ("lightcurves", plot_lightcurves),
        ("spectra", plot_spectra),
        ("reverse_shock", plot_reverse_shock),
        ("ssc", plot_ssc),
        ("pairs", plot_pairs_and_exposures),
        ("sky_single", plot_sky_image_single),
        ("sky_offaxis", plot_sky_image_offaxis),
        ("sky_flux", plot_sky_image_flux_comparison),
    ]
    total = len(plots)
    for idx, (name, fn) in enumerate(plots, start=1):
        print(f"[{idx}/{total}] {name} ...", flush=True)
        print(fn(), flush=True)


if __name__ == "__main__":
    main()
