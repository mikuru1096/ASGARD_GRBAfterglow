from __future__ import annotations

from pathlib import Path
import ssl
import sys
from urllib.request import urlopen, urlretrieve

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_paths import ASGARD_DOC_DIR
from ASGARD import ISM, Model, Observer, PowerLawJet, Radiation, Setups, TophatJet, TwoComponentJet, Wind, units

DOC_ROOT = "https://yihanwangastro.github.io/VegasAfterglow/docs"
OUTPUT_DIR = ASGARD_DOC_DIR / "vegas_afterglow_compare"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def _label(value: float, unit: str) -> str:
    exp = int(np.floor(np.log10(value)))
    base = value / (10.0**exp)
    return fr"${base:.1f} \times 10^{{{exp}}}$ {unit}"


def _progress(i: int, n: int, name: str) -> None:
    width = 20
    done = int(width * i / n)
    bar = f"[{'#' * done}{'-' * (width - done)}] {100 * i / n:4.0f}%"
    print(f"{bar} {i}/{n} {name}")


def _save(fig, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {path.name}")
    return path


def _download(url: str, local: Path) -> Path:
    if local.exists():
        return local
    with urlopen(url, context=ssl._create_unverified_context()) as fp:
        local.write_bytes(fp.read())
    print(f"[down] {local.name}")
    return local


def _model_base(resolution: tuple[float, float, int] = (0.12, 0.25, 8), *, include_reverse: bool = False, include_ssc: bool = False, with_kn: bool = True) -> Model:
    medium = ISM(n_ism=1.0)
    jet = TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0)
    obs = Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0)
    fwd_rad = Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.3, ssc=include_ssc, kn=with_kn)
    kwargs = dict(
        num_gam_e=81,
        num_nu=49,
        num_r=80,
        num_theta=80,
        num_tobs=48,
        electron_adaptive_substeps=False,
    )
    return Model(
        jet=jet,
        medium=medium,
        observer=obs,
        fwd_rad=fwd_rad,
        rvs_rad=Radiation(eps_e=1.0e-2, eps_B=1.0e-4, p=2.4) if include_reverse else None,
        setups=Setups(**kwargs),
        resolutions=resolution,
    )


def _with_note(name: str, note: str) -> Path:
    fig = plt.figure(figsize=(6.2, 2.6))
    plt.text(0.02, 0.6, f"{name}", fontsize=11)
    plt.text(0.02, 0.2, note, fontsize=10)
    plt.axis("off")
    return _save(fig, OUTPUT_DIR / f"asgard_{name}.png")


def _build_basic_lc_spec() -> Path:
    model = _model_base()
    times = np.logspace(2.0, 8.0, 200)
    bands = np.array([1.0e9, 1.0e14, 1.0e17], dtype=float)
    lc = model.flux_density_grid(times, bands).total
    freqs = np.logspace(5.0, 22.0, 200)
    epochs = np.array([1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8], dtype=float)
    spec = model.flux_density_grid(epochs, freqs).total

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6), dpi=200)
    for i, nu in enumerate(bands):
        axes[0].loglog(times, lc[i, :], label=_label(nu, "Hz"))
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Flux Density (erg/cm²/s/Hz)")
    axes[0].legend()

    colors = plt.cm.viridis(np.linspace(0, 1, len(epochs)))
    for i, t in enumerate(epochs):
        axes[1].loglog(freqs, spec[:, i], color=colors[i], label=_label(t, "s"))
    for i, band in enumerate(bands):
        axes[1].axvline(band, ls="--", color=f"C{i}")
    axes[1].set_xlabel("frequency (Hz)")
    axes[1].set_ylabel("flux density (erg/cm²/s/Hz)")
    axes[1].legend(ncol=2)
    return _save(fig, OUTPUT_DIR / "asgard_basic_lc_spec.png")


def _build_basic_bolometric() -> Path:
    model = _model_base()
    times = np.logspace(2.0, 8.0, 100)
    bat = model.flux(times, 3.6e18, 3.6e19, num_points=20)
    opt = model.flux(times, 4.6e14, 5.6e14, num_points=20)
    fig = plt.figure(figsize=(8, 6), dpi=200)
    plt.loglog(times, bat.total, label="Swift/BAT (15-150 keV)", linewidth=2)
    plt.loglog(times, opt.total, label="V-band optical", linewidth=2)
    plt.xlabel("Time [s]")
    plt.ylabel("Integrated Flux [erg/cm²/s]")
    plt.legend()
    plt.title("Broadband Light Curves")
    return _save(fig, OUTPUT_DIR / "asgard_basic_bolometric.png")


def _build_reverse_shock_lc() -> Path:
    model = _model_base(include_reverse=True)
    times = np.logspace(2.0, 8.0, 200)
    bands = np.array([1.0e9, 1.0e14, 1.0e17], dtype=float)
    res = model.flux_density_grid(times, bands)
    fig = plt.figure(figsize=(4.8, 3.6), dpi=200)
    for i, nu in enumerate(bands):
        plt.loglog(times, res.fwd.sync[i, :], label=f"{_label(nu, 'Hz')} (fwd)")
        plt.loglog(times, res.rvs.sync[i, :], ls="--", label=f"{_label(nu, 'Hz')} (rvs)")
    plt.xlabel("Time (s)")
    plt.ylabel("Flux Density (erg/cm²/s/Hz)")
    plt.legend(ncol=2, fontsize=8)
    plt.title("Forward/reverse synchrotron light curves")
    return _save(fig, OUTPUT_DIR / "asgard_reverse_shock_lc.png")


def _build_ssc_lc() -> Path:
    model = _model_base(include_ssc=True)
    times = np.logspace(2.0, 8.0, 200)
    bands = np.array([1.0e9, 1.0e14, 1.0e17], dtype=float)
    res = model.flux_density_grid(times, bands)
    fig = plt.figure(figsize=(4.8, 3.6), dpi=200)
    for i, nu in enumerate(bands):
        plt.loglog(times, res.fwd.sync[i, :], label=f"{_label(nu, 'Hz')} (sync)")
        plt.loglog(times, res.fwd.ssc[i, :], ls="--", label=f"{_label(nu, 'Hz')} (SSC)")
    plt.xlabel("Time (s)")
    plt.ylabel("Flux Density (erg/cm²/s/Hz)")
    plt.legend(ncol=2, fontsize=8)
    plt.title("Synchrotron vs SSC")
    return _save(fig, OUTPUT_DIR / "asgard_ssc_lc.png")


def _build_shock_quantities() -> Path:
    model = _model_base(resolution=(0.1, 0.25, 10))
    details = model.details(t_min=1.0, t_max=1.0e8)
    attrs = ["Gamma", "B_comv", "N_p", "radius", "nu_m", "nu_c", "nu_a"]
    t = details.fwd.t_obs
    if t.ndim != 1:
        return _with_note("shock_quantities", "details field shape differs from expected 1D observer-time form")
    n = len(attrs)
    fig, axes = plt.subplots(1, n, figsize=(4.0 * n, 3.0), dpi=200)
    for ax, attr in zip(axes, attrs):
        if hasattr(details.fwd, attr):
            ax.loglog(t, getattr(details.fwd, attr))
            ax.set_xlabel("t_obs [s]")
            ax.set_ylabel(attr)
    return _save(fig, OUTPUT_DIR / "asgard_shock_quantities.png")


def _build_photon_quantities() -> Path:
    model = _model_base(resolution=(0.1, 0.25, 10))
    details = model.details(t_min=1.0, t_max=1.0e8)
    t = details.fwd.t_obs
    fig = plt.figure(figsize=(11.0, 3.6), dpi=200)
    for i, attr in enumerate(["nu_a", "nu_m", "nu_c"]):
        ax = plt.subplot(1, 3, i + 1)
        ax.loglog(t, getattr(details.fwd, attr))
        ax.set_xlabel("t_obs [s]")
        ax.set_ylabel(fr"${attr}$")
    return _save(fig, OUTPUT_DIR / "asgard_photon_quantities.png")


def _build_doppler() -> Path:
    model = _model_base(resolution=(0.1, 0.25, 10))
    details = model.details(t_min=1.0, t_max=1.0e8)
    return _with_note("doppler", "Current asgard details API no longer exposes per-patch Doppler grids.")


def _build_eat() -> Path:
    model = _model_base(resolution=(0.1, 0.25, 10))
    details = model.details(t_min=1.0, t_max=1.0e8)
    if details.fwd.t_obs.ndim != 1:
        return _with_note("EAT", "Current asgard details API no longer exposes per-patch EAT grids.")
    ratio = np.ones_like(details.fwd.t_obs)
    fig = plt.figure(figsize=(6.0, 4.0), dpi=200)
    plt.plot(details.fwd.t_obs, ratio, marker=".")
    plt.xscale("log")
    plt.xlabel("t_obs [s]")
    plt.ylabel("ratio")
    plt.title("EAT check placeholder")
    return _save(fig, OUTPUT_DIR / "asgard_EAT.png")


def _build_electron_quantities() -> Path:
    return _with_note("electron_quantities", "Current asgard details API does not expose gamma_a/m/c fields used by this example.")


def _build_sky_single() -> Path:
    model = _model_base()
    img = model.sky_image([1.0e6], 1.0e9, 500.0 * units.uas, npixel=64)
    fig, ax = plt.subplots(dpi=100)
    extent = img.extent / units.uas
    pcm = ax.imshow(img.image[0].T, origin="lower", extent=extent, cmap="inferno")
    ax.set_xlabel(r"$\Delta x$ ($\mu$as)")
    ax.set_ylabel(r"$\Delta y$ ($\mu$as)")
    ax.set_title(r"$t_{obs}=10^6\,$s, $\nu=1\,$GHz")
    fig.colorbar(pcm, ax=ax, label="Surface brightness (erg/cm²/s/Hz/sr)")
    return _save(fig, OUTPUT_DIR / "asgard_sky_image_single.png")


def _build_sky_offaxis() -> Path:
    model_offaxis = _model_base()
    model_offaxis.observer.theta_obs = 0.4
    _ = model_offaxis.details()
    times = np.logspace(5.0, 8.0, 30)
    imgs = model_offaxis.sky_image(times, nu_obs=1.0e9, fov=5000.0 * units.uas, npixel=48)
    idx = [0, len(times) // 2, -1]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.6), dpi=180)
    for ax, i in zip(axes, idx):
        extent = imgs.extent / units.uas
        pcm = ax.imshow(imgs.image[i].T, origin="lower", extent=extent, cmap="inferno")
        ax.set_title(f"t={times[i]:.1e}s")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(pcm, ax=ax, label="Surface brightness (erg/cm²/s/Hz/sr)")
    plt.suptitle("Off-axis sky image sequence")
    return _save(fig, OUTPUT_DIR / "asgard_sky_image_offaxis.png")


def _build_sky_flux_comparison() -> Path:
    model = _model_base()
    t_obs = np.logspace(3.0, 8.0, 30)
    nu_obs = 1.0e9
    img = model.sky_image(t_obs, nu_obs=nu_obs, fov=2000.0 * units.uas, npixel=64)
    flux_from_image = img.image.sum(axis=(1, 2)) * img.pixel_solid_angle
    flux_direct = model.flux_density_grid(t_obs, np.array([nu_obs])).total[0, :]
    fig, axes = plt.subplots(2, 1, figsize=(6.0, 6.0), dpi=180, sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    axes[0].loglog(t_obs, flux_direct, "k-", label="flux_density_grid")
    axes[0].plot(t_obs, flux_from_image, "o", ms=3.5, color="C1", label="sky_image (integrated)")
    axes[0].set_ylabel("Flux density (erg/cm²/s/Hz)")
    axes[0].legend()
    ratio = flux_from_image / flux_direct
    axes[1].semilogx(t_obs, ratio, "o-", ms=3.5, color="C1")
    axes[1].axhline(1.0, color="k", ls="--", lw=0.8)
    axes[1].set_xlabel("Observer time (s)")
    axes[1].set_ylabel("image / direct")
    return _save(fig, OUTPUT_DIR / "asgard_sky_image_flux_comparison.png")


def _build_introspection_jet() -> Path:
    model = Model(
        jet=PowerLawJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0, k_e=2.0, k_g=1.5),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.3),
        setups=Setups(num_tobs=48, num_nu=49, num_gam_e=81, num_r=80, num_theta=80),
        resolutions=(0.12, 0.25, 8),
    )
    theta = np.linspace(0.0, 0.5, 100)
    eiso = np.array([model.jet_E_iso(0.0, t) for t in theta], dtype=float)
    g0 = np.array([model.jet_Gamma0(0.0, t) for t in theta], dtype=float)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.0, 4.0), dpi=180)
    ax1.semilogy(np.degrees(theta), eiso)
    ax1.set_xlabel("Polar Angle [degrees]")
    ax1.set_ylabel(r"$E_{\rm iso}$ [erg]")
    ax1.set_title("Jet Energy Profile")
    ax1.grid(True, alpha=0.3)
    ax2.semilogy(np.degrees(theta), g0)
    ax2.set_xlabel("Polar Angle [degrees]")
    ax2.set_ylabel(r"$\Gamma_0$")
    ax2.set_title("Jet Lorentz Factor Profile")
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    return _save(fig, OUTPUT_DIR / "asgard_introspection_jet.png")


def _build_introspection_medium() -> Path:
    medium = Wind(A_star=1.0, n_ism=0.1, n0=1.0e3, k=2.0)
    model = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0),
        medium=medium,
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.3),
    )
    phi = 0.0
    theta = 0.1
    r = np.logspace(15.0, 20.0, 100)
    rho = np.array([model.medium_rho(phi, theta, rr) for rr in r], dtype=float)
    n = rho / 1.67e-24
    fig = plt.figure(figsize=(8.0, 6.0), dpi=180)
    plt.loglog(r, n)
    plt.xlabel(r"Radius [cm]")
    plt.ylabel(r"Number Density [cm$^{-3}$]")
    plt.title("Medium Density Profile")
    plt.grid(True, alpha=0.3)
    plt.axhline(1.0e3, color="red", linestyle="--", alpha=0.7, label="Inner constant density")
    plt.axhline(0.1, color="blue", linestyle="--", alpha=0.7, label="Outer ISM density")
    plt.legend()
    return _save(fig, OUTPUT_DIR / "asgard_introspection_medium.png")


def _build_introspection_twocomp() -> Path:
    model = Model(
        jet=TwoComponentJet(
            theta_c=0.05,
            E_iso_c=1.0e53,
            Gamma0_c=300.0,
            theta_w=0.15,
            E_iso_w=1.0e52,
            Gamma0_w=100.0,
        ),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=1.0e-1, eps_B=1.0e-3, p=2.3),
    )
    theta = np.linspace(0.0, 0.3, 200)
    eiso = np.array([model.jet_E_iso(0.0, t) for t in theta], dtype=float)
    g0 = np.array([model.jet_Gamma0(0.0, t) for t in theta], dtype=float)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.0, 8.0), dpi=180)
    ax1.semilogy(np.degrees(theta), eiso)
    ax1.axvline(np.degrees(0.05), color="red", ls="--", alpha=0.7, label="Core boundary")
    ax1.axvline(np.degrees(0.15), color="blue", ls="--", alpha=0.7, label="Wide component boundary")
    ax1.set_ylabel(r"$E_{\rm iso}$ [erg]")
    ax1.set_title("Two-Component Jet: Energy Profile")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax2.semilogy(np.degrees(theta), g0)
    ax2.axvline(np.degrees(0.05), color="red", ls="--", alpha=0.7, label="Core boundary")
    ax2.axvline(np.degrees(0.15), color="blue", ls="--", alpha=0.7, label="Wide component boundary")
    ax2.set_xlabel("Polar Angle [degrees]")
    ax2.set_ylabel(r"$\Gamma_0$")
    ax2.set_title("Two-Component Jet: Lorentz Factor Profile")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    return _save(fig, OUTPUT_DIR / "asgard_introspection_twocomp.png")


def _compare(name: str, ref_name: str, local: Path) -> Path:
    if not local.exists():
        raise FileNotFoundError(f"local image not found: {local}")
    ref_url = f"{DOC_ROOT}/_static/images/{ref_name}"
    ref = _download(ref_url, OUTPUT_DIR / f"ref_{ref_name.replace('/', '_')}")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    ref_img = plt.imread(ref)
    local_img = plt.imread(local)
    axes[0].imshow(ref_img)
    axes[0].set_axis_off()
    axes[0].set_title(f"VegasAfterglow {name}")
    axes[1].imshow(local_img)
    axes[1].set_axis_off()
    axes[1].set_title("ASGARD")
    return _save(fig, OUTPUT_DIR / f"compare_{name}.png")


def main() -> None:
    builders = [
        ("basic_lc_spec", "basic_lc_spec.png", _build_basic_lc_spec),
        ("basic_bolometric", "basic_bolometric.png", _build_basic_bolometric),
        ("reverse_shock_lc", "reverse_shock_lc.png", _build_reverse_shock_lc),
        ("ssc_lc", "ssc_lc.png", _build_ssc_lc),
        ("shock_quantities", "shock_quantities.png", _build_shock_quantities),
        ("electron_quantities", "electron_quantities.png", _build_electron_quantities),
        ("photon_quantities", "photon_quantities.png", _build_photon_quantities),
        ("doppler", "doppler.png", _build_doppler),
        ("EAT", "EAT.png", _build_eat),
        ("sky_image_single", "sky_image_single.png", _build_sky_single),
        ("sky_image_offaxis", "sky_image_offaxis.png", _build_sky_offaxis),
        ("sky_image_flux_comparison", "sky_image_flux_comparison.png", _build_sky_flux_comparison),
        ("introspection_jet", "introspection_jet.png", _build_introspection_jet),
        ("introspection_medium", "introspection_medium.png", _build_introspection_medium),
        ("introspection_twocomp", "introspection_twocomp.png", _build_introspection_twocomp),
    ]

    outputs = {}
    for i, (name, ref_name, builder) in enumerate(builders, start=1):
        _progress(i, len(builders), name)
        try:
            local = builder()
        except Exception as exc:
            local = _with_note(name, f"{type(exc).__name__}: {exc}")
        outputs[name] = (ref_name, local)

    compare_dir = OUTPUT_DIR
    for name, (ref_name, local) in outputs.items():
        _compare(name, ref_name, local)

    print(f"output: {OUTPUT_DIR}")
    print("generated:", len(outputs))
    for name, (ref_name, local) in outputs.items():
        print(f"{name}: {local}")


if __name__ == "__main__":
    main()
