from __future__ import annotations

from pathlib import Path
from functools import lru_cache
import time
from typing import Callable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
import sys

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_paths import ASGARD_DOC_DIR
from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet, Wind, units
from ASGARD import PowerLawJet as ASGARD_PowerLawJet
from ASGARD import TwoComponentJet
from src import constants
from VegasAfterglow import ISM as VegasISM
from VegasAfterglow import Model as VegasModel
from VegasAfterglow import Observer as VegasObserver
from VegasAfterglow import PowerLawJet as VegasPowerLawJet
from VegasAfterglow import Radiation as VegasRadiation
from VegasAfterglow import TophatJet as VegasTophatJet
from VegasAfterglow import Wind as VegasWind
from VegasAfterglow import TwoComponentJet as VegasTwoComponentJet
from VegasAfterglow import units as vegas_units


OUTPUT_DIR = ASGARD_DOC_DIR / "vegas_afterglow_compare"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MEDIUM_ISM_N = 1.0
MEDIUM_WIND_ASTAR = 1.0
MEDIUM_WIND_NISM = 0.1
MEDIUM_WIND_N0 = 1.0e3
MEDIUM_WIND_K = 2.0
PROTON_MASS_G = constants.para_m_p
ASGARD_CHARINT_NUM_GAM_E = 41
BASE_NUM_R = 80
REVERSE_NUM_R = 120
REVERSE_DURATION_S = 10.0
STRONG_IC_EPS_E = 3.0e-1
STRONG_IC_EPS_B = 1.0e-5
STRONG_IC_XI = 1.0e-1
FWD_P = 2.3
RVS_P = 2.4
BASIC_TIMES = np.logspace(0.0, 8.0, 160)
BASIC_BANDS = np.array([1.0e9, 1.0e14, 1.0e17], dtype=float)
BASIC_FREQS = np.logspace(5.0, 29.0, 160)
BASIC_EPOCHS = np.array([1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8], dtype=float)
OFFAXIS_SKY_TIMES = np.logspace(5.0, 8.0, 3)
SINGLE_SKY_TIME = np.array([1.0e6], dtype=float)
SINGLE_SKY_NPIXEL = 48
FLUX_SKY_TIMES = np.logspace(0.0, 8.0, 20)
FLUX_SKY_NPIXEL = 48
SPEED_TIMES = np.logspace(2.0, 8.0, 16)
SPEED_BANDS = np.array([1.0e9, 1.0e14, 1.0e17], dtype=float)
SPEED_FREQS = np.logspace(9.0, 18.0, 16)
SPEED_EXPO = np.logspace(2.0, 6.0, 8)
SPEED_SPEC_TIMES = np.array([1.0e3, 1.0e5, 1.0e7], dtype=float)
SPEED_SPEC_FREQS = np.logspace(9.0, 22.0, 24)
SPEED_SKY_NPIXEL = 20
SPECTRUM_COMPARE_FREQS = np.logspace(8.0, 29.0, 240)
SPECTRUM_COMPARE_NUM_NU = 81


def _save(fig, path: Path) -> Path:
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {path.name}")
    return path


def _label(value: float, unit: str) -> str:
    exp = int(np.floor(np.log10(value)))
    base = value / (10.0**exp)
    return fr"${base:.1f}\times 10^{{{exp}}}$ {unit}"


def _run_timed(name: str, fn: Callable[[], object]) -> tuple[float, object]:
    t0 = time.perf_counter()
    out = fn()
    return time.perf_counter() - t0, out


def _collapse_curve(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    while arr.ndim > 1:
        arr = np.nanmean(arr, axis=0)
    return arr


def _reference_series(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    while arr.ndim > 1:
        arr = arr[(0,) * (arr.ndim - 1)]
    return arr


def _calibrate_sky_image_to_asgard_basis(image: np.ndarray, backend: str) -> np.ndarray:
    if backend.lower() == "vegasafterglow":
        return np.rot90(np.asarray(image, dtype=float), 3)
    return np.asarray(image, dtype=float)


def _shared_positive_log_norm(*images: np.ndarray) -> LogNorm:
    positive = [np.asarray(img, dtype=float)[np.asarray(img, dtype=float) > 0.0] for img in images]
    positive = [arr for arr in positive if arr.size > 0]
    if not positive:
        return LogNorm()
    stacked = np.concatenate(positive)
    return LogNorm(vmin=float(np.min(stacked)), vmax=float(np.max(stacked)))


def _safe_log_interp(x_ref: np.ndarray, x_src: np.ndarray, y_src: np.ndarray) -> np.ndarray:
    x_ref = np.asarray(x_ref, dtype=float)
    x_src = np.asarray(x_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    mask = np.isfinite(x_src) & np.isfinite(y_src) & (x_src > 0.0) & (y_src > 0.0)
    x_src = x_src[mask]
    y_src = y_src[mask]
    if x_src.size < 2:
        return np.full_like(x_ref, np.nan, dtype=float)
    order = np.argsort(x_src)
    x_src = x_src[order]
    y_src = y_src[order]
    lx_src = np.log10(x_src)
    ly_src = np.log10(y_src)
    lx_ref = np.log10(x_ref)
    ly_ref = np.interp(lx_ref, lx_src, ly_src, left=np.nan, right=np.nan)
    return np.where(np.isfinite(ly_ref), 10.0**ly_ref, np.nan)


def _to_lab_frequency_frame(freq_hz: np.ndarray, doppler: np.ndarray, z: float) -> np.ndarray:
    freq_hz = np.asarray(freq_hz, dtype=float)
    doppler = np.asarray(doppler, dtype=float)
    return freq_hz * doppler / (1.0 + float(z))


@lru_cache(maxsize=None)
def _build_asgard_model(
    include_reverse: bool = False,
    include_ssc: bool = False,
    theta_obs: float = 0.0,
    num_nu: int = 49,
) -> Model:
    medium = ISM(n_ism=1.0)
    jet = TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0, duration=REVERSE_DURATION_S)
    observer = Observer(lumi_dist=1.0e26, z=0.1, theta_obs=theta_obs)
    fwd_rad = Radiation(
        eps_e=STRONG_IC_EPS_E,
        eps_B=STRONG_IC_EPS_B,
        p=FWD_P,
        xi_N=STRONG_IC_XI,
        ssc=include_ssc,
        kn=True,
    )
    kwargs = dict(
        num_gam_e=ASGARD_CHARINT_NUM_GAM_E,
        num_nu=num_nu,
        num_r=REVERSE_NUM_R if include_reverse else BASE_NUM_R,
        num_theta=80,
        num_tobs=48,
        electron_adaptive_substeps=False,
        electron_solver="charint",
    )
    return Model(
        jet=jet,
        medium=medium,
        observer=observer,
        fwd_rad=fwd_rad,
        rvs_rad=Radiation(
            eps_e=STRONG_IC_EPS_E,
            eps_B=STRONG_IC_EPS_B,
            p=RVS_P,
            xi_N=STRONG_IC_XI,
            ssc=include_ssc,
            kn=True,
        ) if include_reverse else None,
        setups=Setups(reverse_delta_t_s=REVERSE_DURATION_S, **kwargs),
    )


@lru_cache(maxsize=None)
def _build_vegas_model(
    include_reverse: bool = False,
    include_ssc: bool = False,
    theta_obs: float = 0.0,
) -> VegasModel:
    medium = VegasISM(n_ism=1.0)
    jet = VegasTophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0, duration=REVERSE_DURATION_S)
    observer = VegasObserver(lumi_dist=1.0e26, z=0.1, theta_obs=theta_obs)
    fwd_rad = VegasRadiation(
        eps_e=STRONG_IC_EPS_E,
        eps_B=STRONG_IC_EPS_B,
        p=FWD_P,
        xi_e=STRONG_IC_XI,
        ssc=include_ssc,
        kn=True,
    )
    return VegasModel(
        jet=jet,
        medium=medium,
        observer=observer,
        fwd_rad=fwd_rad,
        rvs_rad=VegasRadiation(
            eps_e=STRONG_IC_EPS_E,
            eps_B=STRONG_IC_EPS_B,
            p=RVS_P,
            xi_e=STRONG_IC_XI,
            ssc=include_ssc,
            kn=True,
        ) if include_reverse else None,
        resolutions=(0.12, 0.25, 8),
    )


@lru_cache(maxsize=None)
def _cached_details_pair(
    include_reverse: bool = False,
    include_ssc: bool = False,
    theta_obs: float = 0.0,
    t_min: float = 1.0,
    t_max: float = 1.0e8,
) -> tuple[object, object]:
    model_asgard = _build_asgard_model(include_reverse=include_reverse, include_ssc=include_ssc, theta_obs=theta_obs)
    model_vegas = _build_vegas_model(include_reverse=include_reverse, include_ssc=include_ssc, theta_obs=theta_obs)
    return model_asgard.details(t_min=t_min, t_max=t_max), model_vegas.details(t_min=t_min, t_max=t_max)


def _build_basic_lc_spec() -> Path:
    model_asgard = _build_asgard_model()
    model_vegas = _build_vegas_model()
    times = BASIC_TIMES
    bands = BASIC_BANDS
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.0), dpi=200)

    lc_as = model_asgard.flux_density_grid(times, bands).total
    lc_vg = model_vegas.flux_density_grid(times, bands).total
    colors = plt.cm.tab10(np.linspace(0, 1, len(bands)))
    for i, nu in enumerate(bands):
        label = _label(nu, "Hz")
        axes[0].loglog(times, lc_as[i, :], color=colors[i], ls="-", lw=1.8, label=f"ASGARD {label}")
        axes[0].loglog(times, lc_vg[i, :], color=colors[i], ls="--", lw=1.4, label=f"Vegas {label}")

    freqs = BASIC_FREQS
    epochs = BASIC_EPOCHS
    sed_as = model_asgard.flux_density_grid(epochs, freqs).total * freqs[:, None]
    sed_vg = model_vegas.flux_density_grid(epochs, freqs).total * freqs[:, None]
    spec_peak = float(max(np.nanmax(sed_as), np.nanmax(sed_vg)))
    colors2 = plt.cm.viridis(np.linspace(0, 1, len(epochs)))
    for i, t in enumerate(epochs):
        label = _label(t, "s")
        axes[1].loglog(freqs, sed_as[:, i], color=colors2[i], ls="-", lw=1.6, label=f"ASGARD t={label}")
        axes[1].loglog(freqs, sed_vg[:, i], color=colors2[i], ls="--", lw=1.2, label=f"Vegas t={label}")
    for idx, band in enumerate(bands):
        axes[1].axvline(band, ls="--", alpha=0.4, color=f"C{idx}")

    axes[0].set_xlabel("Time [s]")
    axes[0].set_ylabel(r"Flux Density (erg/cm$^2$/s/Hz)")
    axes[0].set_title("Light Curves")
    axes[1].set_xlabel("Frequency [Hz]")
    axes[1].set_ylabel(r"$\nu F_\nu$ (erg/cm$^2$/s)")
    axes[1].set_title("SED")
    if spec_peak > 0.0:
        axes[1].set_ylim(bottom=spec_peak * 1.0e-10)
    axes[0].legend(fontsize=7, ncol=1)
    axes[1].legend(fontsize=6, ncol=2)
    return _save(fig, OUTPUT_DIR / "compare_basic_lc_spec.png")


def _build_basic_bolometric() -> Path:
    model_asgard = _build_asgard_model()
    model_vegas = _build_vegas_model()
    times = np.logspace(0.0, 8.0, 80)
    bat_as = model_asgard.flux(times, 3.6e18, 3.6e19, 20).total
    opt_as = model_asgard.flux(times, 4.6e14, 5.6e14, 20).total
    bat_vg = model_vegas.flux(times, 3.6e18, 3.6e19, 20).total
    opt_vg = model_vegas.flux(times, 4.6e14, 5.6e14, 20).total

    fig = plt.figure(figsize=(8.4, 4.8), dpi=200)
    ax = plt.gca()
    ax.loglog(times, bat_as, label="ASGARD Swift/BAT", lw=2, color="C0")
    ax.loglog(times, bat_vg, label="VegasAfterglow Swift/BAT", lw=2, ls="--", color="C0")
    ax.loglog(times, opt_as, label="ASGARD V-band", lw=2, color="C1")
    ax.loglog(times, opt_vg, label="VegasAfterglow V-band", lw=2, ls="--", color="C1")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"Integrated Flux (erg/cm$^2$/s)")
    ax.set_title("Broadband Light Curves")
    ax.legend(fontsize=8)
    return _save(fig, OUTPUT_DIR / "compare_basic_bolometric.png")


def _build_reverse_shock_lc() -> Path:
    model_asgard = _build_asgard_model(include_reverse=True, include_ssc=True)
    model_vegas = _build_vegas_model(include_reverse=True, include_ssc=True)
    times = BASIC_TIMES
    bands = BASIC_BANDS

    as_res = model_asgard.flux_density_grid(times, bands)
    vg_res = model_vegas.flux_density_grid(times, bands)

    fig, ax = plt.subplots(1, 1, figsize=(7.4, 5.2), dpi=200)
    all_series: list[np.ndarray] = []
    for i, nu in enumerate(bands):
        label = _label(nu, "Hz")
        fwd_a = np.asarray(as_res.fwd.sync[i, :], dtype=float)
        fwd_v = np.asarray(vg_res.fwd.sync[i, :], dtype=float)
        rvs_a = np.asarray(as_res.rvs.sync[i, :], dtype=float)
        rvs_v = np.asarray(vg_res.rvs.sync[i, :], dtype=float)
        ax.loglog(times, fwd_a, color=f"C{i}", lw=1.8, label=f"ASGARD fwd {label}")
        ax.loglog(times, fwd_v, color=f"C{i}", ls="--", lw=1.3, label=f"Vegas fwd {label}")
        ax.loglog(times, rvs_a, color=f"C{i}", ls="-.", lw=1.8, alpha=0.8, label=f"ASGARD rvs {label}")
        ax.loglog(times, rvs_v, color=f"C{i}", ls=":", lw=1.6, alpha=0.8, label=f"Vegas rvs {label}")
        all_series.extend([fwd_a, fwd_v, rvs_a, rvs_v])

    peak_flux = 0.0
    for s in all_series:
        arr = np.asarray(s, dtype=float)
        pos = arr[np.isfinite(arr) & (arr > 0.0)]
        if pos.size > 0:
            peak_flux = max(peak_flux, float(np.nanmax(pos)))
    if peak_flux > 0.0:
        ax.set_ylim(peak_flux * 1.0e-10, peak_flux * 1.25)

    ax.set_xlabel("Observer time [s]")
    ax.set_ylabel(r"Flux Density (erg/cm$^2$/s/Hz)")
    ax.set_title("Forward + reverse shock")
    ax.legend(fontsize=6.2, ncol=2)
    return _save(fig, OUTPUT_DIR / "compare_reverse_shock_lc.png")


def _build_ssc_lc() -> Path:
    model_asgard = _build_asgard_model(include_ssc=True)
    model_vegas = _build_vegas_model(include_ssc=True)
    times = BASIC_TIMES
    bands = BASIC_BANDS

    as_res = model_asgard.flux_density_grid(times, bands)
    vg_res = model_vegas.flux_density_grid(times, bands)

    fig = plt.figure(figsize=(6.0, 4.2), dpi=200)
    ax = plt.gca()
    for i, nu in enumerate(bands):
        label = _label(nu, "Hz")
        ax.loglog(times, as_res.fwd.sync[i, :], color=f"C{i}", lw=1.8, label=f"ASGARD sync {label}")
        ax.loglog(times, as_res.fwd.ssc[i, :], ls="--", color=f"C{i}", lw=1.4, label=f"ASGARD SSC {label}")
        ax.loglog(times, vg_res.fwd.sync[i, :], color=f"C{i}", ls=":", lw=1.2, label=f"Vegas sync {label}")
        ax.loglog(times, vg_res.fwd.ssc[i, :], color=f"C{i}", ls="-.", lw=1.0, label=f"Vegas SSC {label}")

    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"Flux Density (erg/cm$^2$/s/Hz)")
    ax.set_title("Synchrotron and SSC")
    ax.legend(fontsize=6.8, ncol=2)
    return _save(fig, OUTPUT_DIR / "compare_ssc_lc.png")


def _plot_two_line_sets(ax, x: np.ndarray, lines: list[tuple[np.ndarray, str, str]], x_label: str) -> None:
    x = np.asarray(x, dtype=float)
    for y, color, label in lines:
        y = np.asarray(y, dtype=float)
        mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
        if not mask.any():
            continue
        ax.loglog(x[mask], y[mask], color=color, label=label)
    ax.set_xlabel(x_label)


def _sample_theta_curve(model, attribute: str, theta: np.ndarray) -> np.ndarray:
    theta = np.asarray(theta, dtype=float)
    fn = getattr(model, attribute)
    try:
        values = np.asarray(fn(0.0, theta), dtype=float)
        if values.shape == theta.shape:
            return values
    except Exception:
        pass
    values = [fn(0.0, float(t)) for t in theta]
    return np.asarray(values, dtype=float)


def _with_note(name: str, note: str) -> Path:
    fig = plt.figure(figsize=(6.2, 2.6))
    plt.text(0.05, 0.55, name, fontsize=11)
    plt.text(0.05, 0.2, note, fontsize=10)
    plt.axis("off")
    return _save(fig, OUTPUT_DIR / f"compare_{name}.png")


def _build_shock_quantities() -> Path:
    da, dv = _cached_details_pair(False, False, 0.0, 1.0, 1.0e8)
    vegas_z = 0.1

    attrs = ["Gamma", "B_comv", "N_p", "nu_m", "nu_c", "nu_a"]
    t_ref = _reference_series(da.fwd.t_obs)
    if t_ref.ndim != 1:
        return _with_note("shock_quantities", "ASGARD shock details are not 1D on this baseline.")

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.0), dpi=200)
    axes = axes.ravel()
    for i, attr in enumerate(attrs):
        ax = axes[i]
        if hasattr(da.fwd, attr) and hasattr(dv.fwd, attr):
            ya = _reference_series(getattr(da.fwd, attr))
            yv = _reference_series(getattr(dv.fwd, attr))
            if attr == "N_p":
                yv = yv * (4.0 * np.pi)
            if attr in {"nu_m", "nu_c"}:
                doppler_v = _reference_series(dv.fwd.Doppler)
                yv = _to_lab_frequency_frame(yv, doppler_v, vegas_z)
            tv = _reference_series(dv.fwd.t_obs)
            yv = _safe_log_interp(t_ref, tv, yv)
            _plot_two_line_sets(ax, t_ref, [(ya, "C0", f"ASGARD {attr}"), (yv, "C1", f"Vegas {attr}")], "t_obs [s]")
            ax.set_ylabel(attr)
            ax.set_xlabel("t_obs [s]")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.xaxis.set_major_formatter(LogFormatter())
        else:
            ax.text(0.1, 0.55, f"{attr}: not available in both", fontsize=9)
        ax.set_title(attr)
    plt.tight_layout()
    return _save(fig, OUTPUT_DIR / "compare_shock_quantities.png")


def _build_photon_quantities() -> Path:
    da, dv = _cached_details_pair(False, False, 0.0, 1.0, 1.0e8)
    vegas_z = 0.1
    attrs = ["nu_a", "nu_m", "nu_c"]
    t_ref = _reference_series(da.fwd.t_obs)
    if t_ref.ndim != 1:
        return _with_note("photon_quantities", "ASGARD photon quantities are not 1D on this baseline.")

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.6), dpi=200)
    for i, attr in enumerate(attrs):
        ax = axes[i]
        ya = _reference_series(getattr(da.fwd, attr))
        yv_raw = _reference_series(getattr(dv.fwd, attr))
        if attr in {"nu_m", "nu_c"}:
            doppler_v = _reference_series(dv.fwd.Doppler)
            yv_raw = _to_lab_frequency_frame(yv_raw, doppler_v, vegas_z)
        tv = _reference_series(dv.fwd.t_obs)
        yv = _safe_log_interp(t_ref, tv, yv_raw)
        _plot_two_line_sets(ax, t_ref, [(ya, "C0", f"ASGARD {attr}"), (yv, "C1", f"Vegas {attr}")], "t_obs [s]")
        ax.set_xlabel(r"t_obs [s]")
        ax.set_ylabel(attr)
        ax.set_title(attr)
        ax.set_xscale("log")
        ax.set_yscale("log")
    axes[0].legend(fontsize=7)
    plt.tight_layout()
    return _save(fig, OUTPUT_DIR / "compare_photon_quantities.png")


def _build_sky_single() -> Path:
    model_a = _build_asgard_model()
    model_v = _build_vegas_model()
    sky_a = model_a.sky_image(SINGLE_SKY_TIME, nu_obs=1.0e9, fov=500.0 * units.uas, npixel=SINGLE_SKY_NPIXEL)
    sky_v = model_v.sky_image(SINGLE_SKY_TIME, nu_obs=1.0e9, fov=500.0 * vegas_units.uas, npixel=SINGLE_SKY_NPIXEL)
    img_v = _calibrate_sky_image_to_asgard_basis(sky_v.image[0], "VegasAfterglow")
    img_a = _calibrate_sky_image_to_asgard_basis(sky_a.image[0], "ASGARD")
    norm = _shared_positive_log_norm(img_v, img_a)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), dpi=180)
    for ax, sky, title, cmap, unit in [
        (axes[0], sky_v, "VegasAfterglow", "magma", vegas_units.uas),
        (axes[1], sky_a, "ASGARD", "inferno", units.uas),
    ]:
        img = img_v if title == "VegasAfterglow" else img_a
        extent = sky.extent / unit
        im = ax.imshow(np.where(img > 0, img, np.nan).T, origin="lower", extent=extent, cmap="magma", norm=norm)
        ax.set_title(f"{title}: t=10^6 s, nu=1 GHz")
        ax.set_xlabel(r"$\Delta x$ ($\mu$as)")
        ax.set_ylabel(r"$\Delta y$ ($\mu$as)")
        fig.colorbar(im, ax=ax, shrink=0.9, label="Surface brightness")
    return _save(fig, OUTPUT_DIR / "compare_sky_image_single.png")


def _build_sky_offaxis() -> Path:
    model_a = _build_asgard_model(theta_obs=0.4)
    model_v = _build_vegas_model(theta_obs=0.4)
    times = OFFAXIS_SKY_TIMES
    img_a = model_a.sky_image(times, nu_obs=1.0e9, fov=5000.0 * units.uas, npixel=40)
    img_v = model_v.sky_image(times, nu_obs=1.0e9, fov=5000.0 * vegas_units.uas, npixel=40)
    idx = -1
    img_v_cal = _calibrate_sky_image_to_asgard_basis(img_v.image[idx], "VegasAfterglow")
    img_a_cal = _calibrate_sky_image_to_asgard_basis(img_a.image[idx], "ASGARD")
    norm = _shared_positive_log_norm(img_v_cal, img_a_cal)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), dpi=180)
    extent_v = img_v.extent / vegas_units.uas
    extent_a = img_a.extent / units.uas
    title = f"{times[idx]:.1e} s"
    for i, (ax, image, extent, name, cmap) in enumerate(
        [
            (axes[0], img_v_cal.T, extent_v, "VegasAfterglow", "magma"),
            (axes[1], img_a_cal.T, extent_a, "ASGARD", "inferno"),
        ]
    ):
        im = ax.imshow(np.where(image > 0, image, np.nan), origin="lower", extent=extent, cmap="magma", norm=norm)
        ax.set_title(f"{name}, {title}")
        ax.set_xlabel(r"$\Delta x$ ($\mu$as)")
        if i == 0:
            ax.set_ylabel(r"$\Delta y$ ($\mu$as)")
            cax = axes[0]
            fig.colorbar(im, ax=cax, shrink=0.9, label="Surface brightness")
        else:
            ax.set_xticks([])
            ax.set_yticks([])
    return _save(fig, OUTPUT_DIR / "compare_sky_image_offaxis.png")


def _build_sky_flux_comparison() -> Path:
    model_a = _build_asgard_model()
    model_v = _build_vegas_model()
    times = FLUX_SKY_TIMES
    nu_obs = 1.0e9

    img_a = model_a.sky_image(times, nu_obs=nu_obs, fov=2000.0 * units.uas, npixel=FLUX_SKY_NPIXEL)
    img_v = model_v.sky_image(times, nu_obs=nu_obs, fov=2000.0 * vegas_units.uas, npixel=FLUX_SKY_NPIXEL)
    flux_from_image_a = img_a.image.sum(axis=(1, 2)) * img_a.pixel_solid_angle
    flux_from_image_v = img_v.image.sum(axis=(1, 2)) * img_v.pixel_solid_angle
    direct_a = model_a.flux_density_grid(times, np.array([nu_obs])).total[0, :]
    direct_v = model_v.flux_density_grid(times, np.array([nu_obs])).total[0, :]

    fig, axes = plt.subplots(2, 1, figsize=(6.6, 6.0), dpi=180, sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    axes[0].loglog(times, direct_a, "k-", lw=2.0, label="ASGARD direct")
    axes[0].loglog(times, flux_from_image_a, "o", ms=3.5, color="C1", label="ASGARD image")
    axes[0].loglog(times, direct_v, "k--", lw=1.8, label="Vegas direct")
    axes[0].loglog(times, flux_from_image_v, "x", ms=3.2, color="C3", label="Vegas image")
    axes[0].set_ylabel(r"Flux density (erg/cm$^2$/s/Hz)")
    axes[0].legend(fontsize=8)

    ratio_a = flux_from_image_a / np.maximum(direct_a, np.finfo(float).tiny)
    ratio_v = flux_from_image_v / np.maximum(direct_v, np.finfo(float).tiny)
    axes[1].semilogx(times, ratio_a, "o-", ms=3.5, color="C1", label="ASGARD")
    axes[1].semilogx(times, ratio_v, "x-", ms=3.5, color="C3", label="Vegas")
    axes[1].axhline(1.0, color="k", ls="--", lw=0.8)
    axes[1].set_xlabel("Observer time [s]")
    axes[1].set_ylabel("ratio")
    axes[1].legend(fontsize=8)
    axes[1].set_ylim(0.95, 1.05)

    return _save(fig, OUTPUT_DIR / "compare_sky_image_flux_comparison.png")


def _build_introspection_jet() -> Path:
    model_a = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=STRONG_IC_EPS_E, eps_B=STRONG_IC_EPS_B, p=FWD_P, xi_N=STRONG_IC_XI),
    )
    model_v = VegasModel(
        jet=VegasTophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0),
        medium=VegasISM(n_ism=1.0),
        observer=VegasObserver(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=VegasRadiation(eps_e=STRONG_IC_EPS_E, eps_B=STRONG_IC_EPS_B, p=FWD_P, xi_e=STRONG_IC_XI),
    )

    theta = np.linspace(0.0, 0.5, 100)
    e_iso_a = _sample_theta_curve(model_a, "jet_E_iso", theta)
    g0_a = _sample_theta_curve(model_a, "jet_Gamma0", theta)
    e_iso_v = _sample_theta_curve(model_v, "jet_E_iso", theta)
    g0_v = _sample_theta_curve(model_v, "jet_Gamma0", theta)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.0, 3.8), dpi=180)
    ax1.semilogy(np.degrees(theta), e_iso_a, label="ASGARD")
    ax1.semilogy(np.degrees(theta), e_iso_v, ls="--", label="Vegas")
    ax1.set_xlabel("Polar Angle [deg]")
    ax1.set_ylabel(r"$E_{\mathrm{iso}}$ [erg]")
    ax1.set_title("Jet Energy Profile")
    ax1.legend()
    ax1.grid(alpha=0.3)

    ax2.semilogy(np.degrees(theta), g0_a, label="ASGARD")
    ax2.semilogy(np.degrees(theta), g0_v, ls="--", label="Vegas")
    ax2.set_xlabel("Polar Angle [deg]")
    ax2.set_ylabel(r"$\Gamma_0$")
    ax2.set_title("Jet Initial Lorentz Factor")
    ax2.legend()
    ax2.grid(alpha=0.3)
    return _save(fig, OUTPUT_DIR / "compare_introspection_jet.png")


def _build_introspection_medium() -> Path:
    model_a = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0),
        medium=Wind(MEDIUM_WIND_ASTAR, MEDIUM_WIND_NISM, MEDIUM_WIND_N0, MEDIUM_WIND_K),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=STRONG_IC_EPS_E, eps_B=STRONG_IC_EPS_B, p=FWD_P, xi_N=STRONG_IC_XI),
    )
    model_v = VegasModel(
        jet=VegasTophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0),
        medium=VegasWind(MEDIUM_WIND_ASTAR, MEDIUM_WIND_NISM, MEDIUM_WIND_N0, int(MEDIUM_WIND_K)),
        observer=VegasObserver(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=VegasRadiation(eps_e=STRONG_IC_EPS_E, eps_B=STRONG_IC_EPS_B, p=FWD_P, xi_e=STRONG_IC_XI),
    )

    r = np.logspace(13.0, 20.0, 200)
    rho_a = np.array([model_a.medium.density(0.0, 0.1, float(rr)) for rr in r], dtype=float)
    rho_v = np.asarray(model_v.medium(0.0, 0.1, r), dtype=float)
    if rho_v.shape != r.shape:
        rho_v = rho_v.reshape(r.shape)
    n_a = rho_a
    n_v = rho_v / PROTON_MASS_G

    fig = plt.figure(figsize=(6.8, 4.0), dpi=180)
    ax = plt.gca()
    ax.loglog(r, n_a, label="ASGARD")
    ax.loglog(r, n_v, ls="--", label="Vegas")
    ax.set_xlabel("Radius [cm]")
    ax.set_ylabel(r"Number Density [cm$^{-3}$]")
    ax.set_title("Medium Density")
    ax.grid(alpha=0.3)
    ax.axhline(1.0e3, color="red", ls="--", alpha=0.7, label="inner $n_0$")
    ax.axhline(0.1, color="blue", ls="--", alpha=0.7, label="outer ISM")
    ax.legend(fontsize=8)
    return _save(fig, OUTPUT_DIR / "compare_introspection_medium.png")


def _build_introspection_twocomp() -> Path:
    model_a = Model(
        jet=TwoComponentJet(0.05, 1.0e53, 300.0, 0.15, 1.0e52, 100.0),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=STRONG_IC_EPS_E, eps_B=STRONG_IC_EPS_B, p=FWD_P, xi_N=STRONG_IC_XI),
    )
    model_v = VegasModel(
        jet=VegasTwoComponentJet(0.05, 1.0e53, 300.0, 0.15, 1.0e52, 100.0),
        medium=VegasISM(n_ism=1.0),
        observer=VegasObserver(lumi_dist=1.0e26, z=0.1, theta_obs=0.0),
        fwd_rad=VegasRadiation(eps_e=STRONG_IC_EPS_E, eps_B=STRONG_IC_EPS_B, p=FWD_P, xi_e=STRONG_IC_XI),
    )

    theta = np.linspace(0.0, 0.3, 200)
    e_a = _sample_theta_curve(model_a, "jet_E_iso", theta)
    g_a = _sample_theta_curve(model_a, "jet_Gamma0", theta)
    e_v = _sample_theta_curve(model_v, "jet_E_iso", theta)
    g_v = _sample_theta_curve(model_v, "jet_Gamma0", theta)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.5, 7.0), dpi=180)
    for ax, y_a, y_v, ylab in [
        (ax1, e_a, e_v, r"$E_{\mathrm{iso}}$ [erg]"),
        (ax2, g_a, g_v, r"$\Gamma_0$"),
    ]:
        ax.semilogy(np.degrees(theta), y_a, label="ASGARD")
        ax.semilogy(np.degrees(theta), y_v, ls="--", label="Vegas")
        ax.axvline(np.degrees(0.05), color="red", ls="--", alpha=0.7, label="Core")
        ax.axvline(np.degrees(0.15), color="blue", ls="--", alpha=0.7, label="Wide")
        ax.set_xlabel("Polar Angle [deg]")
        ax.set_ylabel(ylab)
        ax.set_title("two-component profile")
        ax.legend()
        ax.grid(alpha=0.3)
    return _save(fig, OUTPUT_DIR / "compare_introspection_twocomp.png")


def _build_speed_compare() -> Path:
    times = SPEED_TIMES
    bands = SPEED_BANDS
    freqs = SPEED_FREQS
    expo = SPEED_EXPO
    pairs = np.full(expo.shape, 1.0e14)
    asgard_model = _build_asgard_model()
    vegas_model = _build_vegas_model()

    def asgard_quickstart():
        return asgard_model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total

    def asgard_lc():
        return asgard_model.flux_density_grid(times, bands).total

    def asgard_spec():
        return asgard_model.flux_density_grid(SPEED_SPEC_TIMES, SPEED_SPEC_FREQS).total

    def asgard_pairs():
        return asgard_model.flux_density(times, freqs).total

    def asgard_expo():
        return asgard_model.flux_density_exposures(expo, pairs, 0.2 * expo, num_subsamples=4).total

    def asgard_details():
        return asgard_model.details(t_min=1.0e2, t_max=1.0e6)

    def asgard_sky():
        return asgard_model.sky_image(SINGLE_SKY_TIME, nu_obs=1.0e9, fov=500.0 * units.uas, npixel=SPEED_SKY_NPIXEL).image

    def vegas_quickstart():
        return vegas_model.flux_density(np.array([1.0e4]), np.array([1.0e14]))

    def vegas_lc():
        return vegas_model.flux_density_grid(times, bands)

    def vegas_spec():
        return vegas_model.flux_density_grid(SPEED_SPEC_TIMES, SPEED_SPEC_FREQS)

    def vegas_pairs():
        return vegas_model.flux_density(times, freqs)

    def vegas_expo():
        return vegas_model.flux_density_exposures(expo, pairs, 0.2 * expo, num_points=4)

    def vegas_details():
        return vegas_model.details(1.0e2, 1.0e6)

    def vegas_sky():
        return vegas_model.sky_image(SINGLE_SKY_TIME, nu_obs=1.0e9, fov=500.0 * vegas_units.uas, npixel=SPEED_SKY_NPIXEL).image

    asgard_cases = [
        ("quickstart", asgard_quickstart),
        ("lightcurve", asgard_lc),
        ("spectrum", asgard_spec),
        ("pair", asgard_pairs),
        ("exposure", asgard_expo),
        ("details", asgard_details),
        ("sky_image", asgard_sky),
    ]
    vegas_cases = [
        ("quickstart", vegas_quickstart),
        ("lightcurve", vegas_lc),
        ("spectrum", vegas_spec),
        ("pair", vegas_pairs),
        ("exposure", vegas_expo),
        ("details", vegas_details),
        ("sky_image", vegas_sky),
    ]

    asgard_times: list[float] = []
    vegas_times: list[float] = []
    for (name_a, fn_a), (name_v, fn_v) in zip(asgard_cases, vegas_cases):
        if name_a != name_v:
            continue
        t_a, _ = _run_timed(name_a, fn_a)
        t_v, _ = _run_timed(name_v, fn_v)
        asgard_times.append(t_a)
        vegas_times.append(t_v)

    x = np.arange(len(asgard_cases), dtype=float)
    labels = [name for name, _ in asgard_cases]
    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(11.0, 7.0), dpi=180, sharex=True)
    ax_top.bar(x, asgard_times, width=0.45, color="C0", alpha=0.8)
    ax_bottom.bar(x, vegas_times, width=0.45, color="C1", alpha=0.8)
    ax_top.set_ylabel("ASGARD (s)")
    ax_bottom.set_ylabel("VegasAfterglow (s)")
    ax_bottom.set_xlabel("Benchmark case")
    ax_top.set_title("Speed test: ASGARD (top) vs VegasAfterglow (bottom)")
    ax_top.set_xticks(x)
    ax_top.set_xticklabels(labels, rotation=20, ha="right")
    ax_bottom.set_xticks(x)
    ax_bottom.set_xticklabels(labels, rotation=20, ha="right")
    for ax in (ax_top, ax_bottom):
        ax.grid(axis="y", alpha=0.25)
    return _save(fig, OUTPUT_DIR / "compare_speed_profile.png")


def _build_spectrum_compare(*, quantity: str = "sed") -> Path:
    model_asgard = _build_asgard_model(include_ssc=True, num_nu=SPECTRUM_COMPARE_NUM_NU)
    model_vegas = _build_vegas_model(include_ssc=True)
    times = SPEED_SPEC_TIMES
    freqs = SPECTRUM_COMPARE_FREQS

    kind = quantity.lower()
    if kind == "sed":
        kind = "nufnu"
    if kind not in ("nufnu", "fnu"):
        raise ValueError("quantity must be 'sed', 'nufnu', or 'fnu'.")

    as_res = model_asgard.flux_density_grid(times, freqs).total
    vg_res = model_vegas.flux_density_grid(times, freqs).total
    if kind == "nufnu":
        as_res = as_res * freqs[:, None]
        vg_res = vg_res * freqs[:, None]

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0), dpi=200, sharey=True)
    for i, t_obs in enumerate(times):
        ax = axes[i]
        ax.loglog(freqs, np.asarray(as_res[:, i], dtype=float), color="C0", lw=1.8, label="ASGARD")
        ax.loglog(freqs, np.asarray(vg_res[:, i], dtype=float), color="C1", ls="--", lw=1.5, label="Vegas")
        ax.set_title(_label(t_obs, "s"))
        ax.set_xlabel("Frequency [Hz]")
        ax.grid(alpha=0.25, which="both")
    axes[0].set_ylabel(r"$\nu F_\nu$ [erg s$^{-1}$ cm$^{-2}$]" if kind == "nufnu" else r"$F_\nu$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]")
    axes[0].legend(fontsize=8, loc="best")
    fig.suptitle("SED comparison" if kind == "nufnu" else "Spectrum comparison", y=1.02)
    return _save(fig, OUTPUT_DIR / "compare_spectrum.png")


def main(*, spectrum_quantity: str = "sed") -> None:
    builders = [
        ("basic_lc_spec", _build_basic_lc_spec),
        ("basic_bolometric", _build_basic_bolometric),
        ("reverse_shock_lc", _build_reverse_shock_lc),
        ("ssc_lc", _build_ssc_lc),
        ("spectrum_compare", lambda: _build_spectrum_compare(quantity=spectrum_quantity)),
        ("shock_quantities", _build_shock_quantities),
        ("photon_quantities", _build_photon_quantities),
        ("sky_image_single", _build_sky_single),
        ("sky_image_offaxis", _build_sky_offaxis),
        ("sky_image_flux_comparison", _build_sky_flux_comparison),
        ("introspection_jet", _build_introspection_jet),
        ("introspection_medium", _build_introspection_medium),
        ("introspection_twocomp", _build_introspection_twocomp),
    ]

    generated: list[Path] = []
    for i, (name, builder) in enumerate(builders, start=1):
        try:
            print(f"[ {i}/{len(builders)} ] build {name}")
            generated.append(builder())
        except Exception as exc:
            generated.append(_with_note(name, f"{type(exc).__name__}: {exc}"))

    generated.append(_build_speed_compare())

    print(f"output: {OUTPUT_DIR}")
    for path in generated:
        print(path)


if __name__ == "__main__":
    main()
