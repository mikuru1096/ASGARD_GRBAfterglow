from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from functools import lru_cache
import time
from typing import Callable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
import sys

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.asgard_paths import ASGARD_DOC_DIR
from asgard_core.asgard_state import project_flux_grid
from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet, Wind, units
from ASGARD import PowerLawJet as ASGARD_PowerLawJet
from ASGARD import TwoComponentJet
from ASGARD.api_observe import _build_fit_config_for_patch, _solve_patch_state
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

# ============================================================================
# Unified Color Scheme for Scientific Aesthetics
# ============================================================================
ASGARD_COLOR = "#1f77b4"  # Professional blue
VEGAS_COLOR = "#ff7f0e"   # Professional orange
ASGARD_ALPHA = 0.8
VEGAS_ALPHA = 0.7
GRID_ALPHA = 0.25
GRID_STYLE = {"alpha": GRID_ALPHA, "linestyle": ":", "linewidth": 0.5}
HADRONIC_TOTAL_LS = "-."
HADRONIC_ONLY_LS = (0, (1.6, 1.2))
HADRONIC_SOLVER = "am3_1d"
HADRONIC_PGAMMA_SCHEME = "hummer_2010_response"


@dataclass(frozen=True)
class BenchmarkScenario:
    output_subdir: str
    n_ism: float
    theta_c: float
    e_iso: float
    gamma0: float
    lumi_dist: float
    z: float
    eps_e: float
    eps_b: float
    xi_n: float
    epsilon_p: float


BENCHMARK_SCENARIOS: dict[str, BenchmarkScenario] = {
    "baseline": BenchmarkScenario(
        output_subdir="vegas_afterglow_compare",
        n_ism=1.0,
        theta_c=0.1,
        e_iso=1.0e52,
        gamma0=300.0,
        lumi_dist=1.0e26,
        z=0.1,
        eps_e=3.0e-1,
        eps_b=1.0e-5,
        xi_n=1.0e-1,
        epsilon_p=0.2,
    ),
    "hadronic_dominated": BenchmarkScenario(
        output_subdir="vegas_afterglow_compare_hadronic_dominated",
        n_ism=30.0,
        theta_c=0.08,
        e_iso=10.0**54.5,
        gamma0=300.0,
        lumi_dist=1.0e26,
        z=0.1,
        eps_e=1.0e-2,
        eps_b=1.0e-2,
        xi_n=3.0e-2,
        epsilon_p=0.8,
    ),
}
BENCHMARK_SCENARIO = "baseline"

# Plot style settings
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 180,
    'savefig.dpi': 200,
    'axes.grid': True,
    'grid.alpha': GRID_ALPHA,
    'grid.linestyle': ':',
    'grid.linewidth': 0.5,
})

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
STRONG_IC_EPS_E = BENCHMARK_SCENARIOS["baseline"].eps_e
STRONG_IC_EPS_B = BENCHMARK_SCENARIOS["baseline"].eps_b
STRONG_IC_XI = BENCHMARK_SCENARIOS["baseline"].xi_n
DECAY_ALPHA_T = -0.5
DECAY_T0_S = 1.0e2
DECAY_EPSB_FLOOR_FACTOR = 1.0e-3
FWD_P = 2.3
RVS_P = 2.4
BASIC_TIMES = np.logspace(0.0, 8.0, 160)
BASIC_BANDS = np.array([1.0e9, 1.0e14, 1.0e17], dtype=float)
BASIC_FREQS = np.logspace(5.0, 29.0, 160)
BASIC_EPOCHS = np.array([1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8], dtype=float)
MAGNETIC_DECAY_BANDS = np.array([1.0e9, 1.0e12, 1.0e14, 1.0e17, 1.0e20, 1.0e23, 1.0e16 * constants.para_ev2hz], dtype=float)
MAGNETIC_DECAY_FREQS = np.logspace(8.0, np.log10(1.0e16 * constants.para_ev2hz), 320)
MAGNETIC_DECAY_EPOCHS = np.array([1.0e2, 1.0e4, 1.0e6, 1.0e8], dtype=float)
MAGNETIC_DECAY_NUM_NU = 151
MAGNETIC_DECAY_NUM_GAM_E = 121
MAGNETIC_DECAY_NUM_CHI = 16
MAGNETIC_DECAY_LC_SHIFT_DEX = 1.0
OFFAXIS_SKY_TIMES = np.logspace(5.0, 8.0, 3)
SINGLE_SKY_TIME = np.array([1.0e6], dtype=float)
SINGLE_SKY_NPIXEL = 48
FLUX_SKY_TIMES = np.logspace(0.0, 8.0, 20)
FLUX_SKY_NPIXEL = 48
SKY_CENTROID_TIMES = np.logspace(5.0, 8.0, 12)
SPEED_TIMES = np.logspace(2.0, 8.0, 16)
SPEED_BANDS = np.array([1.0e9, 1.0e14, 1.0e17], dtype=float)
SPEED_FREQS = np.logspace(9.0, 18.0, 16)
SPEED_EXPO = np.logspace(2.0, 6.0, 8)
SPEED_SPEC_TIMES = np.array([1.0e3, 1.0e5, 1.0e7], dtype=float)
SPEED_SPEC_FREQS = np.logspace(9.0, 22.0, 24)
SPEED_SKY_NPIXEL = 20
SPECTRUM_COMPARE_FREQS = np.logspace(8.0, 29.0, 240)
SPECTRUM_COMPARE_NUM_NU = 81
ELECTRON_COMPARE_TIMES = np.array([1.0e2, 3.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7], dtype=float)
ELECTRON_COMPARE_SOLVERS = ("fullhide_1d", "fullhide_2d", "slc1_1d", "charint_1d", "charint_2d")
ELECTRON_COMPARE_NUM_GAM_E = {"fullhide_1d": 81, "fullhide_2d": 81, "slc1_1d": 81, "charint_1d": 41, "charint_2d": 41}
ELECTRON_COMPARE_NUM_CHI = {"fullhide_1d": None, "fullhide_2d": 8, "slc1_1d": None, "charint_1d": None, "charint_2d": 8}
ELECTRON_COMPARE_LINESTYLES = {"fullhide_1d": "-", "fullhide_2d": ":", "slc1_1d": "--", "charint_1d": "-.", "charint_2d": (0, (5, 1.4))}


def _benchmark_scenario(name: str | None = None) -> BenchmarkScenario:
    scenario_name = BENCHMARK_SCENARIO if name is None else str(name)
    try:
        return BENCHMARK_SCENARIOS[scenario_name]
    except KeyError as exc:
        known = ", ".join(sorted(BENCHMARK_SCENARIOS))
        raise ValueError(f"unknown benchmark scenario {scenario_name!r}; expected one of: {known}") from exc


def _set_benchmark_context(scenario: str) -> None:
    global BENCHMARK_SCENARIO, OUTPUT_DIR
    _benchmark_scenario(scenario)
    BENCHMARK_SCENARIO = scenario
    OUTPUT_DIR = ASGARD_DOC_DIR / BENCHMARK_SCENARIOS[scenario].output_subdir
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    _build_asgard_model.cache_clear()
    _build_vegas_model.cache_clear()
    _cached_details_pair.cache_clear()


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
    """
    Calibrate sky image to ASGARD coordinate basis.

    Current off-axis comparison baseline shows that the VegasAfterglow image
    basis is rotated relative to the ASGARD basis. A 90-degree clockwise
    rotation is required before comparing image morphology or center-of-light.

    Note: FOV semantics still differ between codes:
    - ASGARD: FOV is a lower bound and may auto-expand to include emission
    - VegasAfterglow: FOV is the fixed rendered window
    """
    arr = np.asarray(image, dtype=float)
    if backend == "VegasAfterglow":
        return np.rot90(arr, 3)
    return arr


def _shared_positive_log_norm(*images: np.ndarray) -> LogNorm:
    positive = [np.asarray(img, dtype=float)[np.asarray(img, dtype=float) > 0.0] for img in images]
    positive = [arr for arr in positive if arr.size > 0]
    if not positive:
        return LogNorm()
    stacked = np.concatenate(positive)
    return LogNorm(vmin=float(np.min(stacked)), vmax=float(np.max(stacked)))


def _centroid_from_image_stack(image_stack: np.ndarray, extent: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    arr = np.asarray(image_stack, dtype=float)
    x_axis = np.linspace(float(extent[0]), float(extent[1]), arr.shape[1])
    y_axis = np.linspace(float(extent[2]), float(extent[3]), arr.shape[2])
    x_weights = arr.sum(axis=2)
    y_weights = arr.sum(axis=1)
    total = arr.sum(axis=(1, 2))
    x_centroid = np.zeros(arr.shape[0], dtype=float)
    y_centroid = np.zeros(arr.shape[0], dtype=float)
    valid = total > 0.0
    if np.any(valid):
        x_centroid[valid] = np.sum(x_weights[valid, :] * x_axis[None, :], axis=1) / total[valid]
        y_centroid[valid] = np.sum(y_weights[valid, :] * y_axis[None, :], axis=1) / total[valid]
    return x_centroid, y_centroid


def _vegas_sky_compare_payload(sky_image) -> dict[str, np.ndarray]:
    """
    Prepare VegasAfterglow sky-image outputs for direct comparison against ASGARD.

    This helper is intentionally comparison-only: it applies the current
    rotation convention, computes image-integrated flux, and derives centroids
    in the ASGARD display basis.
    """
    image = np.asarray(sky_image.image, dtype=float)
    image_cal = np.asarray(
        [_calibrate_sky_image_to_asgard_basis(frame, "VegasAfterglow") for frame in image],
        dtype=float,
    )
    x_centroid, y_centroid = _centroid_from_image_stack(image_cal, sky_image.extent)
    rendered_flux = image.sum(axis=(1, 2)) * float(sky_image.pixel_solid_angle)
    return {
        "image_cal": image_cal,
        "extent_uas": np.asarray(sky_image.extent, dtype=float) / vegas_units.uas,
        "rendered_flux": rendered_flux,
        "x_centroid_uas": x_centroid / vegas_units.uas,
        "y_centroid_uas": y_centroid / vegas_units.uas,
    }


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


def _asgard_radiation(
    params: BenchmarkScenario,
    *,
    p: float,
    include_ssc: bool,
    include_hadronic: bool,
    magnetic_decay_alpha_t: float,
    magnetic_decay_t0_s: float,
    epsilon_b_floor: float | None,
) -> Radiation:
    return Radiation(
        eps_e=params.eps_e,
        eps_B=params.eps_b,
        epsilon_b_floor=epsilon_b_floor,
        magnetic_decay_alpha_t=magnetic_decay_alpha_t,
        magnetic_decay_t0_s=magnetic_decay_t0_s,
        p=p,
        xi_N=params.xi_n,
        ssc=include_ssc,
        kn=True,
        epsilon_p=params.epsilon_p if include_hadronic else 0.0,
        proton_synch=bool(include_hadronic),
        pg=bool(include_hadronic),
        bethe_heitler=bool(include_hadronic),
        hadronic_inverse_compton=bool(include_hadronic),
        pp=bool(include_hadronic),
        neutrino=False,
        pgamma_scheme=HADRONIC_PGAMMA_SCHEME if include_hadronic else "disabled",
    )


def _vegas_radiation(params: BenchmarkScenario, *, p: float, include_ssc: bool) -> VegasRadiation:
    return VegasRadiation(
        eps_e=params.eps_e,
        eps_B=params.eps_b,
        p=p,
        xi_e=params.xi_n,
        ssc=include_ssc,
        kn=True,
    )


@lru_cache(maxsize=None)
def _build_asgard_model(
    include_reverse: bool = False,
    include_ssc: bool = False,
    theta_obs: float = 0.0,
    num_nu: int = 49,
    electron_solver: str = "charint_1d",
    num_gam_e: int | None = None,
    num_chi: int | None = None,
    magnetic_decay_alpha_t: float = 0.0,
    magnetic_decay_t0_s: float = 1.0,
    epsilon_b_floor: float | None = None,
    include_hadronic: bool = False,
    scenario: str | None = None,
) -> Model:
    params = _benchmark_scenario(scenario)
    medium = ISM(n_ism=params.n_ism)
    jet = TophatJet(
        theta_c=params.theta_c,
        E_iso=params.e_iso,
        Gamma0=params.gamma0,
        duration=REVERSE_DURATION_S,
    )
    observer = Observer(lumi_dist=params.lumi_dist, z=params.z, theta_obs=theta_obs)
    fwd_rad = _asgard_radiation(
        params,
        p=FWD_P,
        include_ssc=include_ssc,
        include_hadronic=include_hadronic,
        magnetic_decay_alpha_t=magnetic_decay_alpha_t,
        magnetic_decay_t0_s=magnetic_decay_t0_s,
        epsilon_b_floor=epsilon_b_floor,
    )
    kwargs = dict(
        num_gam_e=ASGARD_CHARINT_NUM_GAM_E if num_gam_e is None else int(num_gam_e),
        num_nu=num_nu,
        num_r=REVERSE_NUM_R if include_reverse else BASE_NUM_R,
        num_theta=80,
        num_tobs=48,
        electron_adaptive_substeps=False,
        electron_solver=electron_solver,
        num_chi=num_chi,
        hadronic_enabled=bool(include_hadronic),
        hadronic_solver=HADRONIC_SOLVER if include_hadronic else "legacy_1d",
        num_gam_p=81 if include_hadronic else 161,
        num_nu_nu=81 if include_hadronic else 121,
        pgamma_scheme=HADRONIC_PGAMMA_SCHEME if include_hadronic else "disabled",
    )
    return Model(
        jet=jet,
        medium=medium,
        observer=observer,
        fwd_rad=fwd_rad,
        rvs_rad=_asgard_radiation(
            params,
            p=RVS_P,
            include_ssc=include_ssc,
            include_hadronic=include_hadronic,
            magnetic_decay_alpha_t=magnetic_decay_alpha_t,
            magnetic_decay_t0_s=magnetic_decay_t0_s,
            epsilon_b_floor=epsilon_b_floor,
        ) if include_reverse else None,
        setups=Setups(reverse_delta_t_s=REVERSE_DURATION_S, **kwargs),
    )


def _project_asgard_hadronic_flux(model: Model, times_s: np.ndarray, frequencies_hz: np.ndarray) -> np.ndarray:
    config = _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )
    state = _solve_patch_state(model, config, np.asarray(times_s, dtype=float), np.asarray(frequencies_hz, dtype=float))
    observed = project_flux_grid(state, np.asarray(times_s, dtype=float), np.asarray(frequencies_hz, dtype=float))
    hadronic = observed.components.get("fwd_hadronic")
    if hadronic is None:
        return np.zeros((np.asarray(frequencies_hz).size, np.asarray(times_s).size), dtype=float)
    return np.asarray(hadronic, dtype=float)


@lru_cache(maxsize=None)
def _build_vegas_model(
    include_reverse: bool = False,
    include_ssc: bool = False,
    theta_obs: float = 0.0,
    scenario: str | None = None,
) -> VegasModel:
    params = _benchmark_scenario(scenario)
    medium = VegasISM(n_ism=params.n_ism)
    jet = VegasTophatJet(
        theta_c=params.theta_c,
        E_iso=params.e_iso,
        Gamma0=params.gamma0,
        duration=REVERSE_DURATION_S,
    )
    observer = VegasObserver(lumi_dist=params.lumi_dist, z=params.z, theta_obs=theta_obs)
    fwd_rad = _vegas_radiation(params, p=FWD_P, include_ssc=include_ssc)
    return VegasModel(
        jet=jet,
        medium=medium,
        observer=observer,
        fwd_rad=fwd_rad,
        rvs_rad=_vegas_radiation(params, p=RVS_P, include_ssc=include_ssc) if include_reverse else None,
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
    model_asgard_had = _build_asgard_model(include_hadronic=True)
    model_vegas = _build_vegas_model()
    times = BASIC_TIMES
    bands = BASIC_BANDS
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.0), dpi=200)

    lc_as = model_asgard.flux_density_grid(times, bands).total
    lc_as_had = model_asgard_had.flux_density_grid(times, bands).total
    lc_vg = model_vegas.flux_density_grid(times, bands).total
    colors = plt.cm.tab10(np.linspace(0, 1, len(bands)))
    for i, nu in enumerate(bands):
        label = _label(nu, "Hz")
        axes[0].loglog(times, lc_as[i, :], color=colors[i], ls="-", lw=1.8, alpha=ASGARD_ALPHA, label=f"ASGARD {label}")
        axes[0].loglog(times, lc_as_had[i, :], color=colors[i], ls=HADRONIC_TOTAL_LS, lw=1.5, alpha=0.9, label=f"ASGARD+had {label}")
        axes[0].loglog(times, lc_vg[i, :], color=colors[i], ls="--", lw=1.4, alpha=VEGAS_ALPHA, label=f"Vegas {label}")

    freqs = BASIC_FREQS
    epochs = BASIC_EPOCHS
    sed_as = model_asgard.flux_density_grid(epochs, freqs).total * freqs[:, None]
    sed_as_had_total = model_asgard_had.flux_density_grid(epochs, freqs).total * freqs[:, None]
    sed_as_had_only = _project_asgard_hadronic_flux(model_asgard_had, epochs, freqs) * freqs[:, None]
    sed_vg = model_vegas.flux_density_grid(epochs, freqs).total * freqs[:, None]
    spec_peak = float(max(np.nanmax(sed_as), np.nanmax(sed_as_had_total), np.nanmax(sed_as_had_only), np.nanmax(sed_vg)))
    colors2 = plt.cm.viridis(np.linspace(0, 1, len(epochs)))
    for i, t in enumerate(epochs):
        label = _label(t, "s")
        axes[1].loglog(freqs, sed_as[:, i], color=colors2[i], ls="-", lw=1.6, alpha=ASGARD_ALPHA, label=f"ASGARD t={label}")
        axes[1].loglog(freqs, sed_as_had_total[:, i], color=colors2[i], ls=HADRONIC_TOTAL_LS, lw=1.35, alpha=0.95, label=f"ASGARD+had t={label}")
        axes[1].loglog(freqs, np.maximum(sed_as_had_only[:, i], 1.0e-300), color=colors2[i], ls=HADRONIC_ONLY_LS, lw=1.15, alpha=0.95, label=f"had-only t={label}")
        axes[1].loglog(freqs, sed_vg[:, i], color=colors2[i], ls="--", lw=1.2, alpha=VEGAS_ALPHA, label=f"Vegas t={label}")
    for idx, band in enumerate(bands):
        axes[1].axvline(band, ls="--", alpha=0.4, color=f"C{idx}")

    axes[0].set_xlabel("Time [s]")
    axes[0].set_ylabel(r"Flux Density (erg/cm$^2$/s/Hz)")
    axes[0].set_title("Light Curves")
    axes[1].set_xlabel("Frequency [Hz]")
    axes[1].set_ylabel(r"$\nu F_\nu$ (erg/cm$^2$/s)")
    axes[1].set_title("SED")
    if spec_peak > 0.0:
        axes[1].set_ylim(bottom=spec_peak * 1.0e-17)
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
    ax.loglog(times, bat_as, label="ASGARD Swift/BAT", lw=2, color="C0", alpha=ASGARD_ALPHA)
    ax.loglog(times, bat_vg, label="VegasAfterglow Swift/BAT", lw=2, ls="--", color="C0", alpha=VEGAS_ALPHA)
    ax.loglog(times, opt_as, label="ASGARD V-band", lw=2, color="C1", alpha=ASGARD_ALPHA)
    ax.loglog(times, opt_vg, label="VegasAfterglow V-band", lw=2, ls="--", color="C1", alpha=VEGAS_ALPHA)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"Integrated Flux (erg/cm$^2$/s)")
    ax.set_title("Broadband Light Curves")
    ax.legend(fontsize=8)
    ax.grid(**GRID_STYLE)
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
        rvs_a = np.asarray(as_res.rev.sync[i, :], dtype=float)
        rvs_v = np.asarray(vg_res.rvs.sync[i, :], dtype=float)
        ax.loglog(times, fwd_a, color=f"C{i}", lw=1.8, alpha=ASGARD_ALPHA, label=f"ASGARD fwd {label}")
        ax.loglog(times, fwd_v, color=f"C{i}", ls="--", lw=1.3, alpha=VEGAS_ALPHA, label=f"Vegas fwd {label}")
        ax.loglog(times, rvs_a, color=f"C{i}", ls="-.", lw=1.8, alpha=0.7, label=f"ASGARD rvs {label}")
        ax.loglog(times, rvs_v, color=f"C{i}", ls=":", lw=1.6, alpha=0.6, label=f"Vegas rvs {label}")
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
    ax.grid(**GRID_STYLE)
    return _save(fig, OUTPUT_DIR / "compare_reverse_shock_lc.png")


def _single_patch_state(model: Model, times_s: np.ndarray, frequencies_hz: np.ndarray):
    config = _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )
    return _solve_patch_state(model, config, times_s, frequencies_hz, solve_reference_times_s=times_s)


def _plot_positive_series(ax, x: np.ndarray, y: np.ndarray, *, label: str, color: str, linestyle: str = "-") -> None:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr) & (x_arr > 0.0) & (y_arr > 0.0)
    if np.any(mask):
        ax.loglog(x_arr[mask], y_arr[mask], color=color, ls=linestyle, lw=1.45, label=label)


def _build_reverse_shock_thermal_benchmark() -> Path:
    model_asgard = _build_asgard_model(include_reverse=True, include_ssc=True)
    model_vegas = _build_vegas_model(include_reverse=True, include_ssc=True)
    detail_asgard = model_asgard.details(t_min=1.0, t_max=1.0e8)
    detail_vegas = model_vegas.details(t_min=1.0, t_max=1.0e8)
    state = _single_patch_state(model_asgard, BASIC_TIMES, BASIC_BANDS)
    if state.dynamics.reverse_shock is None:
        raise RuntimeError("reverse-shock dynamics are required for the RS thermal benchmark.")

    t_ref = _reference_series(detail_asgard.rev.t_obs)
    t_vegas = _reference_series(detail_vegas.rvs.t_obs)
    attrs = ("Gamma", "B_comv", "N_p", "nu_m", "nu_c", "nu_a")
    vegas_z = _benchmark_scenario().z

    fig, axes = plt.subplots(3, 3, figsize=(12.8, 9.2), dpi=200)
    ratio_axes = axes[:2, :].ravel()
    for ax, attr in zip(ratio_axes, attrs):
        asgard_y = _reference_series(getattr(detail_asgard.rev, attr))
        vegas_y = _reference_series(getattr(detail_vegas.rvs, attr))
        if attr == "N_p":
            vegas_y = vegas_y * (4.0 * np.pi)
        if attr in {"nu_m", "nu_c"}:
            vegas_y = _to_lab_frequency_frame(vegas_y, _reference_series(detail_vegas.rvs.Doppler), vegas_z)
        vegas_interp = _safe_log_interp(t_ref, t_vegas, vegas_y)
        ratio = np.divide(
            asgard_y,
            vegas_interp,
            out=np.full_like(asgard_y, np.nan, dtype=float),
            where=(asgard_y > 0.0) & (vegas_interp > 0.0),
        )
        _plot_positive_series(ax, t_ref, ratio, label="ASGARD / Vegas", color=ASGARD_COLOR)
        ax.axhline(1.0, color="black", ls=":", lw=0.9)
        ax.set_title(f"RS {attr} ratio")
        ax.set_xlabel("t_obs [s]")
        ax.set_ylabel("ratio")
        ax.grid(**GRID_STYLE)

    rs = state.dynamics.reverse_shock
    t_dyn = np.asarray(state.dynamics.r_tobs, dtype=float)
    radius_dyn = np.asarray(state.dynamics.radius, dtype=float)
    pre_cross = radius_dyn <= float(rs.r_cross)
    t_cross_obs = float(t_dyn[int(np.argmin(np.abs(radius_dyn - float(rs.r_cross))))])
    e3 = np.asarray(rs.internal_energy_erg, dtype=float) / np.asarray(rs.comoving_volume_cm3, dtype=float)
    u3_per_m3 = np.asarray(rs.internal_energy_erg, dtype=float) / np.asarray(rs.swept_mass_g, dtype=float)
    gamma34_injection = np.where(pre_cross, np.asarray(rs.gamma34, dtype=float), np.nan)

    _plot_positive_series(axes[2, 0], t_dyn, gamma34_injection, label=r"$\gamma_{34}$", color=ASGARD_COLOR)
    axes[2, 0].set_title(r"ASGARD pre-crossing injection $\gamma_{34}$")
    axes[2, 0].set_xlabel("t_obs [s]")
    axes[2, 0].set_ylabel(r"$\gamma_{34}$")

    _plot_positive_series(axes[2, 1], t_dyn, e3, label=r"$U_3/V_3$", color=ASGARD_COLOR)
    axes[2, 1].set_title(r"ASGARD region-3 thermal density")
    axes[2, 1].set_xlabel("t_obs [s]")
    axes[2, 1].set_ylabel(r"$U_3/V_3$ [erg cm$^{-3}$]")

    _plot_positive_series(axes[2, 2], t_dyn, u3_per_m3, label=r"$U_3/M_3$", color=ASGARD_COLOR)
    axes[2, 2].set_title(r"ASGARD post-crossing thermal scale")
    axes[2, 2].set_xlabel("t_obs [s]")
    axes[2, 2].set_ylabel(r"$U_3/M_3$ [erg g$^{-1}$]")

    for ax in axes[2, :]:
        ax.axvline(t_cross_obs, color="black", ls=":", lw=0.9, label="crossing")
        ax.grid(**GRID_STYLE)
        ax.legend(fontsize=8)
    fig.suptitle("Reverse-shock benchmark: Vegas comparison and ASGARD thermal closure", y=0.995)
    fig.tight_layout()
    return _save(fig, OUTPUT_DIR / "compare_reverse_shock_thermal_benchmark.png")


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
        ax.loglog(times, as_res.fwd.sync[i, :], color=f"C{i}", lw=1.8, alpha=ASGARD_ALPHA, label=f"ASGARD sync {label}")
        ax.loglog(times, as_res.fwd.ssc[i, :], ls="--", color=f"C{i}", lw=1.4, alpha=ASGARD_ALPHA, label=f"ASGARD SSC {label}")
        ax.loglog(times, vg_res.fwd.sync[i, :], color=f"C{i}", ls=":", lw=1.2, alpha=VEGAS_ALPHA, label=f"Vegas sync {label}")
        ax.loglog(times, vg_res.fwd.ssc[i, :], color=f"C{i}", ls="-.", lw=1.0, alpha=VEGAS_ALPHA, label=f"Vegas SSC {label}")

    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"Flux Density (erg/cm$^2$/s/Hz)")
    ax.set_title("Synchrotron and SSC")
    ax.legend(fontsize=6.8, ncol=2)
    ax.grid(**GRID_STYLE)
    return _save(fig, OUTPUT_DIR / "compare_ssc_lc.png")


def _plot_two_line_sets(ax, x: np.ndarray, lines: list[tuple[np.ndarray, str, str]], x_label: str) -> None:
    x = np.asarray(x, dtype=float)
    for y, color, label in lines:
        y = np.asarray(y, dtype=float)
        mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
        if not mask.any():
            continue
        ax.loglog(x[mask], y[mask], color=color, label=label, alpha=0.8)
    ax.set_xlabel(x_label)
    ax.grid(**GRID_STYLE)


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


def _build_shock_quantities() -> Path:
    da, dv = _cached_details_pair(False, False, 0.0, 1.0, 1.0e8)
    vegas_z = 0.1

    attrs = ["Gamma", "B_comv", "N_p", "nu_m", "nu_c", "nu_a"]
    t_ref = _reference_series(da.fwd.t_obs)
    if t_ref.ndim != 1:
        raise RuntimeError("ASGARD shock details are not 1D; compare_shock_quantities requires a 1D characteristic-time track.")

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.0), dpi=200)
    axes = axes.ravel()
    for i, attr in enumerate(attrs):
        ax = axes[i]
        if not (hasattr(da.fwd, attr) and hasattr(dv.fwd, attr)):
            raise RuntimeError(f"{attr} is not available in both ASGARD and VegasAfterglow shock details.")
        ya = _reference_series(getattr(da.fwd, attr))
        yv = _reference_series(getattr(dv.fwd, attr))
        if attr == "N_p":
            yv = yv * (4.0 * np.pi)
        if attr in {"nu_m", "nu_c"}:
            doppler_v = _reference_series(dv.fwd.Doppler)
            yv = _to_lab_frequency_frame(yv, doppler_v, vegas_z)
        tv = _reference_series(dv.fwd.t_obs)
        yv = _safe_log_interp(t_ref, tv, yv)
        _plot_two_line_sets(ax, t_ref, [(ya, ASGARD_COLOR, f"ASGARD {attr}"), (yv, VEGAS_COLOR, f"Vegas {attr}")], "t_obs [s]")
        ax.set_ylabel(attr)
        ax.set_xlabel("t_obs [s]")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.xaxis.set_major_formatter(LogFormatter())
        ax.set_title(attr)
    plt.tight_layout()
    return _save(fig, OUTPUT_DIR / "compare_shock_quantities.png")


def _build_photon_quantities() -> Path:
    da, dv = _cached_details_pair(False, False, 0.0, 1.0, 1.0e8)
    vegas_z = 0.1
    attrs = ["nu_a", "nu_m", "nu_c"]
    t_ref = _reference_series(da.fwd.t_obs)
    if t_ref.ndim != 1:
        raise RuntimeError("ASGARD photon quantities are not 1D; compare_photon_quantities requires a 1D characteristic-time track.")

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
        _plot_two_line_sets(ax, t_ref, [(ya, ASGARD_COLOR, f"ASGARD {attr}"), (yv, VEGAS_COLOR, f"Vegas {attr}")], "t_obs [s]")
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
    vegas_payload = _vegas_sky_compare_payload(sky_v)
    img_v = vegas_payload["image_cal"][0]
    img_a = _calibrate_sky_image_to_asgard_basis(sky_a.image[0], "ASGARD")
    norm = _shared_positive_log_norm(img_v, img_a)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), dpi=180)
    for ax, sky, title, cmap, unit in [
        (axes[0], sky_v, "VegasAfterglow", "magma", vegas_units.uas),
        (axes[1], sky_a, "ASGARD", "inferno", units.uas),
    ]:
        img = img_v if title == "VegasAfterglow" else img_a
        extent = vegas_payload["extent_uas"] if title == "VegasAfterglow" else sky.extent / unit
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
    vegas_payload = _vegas_sky_compare_payload(img_v)
    idx = -1
    img_v_cal = vegas_payload["image_cal"][idx]
    img_a_cal = _calibrate_sky_image_to_asgard_basis(img_a.image[idx], "ASGARD")
    norm = _shared_positive_log_norm(img_v_cal, img_a_cal)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), dpi=180)
    extent_v = vegas_payload["extent_uas"]
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
    vegas_payload = _vegas_sky_compare_payload(img_v)
    flux_from_image_a = np.asarray(img_a.rendered_flux, dtype=float)
    flux_from_image_v = vegas_payload["rendered_flux"]
    direct_a = np.asarray(img_a.direct_flux, dtype=float)
    direct_v = model_v.flux_density_grid(times, np.array([nu_obs])).total[0, :]

    fig, axes = plt.subplots(2, 1, figsize=(6.6, 6.0), dpi=180, sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    axes[0].loglog(times, direct_a, "k-", lw=2.0, label="ASGARD direct")
    axes[0].loglog(times, flux_from_image_a, "o", ms=3.5, color="C1", label="ASGARD image")
    axes[0].loglog(times, direct_v, "k--", lw=1.8, label="Vegas direct")
    axes[0].loglog(times, flux_from_image_v, "x", ms=3.2, color="C3", label="Vegas image")
    axes[0].set_ylabel(r"Flux density (erg/cm$^2$/s/Hz)")
    axes[0].legend(fontsize=8)

    ratio_a = np.asarray(img_a.flux_ratio, dtype=float)
    ratio_v = flux_from_image_v / np.maximum(direct_v, np.finfo(float).tiny)
    axes[1].semilogx(times, ratio_a, "o-", ms=3.5, color="C1", label="ASGARD")
    axes[1].semilogx(times, ratio_v, "x-", ms=3.5, color="C3", label="Vegas")
    axes[1].axhline(1.0, color="k", ls="--", lw=0.8)
    axes[1].set_xlabel("Observer time [s]")
    axes[1].set_ylabel("ratio")
    axes[1].legend(fontsize=8)
    axes[1].set_ylim(0.95, 1.05)

    return _save(fig, OUTPUT_DIR / "compare_sky_image_flux_comparison.png")


def _build_sky_centroid_comparison() -> Path:
    model_a = _build_asgard_model(theta_obs=0.4)
    model_v = _build_vegas_model(theta_obs=0.4)
    times = SKY_CENTROID_TIMES
    nu_obs = 1.0e9

    img_a = model_a.sky_image(times, nu_obs=nu_obs, fov=5000.0 * units.uas, npixel=40)
    img_v = model_v.sky_image(times, nu_obs=nu_obs, fov=5000.0 * vegas_units.uas, npixel=40)
    vegas_payload = _vegas_sky_compare_payload(img_v)

    centroid_a = np.asarray(img_a.centroid, dtype=float) / units.uas
    x_a = centroid_a[:, 0]
    y_a = centroid_a[:, 1]
    x_v = vegas_payload["x_centroid_uas"]
    y_v = vegas_payload["y_centroid_uas"]

    fig, axes = plt.subplots(2, 1, figsize=(6.8, 6.0), dpi=180, sharex=True, gridspec_kw={"hspace": 0.08})
    axes[0].semilogx(times, x_a, color=ASGARD_COLOR, lw=1.6, alpha=ASGARD_ALPHA, label="ASGARD")
    axes[0].semilogx(times, x_v, color=VEGAS_COLOR, ls="--", lw=1.4, alpha=VEGAS_ALPHA, label="Vegas")
    axes[0].set_ylabel(r"$x_c$ ($\mu$as)")
    axes[0].legend(fontsize=8)
    axes[0].grid(**GRID_STYLE)

    axes[1].semilogx(times, y_a, color=ASGARD_COLOR, lw=1.6, alpha=ASGARD_ALPHA, label="ASGARD")
    axes[1].semilogx(times, y_v, color=VEGAS_COLOR, ls="--", lw=1.4, alpha=VEGAS_ALPHA, label="Vegas")
    axes[1].set_xlabel("Observer time [s]")
    axes[1].set_ylabel(r"$y_c$ ($\mu$as)")
    axes[1].grid(**GRID_STYLE)

    fig.suptitle("Sky-image centroid comparison", y=0.98)
    return _save(fig, OUTPUT_DIR / "compare_sky_image_centroid.png")


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

    def asgard_cases(model: Model) -> list[tuple[str, Callable[[], object]]]:
        return [
            ("quickstart", lambda: model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total),
            ("lightcurve", lambda: model.flux_density_grid(times, bands).total),
            ("spectrum", lambda: model.flux_density_grid(SPEED_SPEC_TIMES, SPEED_SPEC_FREQS).total),
            ("pair", lambda: model.flux_density(times, freqs).total),
            ("exposure", lambda: model.flux_density_exposures(expo, pairs, 0.2 * expo, num_subsamples=4).total),
            ("details", lambda: model.details(t_min=1.0e2, t_max=1.0e6)),
            ("sky_image", lambda: model.sky_image(SINGLE_SKY_TIME, nu_obs=1.0e9, fov=500.0 * units.uas,
                                                   npixel=SPEED_SKY_NPIXEL).image),
        ]

    def make_vegas_cases(model: VegasModel) -> list[tuple[str, Callable[[], object]]]:
        def vegas_quickstart():
            return model.flux_density(np.array([1.0e4]), np.array([1.0e14]))

        def vegas_lc():
            return model.flux_density_grid(times, bands)

        def vegas_spec():
            return model.flux_density_grid(SPEED_SPEC_TIMES, SPEED_SPEC_FREQS)

        def vegas_pairs():
            return model.flux_density(times, freqs)

        def vegas_expo():
            return model.flux_density_exposures(expo, pairs, 0.2 * expo, num_points=4)

        def vegas_details():
            return model.details(1.0e2, 1.0e6)

        def vegas_sky():
            return model.sky_image(SINGLE_SKY_TIME, nu_obs=1.0e9, fov=500.0 * vegas_units.uas,
                                   npixel=SPEED_SKY_NPIXEL).image

        return [
            ("quickstart", vegas_quickstart),
            ("lightcurve", vegas_lc),
            ("spectrum", vegas_spec),
            ("pair", vegas_pairs),
            ("exposure", vegas_expo),
            ("details", vegas_details),
            ("sky_image", vegas_sky),
        ]

    builders = [
        (
            "ASGARD fullhide_1d",
            lambda: _build_asgard_model(electron_solver="fullhide_1d", num_gam_e=81),
            asgard_cases,
        ),
        (
            "ASGARD charint_1d",
            lambda: _build_asgard_model(electron_solver="charint_1d"),
            asgard_cases,
        ),
        ("VegasAfterglow", _build_vegas_model, make_vegas_cases),
    ]

    labels = ["quickstart", "lightcurve", "spectrum", "pair", "exposure", "details", "sky_image"]
    cold_times = np.zeros((len(builders), len(labels)), dtype=float)
    cached_times = np.zeros_like(cold_times)
    for i_backend, (_backend_label, build_model, make_cases) in enumerate(builders):
        for i_case, label in enumerate(labels):
            _build_asgard_model.cache_clear()
            _build_vegas_model.cache_clear()
            model = build_model()
            cases = dict(make_cases(model))
            cold_times[i_backend, i_case], _ = _run_timed(label, cases[label])

            _build_asgard_model.cache_clear()
            _build_vegas_model.cache_clear()
            model = build_model()
            cases = dict(make_cases(model))
            cases[label]()
            cached_times[i_backend, i_case], _ = _run_timed(label, cases[label])

    x = np.arange(len(labels), dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.6), dpi=180, sharey=True)
    width = 0.25
    offsets = (-width, 0.0, width)
    colors = ("C0", "C2", "C1")
    for ax, matrix, title in zip(axes, (cold_times, cached_times), ("Cold solve", "Cached query")):
        for i_backend, (backend_label, _build_model, _make_cases) in enumerate(builders):
            ax.bar(x + offsets[i_backend], matrix[i_backend], width=width, color=colors[i_backend], alpha=0.85,
                   label=backend_label)
        ax.set_yscale("log")
        ax.set_title(title)
        ax.set_xlabel("Benchmark case")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=24, ha="right")
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Runtime (s)")
    axes[0].legend(frameon=False, loc="upper right")
    fig.suptitle("Speed profile: explicit first solve and warmed-cache query")
    fig.tight_layout()
    return _save(fig, OUTPUT_DIR / "compare_speed_profile.png")


def _build_spectrum_compare(*, quantity: str = "sed") -> Path:
    model_asgard = _build_asgard_model(include_ssc=True, num_nu=SPECTRUM_COMPARE_NUM_NU)
    model_asgard_had = _build_asgard_model(include_ssc=True, num_nu=SPECTRUM_COMPARE_NUM_NU, include_hadronic=True)
    model_vegas = _build_vegas_model(include_ssc=True)
    times = SPEED_SPEC_TIMES
    freqs = SPECTRUM_COMPARE_FREQS

    kind = quantity.lower()
    if kind == "sed":
        kind = "nufnu"
    if kind not in ("nufnu", "fnu"):
        raise ValueError("quantity must be 'sed', 'nufnu', or 'fnu'.")

    as_res = model_asgard.flux_density_grid(times, freqs).total
    as_had_total = model_asgard_had.flux_density_grid(times, freqs).total
    as_had_only = _project_asgard_hadronic_flux(model_asgard_had, times, freqs)
    vg_res = model_vegas.flux_density_grid(times, freqs).total
    if kind == "nufnu":
        as_res = as_res * freqs[:, None]
        as_had_total = as_had_total * freqs[:, None]
        as_had_only = as_had_only * freqs[:, None]
        vg_res = vg_res * freqs[:, None]

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.0), dpi=200, sharey=True)
    for i, t_obs in enumerate(times):
        ax = axes[i]
        ax.loglog(freqs, np.asarray(as_res[:, i], dtype=float), color=ASGARD_COLOR, lw=1.8, alpha=ASGARD_ALPHA, label="ASGARD")
        ax.loglog(freqs, np.asarray(as_had_total[:, i], dtype=float), color=ASGARD_COLOR, ls=HADRONIC_TOTAL_LS, lw=1.45, alpha=0.95, label="ASGARD+had")
        ax.loglog(freqs, np.maximum(np.asarray(as_had_only[:, i], dtype=float), 1.0e-300), color=ASGARD_COLOR, ls=HADRONIC_ONLY_LS, lw=1.15, alpha=0.95, label="had-only")
        ax.loglog(freqs, np.asarray(vg_res[:, i], dtype=float), color=VEGAS_COLOR, ls="--", lw=1.5, alpha=VEGAS_ALPHA, label="Vegas")
        ax.set_title(_label(t_obs, "s"))
        ax.set_xlabel("Frequency [Hz]")
        ax.grid(**GRID_STYLE)
    axes[0].set_ylabel(r"$\nu F_\nu$ [erg s$^{-1}$ cm$^{-2}$]" if kind == "nufnu" else r"$F_\nu$ [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]")
    axes[0].legend(fontsize=8, loc="best")
    fig.suptitle("SED comparison" if kind == "nufnu" else "Spectrum comparison", y=1.02)
    return _save(fig, OUTPUT_DIR / "compare_spectrum.png")


def _build_solver_sed_compare() -> Path:
    times = SPEED_SPEC_TIMES
    freqs_hz = np.logspace(np.log10(SPECTRUM_COMPARE_FREQS[0]), np.log10(1.0e16 * constants.para_ev2hz), 280)
    energy_ev = freqs_hz / constants.para_ev2hz

    vegas_sed = np.asarray(_build_vegas_model(include_ssc=True).flux_density_grid(times, freqs_hz).total, dtype=float)
    vegas_sed = vegas_sed * freqs_hz[:, None]

    asgard_sed: dict[str, np.ndarray] = {}
    for solver in ELECTRON_COMPARE_SOLVERS:
        model_asgard = _build_asgard_model(
            include_ssc=True,
            num_nu=SPECTRUM_COMPARE_NUM_NU,
            electron_solver=solver,
            num_gam_e=ELECTRON_COMPARE_NUM_GAM_E[solver],
            num_chi=ELECTRON_COMPARE_NUM_CHI[solver],
        )
        sed = np.asarray(model_asgard.flux_density_grid(times, freqs_hz).total, dtype=float)
        asgard_sed[solver] = sed * freqs_hz[:, None]

    spec_peak = float(max(np.nanmax(vegas_sed), max(np.nanmax(sed) for sed in asgard_sed.values())))
    late_spec_peak = float(
        max(
            np.nanmax(vegas_sed[:, -1]),
            max(np.nanmax(sed[:, -1]) for sed in asgard_sed.values()),
        )
    )
    colors = plt.cm.viridis(np.linspace(0.0, 1.0, len(times)))
    fig, ax = plt.subplots(1, 1, figsize=(7.6, 5.2), dpi=200)
    time_handles: list[Line2D] = []
    solver_handles: list[Line2D] = []

    for i_time, t_obs in enumerate(times):
        label = _label(t_obs, "s")
        time_handles.append(Line2D([0], [0], color=colors[i_time], lw=1.8, ls="-", alpha=0.9, label=f"t={label}"))
        ax.loglog(
            energy_ev,
            np.asarray(vegas_sed[:, i_time], dtype=float),
            color=colors[i_time],
            ls="--",
            lw=1.2,
            alpha=VEGAS_ALPHA,
        )
        for solver in ELECTRON_COMPARE_SOLVERS:
            ax.loglog(
                energy_ev,
                np.asarray(asgard_sed[solver][:, i_time], dtype=float),
                color=colors[i_time],
                ls=ELECTRON_COMPARE_LINESTYLES[solver],
                lw=1.6,
                alpha=ASGARD_ALPHA,
            )

    for solver in ELECTRON_COMPARE_SOLVERS:
        solver_handles.append(
            Line2D([0], [0], color="black", lw=1.6, ls=ELECTRON_COMPARE_LINESTYLES[solver], label=f"ASGARD {solver}")
        )
    solver_handles.append(
        Line2D([0], [0], color="black", lw=1.2, ls="--", label="VegasAfterglow")
    )

    ax.set_xlabel("Energy [eV]")
    ax.set_ylabel(r"$\nu F_\nu$ (erg/cm$^2$/s)")
    ax.set_title("SED Comparison")
    ax.grid(**GRID_STYLE)
    if spec_peak > 0.0:
        ax.set_ylim(late_spec_peak * 1.0e-5, spec_peak * 2.0)
    ax.set_xlim(float(energy_ev[0]), 1.0e16)
    legend_times = ax.legend(
        handles=time_handles,
        fontsize=7,
        ncol=1,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0,
        title="Epochs",
        frameon=False,
    )
    ax.add_artist(legend_times)
    ax.legend(
        handles=solver_handles,
        fontsize=7,
        ncol=1,
        loc="upper left",
        bbox_to_anchor=(1.01, 0.55),
        borderaxespad=0.0,
        title="Curves",
        frameon=False,
    )
    return _save(fig, OUTPUT_DIR / "compare_sed_three_solvers_vs_vegas.png")


def _build_electron_spectrum_compare(*, magnetic_decay: bool = False) -> Path:
    """Plot ASGARD electron-spectrum evolution for multiple solvers and times."""
    n_times = ELECTRON_COMPARE_TIMES.size
    n_cols = 4
    n_rows = int(np.ceil(n_times / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14.2, 7.2), dpi=200, sharex=True, sharey=True)
    axes_flat = np.asarray(axes, dtype=object).ravel()
    solver_colors = {
        "fullhide_1d": "#111111",
        "fullhide_2d": "#0072B2",
        "slc1_1d": "#D55E00",
        "charint_1d": "#009E73",
        "charint_2d": "#D55E00",
    }
    solver_handles: list[Line2D] = []
    for solver in ELECTRON_COMPARE_SOLVERS:
        solver_handles.append(
            Line2D(
                [0],
                [0],
                color=solver_colors[solver],
                lw=1.4,
                ls=ELECTRON_COMPARE_LINESTYLES[solver],
                label=solver,
            )
        )

    for solver in ELECTRON_COMPARE_SOLVERS:
        decay_kwargs = {}
        if magnetic_decay:
            decay_kwargs = dict(
                magnetic_decay_alpha_t=DECAY_ALPHA_T,
                magnetic_decay_t0_s=DECAY_T0_S,
                epsilon_b_floor=STRONG_IC_EPS_B * DECAY_EPSB_FLOOR_FACTOR,
            )
        model_asgard = _build_asgard_model(
            electron_solver=solver,
            num_gam_e=ELECTRON_COMPARE_NUM_GAM_E[solver],
            num_chi=ELECTRON_COMPARE_NUM_CHI[solver],
            **decay_kwargs,
        )
        details = model_asgard.details(float(ELECTRON_COMPARE_TIMES[0]), float(ELECTRON_COMPARE_TIMES[-1]))
        if details.fwd.gamma_e is None or details.fwd.dN_dgamma_e is None:
            continue
        gamma_e_asgard = np.asarray(details.fwd.gamma_e, dtype=float)
        dnde_asgard = np.asarray(details.fwd.dN_dgamma_e, dtype=float)
        characteristic_times = np.asarray(details.fwd.t_obs, dtype=float)
        for i_panel, t_obs in enumerate(ELECTRON_COMPARE_TIMES):
            ax = axes_flat[i_panel]
            i_time = int(np.argmin(np.abs(characteristic_times - t_obs)))
            asgard_slice = dnde_asgard[:, i_time]
            mask_a = (
                np.isfinite(gamma_e_asgard)
                & np.isfinite(asgard_slice)
                & (gamma_e_asgard > 0.0)
                & (asgard_slice > 0.0)
            )
            if np.any(mask_a):
                ax.loglog(
                    gamma_e_asgard[mask_a],
                    asgard_slice[mask_a],
                    color=solver_colors[solver],
                    lw=1.05,
                    alpha=0.9,
                    ls=ELECTRON_COMPARE_LINESTYLES[solver],
                )
                ax.set_title(fr"$t_{{\rm obs}}\approx {t_obs:.1e}\,\rm s$" "\n" fr"shell $t={characteristic_times[i_time]:.2e}\,\rm s$")

    for ax in axes_flat[n_times:]:
        ax.set_visible(False)
    for ax in axes_flat[:n_times]:
        ax.set_ylim(bottom=1.0e15)
        ax.grid(**GRID_STYLE)
    for ax in axes_flat[(n_rows - 1) * n_cols:n_times]:
        ax.set_xlabel(r"Electron Lorentz factor $\gamma_e$")
    for ax in axes_flat[::n_cols]:
        ax.set_ylabel(r"$dN_e/d\gamma_e$")
    axes_flat[0].legend(handles=solver_handles, fontsize=8.0, loc="lower left", title="Solver")
    if magnetic_decay:
        fig.suptitle(
            rf"ASGARD Electron-Spectrum Evolution With Decaying Magnetic Field "
            rf"($\alpha_t={DECAY_ALPHA_T:.1f}$)",
            y=0.995,
        )
    else:
        fig.suptitle("ASGARD Electron-Spectrum Evolution", y=0.995)
    fig.tight_layout()

    filename = "compare_electron_spectrum_magnetic_decay.png" if magnetic_decay else "compare_electron_spectrum.png"
    return _save(fig, OUTPUT_DIR / filename)


def _build_magnetic_decay_compare() -> Path:
    def magnetic_decay_pair(electron_solver: str, num_chi: int | None = None) -> tuple[Model, Model]:
        base = _build_asgard_model(
            include_ssc=True,
            electron_solver=electron_solver,
            num_nu=MAGNETIC_DECAY_NUM_NU,
            num_gam_e=MAGNETIC_DECAY_NUM_GAM_E,
            num_chi=num_chi,
        )
        decay = _build_asgard_model(
            include_ssc=True,
            electron_solver=electron_solver,
            num_nu=MAGNETIC_DECAY_NUM_NU,
            num_gam_e=MAGNETIC_DECAY_NUM_GAM_E,
            num_chi=num_chi,
            magnetic_decay_alpha_t=DECAY_ALPHA_T,
            magnetic_decay_t0_s=DECAY_T0_S,
            epsilon_b_floor=STRONG_IC_EPS_B * DECAY_EPSB_FLOOR_FACTOR,
        )
        return base, decay

    model_base, model_decay = magnetic_decay_pair("fullhide_2d", MAGNETIC_DECAY_NUM_CHI)
    model_base_1d, model_decay_1d = magnetic_decay_pair("fullhide_1d")
    times = BASIC_TIMES
    bands = MAGNETIC_DECAY_BANDS
    spectrum_freqs = MAGNETIC_DECAY_FREQS
    spectrum_epochs = MAGNETIC_DECAY_EPOCHS
    details_base = model_base.details(1.0e2, 1.0e8)
    details_decay = model_decay.details(1.0e2, 1.0e8)
    lc_base_fnu = np.asarray(model_base.flux_density_grid(times, bands).total, dtype=float)
    lc_decay_fnu = np.asarray(model_decay.flux_density_grid(times, bands).total, dtype=float)
    lc_base = bands[:, None] * lc_base_fnu
    lc_decay = bands[:, None] * lc_decay_fnu

    fig, axes = plt.subplots(2, 2, figsize=(12.6, 7.0), dpi=200, sharex="col")
    colors = plt.cm.tab10(np.linspace(0, 1, bands.size))
    order_desc = np.argsort(bands)[::-1]
    lc_peak = 0.0
    for rank, i in enumerate(order_desc):
        nu = bands[i]
        shift = 10.0 ** (-MAGNETIC_DECAY_LC_SHIFT_DEX * rank)
        label = _label(nu, "Hz")
        shift_tag = f" x1e-{MAGNETIC_DECAY_LC_SHIFT_DEX * rank:.0f}" if rank > 0 else ""
        base_shifted = lc_base[i, :] * shift
        decay_shifted = lc_decay[i, :] * shift
        lc_peak = max(lc_peak, float(np.nanmax(base_shifted)), float(np.nanmax(decay_shifted)))
        axes[0, 0].loglog(times, base_shifted, color=colors[i], lw=1.8, alpha=ASGARD_ALPHA, label=f"baseline {label}{shift_tag}")
        axes[0, 0].loglog(times, decay_shifted, color=colors[i], lw=1.4, ls="--", alpha=0.85, label=f"decay {label}{shift_tag}")
        band_floor = max(1.0e-99, float(np.nanmax(lc_base[i, :])) * 1.0e-4)
        ratio = np.divide(
            lc_decay[i, :],
            lc_base[i, :],
            out=np.full_like(times, np.nan, dtype=float),
            where=lc_base[i, :] > band_floor,
        )
        mask = np.isfinite(ratio) & (ratio > 0.0)
        if np.any(mask):
            axes[1, 0].semilogx(times[mask], ratio[mask], color=colors[i], lw=1.5, label=label)
    axes[0, 0].set_ylabel(r"Energy Flux (erg/cm$^2$/s)")
    if np.isfinite(lc_peak) and lc_peak > 0.0:
        axes[0, 0].set_ylim(bottom=lc_peak * 1.0e-20)
    axes[0, 0].set_title("ASGARD 2D Magnetic-Decay Light Curves")
    axes[0, 0].grid(**GRID_STYLE)
    axes[0, 0].legend(fontsize=5.8, ncol=2)
    axes[1, 0].axhline(1.0, color="k", ls=":", lw=1.0)
    axes[1, 0].set_xlabel("Time [s]")
    axes[1, 0].set_ylabel("decay / baseline")
    axes[1, 0].grid(**GRID_STYLE)
    axes[1, 0].legend(fontsize=7, ncol=1)

    time_base = np.asarray(details_base.fwd.t_obs, dtype=float)
    time_decay = np.asarray(details_decay.fwd.t_obs, dtype=float)
    t_min_plot = float(np.min(times))
    t_max_plot = float(np.max(times))
    freq_sets = [
        (np.asarray(details_base.fwd.nu_m, dtype=float), np.asarray(details_decay.fwd.nu_m, dtype=float), r"$\nu_m$", "C0"),
        (np.asarray(details_base.fwd.nu_c, dtype=float), np.asarray(details_decay.fwd.nu_c, dtype=float), r"$\nu_c$", "C1"),
        (np.asarray(details_base.fwd.nu_a, dtype=float), np.asarray(details_decay.fwd.nu_a, dtype=float), r"$\nu_a$", "C2"),
    ]
    for base_arr, decay_arr, label, color in freq_sets:
        for arr, tt, ls, alpha, prefix in [
            (base_arr, time_base, "-", ASGARD_ALPHA, "baseline"),
            (decay_arr, time_decay, "--", 0.85, "decay"),
        ]:
            mask = (
                np.isfinite(tt)
                & np.isfinite(arr)
                & (tt > 0.0)
                & (arr > 0.0)
                & (tt >= t_min_plot)
                & (tt <= t_max_plot)
            )
            if np.any(mask):
                axes[0, 1].loglog(tt[mask], arr[mask], color=color, lw=1.5, ls=ls, alpha=alpha, label=f"{prefix} {label}")
        base_mask = np.isfinite(time_base) & np.isfinite(base_arr) & (time_base >= t_min_plot) & (time_base <= t_max_plot)
        decay_mask = np.isfinite(time_decay) & np.isfinite(decay_arr) & (time_decay >= t_min_plot) & (time_decay <= t_max_plot)
        if np.any(base_mask) and np.any(decay_mask):
            time_base_clip = time_base[base_mask]
            base_clip = base_arr[base_mask]
            decay_interp = _safe_log_interp(time_base_clip, time_decay[decay_mask], decay_arr[decay_mask])
            ratio = np.divide(decay_interp, base_clip, out=np.full_like(base_clip, np.nan, dtype=float), where=base_clip > 0.0)
        else:
            time_base_clip = np.array([], dtype=float)
            ratio = np.array([], dtype=float)
        mask = np.isfinite(time_base_clip) & np.isfinite(ratio) & (time_base_clip > 0.0) & (ratio > 0.0)
        if np.any(mask):
            axes[1, 1].semilogx(time_base_clip[mask], ratio[mask], color=color, lw=1.5, label=label)
    axes[0, 1].set_ylabel("Frequency [Hz]")
    axes[0, 1].set_title(
        rf"$\alpha_t={DECAY_ALPHA_T:.1f},\ t_0'={DECAY_T0_S:.1e}\,$s, "
        rf"$\epsilon_{{B,\rm floor}}/\epsilon_{{B,+}}={DECAY_EPSB_FLOOR_FACTOR:.1e}$"
    )
    axes[0, 1].grid(**GRID_STYLE)
    axes[0, 1].legend(fontsize=6.5, ncol=2)
    axes[1, 1].axhline(1.0, color="k", ls=":", lw=1.0)
    axes[1, 1].set_xlabel("Time [s]")
    axes[1, 1].set_ylabel("decay / baseline")
    axes[1, 1].grid(**GRID_STYLE)
    axes[1, 1].legend(fontsize=7)
    plt.tight_layout()
    lc_path = _save(fig, OUTPUT_DIR / "compare_magnetic_decay_2d.png")

    spectrum_sets = []
    for label, base_model, decay_model in [
        ("1D fullhide", model_base_1d, model_decay_1d),
        ("2D fullhide", model_base, model_decay),
    ]:
        spec_base = np.asarray(base_model.flux_density_grid(spectrum_epochs, spectrum_freqs).total, dtype=float) * spectrum_freqs[:, None]
        spec_decay = np.asarray(decay_model.flux_density_grid(spectrum_epochs, spectrum_freqs).total, dtype=float) * spectrum_freqs[:, None]
        spectrum_sets.append((label, spec_base, spec_decay))

    fig, axes = plt.subplots(
        2,
        len(spectrum_sets),
        figsize=(13.6, 7.0),
        dpi=200,
        sharex=True,
        height_ratios=[3.0, 1.25],
    )
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, spectrum_epochs.size))
    energy_ev = spectrum_freqs / constants.para_ev2hz
    spec_peak = float(np.nanmax([np.nanmax(arr) for _, base_arr, decay_arr in spectrum_sets for arr in (base_arr, decay_arr)]))
    for i_col, (solver_label, spec_base, spec_decay) in enumerate(spectrum_sets):
        ax_spec = axes[0, i_col]
        ax_ratio = axes[1, i_col]
        for i_epoch, (color, t_obs) in enumerate(zip(colors, spectrum_epochs)):
            label = fr"$t={t_obs:.0e}\,$s"
            ax_spec.loglog(energy_ev, spec_base[:, i_epoch], color=color, lw=1.7, alpha=ASGARD_ALPHA, label=f"baseline {label}")
            ax_spec.loglog(energy_ev, spec_decay[:, i_epoch], color=color, lw=1.4, ls="--", alpha=0.88, label=f"decay {label}")
            base_floor = max(1.0e-99, float(np.nanmax(spec_base[:, i_epoch])) * 1.0e-3)
            decay_floor = max(1.0e-99, float(np.nanmax(spec_decay[:, i_epoch])) * 1.0e-3)
            valid = (spec_base[:, i_epoch] > base_floor) & (spec_decay[:, i_epoch] > decay_floor)
            if np.any(valid):
                peak_idx = int(np.argmax(spec_base[:, i_epoch]))
                idx_valid = np.where(valid)[0]
                split_points = np.where(np.diff(idx_valid) > 1)[0]
                starts = np.concatenate(([0], split_points + 1))
                ends = np.concatenate((split_points + 1, [idx_valid.size]))
                keep = np.zeros_like(valid)
                for i_seg in range(starts.size):
                    seg_idx = idx_valid[starts[i_seg]:ends[i_seg]]
                    if seg_idx.size == 0:
                        continue
                    if seg_idx[0] <= peak_idx <= seg_idx[-1]:
                        keep[seg_idx] = True
                        break
                if not np.any(keep):
                    keep[valid] = True
            else:
                keep = valid
            ratio = np.divide(
                spec_decay[:, i_epoch],
                spec_base[:, i_epoch],
                out=np.full_like(spectrum_freqs, np.nan, dtype=float),
                where=keep,
            )
            mask = np.isfinite(ratio) & (ratio > 0.0)
            if np.any(mask):
                ax_ratio.semilogx(energy_ev[mask], ratio[mask], color=color, lw=1.4, label=label)
        ax_spec.set_title(f"{solver_label}, {MAGNETIC_DECAY_NUM_NU} frequency bins")
        if np.isfinite(spec_peak) and spec_peak > 0.0:
            ax_spec.set_ylim(bottom=spec_peak * 1.0e-15)
        ax_spec.grid(**GRID_STYLE)
        ax_spec.legend(fontsize=5.8, ncol=2)
        ax_ratio.axhline(1.0, color="k", ls=":", lw=1.0)
        ax_ratio.set_xlabel("Photon energy [eV]")
        ax_ratio.grid(**GRID_STYLE)
        ax_ratio.legend(fontsize=6.5, ncol=2)
        ax_spec.set_xlim(energy_ev[0], energy_ev[-1])
    axes[0, 0].set_ylabel(r"$\nu F_\nu$ (erg/cm$^2$/s)")
    axes[1, 0].set_ylabel("decay / baseline")
    fig.suptitle("ASGARD Magnetic-Decay Broadband Spectra", y=0.995)
    plt.tight_layout()
    _save(fig, OUTPUT_DIR / "compare_magnetic_decay_2d_broadband_spectrum.png")
    return lc_path


def _benchmark_builders(spectrum_quantity: str) -> list[tuple[str, Callable[[], Path]]]:
    return [
        ("basic_lc_spec", _build_basic_lc_spec),
        ("basic_bolometric", _build_basic_bolometric),
        ("reverse_shock_lc", _build_reverse_shock_lc),
        ("reverse_shock_thermal", _build_reverse_shock_thermal_benchmark),
        ("ssc_lc", _build_ssc_lc),
        ("spectrum_compare", lambda: _build_spectrum_compare(quantity=spectrum_quantity)),
        ("solver_sed_compare", _build_solver_sed_compare),
        ("electron_spectrum", _build_electron_spectrum_compare),
        ("magnetic_decay_2d", _build_magnetic_decay_compare),
        ("shock_quantities", _build_shock_quantities),
        ("photon_quantities", _build_photon_quantities),
        ("sky_image_single", _build_sky_single),
        ("sky_image_offaxis", _build_sky_offaxis),
        ("sky_image_flux_comparison", _build_sky_flux_comparison),
        ("sky_image_centroid", _build_sky_centroid_comparison),
        ("introspection_jet", _build_introspection_jet),
        ("introspection_medium", _build_introspection_medium),
        ("introspection_twocomp", _build_introspection_twocomp),
        ("speed_compare", _build_speed_compare),
    ]


def main(*, spectrum_quantity: str = "sed", scenario: str = "baseline", only: tuple[str, ...] | None = None) -> None:
    _set_benchmark_context(scenario)
    params = _benchmark_scenario()
    print(
        "[scenario] "
        f"{scenario}: E_iso={params.e_iso:.3e}, n_ism={params.n_ism:.3e}, "
        f"eps_e={params.eps_e:.3e}, eps_B={params.eps_b:.3e}, "
        f"xi_N={params.xi_n:.3e}, epsilon_p={params.epsilon_p:.3e}"
    )
    builders = _benchmark_builders(spectrum_quantity)
    if only is not None:
        builder_map = dict(builders)
        unknown = [name for name in only if name not in builder_map]
        if unknown:
            known = ", ".join(name for name, _ in builders)
            raise ValueError(f"unknown benchmark builder(s): {', '.join(unknown)}; expected one of: {known}")
        builders = [(name, builder_map[name]) for name in only]

    generated: list[Path] = []
    for i, (name, builder) in enumerate(builders, start=1):
        print(f"[ {i}/{len(builders)} ] build {name}")
        generated.append(builder())

    print(f"output: {OUTPUT_DIR}")
    for path in generated:
        print(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate ASGARD/VegasAfterglow comparison figures.")
    parser.add_argument(
        "--scenario",
        choices=sorted(BENCHMARK_SCENARIOS),
        default="baseline",
        help="Benchmark parameter set.",
    )
    parser.add_argument(
        "--spectrum-quantity",
        choices=("sed", "flux_density"),
        default="sed",
        help="Quantity used in the three-solver spectrum comparison.",
    )
    parser.add_argument(
        "--only",
        nargs="+",
        default=None,
        help="Run selected figure builders only, e.g. reverse_shock_lc reverse_shock_thermal.",
    )
    args = parser.parse_args()
    main(
        spectrum_quantity=args.spectrum_quantity,
        scenario=args.scenario,
        only=None if args.only is None else tuple(args.only),
    )
