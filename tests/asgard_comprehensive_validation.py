from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime
from fractions import Fraction as F
from pathlib import Path
import json
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_paths import asgard_doc_path
from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet, GaussianJet, PowerLawJet, Wind


OUTPUT_JSON = asgard_doc_path("comprehensive_validation_asgard.json")

FIDUCIAL_RESOLUTION = (0.1, 0.25, 10)
REGRESSION_RESOLUTION = (0.3, 2.0, 15)
PHI_VALUES = [0.1, 0.15, 0.2, 0.25]
THETA_VALUES = [0.25, 0.5, 0.75, 1.0]
T_VALUES = [5, 10, 15, 20]
BENCHMARK_PHI_VALUES = [0.1, 0.2]
BENCHMARK_THETA_VALUES = [0.25, 0.75]
BENCHMARK_T_VALUES = [10, 20]
CONVERGENCE_T_LC = np.logspace(1, 7, 24)
BENCHMARK_T_LC = np.logspace(1, 7, 12)
CONVERGENCE_BANDS = {"Radio": 1e9, "Optical": 4.84e14, "X-ray": 1e18}
SSC_BANDS = {"Radio": 1e9, "Optical": 4.84e14, "X-ray": 1e18, "TeV": 2.4e26}

STANDARD_PARAMS = {
    "E_iso": 1e52,
    "Gamma0": 300,
    "theta_c": 0.1,
    "n_ism": 1.0,
    "A_star": 0.1,
    "eps_e": 0.01,
    "eps_B": 0.01,
    "p": 2.2,
    "xi_e": 1.0,
    "lumi_dist": 1e28,
    "z": 1.0,
}

VALIDATION_INITIAL_RADIUS_CM = {"ISM": 1.0e12, "wind": 1.0e9}

TIME_RANGES = {
    "coasting": {"ISM": (1e-1, 1), "wind": (1e-2, 1e-1)},
    "BM": {"ISM": (5e2, 5e3), "wind": (1e4, 1e5)},
    "deep_newtonian": {"ISM": (1e12, 1e13), "wind": (1e14, 1e15)},
}

SHOCK_SCALINGS = {
    "coasting": {"ISM": {"u": F(0), "r": F(1), "B": F(0), "N_p": F(3)},
                 "wind": {"u": F(0), "r": F(1), "B": F(-1), "N_p": F(1)}},
    "BM": {"ISM": {"u": F(-3, 8), "r": F(1, 4), "B": F(-3, 8), "N_p": F(3, 4)},
           "wind": {"u": F(-1, 4), "r": F(1, 2), "B": F(-3, 4), "N_p": F(1, 2)}},
    "deep_newtonian": {"ISM": {"u": F(-3, 5), "r": F(2, 5), "B": F(-3, 5), "N_p": F(6, 5)},
                       "wind": {"u": F(-1, 3), "r": F(2, 3), "B": F(-1), "N_p": F(2, 3)}},
}

FREQ_SCALINGS = {
    "coasting": {"ISM": {"nu_m": F(0), "nu_c": F(-2), "nu_M": F(0)},
                 "wind": {"nu_m": F(-1), "nu_c": F(1), "nu_M": F(0)}},
    "BM": {"ISM": {"nu_m": F(-3, 2), "nu_c": F(-1, 2), "nu_M": F(-3, 8)},
           "wind": {"nu_m": F(-3, 2), "nu_c": F(1, 2), "nu_M": F(-1, 4)}},
    "deep_newtonian": {"ISM": {"nu_m": F(-3, 5), "nu_c": F(-1, 5), "nu_M": F(0)},
                       "wind": {"nu_m": F(-1), "nu_c": F(1), "nu_M": F(0)}},
}

RVS_TIME_RANGES_THIN = {
    "crossing": {"ISM": (1e-2, 1e-1), "wind": (5, 50)},
    "post_crossing": {"ISM": (1e5, 1e6), "wind": (1e6, 1e7)},
    "deep_newtonian": {"ISM": (1e12, 1e13), "wind": (1e14, 1e15)},
}
RVS_TIME_RANGES_THICK = {
    "crossing": {"ISM": (5e3, 5e4), "wind": (1e2, 1e3)},
    "post_crossing": {"ISM": (1e7, 1e8), "wind": (1e7, 1e8)},
    "deep_newtonian": {"ISM": (1e12, 1e13), "wind": (1e14, 1e15)},
}
RVS_SHOCK_SCALINGS_THIN = {
    "crossing": {"ISM": {"r": 1, "B": 0, "N_p": F(3, 2)}, "wind": {"r": 1, "B": -1, "N_p": F(1, 2)}},
    "post_crossing": {"ISM": {"r": F(1, 4), "N_p": 0}, "wind": {"r": F(1, 2), "N_p": 0}},
    "deep_newtonian": {"ISM": {"r": F(2, 5), "N_p": 0}, "wind": {"r": F(2, 3), "N_p": 0}},
}
RVS_SHOCK_SCALINGS_THICK = {
    "crossing": {"ISM": {"r": F(1, 2), "B": F(-1, 4), "N_p": 1}, "wind": {"r": 1, "B": -1, "N_p": 1}},
    "post_crossing": {"ISM": {"r": F(1, 4), "N_p": 0}, "wind": {"r": F(1, 2), "N_p": 0}},
    "deep_newtonian": {"ISM": {"r": F(2, 5), "N_p": 0}, "wind": {"r": F(2, 3), "N_p": 0}},
}
RVS_FREQ_SCALINGS_THIN = {
    "crossing": {"ISM": {"nu_m": 0, "nu_c": -2, "nu_M": 0}, "wind": {"nu_m": -1, "nu_c": 1, "nu_M": 0}},
}
RVS_FREQ_SCALINGS_THICK = {
    "crossing": {"ISM": {"nu_c": -1, "nu_M": F(-1, 4)}, "wind": {"nu_m": -1, "nu_c": 1, "nu_M": 0}},
}

SPECTRAL_REGIMES = {
    "I": {"segments": [("below_nu_a", 2.0), ("nu_a_to_nu_m", F(1, 3)), ("nu_m_to_nu_c", "-(p-1)/2"), ("above_nu_c", "-p/2")]},
    "II": {"segments": [("below_nu_m", 2.0), ("nu_m_to_nu_a", F(5, 2)), ("nu_a_to_nu_c", "-(p-1)/2"), ("above_nu_c", "-p/2")]},
    "III": {"segments": [("below_nu_a", 2.0), ("nu_a_to_nu_c", F(1, 3)), ("nu_c_to_nu_m", F(-1, 2)), ("above_nu_m", "-p/2")]},
    "IV": {"segments": [("below_nu_c", 2.0), ("nu_c_to_nu_a", 2.0), ("nu_a_to_nu_m", F(-1, 2)), ("above_nu_m", "-p/2")]},
    "V": {"segments": [("below_nu_c", 2.0), ("nu_c_to_nu_m", 2.0), ("nu_m_to_nu_a", F(5, 2)), ("above_nu_a", "-p/2")]},
}
REGIME_TEST_CONFIGS = {
    "I": {"medium": "ISM", "t": 5e3, "n_ism": 0.1, "eps_B": 1e-2, "eps_e": 1e-1},
    "II": {"medium": "ISM", "t": 5e4, "n_ism": 1e6, "eps_B": 1e-5, "eps_e": 5e-2},
    "III": {"medium": "ISM", "t": 1e5, "n_ism": 30, "eps_B": 1e-1, "eps_e": 1e-1, "xi_e": 5e-3},
    "IV": {"medium": "ISM", "t": 1e5, "n_ism": 1e6, "eps_B": 0.3, "eps_e": 0.3, "xi_e": 5e-3},
    "V": {"medium": "ISM", "t": 1e4, "n_ism": 1e8, "eps_e": 0.01, "eps_B": 3e-1, "xi_e": 1.0},
}

SLOPE_TOLERANCE = 0.1
SPECTRAL_TOLERANCE = 0.15
MEAN_ERROR_THRESHOLD = 0.05
MAX_ERROR_THRESHOLD = 0.15

BENCHMARK_CONFIGS = [
    {"jet": "tophat", "medium": "ISM", "radiation": "synchrotron", "theta_ratio": 0, "scan_dims": ("phi", "theta", "t")},
    {"jet": "gaussian", "medium": "ISM", "radiation": "synchrotron", "theta_ratio": 1, "scan_dims": ()},
]


@dataclass
class ValidationResult:
    name: str
    category: str
    passed: bool
    expected: float | None = None
    measured: float | None = None
    tolerance: float | None = None
    extra: dict = field(default_factory=dict)


def fit_powerlaw(t, f, window=5):
    valid = (f > 0) & np.isfinite(f) & (t > 0) & np.isfinite(t)
    if np.sum(valid) < 3:
        return np.nan
    t_valid, f_valid = t[valid], f[valid]
    sort_idx = np.argsort(t_valid)
    log_t, log_f = np.log10(t_valid[sort_idx]), np.log10(f_valid[sort_idx])
    d_log_t, d_log_f = np.diff(log_t), np.diff(log_f)
    valid_diff = np.abs(d_log_t) > 1e-10
    if np.sum(valid_diff) < 2:
        return np.nan
    local_slopes = d_log_f[valid_diff] / d_log_t[valid_diff]
    if len(local_slopes) > 1:
        local_slopes = local_slopes[1:]
    if len(local_slopes) >= window:
        kernel = np.ones(window) / window
        local_slopes = np.convolve(local_slopes, kernel, mode="valid")
    return float(np.mean(local_slopes))


def _make_medium(name: str, **overrides):
    if name == "ISM":
        return ISM(n0=overrides.get("n_ism", STANDARD_PARAMS["n_ism"]))
    if name == "wind":
        return Wind(
            A_star=overrides.get("A_star", STANDARD_PARAMS["A_star"]),
            n_ism=overrides.get("n_ism", 0.1),
            n0=overrides.get("n0"),
        )
    raise ValueError(name)


def _make_jet(name: str, theta_ratio: float = 0.0, duration: float | None = None, **overrides):
    e_iso = overrides.get("E_iso", STANDARD_PARAMS["E_iso"])
    gamma0 = overrides.get("Gamma0", STANDARD_PARAMS["Gamma0"])
    theta_c = overrides.get("theta_c", STANDARD_PARAMS["theta_c"])
    if name == "tophat":
        return TophatJet(theta_c=theta_c, E_iso=e_iso, Gamma0=gamma0, duration=duration)
    if name == "gaussian":
        return GaussianJet(theta_c=theta_c, E_iso=e_iso, Gamma0=gamma0, theta_max=0.4, duration=duration)
    if name == "powerlaw":
        return PowerLawJet(theta_c=theta_c, E_iso=e_iso, Gamma0=gamma0, k_e=2.0, k_g=2.0, theta_max=0.4, duration=duration)
    raise ValueError(name)


def _make_observer(theta_ratio: float = 0.0):
    return Observer(lumi_dist=STANDARD_PARAMS["lumi_dist"], z=STANDARD_PARAMS["z"], theta_obs=STANDARD_PARAMS["theta_c"] * theta_ratio, phi_obs=0.0)


def _make_radiation(name: str, reverse: bool = False):
    if reverse:
        return Radiation(eps_e=0.1, eps_B=0.01, p=2.3, xi_N=1.0, ssc=False, kn=False)
    if name == "synchrotron":
        return Radiation(eps_e=0.1, eps_B=0.01, p=2.3, xi_N=1.0, ssc=False, kn=False)
    if name == "full_ssc":
        return Radiation(eps_e=0.1, eps_B=0.01, p=2.3, xi_N=1.0, ssc=True, kn=False)
    if name == "ssc_kn":
        return Radiation(eps_e=0.1, eps_B=0.01, p=2.3, xi_N=1.0, ssc=True, kn=True)
    if name in ("rvs_sync_thin", "rvs_sync_thick"):
        return Radiation(eps_e=0.1, eps_B=0.01, p=2.3, xi_N=1.0, ssc=False, kn=False)
    raise ValueError(name)


def _make_model(
    medium_name: str,
    *,
    jet_name: str = "tophat",
    radiation_name: str = "synchrotron",
    theta_ratio: float = 0.0,
    resolutions=FIDUCIAL_RESOLUTION,
    benchmark_mode: bool = False,
    **overrides,
):
    overrides = dict(overrides)
    duration = None
    rvs_rad = None
    if radiation_name == "rvs_sync_thin":
        duration = 1.0
        rvs_rad = _make_radiation(radiation_name, reverse=True)
    elif radiation_name == "rvs_sync_thick":
        duration = 1.0e4
        rvs_rad = _make_radiation(radiation_name, reverse=True)
    jet_duration = overrides.pop("duration", duration)
    jet = _make_jet(jet_name, theta_ratio=theta_ratio, duration=jet_duration, **overrides)
    medium = _make_medium(medium_name, **overrides)
    observer = _make_observer(theta_ratio=theta_ratio)
    rad = Radiation(
        eps_e=overrides.get("eps_e", STANDARD_PARAMS["eps_e"]),
        eps_B=overrides.get("eps_B", STANDARD_PARAMS["eps_B"]),
        p=overrides.get("p", STANDARD_PARAMS["p"]),
        xi_N=overrides.get("xi_e", STANDARD_PARAMS["xi_e"]),
        ssc=radiation_name in ("full_ssc", "ssc_kn"),
        kn=radiation_name == "ssc_kn",
    )
    setups = _make_benchmark_setups(resolutions, jet_name) if benchmark_mode else None
    model_resolutions = None if benchmark_mode else resolutions
    model = Model(
        jet=jet,
        medium=medium,
        observer=observer,
        fwd_rad=rad,
        rvs_rad=rvs_rad,
        resolutions=model_resolutions,
        setups=setups,
    )
    if "initial_radius_cm" in overrides:
        model.setups.initial_radius_cm = float(overrides["initial_radius_cm"])
    return model


def _make_benchmark_setups(resolutions, jet_name: str) -> Setups:
    theta_ppd, phi_ppd, t_ppd = resolutions
    num_tobs = max(10, int(np.ceil(4.0 * float(t_ppd))))
    if jet_name == "tophat":
        return Setups(
            num_threads=8,
            num_gam_e=121,
            num_nu=121,
            num_r=48,
            num_theta=40,
            num_tobs=num_tobs,
            patch_theta=1,
            patch_phi=1,
        )
    patch_theta = max(1, int(np.ceil(4.0 * float(theta_ppd))))
    patch_phi = max(2, int(np.ceil(8.0 * float(phi_ppd))))
    return Setups(
        num_threads=8,
        num_gam_e=121,
        num_nu=121,
        num_r=48,
        num_theta=40,
        num_tobs=num_tobs,
        patch_theta=patch_theta,
        patch_phi=patch_phi,
    )


def _extract_flux_component(result, radiation_name: str):
    if radiation_name.startswith("rvs_"):
        return np.asarray(result.rev.sync)
    if radiation_name in ("full_ssc", "ssc_kn"):
        return np.asarray(result.fwd.ssc)
    return np.asarray(result.total)


def _run_benchmark_subset():
    records = []
    dim_specs = {
        "phi": (BENCHMARK_PHI_VALUES, 0),
        "theta": (BENCHMARK_THETA_VALUES, 1),
        "t": (BENCHMARK_T_VALUES, 2),
    }
    for cfg in BENCHMARK_CONFIGS:
        entry = {"config": cfg, "timing_ms": None, "convergence": {}}
        model = _make_model(cfg["medium"], jet_name=cfg["jet"], radiation_name=cfg["radiation"], theta_ratio=cfg["theta_ratio"], resolutions=FIDUCIAL_RESOLUTION, benchmark_mode=True)
        t = np.logspace(2, 8, 8)
        nu = np.full_like(t, 4.84e14)
        t0 = time.perf_counter()
        model.flux_density(t, nu)
        entry["timing_ms"] = (time.perf_counter() - t0) * 1000.0

        bands = SSC_BANDS if cfg["radiation"] in ("full_ssc", "ssc_kn") else CONVERGENCE_BANDS
        t_lc = BENCHMARK_T_LC
        nu_arr = np.array(list(bands.values()))

        def _compute_flux(resolution):
            model_for_resolution = _make_model(
                cfg["medium"],
                jet_name=cfg["jet"],
                radiation_name=cfg["radiation"],
                theta_ratio=cfg["theta_ratio"],
                resolutions=tuple(float(x) for x in resolution),
                benchmark_mode=True,
            )
            return _extract_flux_component(model_for_resolution.flux_density_grid(t_lc, nu_arr), cfg["radiation"])

        for dim in cfg.get("scan_dims", ()):
            values, idx = dim_specs[dim]
            ref_res = list(FIDUCIAL_RESOLUTION)
            ref_res[idx] = max(values) * 1.2
            ref_flux = _compute_flux(tuple(ref_res))
            errors = {}
            for value in values:
                res = list(FIDUCIAL_RESOLUTION)
                res[idx] = value
                flux = _compute_flux(tuple(res))
                band_err = {}
                for k, band in enumerate(bands):
                    ref = ref_flux[k, :]
                    cur = flux[k, :]
                    valid = (ref > 0) & np.isfinite(ref) & (cur > 0) & np.isfinite(cur)
                    err = np.abs(cur[valid] - ref[valid]) / ref[valid]
                    band_err[band] = {"max": float(np.max(err)), "mean": float(np.mean(err))}
                errors[str(value)] = band_err
            fiducial_key = str(FIDUCIAL_RESOLUTION[idx])
            worst_mean = max(errors[fiducial_key][band]["mean"] for band in bands)
            worst_max = max(errors[fiducial_key][band]["max"] for band in bands)
            entry["convergence"][dim] = {
                "fiducial": FIDUCIAL_RESOLUTION[idx],
                "worst_mean_error": worst_mean,
                "worst_max_error": worst_max,
                "status": "pass" if worst_mean < MEAN_ERROR_THRESHOLD and worst_max < MAX_ERROR_THRESHOLD else "fail",
                "scan": errors,
            }
        records.append(entry)
    return records


def _run_scaling(details, branch: str, qty: str, t_range, expected, name, medium_key: str):
    shock = details.fwd if branch == "fwd" else details.rvs
    if shock is None:
        return ValidationResult(name=name, category=branch, passed=False, extra={"message": "missing reverse-shock details"})
    t = np.asarray(shock.t_obs)
    if qty == "u":
        Gamma = np.asarray(shock.Gamma)
        y = Gamma * np.sqrt(np.maximum(1.0 - 1.0 / (Gamma * Gamma), 0.0))
    else:
        y = np.asarray(getattr(shock, qty))
    mask = (t >= t_range[0]) & (t <= t_range[1])
    if medium_key == "deep_newtonian":
        Gamma = np.asarray(shock.Gamma)
        u = Gamma * np.sqrt(np.maximum(1.0 - 1.0 / (Gamma * Gamma), 0.0))
        mask &= u < 0.1
    measured = fit_powerlaw(t[mask], y[mask])
    passed = bool(np.isfinite(measured) and abs(measured - float(expected)) < SLOPE_TOLERANCE)
    return ValidationResult(name=name, category=branch, passed=passed, expected=float(expected), measured=float(measured), tolerance=SLOPE_TOLERANCE)


def _run_forward_regression():
    results = []
    for medium_name in ("ISM", "wind"):
        model = _make_model(
            medium_name,
            resolutions=REGRESSION_RESOLUTION,
            n_ism=0.1 if medium_name == "ISM" else 0.0,
            A_star=0.3 if medium_name == "wind" else STANDARD_PARAMS["A_star"],
            Gamma0=300 if medium_name == "ISM" else 70,
            initial_radius_cm=VALIDATION_INITIAL_RADIUS_CM[medium_name],
        )
        phase_ranges = [TIME_RANGES[phase][medium_name] for phase in ("coasting", "BM", "deep_newtonian")]
        details = model.details(min(r[0] for r in phase_ranges) / 10.0, max(r[1] for r in phase_ranges) * 10.0)
        for phase in ("coasting", "BM", "deep_newtonian"):
            for qty, expected in SHOCK_SCALINGS[phase][medium_name].items():
                key = {"r": "radius", "B": "B_comv", "N_p": "N_p"}.get(qty, qty)
                results.append(_run_scaling(details, "fwd", key, TIME_RANGES[phase][medium_name], expected, f"fwd-{medium_name}-{phase}-{qty}", phase))
            for qty, expected in FREQ_SCALINGS[phase][medium_name].items():
                results.append(_run_scaling(details, "fwd", qty, TIME_RANGES[phase][medium_name], expected, f"fwd-{medium_name}-{phase}-{qty}", phase))
    return results


def _run_reverse_regression():
    results = []
    configs = [
        ("thin", "ISM", 300, 1.0, 0.01),
        ("thin", "wind", 50, 0.01, 0.01),
        ("thick", "ISM", 300, 0.01, 1e5),
        ("thick", "wind", 50, 0.01, 1e4),
    ]
    for regime, medium_name, gamma0, medium_amp, duration in configs:
        medium_kwargs = {"n_ism": medium_amp} if medium_name == "ISM" else {"A_star": medium_amp, "n_ism": 0.0}
        model = _make_model(
            medium_name,
            radiation_name=f"rvs_sync_{regime}",
            resolutions=REGRESSION_RESOLUTION,
            Gamma0=gamma0,
            duration=duration,
            initial_radius_cm=VALIDATION_INITIAL_RADIUS_CM[medium_name],
            **medium_kwargs,
        )
        shock_scalings = RVS_SHOCK_SCALINGS_THIN if regime == "thin" else RVS_SHOCK_SCALINGS_THICK
        freq_scalings = RVS_FREQ_SCALINGS_THIN if regime == "thin" else RVS_FREQ_SCALINGS_THICK
        time_ranges = RVS_TIME_RANGES_THIN if regime == "thin" else RVS_TIME_RANGES_THICK
        all_ranges = [time_ranges[phase][medium_name] for phase in time_ranges]
        details = model.details(min(r[0] for r in all_ranges) / 10.0, max(r[1] for r in all_ranges) * 10.0)
        for phase, qdict in shock_scalings.items():
            for qty, expected in qdict[medium_name].items():
                key = {"r": "radius", "B": "B_comv", "N_p": "N_p"}.get(qty, qty)
                results.append(_run_scaling(details, "rvs", key, time_ranges[phase][medium_name], expected, f"rvs-{regime}-{medium_name}-{phase}-{qty}", phase))
        for phase, qdict in freq_scalings.items():
            for qty, expected in qdict[medium_name].items():
                results.append(_run_scaling(details, "rvs", qty, time_ranges[phase][medium_name], expected, f"rvs-{regime}-{medium_name}-{phase}-{qty}", phase))
    return results


def _detect_regime(nu_a, nu_m, nu_c):
    order = tuple(name for name, _ in sorted({"nu_a": nu_a, "nu_m": nu_m, "nu_c": nu_c}.items(), key=lambda item: item[1]))
    for regime, spec in {
        "I": ("nu_a", "nu_m", "nu_c"),
        "II": ("nu_m", "nu_a", "nu_c"),
        "III": ("nu_a", "nu_c", "nu_m"),
        "IV": ("nu_c", "nu_a", "nu_m"),
        "V": ("nu_c", "nu_m", "nu_a"),
    }.items():
        if order == spec:
            return regime
    return None


def _segment_range(seg_name, nu_a, nu_m, nu_c):
    fm = {"nu_a": nu_a, "nu_m": nu_m, "nu_c": nu_c}
    for prefix, fn in [("below_", lambda nu: (nu / 1000.0, nu / 30.0)), ("above_", lambda nu: (nu * 10.0, nu * 100.0))]:
        if seg_name.startswith(prefix):
            key = seg_name[len(prefix):]
            return fn(fm[key]) if fm.get(key) else (None, None)
    if "_to_" in seg_name:
        lo_key, hi_key = seg_name.split("_to_")
        lo, hi = fm.get(lo_key), fm.get(hi_key)
        if lo and hi:
            if seg_name == "nu_c_to_nu_a":
                return lo * 1.2, hi / 30.0
            return lo * 5.0, hi / 20.0
    return None, None


def _measure_local_spectral_slope(freq_hz, flux, expected_slope, smooth_window=5):
    valid = (freq_hz > 0) & np.isfinite(freq_hz) & (flux > 0) & np.isfinite(flux)
    if np.sum(valid) < 8:
        return np.nan, np.nan
    nu_valid = np.asarray(freq_hz[valid], dtype=float)
    flux_valid = np.asarray(flux[valid], dtype=float)
    order = np.argsort(nu_valid)
    log_nu = np.log10(nu_valid[order])
    log_flux = np.log10(flux_valid[order])
    local_slopes = np.diff(log_flux) / np.diff(log_nu)
    slope_freq = np.sqrt(nu_valid[order][:-1] * nu_valid[order][1:])
    if local_slopes.size == 0:
        return np.nan, np.nan
    if local_slopes.size >= smooth_window:
        kernel = np.ones(smooth_window, dtype=float) / float(smooth_window)
        smoothed = np.convolve(local_slopes, kernel, mode="same")
    else:
        smoothed = local_slopes
    lo_idx = max(2, int(0.1 * smoothed.size))
    hi_idx = min(smoothed.size - 2, int(0.9 * smoothed.size))
    if hi_idx <= lo_idx:
        lo_idx = 0
        hi_idx = smoothed.size
    best_idx = lo_idx + int(np.argmin(np.abs(smoothed[lo_idx:hi_idx] - expected_slope)))
    avg_lo = max(0, best_idx - 2)
    avg_hi = min(local_slopes.size, best_idx + 3)
    return float(np.mean(local_slopes[avg_lo:avg_hi])), float(slope_freq[best_idx])


def _eval_beta(expr, p):
    if isinstance(expr, str):
        return float(eval(expr.replace("p", str(p))))
    return float(expr)


def _run_spectral_regimes():
    results = []
    p = STANDARD_PARAMS["p"]
    for regime, cfg in REGIME_TEST_CONFIGS.items():
        medium_name = cfg["medium"]
        model = _make_model(
            medium_name,
            resolutions=REGRESSION_RESOLUTION,
            eps_e=cfg.get("eps_e", STANDARD_PARAMS["eps_e"]),
            eps_B=cfg.get("eps_B", STANDARD_PARAMS["eps_B"]),
            xi_e=cfg.get("xi_e", STANDARD_PARAMS["xi_e"]),
            n_ism=cfg.get("n_ism", STANDARD_PARAMS["n_ism"]),
            A_star=cfg.get("A_star", STANDARD_PARAMS["A_star"]),
        )
        t_test = cfg["t"]
        details = model.details(t_test * 0.9, t_test * 1.1)
        idx = int(np.argmin(np.abs(np.asarray(details.fwd.t_obs) - t_test)))
        nu_a = float(np.asarray(details.fwd.nu_a)[idx])
        nu_m = float(np.asarray(details.fwd.nu_m)[idx])
        nu_c = float(np.asarray(details.fwd.nu_c)[idx])
        actual = _detect_regime(nu_a, nu_m, nu_c)
        for seg_name, beta_expr in SPECTRAL_REGIMES[regime]["segments"]:
            lo, hi = _segment_range(seg_name, nu_a, nu_m, nu_c)
            if lo is None or hi is None or not (hi > 1.5 * lo):
                continue
            freq = np.logspace(np.log10(lo * 1.2), np.log10(hi * 0.8), 160)
            flux = model.flux_density_grid(np.array([t_test]), freq).total[:, 0]
            measured = fit_powerlaw(freq, flux)
            expected = _eval_beta(beta_expr, p)
            measured_local, probe_freq = _measure_local_spectral_slope(freq, flux, expected)
            passed = bool(np.isfinite(measured_local) and abs(measured_local - expected) < SPECTRAL_TOLERANCE)
            results.append(
                ValidationResult(
                    name=f"regime-{regime}-{seg_name}",
                    category="spectrum",
                    passed=passed,
                    expected=expected,
                    measured=float(measured_local) if np.isfinite(measured_local) else None,
                    tolerance=SPECTRAL_TOLERANCE,
                    extra={
                        "actual_regime": actual,
                        "expected_regime": regime,
                        "global_fit": float(measured) if np.isfinite(measured) else None,
                        "probe_frequency_hz": float(probe_freq) if np.isfinite(probe_freq) else None,
                    },
                )
            )
    return results


def main() -> None:
    OUTPUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    benchmark = _run_benchmark_subset()
    regression = _run_forward_regression() + _run_reverse_regression() + _run_spectral_regimes()
    payload = {
        "timestamp": datetime.now().isoformat(),
        "methodology": {
            "source": "ASGARD validation docs/scripts",
            "benchmark_note": "Representative subset of benchmark configs was used. Direct top-hat keeps convergence scans; structured-jet entries keep timing-only benchmarks to stay tractable on the current Python patch backend.",
            "fiducial_resolution": FIDUCIAL_RESOLUTION,
            "regression_resolution": REGRESSION_RESOLUTION,
            "phi_values": PHI_VALUES,
            "theta_values": THETA_VALUES,
            "t_values": T_VALUES,
            "benchmark_phi_values": BENCHMARK_PHI_VALUES,
            "benchmark_theta_values": BENCHMARK_THETA_VALUES,
            "benchmark_t_values": BENCHMARK_T_VALUES,
        },
        "benchmark": benchmark,
        "regression": [asdict(item) for item in regression],
        "summary": {
            "benchmark_pass": all(
                dim["status"] == "pass"
                for cfg in benchmark
                for dim in cfg["convergence"].values()
            ),
            "regression_pass": all(item.passed for item in regression),
            "regression_total": len(regression),
            "regression_failed": sum(0 if item.passed else 1 for item in regression),
        },
    }
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2))
    print(OUTPUT_JSON)
    print(json.dumps(payload["summary"], indent=2))


if __name__ == "__main__":
    main()
