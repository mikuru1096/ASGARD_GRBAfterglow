#!/usr/bin/env python
from __future__ import annotations

import argparse
import contextlib
import csv
import gc
import io
import json
import math
import os
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
import shutil
import sys
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
_MPLCONFIGDIR = ROOT / "output" / ".matplotlib"
_MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPLCONFIGDIR))

import numpy as np
from astropy import units
from astropy.cosmology import FlatLambdaCDM


if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.asgard_config import FitConfig, ReverseShockConfig
from asgard_core.asgard_state import project_flux_grid, solve_state_from_setup
from asgard_core.asgard_setup import build_simulation_setup
from src import constants


def _load_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


TARGET_GRB = "180720B"
TARGET_GRB_LABEL = "GRB 180720B"
TARGET_GRB_SLUG = "grb180720b"
TARGET_REDSHIFT = 0.654
TARGET_R0_CM = 1.0e9
FIXED_A_STAR = 0.0
T90_S = 48.9

DATA_DIR = ROOT / "data_light_curve"
OUTPUT_ROOT = ROOT / "output"
OUTPUT_REL = Path("output")
METRICS_FILE = OUTPUT_ROOT / "grb180720b_model_comparison_metrics.csv"

KEV_TO_HZ = 2.417989242e17
GEV_TO_HZ = 1.0e6 * KEV_TO_HZ
TEV_TO_HZ = 1.0e9 * KEV_TO_HZ
MJY_CGS = 1.0e-26
AB_ZEROPOINT = 48.6
ANGSTROM_TO_KEV = 12.398419843320026
MODEL_T0_OFFSET_S = 0.0

MW_EBV = 0.037
MW_RV = 3.08
HOST_AV = 0.5
HOST_RV = 2.93
HOST_EBV = HOST_AV / HOST_RV
OPTICAL_R_WAVELENGTH_A = 6580.0
OPTICAL_R_LABEL = "Optical R"

XRT_WIDE_BAND_SAMPLES = 5
DEFAULT_WIDE_BAND_SAMPLES = 8
VHE_EBL_MODEL = "saldana21.txt"
VHE_EBL_TABLE = ROOT / "src" / "Radiation" / "EBL" / VHE_EBL_MODEL
XRT_NH_GAL_1E22 = 3.92e20 / 1.0e22
XRT_NH_HOST_1E22 = {
    "WT": 3.71e21 / 1.0e22,
    "PC": 3.4e21 / 1.0e22,
}

MODEL_TO_SOLVER = {
    "fullhide_1d": "fullhide_1d",
    "fullhide_1d_pic": "fullhide_1d_pic",
    "fullhide_2d_pic": "fullhide_2d_pic",
}

_LIKE_PROGRESS: dict[int, dict[str, Any]] = {}

PARAMETER_NAMES_PHYSICAL_BASE = [
    "log10_Ekiso",
    "log10_Gamma0",
    "pf",
    "log10_eps_e_f",
    "log10_eps_B_f",
    "log10_xi_e_f",
    "log10_theta_j",
    "pr",
    "log10_eps_e_r",
    "log10_eps_B_r",
    "log10_xi_e_r",
    "log10_dNe",
]
PARAMETER_NAMES_BASE = PARAMETER_NAMES_PHYSICAL_BASE + ["log10_f_sys"]
PARAMETER_NAMES_FULLHIDE = PARAMETER_NAMES_PHYSICAL_BASE + ["log10_Bw_B0", "log10_f_sys"]

ARTICLE_MEDIAN_PARAMS = {
    "log10_Ekiso": 55.56,
    "log10_Gamma0": 1.93,
    "pf": 2.49,
    "log10_eps_e_f": -2.07,
    "log10_eps_B_f": -6.21,
    "log10_xi_e_f": -1.0,
    "log10_theta_j": -1.04,
    "pr": 2.45,
    "log10_eps_e_r": -1.96,
    "log10_eps_B_r": -5.05,
    "log10_xi_e_r": -1.0,
    "log10_dNe": 0.0,
    "log10_Bw_B0": 2.0,
    "log10_f_sys": -1.0,
}

PRIOR_BASE = {
    "log10_Ekiso": (51.0, 56.0),
    "log10_Gamma0": (1.0, 2.5),
    "pf": (2.0, 3.0),
    "log10_eps_e_f": (-3.0, -0.1),
    "log10_eps_B_f": (-7.0, -0.1),
    "log10_xi_e_f": (-3.0, 0.0),
    "log10_theta_j": (-2.0, 0.0),
    "pr": (2.0, 3.0),
    "log10_eps_e_r": (-3.0, -0.1),
    "log10_eps_B_r": (-7.0, -0.1),
    "log10_xi_e_r": (-3.0, 0.0),
    "log10_dNe": (-3.0, 3.0),
    "log10_Bw_B0": (1.0, 3.0),
    "log10_f_sys": (-2.0, -0.3),
}

FIXED_NUMERICAL = {
    "Num_R": 100,
    "Num_theta": 300,
    "Num_phi": 1,
    "Num_gam_e": 101,
    "Num_nu": 200,
    "Num_tobs": 200,
    "index_Y": 2,
    "index_dyn": 3,
    "index_syn_integr": 2,
}

RADIO_SELECTIONS = (
    (1.2, 1.45, 1.3),
    (5.0, 5.79, 5.0),
    (13.5, 16.0, 15.5),
    (85.0, 93.0, 93.0),
)
OPTICAL_BAND_ALIASES = {"UVU": "U", "UVB": "B", "UVV": "V"}


@dataclass(frozen=True)
class ArticleData:
    kind: str
    file: str
    label: str
    time_s: float
    time_start_s: float
    time_stop_s: float
    model_time_s: float
    model_time_start_s: float
    model_time_stop_s: float
    exposure_s: float
    frequency_hz: float
    flux: float
    error: float
    upper_limit: bool = False
    energy_min_hz: float | None = None
    energy_max_hz: float | None = None
    mode: str = ""

    @property
    def valid_for_fit(self) -> bool:
        return (not self.upper_limit) and self.model_time_s > 0.0 and self.flux > 0.0 and np.isfinite(self.flux)


def parameter_names(model_name: str) -> list[str]:
    return PARAMETER_NAMES_FULLHIDE if model_name == "fullhide_2d_pic" else PARAMETER_NAMES_BASE


def prior_bounds(model_name: str) -> list[tuple[float, float]]:
    bounds = dict(PRIOR_BASE)
    if model_name in {"fullhide_1d_pic", "fullhide_2d_pic"}:
        bounds["log10_eps_e_f"] = (-4.0, -0.1)
        bounds["log10_eps_B_f"] = (-4.0, -2.3)
        bounds["log10_eps_e_r"] = (-4.0, -0.1)
        bounds["log10_eps_B_r"] = (-4.0, -2.3)
    return [bounds[name] for name in parameter_names(model_name)]


def prior(cube, ndim, nparams, model_name: str) -> None:
    for i, (lo, hi) in enumerate(prior_bounds(model_name)):
        cube[i] = cube[i] * (hi - lo) + lo


def params_from_vector(vector, model_name: str) -> dict[str, float]:
    return {name: float(value) for name, value in zip(parameter_names(model_name), vector)}


def parameter_vector_from_article(model_name: str) -> np.ndarray:
    return np.asarray([ARTICLE_MEDIAN_PARAMS[name] for name in parameter_names(model_name)], dtype=float)


def f_sys_from_params(params: dict[str, float], default: float = 0.1) -> float:
    return float(10.0 ** float(params.get("log10_f_sys", math.log10(default))))


def _invalid_physical(params: dict[str, float]) -> bool:
    if not all(np.isfinite(value) for value in params.values()):
        return True
    if params["log10_eps_B_f"] >= params["log10_eps_e_f"]:
        return True
    if params["log10_eps_B_r"] >= params["log10_eps_e_r"]:
        return True
    return False


def _runtime_home(base_dir: Path | None = None) -> Path:
    root = OUTPUT_ROOT if base_dir is None else Path(base_dir)
    home = root / ".runtime_home"
    home.mkdir(parents=True, exist_ok=True)
    os.environ["HOME"] = str(home)
    return home


def _prepare_astromodels(base_dir: Path | None = None) -> None:
    home = _runtime_home(base_dir)
    (home / ".astromodels" / "log").mkdir(parents=True, exist_ok=True)
    (home / ".config" / "astromodels").mkdir(parents=True, exist_ok=True)


@lru_cache(maxsize=4)
def _import_astromodels_functions(base_dir: Path | None = None):
    _prepare_astromodels(base_dir)
    with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
        from astromodels.functions import TbAbs, ZDust
    return ZDust, TbAbs


@lru_cache(maxsize=8)
def _cosmo_d_l(z: float) -> float:
    return FlatLambdaCDM(H0=67.8, Om0=0.308).luminosity_distance(z).to(units.cm).value


def _read_header_rows(path: Path) -> tuple[list[str], list[list[str]]]:
    header = None
    rows: list[list[str]] = []
    with path.open() as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            parts = text.split()
            if header is None:
                header = parts
            else:
                rows.append(parts)
    if header is None or not rows:
        raise ValueError(f"No table rows found in {path}")
    return header, rows


def _energy_range_to_hz(emin: float, emax: float, unit: str) -> tuple[float, float]:
    unit = unit.lower()
    if unit == "kev":
        scale = KEV_TO_HZ
    elif unit == "gev":
        scale = GEV_TO_HZ
    elif unit == "tev":
        scale = TEV_TO_HZ
    else:
        raise ValueError(f"Unsupported energy unit: {unit}")
    return float(emin) * scale, float(emax) * scale


def _time_bin(t_s: float, tp_s: float, tm_s: float) -> tuple[float, float, float]:
    start = min(t_s - tm_s, t_s + tp_s)
    stop = max(t_s - tm_s, t_s + tp_s)
    return start, stop, 0.5 * (start + stop)


def _ab_mag_to_fnu_cgs(mag: float, mag_err: float) -> tuple[float, float]:
    flux = 10.0 ** (-0.4 * (mag + AB_ZEROPOINT))
    if np.isfinite(mag_err) and mag_err > 0.0:
        err = flux * np.log(10.0) * 0.4 * mag_err
    else:
        err = 0.1 * flux
    return float(flux), float(abs(err))


def _optical_band_label(raw_band: str) -> str:
    return OPTICAL_BAND_ALIASES.get(raw_band, raw_band)


@lru_cache(maxsize=3)
def _zdust_model(extinction_law: str):
    ZDust, _ = _import_astromodels_functions()
    model = ZDust()
    model.extinction_law = extinction_law
    return model


@lru_cache(maxsize=256)
def _zdust_transmission(extinction_law: str, energy_kev: float, ebv: float, rv: float, redshift: float) -> float:
    model = _zdust_model(extinction_law)
    values = model.evaluate(np.asarray([float(energy_kev)], dtype=float), float(ebv), float(rv), float(redshift))
    return float(np.asarray(values, dtype=float)[0])


def _astromodels_deredden_factor(wavelength_a: float) -> float:
    energy_kev = ANGSTROM_TO_KEV / float(wavelength_a)
    trans_mw = _zdust_transmission("mw", energy_kev, MW_EBV, MW_RV, 0.0)
    trans_host = _zdust_transmission("smc", energy_kev, HOST_EBV, HOST_RV, TARGET_REDSHIFT)
    return 1.0 / (trans_mw * trans_host)


@lru_cache(maxsize=128)
def _xrt_absorption_transmission(mode: str, energy_kev_tuple: tuple[float, ...]) -> tuple[float, ...]:
    _, TbAbs = _import_astromodels_functions()
    energies = np.asarray(energy_kev_tuple, dtype=float)
    mode_key = str(mode).upper()
    if mode_key not in XRT_NH_HOST_1E22:
        raise ValueError(f"Unsupported XRT mode for absorption: {mode}")
    gal_model = TbAbs()
    host_model = TbAbs()
    gal_model.abundance_table = "WILM"
    host_model.abundance_table = "WILM"
    gal = np.asarray(gal_model.evaluate(energies, XRT_NH_GAL_1E22, 0.0), dtype=float)
    host = np.asarray(host_model.evaluate(energies, XRT_NH_HOST_1E22[mode_key], TARGET_REDSHIFT), dtype=float)
    return tuple((gal * host).tolist())


def _xrt_transmission_for_freq(mode: str, frequencies_hz: np.ndarray) -> np.ndarray:
    energy_kev = np.asarray(frequencies_hz, dtype=float) * constants.para_hz2kev
    return np.asarray(_xrt_absorption_transmission(mode.upper(), tuple(float(v) for v in energy_kev)), dtype=float)


@lru_cache(maxsize=128)
def _vhe_ebl_transmission(frequency_hz_tuple: tuple[float, ...]) -> tuple[float, ...]:
    frequencies = np.asarray(frequency_hz_tuple, dtype=float)
    redshifts, energies_hz, tau_values = _load_saldana21_ebl_table()
    tau_at_z = np.asarray([np.interp(TARGET_REDSHIFT, redshifts, tau_values[i, :]) for i in range(tau_values.shape[0])])
    tau_obs = np.interp(frequencies, energies_hz, tau_at_z, left=0.0, right=1.0e30)
    transmission = np.exp(-tau_obs)
    return tuple(transmission.tolist())


def _vhe_ebl_transmission_for_freq(frequencies_hz: np.ndarray) -> np.ndarray:
    return np.asarray(_vhe_ebl_transmission(tuple(float(v) for v in frequencies_hz)), dtype=float)


@lru_cache(maxsize=1)
def _load_saldana21_ebl_table() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = [line for line in VHE_EBL_TABLE.read_text(encoding="utf-8").splitlines() if line.strip()]
    redshifts = np.asarray([float(value) for value in lines[0].split()], dtype=float)
    rows = [[float(value) for value in line.split()] for line in lines[1:]]
    table = np.asarray(rows, dtype=float)
    energies_hz = table[:, 0] * TEV_TO_HZ
    tau_values = table[:, 1 : 1 + redshifts.size]
    return redshifts, energies_hz, tau_values


def _vhe_ebl_metadata() -> dict[str, Any]:
    redshifts, energies_hz, tau_values = _load_saldana21_ebl_table()
    return {
        "model": VHE_EBL_MODEL,
        "table": str(VHE_EBL_TABLE),
        "z": TARGET_REDSHIFT,
        "applied_to": "vhe_data.txt",
        "energy_min_TeV": float(energies_hz[0] / TEV_TO_HZ),
        "energy_max_TeV": float(energies_hz[-1] / TEV_TO_HZ),
        "n_energy": int(energies_hz.size),
        "n_redshift": int(redshifts.size),
        "tau_shape": list(tau_values.shape),
    }


def _precompute_optical_deredden_factors(data: list[list[str]], columns: dict[str, int]) -> dict[tuple[str, float], float]:
    unique_bands: dict[tuple[str, float], float] = {}
    for row in data:
        band = _optical_band_label(row[columns["band"]])
        wavelength_a = float(row[columns["wavelength[A]"]])
        unique_bands[(band, wavelength_a)] = _astromodels_deredden_factor(wavelength_a)
    return unique_bands


def _wide_band_rows(path: Path, energy_unit: str, label: str) -> list[ArticleData]:
    header, data = _read_header_rows(path)
    columns = {name: index for index, name in enumerate(header)}
    rows = []
    for raw in data:
        t_s = float(raw[columns["t[s]"]])
        tp_s = float(raw[columns["tp[s]"]])
        tm_s = float(raw[columns["tm[s]"]])
        start_s, stop_s, _ = _time_bin(t_s, tp_s, tm_s)
        nu_min, nu_max = _energy_range_to_hz(
            float(raw[columns[f"Emin[{energy_unit}]"]]),
            float(raw[columns[f"Emax[{energy_unit}]"]]),
            energy_unit,
        )
        flux = float(raw[columns["Flux[erg/s/cm2]"]])
        fp = float(raw[columns["Fp[erg/s/cm2]"]]) if raw[columns["Fp[erg/s/cm2]"]].lower() != "nan" else np.nan
        fm = float(raw[columns["Fm[erg/s/cm2]"]]) if raw[columns["Fm[erg/s/cm2]"]].lower() != "nan" else np.nan
        upper = bool(flux <= 0.0 or not (np.isfinite(fp) and np.isfinite(fm)))
        err = 0.5 * (abs(fp) + abs(fm)) if np.isfinite(fp) and np.isfinite(fm) else abs(flux) * 0.3
        mode = raw[columns["mode"]].upper() if "mode" in columns else ""
        rows.append(
            ArticleData(
                kind="wide",
                file=path.name,
                label=label if not mode else f"{label} {mode}",
                time_s=t_s,
                time_start_s=start_s,
                time_stop_s=stop_s,
                model_time_s=t_s - MODEL_T0_OFFSET_S,
                model_time_start_s=start_s - MODEL_T0_OFFSET_S,
                model_time_stop_s=stop_s - MODEL_T0_OFFSET_S,
                exposure_s=tp_s + tm_s,
                frequency_hz=float(np.sqrt(nu_min * nu_max)),
                energy_min_hz=nu_min,
                energy_max_hz=nu_max,
                flux=flux,
                error=float(err),
                upper_limit=upper,
                mode=mode,
            )
        )
    return rows


def _load_optical(path: Path) -> list[ArticleData]:
    header, data = _read_header_rows(path)
    columns = {name: index for index, name in enumerate(header)}
    r_frequency_hz = 2.99792458e18 / OPTICAL_R_WAVELENGTH_A
    r_correction = _astromodels_deredden_factor(OPTICAL_R_WAVELENGTH_A)
    rows = []
    for raw in data:
        t_s = float(raw[columns["time[s]"]])
        mag = float(raw[columns["mag"]])
        mag_err_text = raw[columns["mag_err"]].lower()
        mag_err = float(mag_err_text) if mag_err_text != "nan" else np.nan
        flux, err = _ab_mag_to_fnu_cgs(mag, mag_err)
        rows.append(
            ArticleData(
                kind="fnu",
                file=path.name,
                label=OPTICAL_R_LABEL,
                time_s=t_s,
                time_start_s=t_s,
                time_stop_s=t_s,
                model_time_s=t_s - MODEL_T0_OFFSET_S,
                model_time_start_s=t_s - MODEL_T0_OFFSET_S,
                model_time_stop_s=t_s - MODEL_T0_OFFSET_S,
                exposure_s=0.0,
                frequency_hz=float(r_frequency_hz),
                flux=float(flux * r_correction),
                error=float(err * r_correction),
            )
        )
    return rows


def _radio_nominal_frequency_ghz(nu_ghz: float) -> float | None:
    for low, high, nominal in RADIO_SELECTIONS:
        if low <= nu_ghz <= high:
            return nominal
    return None


def _load_radio(path: Path) -> list[ArticleData]:
    header, data = _read_header_rows(path)
    columns = {name: index for index, name in enumerate(header)}
    rows = []
    for raw in data:
        nu_ghz = float(raw[columns["nu[GHz]"]])
        nominal_ghz = _radio_nominal_frequency_ghz(nu_ghz)
        if nominal_ghz is None:
            continue
        t_s = float(raw[columns["t[d]"]]) * 86400.0
        flux = float(raw[columns["Fnu[mJy]"]]) * MJY_CGS
        err = abs(float(raw[columns["err[mJy]"]])) * MJY_CGS
        rows.append(
            ArticleData(
                kind="fnu",
                file=path.name,
                label=f"Radio {nominal_ghz:g} GHz",
                time_s=t_s,
                time_start_s=t_s,
                time_stop_s=t_s,
                model_time_s=t_s - MODEL_T0_OFFSET_S,
                model_time_start_s=t_s - MODEL_T0_OFFSET_S,
                model_time_stop_s=t_s - MODEL_T0_OFFSET_S,
                exposure_s=0.0,
                frequency_hz=nominal_ghz * 1.0e9,
                flux=flux,
                error=err,
            )
        )
    return rows


@lru_cache(maxsize=1)
def load_grb180720b_data() -> tuple[ArticleData, ...]:
    if not DATA_DIR.exists():
        raise FileNotFoundError(f"Data directory does not exist: {DATA_DIR}")
    rows: list[ArticleData] = []
    rows.extend(_wide_band_rows(DATA_DIR / "XRT.txt", "kev", "XRT 0.3-10 keV"))
    rows.extend(_wide_band_rows(DATA_DIR / "LAT.txt", "gev", "LAT 0.1-1 GeV"))
    rows.extend(_wide_band_rows(DATA_DIR / "vhe_data.txt", "tev", "HESS 0.1-0.4 TeV"))
    rows.extend(_load_optical(DATA_DIR / "OPT.txt"))
    rows.extend(_load_radio(DATA_DIR / "radio.txt"))
    if not any(row.file == "XRT.txt" and row.mode in {"WT", "PC"} for row in rows):
        raise RuntimeError("XRT.txt must include WT/PC mode labels for absorption handling.")
    return tuple(rows)


def _wide_energy_sample_count(row: ArticleData, default_samples: int = DEFAULT_WIDE_BAND_SAMPLES) -> int:
    if row.file == "XRT.txt":
        return XRT_WIDE_BAND_SAMPLES
    return int(default_samples)


def _lookup(values: np.ndarray) -> dict[float, int]:
    return {float(value): index for index, value in enumerate(values)}


def _trapz_weights(x: np.ndarray) -> np.ndarray:
    values = np.asarray(x, dtype=float)
    weights = np.zeros_like(values)
    if values.size < 2:
        return weights
    weights[0] = 0.5 * (values[1] - values[0])
    weights[-1] = 0.5 * (values[-1] - values[-2])
    if values.size > 2:
        weights[1:-1] = 0.5 * (values[2:] - values[:-2])
    return weights


def build_observation_grid(rows: tuple[ArticleData, ...] | None = None) -> dict[str, Any]:
    rows = rows or load_grb180720b_data()
    frequencies = []
    all_times = []
    for row in rows:
        if row.kind == "wide":
            n_samples = _wide_energy_sample_count(row)
            frequencies.extend(np.geomspace(row.energy_min_hz, row.energy_max_hz, n_samples).tolist())
        else:
            frequencies.append(row.frequency_hz)
        all_times.append(max(row.model_time_s, 1.0e-6))
    unique_frequencies = np.unique(np.asarray(frequencies, dtype=float))
    unique_times = np.unique(np.asarray(all_times, dtype=float))
    sample_specs = _build_sample_specs(rows, unique_frequencies, unique_times)
    return {
        "rows": rows,
        "frequencies_hz": unique_frequencies,
        "model_times_s": unique_times,
        "sample_specs_by_row": sample_specs,
        "fit_data": _build_fit_data(rows, sample_specs),
    }


def _build_sample_specs(rows: tuple[ArticleData, ...], frequencies: np.ndarray, times: np.ndarray) -> list[dict[str, Any]]:
    freq_lookup = _lookup(frequencies)
    time_lookup = _lookup(times)
    specs = []
    for row in rows:
        time_index = time_lookup[float(max(row.model_time_s, 1.0e-6))]
        if row.kind == "wide":
            n_samples = _wide_energy_sample_count(row)
            band_freq = np.geomspace(row.energy_min_hz, row.energy_max_hz, n_samples)
            freq_indices = [freq_lookup[float(freq)] for freq in band_freq]
            transmission = np.ones_like(band_freq)
            if row.file == "XRT.txt":
                transmission = _xrt_transmission_for_freq(row.mode, band_freq)
            elif row.file == "vhe_data.txt":
                transmission = _vhe_ebl_transmission_for_freq(band_freq)
            specs.append(
                {
                    "kind": "wide",
                    "time_index": int(time_index),
                    "freq_indices": freq_indices,
                    "band_freq": [float(v) for v in band_freq],
                    "transmission": [float(v) for v in transmission],
                    "mode": row.mode,
                }
            )
        else:
            specs.append(
                {
                    "kind": "fnu",
                    "time_index": int(time_index),
                    "freq_indices": [freq_lookup[float(row.frequency_hz)]],
                    "band_freq": [],
                    "transmission": [],
                    "mode": row.mode,
                }
            )
    return specs


def _build_fit_data(rows: tuple[ArticleData, ...], specs: list[dict[str, Any]]) -> dict[str, Any]:
    fnu_freq_indices = []
    fnu_time_indices = []
    fnu_flux = []
    fnu_error = []
    wide_groups: dict[tuple[int, str], dict[str, list[Any]]] = {}
    for row, spec in zip(rows, specs):
        if not row.valid_for_fit:
            continue
        if spec["kind"] == "wide":
            band_freq = np.asarray(spec["band_freq"], dtype=float)
            transmission = np.asarray(spec["transmission"], dtype=float)
            weights = _trapz_weights(np.log(band_freq)) * band_freq * transmission
            key = (len(spec["freq_indices"]), row.file)
            group = wide_groups.setdefault(key, {"freq_indices": [], "time_indices": [], "weights": [], "flux": [], "error": []})
            group["freq_indices"].append([int(v) for v in spec["freq_indices"]])
            group["time_indices"].append(int(spec["time_index"]))
            group["weights"].append(weights)
            group["flux"].append(row.flux)
            group["error"].append(abs(row.error))
        else:
            fnu_freq_indices.append(int(spec["freq_indices"][0]))
            fnu_time_indices.append(int(spec["time_index"]))
            fnu_flux.append(row.flux)
            fnu_error.append(abs(row.error))
    wide_group_arrays = []
    obs_flux_parts = []
    obs_error_parts = []
    if fnu_flux:
        obs_flux_parts.append(np.asarray(fnu_flux, dtype=float))
        obs_error_parts.append(np.asarray(fnu_error, dtype=float))
    for key in sorted(wide_groups):
        group = wide_groups[key]
        group_flux = np.asarray(group["flux"], dtype=float)
        group_error = np.asarray(group["error"], dtype=float)
        obs_flux_parts.append(group_flux)
        obs_error_parts.append(group_error)
        wide_group_arrays.append(
            {
                "freq_indices": np.asarray(group["freq_indices"], dtype=np.intp),
                "time_indices": np.asarray(group["time_indices"], dtype=np.intp),
                "weights": np.asarray(group["weights"], dtype=float),
                "flux": group_flux,
                "error": group_error,
            }
        )
    obs_flux = np.concatenate(obs_flux_parts) if obs_flux_parts else np.asarray([], dtype=float)
    obs_error = np.concatenate(obs_error_parts) if obs_error_parts else np.asarray([], dtype=float)
    return {
        "fnu_freq_indices": np.asarray(fnu_freq_indices, dtype=np.intp),
        "fnu_time_indices": np.asarray(fnu_time_indices, dtype=np.intp),
        "fnu_flux": np.asarray(fnu_flux, dtype=float),
        "fnu_error": np.asarray(fnu_error, dtype=float),
        "wide_groups": tuple(wide_group_arrays),
        "obs_flux": obs_flux,
        "sigma_obs": obs_error,
    }


@lru_cache(maxsize=1)
def cached_observation_grid() -> dict[str, Any]:
    return build_observation_grid()


def article_params_linear(params: dict[str, float]) -> dict[str, float]:
    return {
        "E_iso": 10.0 ** params["log10_Ekiso"],
        "Eta_0": 10.0 ** params["log10_Gamma0"],
        "p_f": params["pf"],
        "Epsilon_e_f": 10.0 ** params["log10_eps_e_f"],
        "Epsilon_b_f": 10.0 ** params["log10_eps_B_f"],
        "f_e_f": 10.0 ** params["log10_xi_e_f"],
        "theta_j": 10.0 ** params["log10_theta_j"],
        "p_r": params["pr"],
        "Epsilon_e_r": 10.0 ** params["log10_eps_e_r"],
        "Epsilon_b_r": 10.0 ** params["log10_eps_B_r"],
        "f_e_r": 10.0 ** params["log10_xi_e_r"],
        "dNe": 10.0 ** params["log10_dNe"],
        "Bw_factor": 10.0 ** params.get("log10_Bw_B0", 2.0),
    }


def config_from_params(model_name: str, params: dict[str, float], *, num_r: int, num_theta: int) -> FitConfig:
    p = article_params_linear(params)
    is_2d_pic = model_name == "fullhide_2d_pic"
    return FitConfig(
        num_threads=int(os.environ.get("OMP_NUM_THREADS", "20")),
        index_dyn=FIXED_NUMERICAL["index_dyn"],
        index_y=FIXED_NUMERICAL["index_Y"],
        index_syn_integr=FIXED_NUMERICAL["index_syn_integr"],
        electron_solver=MODEL_TO_SOLVER[model_name],
        include_forward_ssc=True,
        electron_adaptive_substeps=False,
        electron_pic_uniform_b=False if is_2d_pic else True,
        electron_pic_eta_acc=1.0,
        electron_pic_kappa_diff_scale=1.0,
        electron_pic_bw_factor=p["Bw_factor"] if is_2d_pic else 100.0,
        num_gam_e=FIXED_NUMERICAL["Num_gam_e"],
        num_nu=FIXED_NUMERICAL["Num_nu"],
        num_r=int(num_r),
        num_theta=int(num_theta),
        num_phi=FIXED_NUMERICAL["Num_phi"],
        num_tobs=FIXED_NUMERICAL["Num_tobs"],
        num_chi=10 if is_2d_pic else None,
        z=TARGET_REDSHIFT,
        eta_0=p["Eta_0"],
        epsilon_e=p["Epsilon_e_f"],
        epsilon_b=p["Epsilon_b_f"],
        p=p["p_f"],
        opening_angle_jet=p["theta_j"],
        theta_v=0.0,
        f_e=p["f_e_f"],
        e_iso=p["E_iso"],
        d_ne=p["dNe"],
        a_star=FIXED_A_STAR,
        r0=TARGET_R0_CM,
        initial_radius_cm=1.0e15,
        t_obs_min_log10=2.0,
        t_obs_max_log10=8.0,
        luminosity_distance_cm_override=_cosmo_d_l(TARGET_REDSHIFT),
        reverse=True,
        reverse_shock=ReverseShockConfig(
            enabled=True,
            delta_t_s=T90_S / (1.0 + TARGET_REDSHIFT),
            epsilon_e=p["Epsilon_e_r"],
            epsilon_b=p["Epsilon_b_r"],
            p=p["p_r"],
            f_e=p["f_e_r"],
            include_ssc=True,
        ),
    )


def _effective_variance(obs_flux: np.ndarray, sigma_obs: np.ndarray, f_sys: float) -> np.ndarray:
    flux_scale = np.abs(np.asarray(obs_flux, dtype=float))
    positive = flux_scale[np.isfinite(flux_scale) & (flux_scale > 0.0)]
    tiny_flux = 1.0e-300 if positive.size == 0 else min(float(positive.min()) * 1.0e-12, 1.0e-300)
    flux_scale = np.maximum(flux_scale, tiny_flux)
    sigma = np.asarray(sigma_obs, dtype=float)
    return sigma * sigma + (float(f_sys) * flux_scale) ** 2


def _collect_fit_arrays(total_grid: np.ndarray, cached_fit_data: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    model_flux_parts = []
    obs_flux_parts = []
    sigma_parts = []
    fnu_freq = cached_fit_data["fnu_freq_indices"]
    if fnu_freq.size:
        fnu_time = cached_fit_data["fnu_time_indices"]
        model_flux_parts.append(total_grid[fnu_freq, fnu_time])
        obs_flux_parts.append(cached_fit_data["fnu_flux"])
        sigma_parts.append(cached_fit_data["fnu_error"])
    for group in cached_fit_data.get("wide_groups", ()):
        wide_freq = group["freq_indices"]
        if wide_freq.size:
            wide_time = group["time_indices"]
            values = total_grid[wide_freq, wide_time[:, None]]
            model_flux_parts.append(np.sum(values * group["weights"], axis=1))
            obs_flux_parts.append(group["flux"])
            sigma_parts.append(group["error"])
    if not model_flux_parts:
        empty = np.asarray([], dtype=float)
        return empty, empty, empty
    return np.concatenate(model_flux_parts), np.concatenate(obs_flux_parts), np.concatenate(sigma_parts)


def model_flux_from_total_grid(total_grid: np.ndarray, cached_fit_data: dict[str, Any]) -> np.ndarray:
    model_flux, _, _ = _collect_fit_arrays(total_grid, cached_fit_data)
    return model_flux


def chi2_from_model_flux(model_flux: np.ndarray, cached_fit_data: dict[str, Any], f_sys: float) -> float:
    obs_flux = cached_fit_data["obs_flux"]
    sigma_obs = cached_fit_data["sigma_obs"]
    if model_flux.size == 0 or model_flux.shape != obs_flux.shape:
        return 1.0e100
    var_eff = _effective_variance(obs_flux, sigma_obs, f_sys)
    if np.any(var_eff <= 0.0) or not np.all(np.isfinite(model_flux)):
        return 1.0e100
    return float(np.sum((model_flux - obs_flux) ** 2 / var_eff))


def loglike_from_model_flux(model_flux: np.ndarray, cached_fit_data: dict[str, Any], f_sys: float) -> float:
    obs_flux = cached_fit_data["obs_flux"]
    sigma_obs = cached_fit_data["sigma_obs"]
    if model_flux.size == 0 or model_flux.shape != obs_flux.shape:
        return -1.0e100
    var_eff = _effective_variance(obs_flux, sigma_obs, f_sys)
    if np.any(var_eff <= 0.0) or not np.all(np.isfinite(model_flux)):
        return -1.0e100
    chi2 = np.sum((model_flux - obs_flux) ** 2 / var_eff)
    return float(-0.5 * (chi2 + np.sum(np.log(2.0 * math.pi * var_eff))))


def _add_timing(timings: dict[str, float] | None, label: str, elapsed: float) -> None:
    if timings is not None:
        timings[label] = timings.get(label, 0.0) + float(elapsed)


def evaluate_model_grid(
    model_name: str,
    params: dict[str, float],
    *,
    num_r: int,
    num_theta: int,
    full_components: bool = False,
    timings: dict[str, float] | None = None,
) -> dict[str, Any]:
    t_step = time.perf_counter()
    grid = cached_observation_grid()
    _add_timing(timings, "fit.cached_observation_grid", time.perf_counter() - t_step)

    t_step = time.perf_counter()
    config = config_from_params(model_name, params, num_r=num_r, num_theta=num_theta)
    _add_timing(timings, "fit.config_from_params", time.perf_counter() - t_step)

    t_step = time.perf_counter()
    setup = build_simulation_setup(config)
    setup.observer_time_s = np.asarray(grid["model_times_s"], dtype=float)
    _add_timing(timings, "fit.build_simulation_setup", time.perf_counter() - t_step)

    t_step = time.perf_counter()
    state = solve_state_from_setup(config, setup, timings=timings, requested_frequencies_hz=grid["frequencies_hz"])
    _add_timing(timings, "fit.solve_state_from_setup", time.perf_counter() - t_step)

    t_step = time.perf_counter()
    observed = project_flux_grid(
        state,
        grid["model_times_s"],
        grid["frequencies_hz"],
        timings=timings,
        mode="full_components" if full_components else "total_only",
    )
    _add_timing(timings, "fit.project_flux_grid", time.perf_counter() - t_step)

    t_step = time.perf_counter()
    components = observed.components
    total = np.asarray(components["total"], dtype=float)
    _add_timing(timings, "fit.extract_total_grid", time.perf_counter() - t_step)

    t_step = time.perf_counter()
    model_fit_flux = model_flux_from_total_grid(total, grid["fit_data"])
    _add_timing(timings, "fit.model_flux_from_total_grid", time.perf_counter() - t_step)

    result = {
        "grid": grid,
        "config": config,
        "frequencies_hz": grid["frequencies_hz"],
        "model_times_s": grid["model_times_s"],
        "total": total,
        "model_fit_flux": model_fit_flux,
    }
    if full_components:
        fwd = np.zeros_like(total)
        for key in ("fwd_sync", "fwd_ssc", "fwd_hadronic", "fwd_hadronic_bethe_heitler", "fwd_hadronic_inverse_compton", "fwd_hadronic_pair_production"):
            if components.get(key) is not None:
                fwd += np.asarray(components[key], dtype=float)
        rev = np.zeros_like(total)
        for key in ("rev_sync", "rev_ssc", "cross_ic"):
            if components.get(key) is not None:
                rev += np.asarray(components[key], dtype=float)
        result.update({"fs": fwd, "rs": rev})
    return result


def fit_chi2(model_name: str, params: dict[str, float], *, num_r: int, num_theta: int) -> float:
    result = evaluate_model_grid(model_name, params, num_r=num_r, num_theta=num_theta)
    return chi2_from_model_flux(result["model_fit_flux"], result["grid"]["fit_data"], f_sys_from_params(params))


def fit_loglike(model_name: str, params: dict[str, float], *, num_r: int, num_theta: int, timings: dict[str, float] | None = None) -> float:
    result = evaluate_model_grid(model_name, params, num_r=num_r, num_theta=num_theta, timings=timings)
    t_step = time.perf_counter()
    value = loglike_from_model_flux(result["model_fit_flux"], result["grid"]["fit_data"], f_sys_from_params(params))
    _add_timing(timings, "fit.loglike_from_model_flux", time.perf_counter() - t_step)
    return value


def _progress_interval_seconds() -> float:
    raw = os.environ.get("ASGARD_LIKE_PROGRESS_INTERVAL", "0").strip()
    if not raw:
        return 0.0
    try:
        return max(0.0, float(raw))
    except ValueError:
        return 0.0


def _record_likelihood_progress(prefix: Path, model_name: str, status: str, elapsed: float, value: float | None = None) -> None:
    interval = _progress_interval_seconds()
    if interval <= 0.0:
        return
    rank = mpi_rank()
    now = time.time()
    state = _LIKE_PROGRESS.setdefault(rank, {"count": 0, "last_write": 0.0})
    state["count"] += 1
    if now - float(state["last_write"]) < interval:
        return
    state["last_write"] = now
    payload = {
        "model": model_name,
        "rank": rank,
        "count": int(state["count"]),
        "status": status,
        "last_elapsed_seconds": float(elapsed),
        "last_loglike": None if value is None else float(value),
        "updated_time": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    progress_dir = prefix.parent / "logs"
    progress_dir.mkdir(parents=True, exist_ok=True)
    tmp_path = progress_dir / f"likelihood_rank_{rank}.json.tmp"
    out_path = progress_dir / f"likelihood_rank_{rank}.json"
    tmp_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp_path.replace(out_path)


def _slow_likelihood_threshold_seconds() -> float:
    raw = os.environ.get("ASGARD_SLOW_LIKE_THRESHOLD", "10").strip()
    if not raw:
        return 0.0
    try:
        return max(0.0, float(raw))
    except ValueError:
        return 10.0


def _record_slow_likelihood(
    prefix: Path,
    model_name: str,
    params: dict[str, float],
    timings: dict[str, float],
    elapsed: float,
    value: float,
) -> None:
    threshold = _slow_likelihood_threshold_seconds()
    if threshold <= 0.0 or elapsed < threshold:
        return
    rank = mpi_rank()
    progress_dir = prefix.parent / "logs"
    progress_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "model": model_name,
        "rank": rank,
        "elapsed_seconds": float(elapsed),
        "loglike": float(value),
        "timings": {name: float(item) for name, item in sorted(timings.items())},
        "params": {name: float(item) for name, item in params.items()},
        "updated_time": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    with (progress_dir / f"slow_likelihood_rank_{rank}.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")


def loglike(cube, ndim, nparams, model_name: str, num_r: int, num_theta: int) -> float:
    prefix = output_prefix(model_name)
    start = time.perf_counter()
    names = parameter_names(model_name)
    if not np.all(np.isfinite(cube[: len(names)])):
        _record_likelihood_progress(prefix, model_name, "nonfinite_cube", time.perf_counter() - start, -1.0e100)
        return -1.0e100
    params = params_from_vector(cube[: len(names)], model_name)
    if _invalid_physical(params):
        _record_likelihood_progress(prefix, model_name, "invalid_physical", time.perf_counter() - start, -1.0e100)
        return -1.0e100
    timings: dict[str, float] = {}
    try:
        value = fit_loglike(model_name, params, num_r=num_r, num_theta=num_theta, timings=timings)
    except Exception as exc:
        _record_likelihood_progress(prefix, model_name, f"exception:{type(exc).__name__}", time.perf_counter() - start, -1.0e100)
        return -1.0e100
    if not np.isfinite(value):
        _record_likelihood_progress(prefix, model_name, "nonfinite_loglike", time.perf_counter() - start, -1.0e100)
        return -1.0e100
    elapsed = time.perf_counter() - start
    _record_likelihood_progress(prefix, model_name, "ok", elapsed, float(value))
    _record_slow_likelihood(prefix, model_name, params, timings, elapsed, float(value))
    return float(value)


def _model_row_flux(row: ArticleData, grid: np.ndarray, frequencies: np.ndarray, times: np.ndarray, spec: dict[str, Any]) -> float:
    time_index = int(spec["time_index"])
    if spec["kind"] == "wide":
        freq_indices = [int(index) for index in spec["freq_indices"]]
        band_freq = np.asarray(spec["band_freq"], dtype=float)
        transmission = np.asarray(spec["transmission"], dtype=float)
        vals = grid[freq_indices, time_index]
        return float(np.sum(vals * _trapz_weights(np.log(band_freq)) * band_freq * transmission))
    freq_index = int(spec["freq_indices"][0])
    return float(grid[freq_index, time_index])


def build_tophat_model(model_name: str, params: dict[str, float], *, num_r: int, num_theta: int) -> dict[str, Any]:
    result = evaluate_model_grid(model_name, params, num_r=num_r, num_theta=num_theta, full_components=True)
    grid = result["grid"]
    model = {
        "data_rows": grid["rows"],
        "sample_specs_by_row": grid["sample_specs_by_row"],
        "fit_data": grid["fit_data"],
        "frequencies_hz": result["frequencies_hz"],
        "model_times_s": result["model_times_s"],
        "total": result["total"],
        "fs": result["fs"],
        "rs": result["rs"],
        "params": params,
        "model_name": model_name,
        "metadata": {
            "grb": TARGET_GRB,
            "geometry": "top_hat",
            "components": "fs_rs_syn_ssc",
            "z": TARGET_REDSHIFT,
            "R0_cm": TARGET_R0_CM,
            "external_medium": "ISM",
            "electron_solver": MODEL_TO_SOLVER[model_name],
            "xrt_absorption": {
                "NH_gal_1e22_cm2": XRT_NH_GAL_1E22,
                "NH_host_1e22_cm2": XRT_NH_HOST_1E22,
                "z_host": TARGET_REDSHIFT,
                "models": "TbAbs(redshift=0) * TbAbs(redshift=z_host)",
                "abundance_table": "WILM",
            },
            "optical_extinction": {
                "models": "ZDust mw + ZDust smc",
                "MW_EBV": MW_EBV,
                "MW_RV": MW_RV,
                "HOST_EBV": HOST_EBV,
                "HOST_RV": HOST_RV,
            },
            "vhe_ebl": _vhe_ebl_metadata(),
        },
    }
    gc.collect()
    return model


def save_model_cache(model: dict[str, Any], output_prefix: Path) -> Path:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    cache_path = output_prefix.with_suffix(".npz")
    arrays: dict[str, Any] = {
        "frequencies_hz": model["frequencies_hz"],
        "model_times_s": model["model_times_s"],
        "rows_json": np.asarray(json.dumps([asdict(row) for row in model["data_rows"]])),
        "sample_specs_json": np.asarray(json.dumps(model["sample_specs_by_row"])),
        "params_json": np.asarray(json.dumps(model.get("params", {}))),
        "metadata_json": np.asarray(json.dumps(model.get("metadata", {}))),
        "model_name": np.asarray(model["model_name"]),
    }
    for key in ("total", "fs", "rs"):
        arrays[key] = model[key]
    np.savez_compressed(cache_path, **arrays)
    return cache_path


def _finite_positive(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float).copy()
    arr[~np.isfinite(arr) | (arr <= 0.0)] = np.nan
    return arr


def _frequency_indices(frequencies: np.ndarray, target_frequencies: np.ndarray) -> list[int]:
    return [int(np.argmin(np.abs(frequencies - freq))) for freq in target_frequencies]


def _model_nufnu_curve(model: dict[str, Any], component: str, freq_hz: float) -> np.ndarray:
    freqs = model["frequencies_hz"]
    i = int(np.argmin(np.abs(freqs - freq_hz)))
    return np.asarray(model[component][i] * freqs[i], dtype=float)


def _model_band_flux_curve(model: dict[str, Any], component: str, row: ArticleData) -> np.ndarray:
    freqs = model["frequencies_hz"]
    band_freq = np.geomspace(row.energy_min_hz, row.energy_max_hz, _wide_energy_sample_count(row))
    indices = _frequency_indices(freqs, band_freq)
    if row.file == "XRT.txt":
        transmission = _xrt_transmission_for_freq(row.mode, band_freq)
    elif row.file == "vhe_data.txt":
        transmission = _vhe_ebl_transmission_for_freq(band_freq)
    else:
        transmission = np.ones_like(band_freq)
    return np.sum(model[component][indices] * _trapz_weights(np.log(band_freq))[:, None] * band_freq[:, None] * transmission[:, None], axis=0)


def _row_plot_y(row: ArticleData) -> float:
    return row.flux if row.kind == "wide" else row.frequency_hz * row.flux


def _stack_factor(rows: list[ArticleData], rank: int, spacing_decades: float = 2.5) -> float:
    y = np.asarray([_row_plot_y(row) for row in rows if not row.upper_limit], dtype=float)
    y = y[np.isfinite(y) & (y > 0.0)]
    median_y = float(np.median(y)) if y.size else 1.0
    exponent = spacing_decades * float(rank) - np.log10(max(median_y, 1.0e-300))
    exponent = round(exponent * 2.0) / 2.0
    return 10.0 ** exponent


def _plot_component_curves(ax, times, curves, factor, color, label) -> None:
    ax.loglog(times, _finite_positive(curves["total"] * factor), color=color, lw=2.0, label=label)
    ax.loglog(times, _finite_positive(curves["fs"] * factor), color=color, ls="--", lw=1.35, alpha=0.95)
    ax.loglog(times, _finite_positive(curves["rs"] * factor), color=color, ls="-.", lw=1.2, alpha=0.85)


def _scatter_rows(ax, rows: list[ArticleData], factor: float, color, marker: str) -> None:
    for row in rows:
        plot_time = max(row.model_time_s, 1.0e-6)
        y = _row_plot_y(row)
        yerr = row.error if row.kind == "wide" else row.frequency_hz * row.error
        if row.upper_limit:
            ax.errorbar(plot_time, y * factor, yerr=max(y, 1.0e-99) * factor * 0.3, uplims=True, fmt="v", color=color, ms=4, alpha=0.9)
        else:
            ax.errorbar(plot_time, y * factor, yerr=abs(yerr * factor), fmt=marker, color=color, ecolor=color, elinewidth=0.7, capsize=1.4, ms=3.1, alpha=0.75)


def plot_fig3_style(model: dict[str, Any], output_prefix: Path) -> Path:
    from collections import defaultdict

    plt = _load_pyplot()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 1, figsize=(16.5, 9.5))
    times = model["model_times_s"]
    rows_by_group: dict[str, list[ArticleData]] = defaultdict(list)
    for row in model["data_rows"]:
        if row.file == "XRT.txt":
            rows_by_group[f"XRT {row.mode}"].append(row)
        elif row.file == "LAT.txt":
            rows_by_group["LAT 0.1-1 GeV"].append(row)
        elif row.file == "vhe_data.txt":
            rows_by_group["HESS 0.1-0.4 TeV"].append(row)
        else:
            rows_by_group[row.label].append(row)
    groups = []
    colors = plt.cm.tab20(np.linspace(0.0, 1.0, max(len(rows_by_group), 1)))
    for i, (label, rows) in enumerate(sorted(rows_by_group.items())):
        ref = rows[0]
        if ref.kind == "wide":
            curves = {component: _model_band_flux_curve(model, component, ref) for component in ("total", "fs", "rs")}
            marker = "^" if ref.file == "XRT.txt" else "o"
        else:
            freq = float(np.median([row.frequency_hz for row in rows]))
            curves = {component: _model_nufnu_curve(model, component, freq) for component in ("total", "fs", "rs")}
            marker = "D" if ref.file == "OPT.txt" else "s"
        groups.append((label, rows, colors[i], marker, curves))
    y_values = []
    x_values = []
    for rank, (label, rows, color, marker, curves) in enumerate(groups):
        factor = _stack_factor(rows, rank)
        _plot_component_curves(ax, times, curves, factor, color, f"{label} x {factor:.1e}")
        _scatter_rows(ax, rows, factor, color, marker)
        y_values.extend([value for value in (curves["total"] * factor) if np.isfinite(value) and value > 0.0])
        x_values.extend([max(row.model_time_s, 1.0e-6) for row in rows])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.35)
    if x_values:
        ax.set_xlim(10.0 ** (np.log10(min(x_values)) - 0.08), 10.0 ** (np.log10(max(x_values)) + 0.08))
    if y_values:
        y_arr = np.asarray(y_values, dtype=float)
        ax.set_ylim(10.0 ** (np.log10(np.nanmin(y_arr)) - 0.25), 10.0 ** (np.log10(np.nanmax(y_arr)) + 0.25))
    ax.text(0.015, 0.035, "solid: total; dashed: FS; dash-dot: RS", transform=ax.transAxes, fontsize=8)
    ax.set_ylabel(r"Scaled flux / $\nu F_\nu$")
    ax.set_xlabel(r"$T-T_0$ [s]")
    ax.legend(fontsize=7.1, ncol=1, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=True)
    ax.set_title(f"{TARGET_GRB_LABEL} top-hat FS+RS light-curve fit")
    fig.tight_layout(rect=(0.0, 0.0, 0.78, 1.0))
    fig.savefig(output_prefix.with_suffix(".pdf"))
    fig.savefig(output_prefix.with_suffix(".png"), dpi=180)
    plt.close(fig)
    return output_prefix.with_suffix(".pdf")


def evidence_from_stats(stats: dict[str, Any]) -> dict[str, Any]:
    return {
        "nested_sampling_logz": stats.get("nested sampling global log-evidence", ""),
        "nested_sampling_logz_err": stats.get("nested sampling global log-evidence error", ""),
        "importance_nested_sampling_logz": stats.get("nested importance sampling global log-evidence", stats.get("global evidence", "")),
        "importance_nested_sampling_logz_err": stats.get("nested importance sampling global log-evidence error", stats.get("global evidence error", "")),
    }


def save_model_metric(model_name: str, chi2_best: float, k: int, chain_prefix: str, result_file: str, plot_file: str, **extra) -> dict[str, Any]:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    n = int(cached_observation_grid()["fit_data"]["obs_flux"].size)
    loglike_best = extra.get("loglike_best")
    bic_best = "" if loglike_best is None else float(k * math.log(n) - 2.0 * float(loglike_best))
    row = {
        "model": model_name,
        "chi2_best": float(chi2_best),
        "loglike_best": "" if loglike_best is None else float(loglike_best),
        "k": int(k),
        "n": n,
        "bic_best": bic_best,
        "chain_prefix": chain_prefix,
        "result_file": result_file,
        "plot_file": plot_file,
    }
    row.update(extra)
    rows = []
    if METRICS_FILE.exists():
        with METRICS_FILE.open(newline="") as handle:
            rows = [item for item in csv.DictReader(handle) if item.get("model") != model_name]
    rows.append(row)
    fieldnames = list(row.keys())
    for item in rows:
        for key in item:
            if key not in fieldnames:
                fieldnames.append(key)
    with METRICS_FILE.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return row


def _posterior_candidates(prefix: Path) -> list[Path]:
    candidates = [
        Path(str(prefix) + "post_equal_weights.dat"),
        Path(str(prefix) + "post_equal_weights"),
        Path(str(prefix) + "post_equal_weights."),
        Path(str(prefix) + "post_equal_wei"),
    ]
    candidates.extend(sorted(prefix.parent.glob(prefix.name + "post*")))
    return candidates


def multinest_resume_enabled() -> bool:
    return os.environ.get("MULTINEST_RESUME", "1").strip().lower() not in {"0", "false", "no", "fresh"}


def output_prefix(model_name: str, output_dir: str | Path | None = None) -> Path:
    default_dir = OUTPUT_REL / f"{TARGET_GRB_SLUG}_chains_{model_name}"
    out_dir = Path(output_dir or os.environ.get("MULTINEST_OUTPUT_DIR", default_dir))
    return out_dir / "grb_"


def guard_resume_dimension(prefix: Path, n_params: int) -> None:
    if not multinest_resume_enabled():
        return
    for candidate in _posterior_candidates(prefix):
        if not candidate.exists() or candidate.stat().st_size == 0:
            continue
        try:
            data = np.loadtxt(candidate, ndmin=2, max_rows=1)
        except Exception:
            continue
        if data.shape[1] != n_params + 1:
            raise RuntimeError(
                f"Existing MultiNest output {candidate} has {data.shape[1] - 1} parameters, "
                f"but this run expects {n_params}. Use MULTINEST_RESUME=0 or a fresh output directory."
            )
        return


def ensure_multinest_expected_files(prefix: Path) -> None:
    aliases = {
        Path(str(prefix) + "stats.dat"): prefix.parent.glob(prefix.name + "stats*"),
        Path(str(prefix) + "post_equal_weights.dat"): prefix.parent.glob(prefix.name + "post*"),
    }
    for expected, candidates in aliases.items():
        if expected.exists() and expected.stat().st_size > 0:
            continue
        for candidate in sorted(candidates):
            if candidate != expected and candidate.exists() and candidate.stat().st_size > 0:
                shutil.copyfile(candidate, expected)
                break


def _git_snapshot() -> dict[str, str]:
    import subprocess

    def run_git(args: list[str]) -> str:
        try:
            return subprocess.check_output(["git", *args], cwd=ROOT, text=True, stderr=subprocess.STDOUT).strip()
        except Exception:
            return "not-a-git-repo"

    return {
        "git_head": run_git(["rev-parse", "HEAD"]),
        "git_status": run_git(["status", "--short", "--branch"]),
        "git_diff_stat": run_git(["diff", "--stat"]),
    }


def write_run_info(model_name: str, prefix: Path, extra: dict[str, Any] | None = None) -> Path:
    log_dir = prefix.parent / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    info = {
        "model_name": model_name,
        "electron_solver": MODEL_TO_SOLVER[model_name],
        "output_prefix": str(prefix),
        "start_time": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS", ""),
        "multinest_resume": os.environ.get("MULTINEST_RESUME", ""),
        "fixed_numerical": FIXED_NUMERICAL,
        "xrt_absorption": {
            "NH_gal_1e22_cm2": XRT_NH_GAL_1E22,
            "NH_host_1e22_cm2": XRT_NH_HOST_1E22,
            "z_host": TARGET_REDSHIFT,
            "models": "TbAbs(redshift=0) * TbAbs(redshift=z_host)",
            "abundance_table": "WILM",
        },
        "vhe_ebl": _vhe_ebl_metadata(),
        **_git_snapshot(),
    }
    if extra:
        info.update(extra)
    path = log_dir / "run_info.txt"
    path.write_text(json.dumps(info, indent=2, sort_keys=True), encoding="utf-8")
    return path


def analyze_and_plot(model_name: str, prefix: Path, num_r: int, num_theta: int, postprocess_article_params: bool = False) -> dict[str, Any]:
    from pymultinest.analyse import Analyzer

    names = parameter_names(model_name)
    ensure_multinest_expected_files(prefix)
    analyzer = Analyzer(n_params=len(names), outputfiles_basename=str(prefix))
    stats = None
    try:
        stats = analyzer.get_stats()
    except Exception as exc:
        if not postprocess_article_params:
            raise
        print(f"Warning: MultiNest stats parsing failed in fast smoke mode: {exc}", flush=True)
    evidence = evidence_from_stats(stats) if stats is not None else {}
    results = []
    if stats is None:
        article = parameter_vector_from_article(model_name)
        results = [[float(value), 0.0, 0.0] for value in article]
        for name, value in zip(names, article):
            print(f"{name:24s} {value:.6f} +0.000000 -0.000000", flush=True)
    else:
        for i, name in enumerate(names):
            marginal = stats["marginals"][i]
            median = float(marginal["median"])
            low, high = marginal["1sigma"]
            results.append([median, float(high - median), float(median - low)])
            print(f"{name:24s} {median:.6f} +{high - median:.6f} -{median - low:.6f}", flush=True)
    model_dir = prefix.parent
    result_file = model_dir / f"results_{model_name}.txt"
    best_file = model_dir / f"results_{model_name}_best.txt"
    best_params = parameter_vector_from_article(model_name) if stats is None else np.asarray(stats["modes"][0]["maximum"], dtype=float)
    median_params = np.asarray([row[0] for row in results], dtype=float)
    if postprocess_article_params:
        best_params = parameter_vector_from_article(model_name)
        median_params = best_params.copy()
    np.savetxt(result_file, np.asarray(results), fmt="%.8f", header="median err_plus err_minus for " + " ".join(names))
    np.savetxt(best_file, best_params[None, :], fmt="%.8f", header=" ".join(names))

    posterior = None
    try:
        posterior = analyzer.get_equal_weighted_posterior()[:, :-1]
    except Exception:
        for candidate in _posterior_candidates(prefix):
            if candidate.exists() and candidate.stat().st_size > 0:
                posterior = np.loadtxt(candidate, ndmin=2)[:, :-1]
                break
    max_corner_samples = int(os.environ.get("CORNER_MAX_SAMPLES", "5000"))
    if posterior is not None and posterior.shape[0] > max_corner_samples:
        rng = np.random.default_rng(130427)
        posterior = posterior[rng.choice(posterior.shape[0], max_corner_samples, replace=False)]
    try:
        plt = _load_pyplot()
        import corner

        if posterior is not None and posterior.shape[0] >= 2:
            truths = np.median(posterior, axis=0)
            fig = corner.corner(
                posterior,
                labels=names,
                truths=truths,
                quantiles=[0.16, 0.5, 0.84],
                plot_datapoints=False,
                plot_density=False,
                plot_contours=True,
                fill_contours=True,
                show_titles=True,
                title_fmt=".2f",
                smooth=2.0,
                smooth1d=2.0,
                truth_color="darkorange",
                color="black",
                bins=22,
                range=[0.99] * len(names),
                max_n_ticks=3,
                levels=(0.393, 0.675, 0.864, 0.955),
            )
            fig.savefig(model_dir / f"corner_{model_name}.pdf", dpi=250, bbox_inches="tight", facecolor="white")
            fig.savefig(model_dir / f"corner_{model_name}.png", dpi=180, bbox_inches="tight", facecolor="white")
            plt.close(fig)
    except Exception as exc:
        print(f"Warning: corner plot skipped: {exc}", flush=True)

    best_dict = params_from_vector(best_params, model_name)
    median_dict = params_from_vector(median_params, model_name)
    chi2_best = fit_chi2(model_name, best_dict, num_r=num_r, num_theta=num_theta)
    chi2_median = fit_chi2(model_name, median_dict, num_r=num_r, num_theta=num_theta)
    loglike_best = fit_loglike(model_name, best_dict, num_r=num_r, num_theta=num_theta)
    loglike_median = fit_loglike(model_name, median_dict, num_r=num_r, num_theta=num_theta)
    model = build_tophat_model(model_name, best_dict, num_r=num_r, num_theta=num_theta)
    plot_prefix = model_dir / f"lightcurve_{model_name}_bestfit"
    save_model_cache(model, plot_prefix)
    plot_file = plot_fig3_style(model, plot_prefix)
    metric = save_model_metric(
        f"{TARGET_GRB_SLUG}_tophat_{model_name}",
        chi2_best,
        len(names),
        str(prefix),
        str(result_file),
        str(plot_file),
        chi2_median=chi2_median,
        best_result_file=str(best_file),
        loglike_best=loglike_best,
        loglike_median=loglike_median,
        **evidence,
    )
    (model_dir / f"metric_{model_name}.json").write_text(json.dumps(metric, indent=2), encoding="utf-8")
    derived = {
        "log10_f_sys_best": float(best_dict["log10_f_sys"]),
        "f_sys_best": f_sys_from_params(best_dict),
        "log10_f_sys_median": float(median_dict["log10_f_sys"]),
        "f_sys_median": f_sys_from_params(median_dict),
    }
    (model_dir / f"results_{model_name}_derived.json").write_text(json.dumps(derived, indent=2), encoding="utf-8")
    write_run_info(model_name, prefix, {"end_time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "metric": metric})
    return {
        "model": model_name,
        "model_dir": str(model_dir),
        "prefix": str(prefix),
        "result_file": str(result_file),
        "best_result_file": str(best_file),
        "plot_file": str(plot_file),
        "metric_file": str(model_dir / f"metric_{model_name}.json"),
        "metric": metric,
    }


def mpi_rank() -> int:
    for name in (
        "OMPI_COMM_WORLD_RANK",
        "PMI_RANK",
        "PMIX_RANK",
        "PMI_ID",
        "PMI_PROCESS_RANK",
        "MPI_RANKID",
        "SLURM_PROCID",
    ):
        value = os.environ.get(name)
        if value is not None:
            return int(value)
    return 0


def env_flag(name: str, default: bool = False) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


def is_main_process() -> bool:
    return mpi_rank() == 0


def run_sampling(
    model_name: str,
    num_r: int,
    num_theta: int,
    *,
    n_live_points: int,
    sampling_efficiency: float,
    evidence_tolerance: float,
    max_iter: int,
    output_dir: str | None = None,
    run_tag: str | None = None,
    fast_smoke_likelihood: bool = False,
    postprocess_article_params: bool = False,
    sampling_only: bool = False,
) -> dict[str, Any] | None:
    import pymultinest

    names = parameter_names(model_name)
    if output_dir is None and run_tag:
        output_dir = OUTPUT_REL / f"{TARGET_GRB_SLUG}_chains_{model_name}" / run_tag
    prefix = output_prefix(model_name, output_dir)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    _prepare_astromodels(prefix.parent)
    init_mpi = env_flag("MULTINEST_INIT_MPI", False)
    use_mpi = env_flag("MULTINEST_USE_MPI", True)
    write_run_info(
        model_name,
        prefix,
        {
            "num_r": num_r,
            "num_theta": num_theta,
            "n_live_points": n_live_points,
            "multinest_init_mpi": init_mpi,
            "multinest_use_mpi": use_mpi,
            "multinest_mpi_np": os.environ.get("MULTINEST_MPI_NP", ""),
            "slurm_ntasks": os.environ.get("SLURM_NTASKS", ""),
            "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK", ""),
            "omp_num_threads": os.environ.get("OMP_NUM_THREADS", ""),
        },
    )
    guard_resume_dimension(prefix, len(names))
    if is_main_process():
        print(f"Running {model_name} ({MODEL_TO_SOLVER[model_name]}) {TARGET_GRB_LABEL}", flush=True)
        print(f"Output prefix: {prefix}", flush=True)
    if fast_smoke_likelihood:
        def loglike_func(cube, ndim, nparams):
            return -0.5 * float(np.sum((np.asarray(cube[: len(names)]) - 0.5) ** 2))
    else:
        loglike_func = lambda cube, ndim, nparams: loglike(cube, ndim, nparams, model_name, num_r, num_theta)
    start_sampling = time.perf_counter()
    pymultinest.run(
        LogLikelihood=loglike_func,
        Prior=lambda cube, ndim, nparams: prior(cube, ndim, nparams, model_name),
        n_dims=len(names),
        n_params=len(names),
        n_live_points=n_live_points,
        sampling_efficiency=sampling_efficiency,
        evidence_tolerance=evidence_tolerance,
        outputfiles_basename=str(prefix),
        resume=multinest_resume_enabled(),
        verbose=is_main_process(),
        importance_nested_sampling=True,
        multimodal=False if sampling_only else True,
        max_modes=1 if sampling_only else 100,
        const_efficiency_mode=False,
        n_iter_before_update=100,
        null_log_evidence=-1e90,
        max_iter=max_iter,
        write_output=True,
        log_zero=-1e100,
        init_MPI=init_mpi,
        use_MPI=use_mpi,
    )
    sampling_seconds = time.perf_counter() - start_sampling
    if is_main_process() and sampling_only:
        result = {
            "model": model_name,
            "electron_solver": MODEL_TO_SOLVER[model_name],
            "prefix": str(prefix),
            "sampling_seconds": float(sampling_seconds),
            "n_live_points": int(n_live_points),
            "max_iter": int(max_iter),
            "num_r": int(num_r),
            "num_theta": int(num_theta),
            "omp_num_threads": os.environ.get("OMP_NUM_THREADS", ""),
        }
        (prefix.parent / f"sampling_time_{model_name}.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        write_run_info(model_name, prefix, {"end_time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "sampling_only": True, "sampling_seconds": sampling_seconds})
        return result
    if is_main_process():
        return analyze_and_plot(model_name, prefix, num_r, num_theta, postprocess_article_params=postprocess_article_params)
    return None


def test_interfaces(model_name: str, num_r: int, num_theta: int) -> dict[str, Any]:
    prefix = output_prefix(model_name, OUTPUT_ROOT / f"{TARGET_GRB_SLUG}_interface_{model_name}")
    _prepare_astromodels(prefix.parent)
    names = parameter_names(model_name)
    params = params_from_vector(parameter_vector_from_article(model_name), model_name)
    grid = cached_observation_grid()
    config = config_from_params(model_name, params, num_r=num_r, num_theta=num_theta)
    return {
        "model": model_name,
        "electron_solver": config.electron_solver,
        "n_params": len(names),
        "invalid_physical": _invalid_physical(params),
        "bounds": {name: bound for name, bound in zip(names, prior_bounds(model_name))},
        "num_fit_points": int(grid["fit_data"]["obs_flux"].size),
        "xrt_modes": sorted({row.mode for row in grid["rows"] if row.file == "XRT.txt"}),
        "fixed_numerical": FIXED_NUMERICAL | {"Num_R": int(num_r), "Num_theta": int(num_theta)},
        "xrt_absorption": {
            "NH_gal_1e22_cm2": XRT_NH_GAL_1E22,
            "NH_host_1e22_cm2": XRT_NH_HOST_1E22,
            "z": TARGET_REDSHIFT,
            "models": "TbAbs(redshift=0) * TbAbs(redshift=z_host)",
            "abundance_table": "WILM",
        },
        "vhe_ebl": _vhe_ebl_metadata(),
    }


def benchmark_fit(model_name: str, num_r: int, num_theta: int, repeat: int) -> dict[str, Any]:
    params = params_from_vector(parameter_vector_from_article(model_name), model_name)
    runs = []
    values = []
    stage_totals: dict[str, float] = {}
    for _ in range(max(1, int(repeat))):
        timings: dict[str, float] = {}
        start = time.perf_counter()
        result = evaluate_model_grid(model_name, params, num_r=num_r, num_theta=num_theta, timings=timings)
        loglike_value = loglike_from_model_flux(result["model_fit_flux"], result["grid"]["fit_data"], f_sys_from_params(params))
        elapsed = time.perf_counter() - start
        runs.append(elapsed)
        values.append(loglike_value)
        for name, value in timings.items():
            stage_totals[name] = stage_totals.get(name, 0.0) + float(value)
    return {
        "model": model_name,
        "num_r": int(num_r),
        "num_theta": int(num_theta),
        "repeat": len(runs),
        "loglike": float(values[-1]),
        "mean_seconds": float(np.mean(runs)),
        "min_seconds": float(np.min(runs)),
        "max_seconds": float(np.max(runs)),
        "stage_mean_seconds": {name: value / len(runs) for name, value in sorted(stage_totals.items())},
    }


def spectrum_test(model_name: str, num_r: int, num_theta: int, output_dir: str | None = None) -> dict[str, Any]:
    params = params_from_vector(parameter_vector_from_article(model_name), model_name)
    out_dir = Path(output_dir) if output_dir else OUTPUT_ROOT / "grb180720b_spectrum_tests"
    out_dir.mkdir(parents=True, exist_ok=True)
    config = config_from_params(model_name, params, num_r=num_r, num_theta=num_theta)
    times_s = np.asarray([100.0, 1000.0, 1.0e4], dtype=float)
    freq_hz = np.logspace(9.0, 27.0, 160)
    setup = build_simulation_setup(config)
    setup.observer_time_s = times_s
    state = solve_state_from_setup(config, setup, requested_frequencies_hz=freq_hz)
    observed = project_flux_grid(state, times_s, freq_hz, mode="total_only")
    fnu = np.asarray(observed.components["total"], dtype=float)
    nufnu = fnu * freq_hz[:, None]
    if not np.all(np.isfinite(fnu)) or not np.any(nufnu > 0.0):
        raise RuntimeError(f"{model_name} spectrum test did not produce finite positive flux.")
    path = out_dir / f"spectrum_{model_name}.npz"
    np.savez(path, times_s=times_s, freq_hz=freq_hz, fnu=fnu, nufnu=nufnu, model=model_name, solver=MODEL_TO_SOLVER[model_name])
    plt = _load_pyplot()
    fig, ax = plt.subplots(figsize=(6.5, 4.8), dpi=180)
    for i, time_s in enumerate(times_s):
        ax.loglog(freq_hz, _finite_positive(nufnu[:, i]), label=f"{time_s:g} s")
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel(r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
    ax.set_title(f"{TARGET_GRB_LABEL} {model_name} test SED")
    ax.grid(True, which="both", ls=":", alpha=0.35)
    ax.legend()
    png = out_dir / f"spectrum_{model_name}.png"
    fig.savefig(png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return {"model": model_name, "solver": MODEL_TO_SOLVER[model_name], "npz": str(path), "png": str(png), "shape": list(fnu.shape)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=f"{TARGET_GRB_LABEL} three-model MultiNest driver.")
    parser.add_argument("model", choices=["fullhide_1d", "fullhide_1d_pic", "fullhide_2d_pic"])
    parser.add_argument("--num-r", type=int, default=int(os.environ.get("ARTICLE_NUM_R", str(FIXED_NUMERICAL["Num_R"]))))
    parser.add_argument("--num-theta", type=int, default=FIXED_NUMERICAL["Num_theta"])
    parser.add_argument("--n-live", type=int, default=int(os.environ.get("MULTINEST_N_LIVE", "500")))
    parser.add_argument("--max-iter", type=int, default=int(os.environ.get("MULTINEST_MAX_ITER", "0")))
    parser.add_argument("--sampling-efficiency", type=float, default=float(os.environ.get("MULTINEST_SAMPLING_EFFICIENCY", "0.3")))
    parser.add_argument("--evidence-tolerance", type=float, default=float(os.environ.get("MULTINEST_EVIDENCE_TOLERANCE", "0.3")))
    parser.add_argument("--output-dir")
    parser.add_argument("--run-tag")
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--fast-smoke-likelihood", action="store_true")
    parser.add_argument("--interface-test", action="store_true")
    parser.add_argument("--benchmark-fit", action="store_true")
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--spectrum-test", action="store_true")
    parser.add_argument("--sampling-only", action="store_true", help="Stop after MultiNest sampling and skip Analyzer, plots, and model postprocess.")
    args = parser.parse_args(argv)

    if args.num_theta != FIXED_NUMERICAL["Num_theta"]:
        raise SystemExit("Num_theta is fixed at 300 for GRB 180720B article-style runs.")
    if args.interface_test:
        print(json.dumps(test_interfaces(args.model, args.num_r, args.num_theta), indent=2))
        return 0
    if args.benchmark_fit:
        print(json.dumps(benchmark_fit(args.model, args.num_r, args.num_theta, args.repeat), indent=2))
        return 0
    if args.spectrum_test:
        print(json.dumps(spectrum_test(args.model, args.num_r, args.num_theta, args.output_dir), indent=2))
        return 0
    if args.smoke_test:
        os.environ.setdefault("MULTINEST_RESUME", "0")
        os.environ.setdefault("CORNER_MAX_SAMPLES", "200")
        args.num_r = min(args.num_r, 20)
        args.n_live = min(args.n_live, max(8, len(parameter_names(args.model)) + 2))
        args.max_iter = args.max_iter if args.max_iter > 0 else 1
        args.sampling_efficiency = max(args.sampling_efficiency, 0.8)
        args.evidence_tolerance = max(args.evidence_tolerance, 10.0)
        args.run_tag = args.run_tag or f"smoke_{args.model}"
        if args.fast_smoke_likelihood:
            args.max_iter = min(args.max_iter, 1)
    result = run_sampling(
        args.model,
        args.num_r,
        args.num_theta,
        n_live_points=args.n_live,
        sampling_efficiency=args.sampling_efficiency,
        evidence_tolerance=args.evidence_tolerance,
        max_iter=args.max_iter,
        output_dir=args.output_dir,
        run_tag=args.run_tag,
        fast_smoke_likelihood=args.fast_smoke_likelihood,
        postprocess_article_params=args.fast_smoke_likelihood,
        sampling_only=args.sampling_only,
    )
    if result is not None:
        print(json.dumps(result, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
