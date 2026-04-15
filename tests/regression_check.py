from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_models import FitConfig
from asgard_component_backend import solve_model_state_from_setup, extract_physical_solution_from_state
from asgard_observables import OUTPUT_BANDS, build_multiband_observer_frequencies, combine_multiband_flux
from asgard_postprocess import compute_observed_flux_matrix, compute_light_curve_redchi
from asgard_presets import build_baseline_config
from asgard_setup import build_simulation_setup


BANDS_TO_CHECK = ("xrt", "optr", "9GHz")


def build_reference_config():
    return build_baseline_config()


def _compute_bands_flux(config: FitConfig):
    setup = build_simulation_setup(config)
    state = solve_model_state_from_setup(config, setup)
    physical = extract_physical_solution_from_state(state)
    num_xrt, band_freqs = build_multiband_observer_frequencies()
    flux_matrix = compute_observed_flux_matrix(setup, physical, band_freqs, config)
    return combine_multiband_flux(flux_matrix, band_freqs, num_xrt)


def extract_bands(bands_flux, band_names: tuple[str, ...]) -> np.ndarray:
    idx = [OUTPUT_BANDS.index(name) for name in band_names]
    return bands_flux[idx]


def write_baseline(path: Path, t_obs_s: np.ndarray, band_flux: np.ndarray, band_names: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, t_obs_s=t_obs_s, bands=np.array(band_names), band_flux=band_flux)


def check_against_baseline(path: Path, threshold: float) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Baseline file not found: {path}")

    data = np.load(path, allow_pickle=False)
    baseline_t = data["t_obs_s"]
    baseline_flux = data["band_flux"]
    baseline_bands = tuple(str(x) for x in data["bands"])

    config = build_reference_config()
    setup = build_simulation_setup(config)
    bands_flux = _compute_bands_flux(config)

    current_flux = extract_bands(bands_flux, baseline_bands)
    for i, band in enumerate(baseline_bands):
        base = baseline_flux[i]
        cur = current_flux[i]

        if np.any(base == 0.0):
            if not np.array_equal(cur[base == 0.0], base[base == 0.0]):
                raise ValueError(f"Baseline contains zero values in {band}, and current values differ on those points.")
            mask = base != 0.0
            rel = np.abs((cur[mask] - base[mask]) / base[mask])
        else:
            rel = np.abs((cur - base) / base)

        max_rel = float(np.max(rel))
        print(f"{band}: max relative deviation = {max_rel:.6f}")
        if max_rel > threshold:
            raise AssertionError(f"{band} exceeds threshold {threshold:.3f}: {max_rel:.6f}")

    print(f"PASS: all checked bands are within threshold {threshold:.3f}.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Regression check for xrt/optr/9GHz light curves.")
    parser.add_argument("--baseline", type=Path, default=Path("tests") / "baseline_lightcurves.npz")
    parser.add_argument("--threshold", type=float, default=0.03)
    parser.add_argument("--write-baseline", action="store_true")
    args = parser.parse_args()

    config = build_reference_config()
    setup = build_simulation_setup(config)
    bands_flux = _compute_bands_flux(config)

    if args.write_baseline:
        band_flux = extract_bands(bands_flux, BANDS_TO_CHECK)
        write_baseline(args.baseline, setup.observer_time_s, band_flux, BANDS_TO_CHECK)
        print(f"Baseline written: {args.baseline}")
        return

    check_against_baseline(args.baseline, args.threshold)


if __name__ == "__main__":
    main()
