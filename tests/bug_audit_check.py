from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass
from pathlib import Path
import json
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_component_backend import build_query_setup, observe_spectra_from_setup, solve_component_spectra_from_setup
from asgard_models import default_num_threads
from asgard_presets import build_baseline_config, build_reverse_demo_config
from asgard_runtime import solve_dynamics, solve_electron
from asgard_setup import build_simulation_setup
from src import constants


MULTI_THREADS = min(8, max(2, default_num_threads()))
DIRECT_TIMES_S = np.logspace(2.0, 6.0, 12)
REVERSE_TIMES_S = np.logspace(1.0, 5.0, 12)
OBS_FREQS_HZ = np.array([9.0e9, 4.84e14, 1.0e18], dtype=float)
THREAD_TOL = 1.0e-7
COMPONENT_TOL = 1.0e-10


@dataclass
class AuditResult:
    name: str
    passed: bool
    extra: dict


def _max_rel(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    scale = np.maximum(np.abs(b), 1.0e-99)
    return float(np.max(np.abs(a - b) / scale))


def _require(name: str, condition: bool, **extra) -> AuditResult:
    return AuditResult(name=name, passed=bool(condition), extra=extra)


def _component_sum(observed: dict[str, np.ndarray | None]) -> np.ndarray:
    total = np.array(observed["fwd_sync"], dtype=float, copy=True)
    total = total + np.array(observed["fwd_ssc"], dtype=float, copy=False)
    if observed["rev_sync"] is not None:
        total = total + np.array(observed["rev_sync"], dtype=float, copy=False)
    if observed["rev_ssc"] is not None:
        total = total + np.array(observed["rev_ssc"], dtype=float, copy=False)
    if observed["cross_ic"] is not None:
        total = total + np.array(observed["cross_ic"], dtype=float, copy=False)
    return total


def _run_component_case(config, times_s: np.ndarray, freqs_hz: np.ndarray):
    setup = build_query_setup(config, times_s, freqs_hz)
    component_spectra = solve_component_spectra_from_setup(config, setup)
    observed = observe_spectra_from_setup(config, component_spectra, setup, freqs_hz)
    return setup, component_spectra, observed


def _ambient_density_scalar(radius_cm: float, config) -> float:
    if config.a_star > 0.0:
        d_ne_wind = config.a_star * 3.0e35 / radius_cm**2
        if d_ne_wind <= config.d_ne / 4.0:
            density = config.d_ne
        else:
            density = d_ne_wind
    else:
        density = config.d_ne * (
            1.0
            + (config.f_jump - 1.0)
            * np.exp(-(np.log10(radius_cm) - np.log10(config.r_tr)) ** 2 / (2.0 * config.f_wide**2))
        )
    if config.a_star > 0.0 and radius_cm < config.r0:
        density = config.a_star * 3.0e35 / config.r0**2
    return float(density)


def _tau_comoving(radius_cm: float, db: float, gam_e: np.ndarray, d_n_gam_e: np.ndarray, nu_hz: float) -> float:
    factor = (3.62 / np.pi) ** 2
    d_n1 = d_n_gam_e / (gam_e * gam_e)
    dd_n = d_n1[:-1] - d_n1[1:]
    gam_e_mean2 = ((gam_e[:-1] + gam_e[1:]) ** 2) / 4.0
    vc = 4.2e6 * gam_e_mean2 * db
    x = nu_hz / vc
    fx = 1.81 * np.exp(-x) / np.sqrt(x ** (-2.0 / 3.0) + factor)
    tau = np.sum(gam_e_mean2 * dd_n * fx)
    return float(1.046e4 * tau * db / (4.0 * np.pi * radius_cm * radius_cm * nu_hz * nu_hz))


def check_thread_consistency() -> list[AuditResult]:
    results: list[AuditResult] = []
    for name, cfg, times in [
        ("forward_sync", build_baseline_config(include_forward_ssc=False, num_r=120, num_nu=121, num_gam_e=121), DIRECT_TIMES_S),
        ("forward_ssc", build_baseline_config(include_forward_ssc=True, num_r=120, num_nu=121, num_gam_e=121), DIRECT_TIMES_S),
        ("reverse_sync", build_reverse_demo_config(include_forward_ssc=False, num_r=120, num_nu=121, num_gam_e=121), REVERSE_TIMES_S),
    ]:
        cfg_serial = cfg
        cfg_serial.num_threads = 1
        _, _, obs_serial = _run_component_case(cfg_serial, times, OBS_FREQS_HZ)

        cfg_parallel = deepcopy(cfg)
        cfg_parallel.num_threads = MULTI_THREADS
        _, _, obs_parallel_a = _run_component_case(cfg_parallel, times, OBS_FREQS_HZ)
        _, _, obs_parallel_b = _run_component_case(cfg_parallel, times, OBS_FREQS_HZ)

        serial_vs_parallel = _max_rel(obs_parallel_a["total"], obs_serial["total"])
        parallel_repeat = _max_rel(obs_parallel_a["total"], obs_parallel_b["total"])
        results.append(
            _require(
                f"thread-{name}",
                serial_vs_parallel <= THREAD_TOL and parallel_repeat <= THREAD_TOL,
                serial_vs_parallel=serial_vs_parallel,
                parallel_repeat=parallel_repeat,
                tolerance=THREAD_TOL,
            )
        )
    return results


def check_nu_a_bounds() -> list[AuditResult]:
    config = build_baseline_config(include_forward_ssc=False, num_r=140, num_nu=121, num_gam_e=121)
    setup = build_simulation_setup(config)
    dynamics = solve_dynamics(setup.boundary, config)
    electron = solve_electron(setup.boundary, dynamics, setup.seed_frequency_hz, config)

    results: list[AuditResult] = []
    sample_indices = np.unique(np.linspace(0, config.num_r - 2, 16, dtype=int))
    for idx in sample_indices:
        radius_cm = float(dynamics.radius[idx])
        gamma_bulk = 0.5 * float(dynamics.r_gamma[idx] + dynamics.r_gamma[idx + 1])
        beta_bulk = np.sqrt(1.0 - gamma_bulk**-2)
        doppler_den = gamma_bulk * (1.0 - beta_bulk) * (1.0 + config.z)
        density = _ambient_density_scalar(radius_cm, config)
        db = 0.39 * np.sqrt(config.epsilon_b * density * (gamma_bulk * (gamma_bulk - 1.0)))
        gam_e_mean2 = ((electron.gam_e[:-1] + electron.gam_e[1:]) ** 2) / 4.0
        v_a_cap = max(1.0e14, min(1.0e30, 10.0 * float(np.max(4.2e6 * gam_e_mean2 * db))))
        nu_a_comoving = float(electron.nu_a[idx] * doppler_den)
        tau_floor = _tau_comoving(radius_cm, db, electron.gam_e, electron.d_n_gam_e[:, idx], 1.0e4)

        lower_bound = 1.0e4 * (1.0 - 1.0e-12)
        passed = lower_bound <= nu_a_comoving <= v_a_cap * (1.0 + 1.0e-9)
        extra = {
            "shell_index": int(idx),
            "nu_a_comoving_hz": nu_a_comoving,
            "nu_a_cap_hz": v_a_cap,
            "tau_floor": tau_floor,
            "nu_a_lower_bound_hz": lower_bound,
        }
        if nu_a_comoving >= 0.999 * v_a_cap and tau_floor > 1.0:
            tau_cap = _tau_comoving(radius_cm, db, electron.gam_e, electron.d_n_gam_e[:, idx], v_a_cap)
            passed = passed and tau_cap <= 1.0
            extra["tau_cap"] = tau_cap
        results.append(_require(f"nu_a-bound-{idx}", passed, **extra))
    return results


def check_monotonic_arrays() -> list[AuditResult]:
    config = build_baseline_config(num_r=120, num_nu=121, num_gam_e=121)
    setup = build_simulation_setup(config)
    dynamics = solve_dynamics(setup.boundary, config)
    electron = solve_electron(setup.boundary, dynamics, setup.seed_frequency_hz, config)

    checks = [
        ("observer-time-grid", np.all(np.diff(setup.observer_time_s) > 0.0), {"min_diff": float(np.min(np.diff(setup.observer_time_s)))}),
        ("seed-frequency-grid", np.all(np.diff(setup.seed_frequency_hz) > 0.0), {"min_diff": float(np.min(np.diff(setup.seed_frequency_hz)))}),
        ("dynamics-r_tobs", np.all(np.diff(dynamics.r_tobs) > 0.0), {"min_diff": float(np.min(np.diff(dynamics.r_tobs)))}),
        ("dynamics-radius", np.all(np.diff(dynamics.radius) > 0.0), {"min_diff": float(np.min(np.diff(dynamics.radius)))}),
        ("electron-gamma-grid", np.all(np.diff(electron.gam_e) > 0.0), {"min_diff": float(np.min(np.diff(electron.gam_e)))}),
    ]
    return [_require(f"monotonic-{name}", cond, **extra) for name, cond, extra in checks]


def check_component_finiteness_and_closure() -> list[AuditResult]:
    results: list[AuditResult] = []
    cases = [
        ("forward_ssc", build_baseline_config(include_forward_ssc=True, num_r=120, num_nu=121, num_gam_e=121), DIRECT_TIMES_S),
        ("reverse_sync", build_reverse_demo_config(include_forward_ssc=False, num_r=120, num_nu=121, num_gam_e=121), REVERSE_TIMES_S),
    ]
    for name, config, times in cases:
        _, component_spectra, observed = _run_component_case(config, times, OBS_FREQS_HZ)
        raw_total = np.array(component_spectra.fwd_sync, dtype=float, copy=True)
        raw_total += np.array(component_spectra.fwd_ssc, dtype=float, copy=False)
        if component_spectra.rev_sync is not None:
            raw_total += np.array(component_spectra.rev_sync, dtype=float, copy=False)
        if component_spectra.rev_ssc is not None:
            raw_total += np.array(component_spectra.rev_ssc, dtype=float, copy=False)
        if component_spectra.cross_ic is not None:
            raw_total += np.array(component_spectra.cross_ic, dtype=float, copy=False)
        results.append(
            _require(
                f"{name}-raw-component-sum",
                _max_rel(component_spectra.total, raw_total) <= COMPONENT_TOL,
                max_relative_error=_max_rel(component_spectra.total, raw_total),
                tolerance=COMPONENT_TOL,
            )
        )
        for label, arr in [
            ("raw-total", component_spectra.total),
            ("raw-fwd-sync", component_spectra.fwd_sync),
            ("raw-fwd-ssc", component_spectra.fwd_ssc),
            ("obs-total", observed["total"]),
            ("obs-fwd-sync", observed["fwd_sync"]),
            ("obs-fwd-ssc", observed["fwd_ssc"]),
        ]:
            if arr is None:
                continue
            arr_np = np.asarray(arr, dtype=float)
            results.append(
                _require(
                    f"{name}-{label}-finite",
                    np.all(np.isfinite(arr_np)) and np.all(arr_np >= 0.0),
                    min_value=float(np.min(arr_np)),
                    max_value=float(np.max(arr_np)),
                )
            )
        summed = _component_sum(observed)
        results.append(
            _require(
                f"{name}-observed-component-sum",
                _max_rel(observed["total"], summed) <= COMPONENT_TOL,
                max_relative_error=_max_rel(observed["total"], summed),
                tolerance=COMPONENT_TOL,
            )
        )
    return results


def main() -> None:
    checks = (
        check_thread_consistency()
        + check_nu_a_bounds()
        + check_monotonic_arrays()
        + check_component_finiteness_and_closure()
    )
    payload = {
        "summary": {
            "total": len(checks),
            "failed": sum(0 if item.passed else 1 for item in checks),
        },
        "checks": [asdict(item) for item in checks],
    }
    print(json.dumps(payload["summary"], indent=2))
    failed = [item for item in checks if not item.passed]
    if failed:
        for item in failed:
            print(json.dumps(asdict(item), ensure_ascii=False))
        raise SystemExit(1)
    print("PASS: bug audit check succeeded.")


if __name__ == "__main__":
    main()
