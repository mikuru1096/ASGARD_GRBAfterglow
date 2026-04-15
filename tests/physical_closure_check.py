from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import json
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet, units
from asgard_component_backend import build_query_setup, observe_spectra_from_setup, solve_component_spectra_from_setup
from asgard_presets import build_baseline_config, build_reverse_demo_config
from asgard_runtime import solve_dynamics, solve_electron, solve_reverse_shock_emission
from asgard_setup import build_simulation_setup
from src import Radiation as RadiationKernel, constants
from tests.asgard_comprehensive_validation import (
    REGIME_TEST_CONFIGS,
    REGRESSION_RESOLUTION,
    SPECTRAL_TOLERANCE,
    _detect_regime,
    _make_model,
    _run_spectral_regimes,
)


DIRECT_TIMES_S = np.logspace(2.0, 6.0, 10)
REVERSE_TIMES_S = np.logspace(1.0, 5.0, 10)
OBS_FREQS_HZ = np.array([9.0e9, 4.84e14, 1.0e18], dtype=float)
PAIR_TIMES_S = np.array([1.0e6, 3.0e6], dtype=float)
PAIR_FREQS_HZ = np.array([1.0e9, 1.0e9], dtype=float)
PAIR_GRID_FREQS_HZ = np.unique(PAIR_FREQS_HZ)

REGIME_TOL = 0.0
DYN_FREQ_TOL = 5.0e-3
RAW_CLOSURE_TOL = 1.0e-10
OBS_CLOSURE_TOL = 1.0e-10
SYNC_STABILITY_TOL = 1.0e-12
SKY_IMAGE_TOL = 6.0e-2
RS_PEAK_TOL_DEX = 1.2


@dataclass
class ClosureResult:
    name: str
    category: str
    passed: bool
    expected: float | None = None
    measured: float | None = None
    tolerance: float | None = None
    extra: dict | None = None


def _max_rel(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    scale = np.maximum(np.abs(b), 1.0e-99)
    return float(np.max(np.abs(a - b) / scale))


def _component_sum(observed: dict[str, np.ndarray | None]) -> np.ndarray:
    total = np.array(observed["fwd_sync"], dtype=float, copy=True)
    total += np.array(observed["fwd_ssc"], dtype=float, copy=False)
    if observed["rev_sync"] is not None:
        total += np.array(observed["rev_sync"], dtype=float, copy=False)
    if observed["rev_ssc"] is not None:
        total += np.array(observed["rev_ssc"], dtype=float, copy=False)
    if observed["cross_ic"] is not None:
        total += np.array(observed["cross_ic"], dtype=float, copy=False)
    return total


def _ambient_density_scalar(radius_cm: float, config) -> float:
    if config.a_star > 0.0:
        d_ne_wind = config.a_star * 3.0e35 / radius_cm**2
        density = config.d_ne if d_ne_wind <= config.d_ne / 4.0 else d_ne_wind
    else:
        density = config.d_ne * (
            1.0
            + (config.f_jump - 1.0)
            * np.exp(-(np.log10(radius_cm) - np.log10(config.r_tr)) ** 2 / (2.0 * config.f_wide**2))
        )
    if config.a_star > 0.0 and radius_cm < config.r0:
        density = config.a_star * 3.0e35 / config.r0**2
    return float(density)


def _extract_pair_from_grid(grid: np.ndarray, times_s: np.ndarray, frequencies_hz: np.ndarray, unique_freqs_hz: np.ndarray) -> np.ndarray:
    inverse = np.searchsorted(unique_freqs_hz, frequencies_hz)
    return np.asarray(grid[inverse, np.arange(times_s.shape[0])], dtype=float)


def _base_model() -> Model:
    return Model(
        medium=ISM(n0=0.1),
        jet=TophatJet(E_iso=1.0e53, Gamma0=100.0, theta_c=0.1),
        observer=Observer(lumi_dist=1.0e28, z=0.4, theta_obs=0.0, phi_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.5, xi_N=0.1, ssc=False, kn=False),
        setups=Setups(
            num_threads=8,
            num_gam_e=121,
            num_nu=121,
            num_r=140,
            num_theta=120,
            num_phi=1,
            num_tobs=80,
            observer_time_min_s=float(DIRECT_TIMES_S.min()),
            observer_time_max_s=float(DIRECT_TIMES_S.max()),
        ),
        resolutions=REGRESSION_RESOLUTION,
    )


def check_regime_ordering_closure() -> list[ClosureResult]:
    results: list[ClosureResult] = []
    for regime, cfg in REGIME_TEST_CONFIGS.items():
        model = _make_model(
            cfg["medium"],
            resolutions=REGRESSION_RESOLUTION,
            eps_e=cfg.get("eps_e", 0.1),
            eps_B=cfg.get("eps_B", 1.0e-3),
            xi_e=cfg.get("xi_e", 0.1),
            n_ism=cfg.get("n_ism", 0.1),
            A_star=cfg.get("A_star", 0.1),
        )
        t_eval = float(cfg["t"])
        details = model.details(t_eval * 0.9, t_eval * 1.1)
        idx = int(np.argmin(np.abs(np.asarray(details.fwd.t_obs) - t_eval)))
        nu_a = float(details.fwd.nu_a[idx])
        nu_m = float(details.fwd.nu_m[idx])
        nu_c = float(details.fwd.nu_c[idx])
        actual = _detect_regime(nu_a, nu_m, nu_c)
        results.append(
            ClosureResult(
                name=f"regime-order-{regime}",
                category="regime-order",
                passed=actual == regime,
                expected=None,
                measured=None,
                tolerance=REGIME_TOL,
                extra={
                    "expected_regime": regime,
                    "actual_regime": actual,
                    "nu_a_hz": nu_a,
                    "nu_m_hz": nu_m,
                    "nu_c_hz": nu_c,
                },
            )
        )
    return results


def check_spectral_slope_closure() -> list[ClosureResult]:
    imported_results = _run_spectral_regimes()
    return [
        ClosureResult(
            name=item.name,
            category="spectrum",
            passed=item.passed,
            expected=item.expected,
            measured=item.measured,
            tolerance=item.tolerance,
            extra=item.extra,
        )
        for item in imported_results
    ]


def check_dynamics_radiation_closure() -> list[ClosureResult]:
    config = build_baseline_config(
        include_forward_ssc=False,
        index_y=0,
        num_r=160,
        num_nu=121,
        num_gam_e=121,
    )
    query_setup = build_query_setup(config, DIRECT_TIMES_S, OBS_FREQS_HZ)
    dynamics = solve_dynamics(query_setup.boundary, config)
    electron = solve_electron(query_setup.boundary, dynamics, query_setup.seed_frequency_hz, config)
    component_spectra = solve_component_spectra_from_setup(config, query_setup)

    sample_indices = np.unique(np.linspace(8, config.num_r - 3, 12, dtype=int))
    nu_m_ref = []
    nu_c_ref = []
    nu_M_ref = []
    nu_m_mod = []
    nu_c_mod = []
    nu_M_mod = []
    for idx in sample_indices:
        gamma = float(0.5 * (dynamics.r_gamma[idx] + dynamics.r_gamma[idx + 1]))
        radius_cm = float(dynamics.radius[idx])
        beta = np.sqrt(1.0 - gamma**-2)
        doppler = 1.0 / (gamma * (1.0 - beta) * (1.0 + config.z))
        density = _ambient_density_scalar(radius_cm, config)
        magnetic_field = 0.39 * np.sqrt(config.epsilon_b * density * gamma * max(gamma - 1.0, 0.0))
        temp_gam = config.epsilon_e / config.f_e * constants.para_m_p_div_m_e * (gamma - 1.0)
        gam_e_m = (config.p - 2.0) / (config.p - 1.0) * temp_gam + 1.0
        gam_e_c = 7.7e8 * (1.0 + config.z) / (gamma * magnetic_field * magnetic_field * dynamics.r_tobs[idx + 1])
        gam_e_M = 3.0 * constants.para_m_energy / np.sqrt(8.0 * magnetic_field * constants.para_e**3)
        nu_m_ref.append(4.2e6 * magnetic_field * gam_e_m**2 * doppler)
        nu_c_ref.append(4.2e6 * magnetic_field * gam_e_c**2 * doppler)
        gamma_nu_M = float(dynamics.r_gamma[idx])
        beta_nu_M = np.sqrt(1.0 - gamma_nu_M**-2)
        doppler_nu_M = 1.0 / (gamma_nu_M * (1.0 - beta_nu_M) * (1.0 + config.z))
        magnetic_field_nu_M = 0.39 * np.sqrt(
            config.epsilon_b * density * gamma_nu_M * max(gamma_nu_M - 1.0, 0.0)
        )
        gam_e_M_nu_M = 3.0 * constants.para_m_energy / np.sqrt(8.0 * magnetic_field_nu_M * constants.para_e**3)
        nu_M_ref.append(4.2e6 * magnetic_field_nu_M * gam_e_M_nu_M**2 * doppler_nu_M)
        nu_m_mod.append(float(electron.nu_m[idx]))
        nu_c_mod.append(float(electron.nu_c[idx]))
        nu_M_mod.append(float(component_spectra.fwd.nu_M[idx]))

    results = [
        ClosureResult(
            name="dynamics-radiation-nu_m",
            category="dynamics-radiation",
            passed=_max_rel(np.array(nu_m_mod), np.array(nu_m_ref)) <= DYN_FREQ_TOL,
            expected=0.0,
            measured=_max_rel(np.array(nu_m_mod), np.array(nu_m_ref)),
            tolerance=DYN_FREQ_TOL,
            extra={"indices": sample_indices.tolist()},
        ),
        ClosureResult(
            name="dynamics-radiation-nu_c",
            category="dynamics-radiation",
            passed=_max_rel(np.array(nu_c_mod), np.array(nu_c_ref)) <= DYN_FREQ_TOL,
            expected=0.0,
            measured=_max_rel(np.array(nu_c_mod), np.array(nu_c_ref)),
            tolerance=DYN_FREQ_TOL,
            extra={"indices": sample_indices.tolist()},
        ),
        ClosureResult(
            name="dynamics-radiation-nu_M",
            category="dynamics-radiation",
            passed=_max_rel(np.array(nu_M_mod), np.array(nu_M_ref)) <= DYN_FREQ_TOL,
            expected=0.0,
            measured=_max_rel(np.array(nu_M_mod), np.array(nu_M_ref)),
            tolerance=DYN_FREQ_TOL,
            extra={"indices": sample_indices.tolist()},
        ),
    ]
    return results


def check_absorption_emission_closure() -> list[ClosureResult]:
    config = build_baseline_config(
        include_forward_ssc=True,
        index_y=0,
        num_r=140,
        num_nu=121,
        num_gam_e=121,
    )
    setup = build_query_setup(config, DIRECT_TIMES_S, OBS_FREQS_HZ)
    dynamics = solve_dynamics(setup.boundary, config)
    electron = solve_electron(setup.boundary, dynamics, setup.seed_frequency_hz, config)
    l_ssc, seed_ssc = RadiationKernel.ssc_spec(
        dynamics.radius,
        electron.gam_e,
        electron.d_n_gam_e,
        setup.seed_frequency_hz,
        electron.seed_syn,
        config.num_threads,
    )
    zero_component = np.zeros_like(electron.seed_syn)
    absorption = RadiationKernel.annihilation(
        dynamics.r_gamma,
        dynamics.radius,
        setup.seed_frequency_hz,
        electron.seed_syn,
        seed_ssc if config.include_forward_ssc else zero_component,
        config.num_threads,
    )
    prefactor = absorption / (4.0 * np.pi * setup.luminosity_distance_cm**2) * (1.0 + config.z)
    component_spectra = solve_component_spectra_from_setup(config, setup)

    results = [
        ClosureResult(
            name="absorption-factor-range",
            category="absorption-emission",
            passed=bool(np.all(np.isfinite(absorption)) and np.all(absorption >= 0.0) and np.all(absorption <= 1.0 + 1.0e-12)),
            expected=1.0,
            measured=float(np.max(absorption)),
            tolerance=1.0e-12,
            extra={"min_absorption": float(np.min(absorption)), "max_absorption": float(np.max(absorption))},
        ),
        ClosureResult(
            name="forward-sync-prefactor-closure",
            category="absorption-emission",
            passed=_max_rel(component_spectra.fwd_sync, electron.l_syn_spec * prefactor) <= RAW_CLOSURE_TOL,
            expected=0.0,
            measured=_max_rel(component_spectra.fwd_sync, electron.l_syn_spec * prefactor),
            tolerance=RAW_CLOSURE_TOL,
        ),
        ClosureResult(
            name="forward-ssc-prefactor-closure",
            category="absorption-emission",
            passed=_max_rel(component_spectra.fwd_ssc, l_ssc * prefactor) <= RAW_CLOSURE_TOL,
            expected=0.0,
            measured=_max_rel(component_spectra.fwd_ssc, l_ssc * prefactor),
            tolerance=RAW_CLOSURE_TOL,
        ),
        ClosureResult(
            name="absorption-total-closure",
            category="absorption-emission",
            passed=_max_rel(component_spectra.total, component_spectra.fwd_sync + component_spectra.fwd_ssc) <= RAW_CLOSURE_TOL,
            expected=0.0,
            measured=_max_rel(component_spectra.total, component_spectra.fwd_sync + component_spectra.fwd_ssc),
            tolerance=RAW_CLOSURE_TOL,
        ),
    ]
    return results


def check_forward_reverse_closure() -> list[ClosureResult]:
    results: list[ClosureResult] = []
    cases = [
        ("fs-only", build_baseline_config(include_forward_ssc=False, num_r=120, num_nu=121, num_gam_e=121), DIRECT_TIMES_S),
        ("fs-rs", build_reverse_demo_config(include_forward_ssc=False, num_r=120, num_nu=121, num_gam_e=121), REVERSE_TIMES_S),
        ("fs-rs-ssc", build_reverse_demo_config(include_forward_ssc=False, num_r=120, num_nu=121, num_gam_e=121), REVERSE_TIMES_S),
        ("fs-cross-ic", build_reverse_demo_config(include_forward_ssc=False, num_r=120, num_nu=121, num_gam_e=121), REVERSE_TIMES_S),
    ]
    cases[2][1].reverse_shock.include_ssc = True
    cases[3][1].reverse_shock.include_cross_zone_ic = True

    observed_map: dict[str, dict[str, np.ndarray | None]] = {}
    raw_map: dict[str, object] = {}
    for name, config, times in cases:
        setup = build_query_setup(config, times, OBS_FREQS_HZ)
        component_spectra = solve_component_spectra_from_setup(config, setup)
        observed = observe_spectra_from_setup(config, component_spectra, setup, OBS_FREQS_HZ)
        observed_map[name] = observed
        raw_map[name] = component_spectra
        results.append(
            ClosureResult(
                name=f"{name}-raw-total-closure",
                category="forward-reverse",
                passed=_max_rel(component_spectra.total, _raw_component_sum(component_spectra)) <= RAW_CLOSURE_TOL,
                expected=0.0,
                measured=_max_rel(component_spectra.total, _raw_component_sum(component_spectra)),
                tolerance=RAW_CLOSURE_TOL,
            )
        )
        results.append(
            ClosureResult(
                name=f"{name}-observed-total-closure",
                category="forward-reverse",
                passed=_max_rel(observed["total"], _component_sum(observed)) <= OBS_CLOSURE_TOL,
                expected=0.0,
                measured=_max_rel(observed["total"], _component_sum(observed)),
                tolerance=OBS_CLOSURE_TOL,
            )
        )

    base = raw_map["fs-rs"]
    rs_ssc = raw_map["fs-rs-ssc"]
    cross_ic = raw_map["fs-cross-ic"]
    sync_checks = []
    for name, config in [
        ("fs-rs", build_reverse_demo_config(include_forward_ssc=False, num_r=120, num_nu=121, num_gam_e=121)),
        ("fs-rs-ssc", build_reverse_demo_config(include_forward_ssc=False, num_r=120, num_nu=121, num_gam_e=121)),
        ("fs-cross-ic", build_reverse_demo_config(include_forward_ssc=False, num_r=120, num_nu=121, num_gam_e=121)),
    ]:
        if name == "fs-rs-ssc":
            config.reverse_shock.include_ssc = True
        if name == "fs-cross-ic":
            config.reverse_shock.include_cross_zone_ic = True
        setup = build_query_setup(config, REVERSE_TIMES_S, OBS_FREQS_HZ)
        dynamics = solve_dynamics(setup.boundary, config)
        electron = solve_electron(setup.boundary, dynamics, setup.seed_frequency_hz, config)
        reverse_emission = solve_reverse_shock_emission(setup.boundary, dynamics, setup.seed_frequency_hz, config)
        sync_checks.append((name, electron.l_syn_spec, None if reverse_emission is None else reverse_emission.l_syn_spec))
    results.extend(
        [
            ClosureResult(
                name="rs-toggle-fwd-sync-stability",
                category="forward-reverse",
                passed=max(
                    _max_rel(sync_checks[0][1], sync_checks[1][1]),
                    _max_rel(sync_checks[0][1], sync_checks[2][1]),
                ) <= SYNC_STABILITY_TOL,
                expected=0.0,
                measured=max(
                    _max_rel(sync_checks[0][1], sync_checks[1][1]),
                    _max_rel(sync_checks[0][1], sync_checks[2][1]),
                ),
                tolerance=SYNC_STABILITY_TOL,
            ),
            ClosureResult(
                name="rs-toggle-rev-sync-stability",
                category="forward-reverse",
                passed=max(
                    _max_rel(sync_checks[0][2], sync_checks[1][2]),
                    _max_rel(sync_checks[0][2], sync_checks[2][2]),
                ) <= SYNC_STABILITY_TOL,
                expected=0.0,
                measured=max(
                    _max_rel(sync_checks[0][2], sync_checks[1][2]),
                    _max_rel(sync_checks[0][2], sync_checks[2][2]),
                ),
                tolerance=SYNC_STABILITY_TOL,
            ),
        ]
    )

    reverse_config = build_reverse_demo_config(include_forward_ssc=False, num_r=140, num_nu=121, num_gam_e=121)
    reverse_setup = build_simulation_setup(reverse_config)
    reverse_dynamics = solve_dynamics(reverse_setup.boundary, reverse_config)
    reverse_emission = solve_reverse_shock_emission(reverse_setup.boundary, reverse_dynamics, reverse_setup.seed_frequency_hz, reverse_config)
    if reverse_emission is None:
        results.append(
            ClosureResult(
                name="reverse-peak-consistency",
                category="forward-reverse",
                passed=False,
                extra={"message": "missing reverse-shock emission"},
            )
        )
        return results

    shell_strength = np.max(reverse_emission.l_syn_spec, axis=0)
    shell_order = np.argsort(shell_strength)[::-1]
    shell_indices: list[int] = []
    for idx in shell_order:
        peak_idx = int(np.argmax(reverse_emission.l_syn_spec[:, idx]))
        if peak_idx == 0 or peak_idx == reverse_emission.l_syn_spec.shape[0] - 1:
            continue
        shell_indices.append(int(idx))
        if len(shell_indices) == 8:
            break
    worst_dex = 0.0
    for idx in shell_indices:
        peak_idx = int(np.argmax(reverse_emission.l_syn_spec[:, idx]))
        peak_nu = float(reverse_setup.seed_frequency_hz[peak_idx])
        refs = np.array(
            [
                reverse_emission.nu_a[idx],
                reverse_emission.nu_m[idx],
                reverse_emission.nu_c[idx],
                reverse_emission.nu_M[idx],
            ],
            dtype=float,
        )
        refs = refs[np.isfinite(refs) & (refs > 0.0)]
        if refs.size == 0:
            continue
        worst_dex = max(worst_dex, float(np.min(np.abs(np.log10(peak_nu / refs)))))
    results.append(
        ClosureResult(
            name="reverse-peak-consistency",
            category="forward-reverse",
            passed=worst_dex <= RS_PEAK_TOL_DEX,
            expected=0.0,
            measured=worst_dex,
            tolerance=RS_PEAK_TOL_DEX,
            extra={"shell_indices": shell_indices},
        )
    )
    return results


def _raw_component_sum(component_spectra) -> np.ndarray:
    total = np.array(component_spectra.fwd_sync, dtype=float, copy=True)
    total += np.array(component_spectra.fwd_ssc, dtype=float, copy=False)
    if component_spectra.rev_sync is not None:
        total += np.array(component_spectra.rev_sync, dtype=float, copy=False)
    if component_spectra.rev_ssc is not None:
        total += np.array(component_spectra.rev_ssc, dtype=float, copy=False)
    if component_spectra.cross_ic is not None:
        total += np.array(component_spectra.cross_ic, dtype=float, copy=False)
    return total


def check_observer_frame_closure() -> list[ClosureResult]:
    model = _base_model()
    pair = []
    pair_from_grid = []
    spectrum = []
    sky = []
    for t_obs, nu_obs in zip(PAIR_TIMES_S, PAIR_FREQS_HZ):
        single_t = np.array([float(t_obs)], dtype=float)
        single_nu = np.array([float(nu_obs)], dtype=float)
        pair.append(float(model.flux_density(single_t, single_nu).total[0]))
        pair_from_grid.append(float(model.flux_density_grid(single_t, single_nu).total[0, 0]))
        spectrum.append(float(model.spectrum(float(t_obs), single_nu)[0]))
        img = model.sky_image(float(t_obs), nu_obs=float(nu_obs), fov=500.0 * units.uas, npixel=64)
        sky.append(float(np.sum(img.image[0]) * img.pixel_solid_angle))
    pair = np.array(pair, dtype=float)
    pair_from_grid = np.array(pair_from_grid, dtype=float)
    spectrum = np.array(spectrum, dtype=float)
    sky = np.array(sky, dtype=float)

    return [
        ClosureResult(
            name="observer-pair-vs-grid",
            category="observer",
            passed=_max_rel(pair, pair_from_grid) <= OBS_CLOSURE_TOL,
            expected=0.0,
            measured=_max_rel(pair, pair_from_grid),
            tolerance=OBS_CLOSURE_TOL,
        ),
        ClosureResult(
            name="observer-pair-vs-spectrum",
            category="observer",
            passed=_max_rel(pair, spectrum) <= OBS_CLOSURE_TOL,
            expected=0.0,
            measured=_max_rel(pair, spectrum),
            tolerance=OBS_CLOSURE_TOL,
        ),
        ClosureResult(
            name="observer-pair-vs-sky-image",
            category="observer",
            passed=_max_rel(pair, sky) <= SKY_IMAGE_TOL,
            expected=0.0,
            measured=_max_rel(pair, sky),
            tolerance=SKY_IMAGE_TOL,
        ),
    ]


def main() -> None:
    checks = (
        check_regime_ordering_closure()
        + check_spectral_slope_closure()
        + check_dynamics_radiation_closure()
        + check_absorption_emission_closure()
        + check_forward_reverse_closure()
        + check_observer_frame_closure()
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
    print("PASS: physical closure check succeeded.")


if __name__ == "__main__":
    main()
