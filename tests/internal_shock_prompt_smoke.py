from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from prompt.eats import EATSGeometry, EATSNumerics, project_branch_flux
from prompt.internal_shock import InternalShockNumerics, InternalShockShell, fast_shock_allowed, simulate_internal_shock
from prompt.radiation import InternalShockMicrophysics, RadiationNumerics, compute_prompt_observed_flux


def _run_case(sigma_slow: float, sigma_fast: float, gamma_fast: float = 600.0):
    slow = InternalShockShell(gamma=200.0, total_energy_iso_erg=4.0e52, duration_s=0.2, sigma=sigma_slow)
    fast = InternalShockShell(gamma=gamma_fast, total_energy_iso_erg=4.0e52, duration_s=0.2, sigma=sigma_fast)
    microphysics = InternalShockMicrophysics(epsilon_e=0.1, epsilon_b=0.01, electron_index_p=2.3)
    solution = simulate_internal_shock(
        slow,
        fast,
        engine_gap_s=0.2,
        redshift=0.5,
        luminosity_distance_cm=1.0e28,
        opening_angle_rad=0.1,
        epsilon_e=microphysics.epsilon_e,
        epsilon_b=microphysics.epsilon_b,
        numerics=InternalShockNumerics(num_branch_steps=24),
    )
    flux = compute_prompt_observed_flux(
        solution,
        observer_frequency_hz=np.logspace(16.0, 24.0, 12),
        observer_time_s=np.linspace(1.0e-4, 2.0, 32),
        microphysics=microphysics,
        radiation_numerics=RadiationNumerics(num_electron_gamma=41, num_photon_frequency=43, num_threads=1),
        eats_numerics=EATSNumerics(num_theta=16, num_phi=1, num_threads=1),
    )
    return solution, flux


def _relative_change(a: float, b: float) -> float:
    return abs(a - b) / abs(a)


def _assert_piecewise_smooth_positive_curve(curve: np.ndarray) -> None:
    scale = float(np.max(curve))
    active = np.flatnonzero(curve > scale * 1.0e-12)
    assert active.size > 0
    breaks = np.where(np.diff(active) > 1)[0] + 1
    for segment in np.split(active, breaks):
        if segment.size < 4:
            continue
        diff = np.diff(curve[segment])
        signs = np.sign(diff[np.abs(diff) > scale * 1.0e-3])
        if signs.size > 1:
            assert np.count_nonzero(signs[1:] * signs[:-1] < 0.0) <= 1


def _assert_prompt_eats_rejects_off_axis_phi_collapse() -> None:
    try:
        project_branch_flux(
            characteristic_time_s=np.array([0.0, 1.0, 2.0], dtype=float),
            gamma=np.full(3, 100.0, dtype=float),
            radius_cm=np.array([1.0e15, 2.0e15, 3.0e15], dtype=float),
            source_flux=np.ones((4, 3), dtype=float),
            seed_frequency_hz=np.logspace(10.0, 13.0, 4),
            observer_frequency_hz=np.array([1.0e11], dtype=float),
            observer_time_s=np.array([0.5], dtype=float),
            geometry=EATSGeometry(redshift=0.5, opening_angle_rad=0.1, viewing_angle_rad=0.02),
            numerics=EATSNumerics(num_theta=4, num_phi=1, num_threads=1),
        )
    except ValueError as exc:
        assert "off-axis EATS projection requires num_phi >= 2" in str(exc)
    else:
        raise AssertionError("prompt EATS accepted off-axis projection with num_phi=1")


def main() -> None:
    _assert_prompt_eats_rejects_off_axis_phi_collapse()
    hydro, hydro_flux = _run_case(0.0, 0.0)
    weak, weak_flux = _run_case(1.0e-4, 1.0e-4)
    magnetized, magnetized_flux = _run_case(0.1, 0.3)
    weak_magnetized, weak_magnetized_flux = _run_case(100.0, 100.0, gamma_fast=201.0)

    assert hydro.fs.valid_shock and hydro.rs.valid_shock
    assert hydro.slow_baryonic_mass_g == hydro.slow_shell.total_energy_iso_erg / (hydro.slow_shell.gamma * 2.99792458e10**2)
    assert _relative_change(hydro.gamma_contact, weak.gamma_contact) < 5.0e-3
    assert _relative_change(hydro.fs.jump.crossing_time_lab_s, weak.fs.jump.crossing_time_lab_s) < 5.0e-2
    assert _relative_change(hydro.rs.jump.crossing_time_lab_s, weak.rs.jump.crossing_time_lab_s) < 5.0e-2
    assert _relative_change(hydro.fs.jump.compression, weak.fs.jump.compression) < 5.0e-2
    assert _relative_change(hydro.rs.jump.specific_internal_energy, weak.rs.jump.specific_internal_energy) < 5.0e-2
    assert magnetized.slow_baryonic_mass_g < hydro.slow_baryonic_mass_g
    assert magnetized.fast_baryonic_mass_g < hydro.fast_baryonic_mass_g
    assert np.mean(magnetized.fs.ordered_b_g) > np.mean(hydro.fs.ordered_b_g)
    assert np.mean(magnetized.rs.ordered_b_g) > np.mean(hydro.rs.ordered_b_g)
    assert magnetized.fs.jump.magnetic_fraction > hydro.fs.jump.magnetic_fraction
    assert magnetized.rs.jump.magnetic_fraction > hydro.rs.jump.magnetic_fraction
    assert not fast_shock_allowed(1.0, 0.1)
    assert not weak_magnetized.fs.valid_shock
    assert np.all(weak_magnetized_flux.fs_sync == 0.0)
    assert np.all(weak_magnetized_flux.fs_ssc == 0.0)

    for flux in (hydro_flux, weak_flux, magnetized_flux):
        for matrix in (flux.fs_sync, flux.fs_ssc, flux.rs_sync, flux.rs_ssc, flux.total):
            assert np.all(np.isfinite(matrix))
            assert np.all(matrix >= 0.0)
        assert np.any(flux.total > 0.0)
        band = np.sum(flux.total, axis=0)
        _assert_piecewise_smooth_positive_curve(band)


if __name__ == "__main__":
    main()
