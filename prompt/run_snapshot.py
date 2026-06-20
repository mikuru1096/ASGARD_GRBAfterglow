from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from prompt.eats import EATSNumerics
from prompt.internal_shock import InternalShockNumerics, InternalShockShell, simulate_internal_shock
from prompt.radiation import InternalShockMicrophysics, RadiationNumerics, compute_prompt_observed_flux


def main() -> None:
    slow = InternalShockShell(gamma=200.0, total_energy_iso_erg=4.0e52, duration_s=0.2, sigma=0.01)
    fast = InternalShockShell(gamma=600.0, total_energy_iso_erg=4.0e52, duration_s=0.2, sigma=0.03)
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
        numerics=InternalShockNumerics(num_branch_steps=64),
    )
    observer_time = np.linspace(1.0e-4, 2.0, 96)
    observer_frequency = np.logspace(16.0, 25.0, 32)
    flux = compute_prompt_observed_flux(
        solution,
        observer_frequency,
        observer_time,
        microphysics=microphysics,
        radiation_numerics=RadiationNumerics(num_threads=1),
        eats_numerics=EATSNumerics(num_theta=64, num_phi=1, num_threads=1),
    )
    output_dir = Path(__file__).resolve().parent / "results"
    output_dir.mkdir(exist_ok=True)
    np.savez(
        output_dir / "snapshot_flux.npz",
        observer_time_s=observer_time,
        observer_frequency_hz=observer_frequency,
        fs_sync=flux.fs_sync,
        fs_ssc=flux.fs_ssc,
        rs_sync=flux.rs_sync,
        rs_ssc=flux.rs_ssc,
        total=flux.total,
    )
    diagnostics = {
        "radius_collision_cm": solution.radius_collision_cm,
        "gamma_contact": solution.gamma_contact,
        "fs_valid_shock": bool(solution.fs.valid_shock),
        "rs_valid_shock": bool(solution.rs.valid_shock),
        "fs_crossing_time_lab_s": solution.fs.jump.crossing_time_lab_s,
        "rs_crossing_time_lab_s": solution.rs.jump.crossing_time_lab_s,
        "sigma_slow": solution.slow_shell.sigma,
        "sigma_fast": solution.fast_shell.sigma,
        "slow_baryonic_mass_g": solution.slow_baryonic_mass_g,
        "fast_baryonic_mass_g": solution.fast_baryonic_mass_g,
        "fs_compression": solution.fs.jump.compression,
        "rs_compression": solution.rs.jump.compression,
        "fs_specific_internal": solution.fs.jump.specific_internal_energy,
        "rs_specific_internal": solution.rs.jump.specific_internal_energy,
        "fs_ordered_b_g_mean": float(np.mean(solution.fs.ordered_b_g)),
        "rs_ordered_b_g_mean": float(np.mean(solution.rs.ordered_b_g)),
        "fs_turbulent_b_g_mean": float(np.mean(solution.fs.turbulent_b_g)),
        "rs_turbulent_b_g_mean": float(np.mean(solution.rs.turbulent_b_g)),
        "fs_total_b_g_mean": float(np.mean(solution.fs.total_b_g)),
        "rs_total_b_g_mean": float(np.mean(solution.rs.total_b_g)),
        "fs_gamma_m_median": float(np.median(flux.fs_radiation.gamma_m)),
        "rs_gamma_m_median": float(np.median(flux.rs_radiation.gamma_m)),
        "fs_gamma_c_median": float(np.median(flux.fs_radiation.gamma_c)),
        "rs_gamma_c_median": float(np.median(flux.rs_radiation.gamma_c)),
        "fs_gamma_max_median": float(np.median(flux.fs_radiation.gamma_max)),
        "rs_gamma_max_median": float(np.median(flux.rs_radiation.gamma_max)),
        "fs_efficiency_median": float(np.median(flux.fs_radiation.efficiency)),
        "rs_efficiency_median": float(np.median(flux.rs_radiation.efficiency)),
        "fs_compactness_median": float(np.median(flux.fs_radiation.compactness)),
        "rs_compactness_median": float(np.median(flux.rs_radiation.compactness)),
        "fs_gamma_gamma_min_absorption": float(np.min(flux.fs_radiation.gamma_gamma_absorption)),
        "rs_gamma_gamma_min_absorption": float(np.min(flux.rs_radiation.gamma_gamma_absorption)),
    }
    (output_dir / "snapshot_diagnostics.json").write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
