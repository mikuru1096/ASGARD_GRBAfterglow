from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

import am3

from asgard_core.hadronic_bethe_heitler import solve_bethe_heitler


def _configure_am3_bethe_heitler_only() -> am3.AM3:
    sim = am3.AM3()
    sim.set_verbosity_level(0)
    sim.set_process_parse_sed(0)
    sim.set_process_hadronic(1)
    sim.set_process_escape(0)
    sim.set_process_expansion(0)
    sim.set_process_adiabatic_cooling(0)
    sim.set_process_electron_syn(0)
    sim.set_process_proton_syn(0)
    sim.set_process_electron_compton(0)
    sim.set_process_proton_compton(0)
    sim.set_process_muon_syn(0)
    sim.set_process_pion_syn(0)
    sim.set_process_muon_compton(0)
    sim.set_process_pion_compton(0)
    sim.set_process_annihilation(0)
    sim.set_process_bethe_heitler(1)
    sim.set_process_pp(0)
    sim.set_process_muon_decay(0)
    sim.set_process_pion_decay(0)
    sim.set_process_photopion(0)
    sim.set_process_photopion_photon_loss(0)
    sim.init_kernels()
    return sim


def main() -> None:
    sim = _configure_am3_bethe_heitler_only()
    proton_energy_gev = np.asarray(sim.get_egrid_had(), dtype=float) * 1.0e-9
    photon_energy_gev = np.asarray(sim.get_egrid_photons(), dtype=float) * 1.0e-9
    electron_energy_gev = np.asarray(sim.get_egrid_lep(), dtype=float) * 1.0e-9

    proton_density_per_gev = 1.0e-18 * (proton_energy_gev / proton_energy_gev[20]) ** -2.2
    photon_density_per_gev = 1.0e-10 * (photon_energy_gev / photon_energy_gev[50]) ** -1.5

    sim.clear_am3()
    sim.set_current_densities_protons(proton_energy_gev * proton_density_per_gev)
    sim.set_current_densities_photons(photon_energy_gev * photon_density_per_gev)
    sim.evolve_step()

    am3_t_bh_s = np.asarray(sim.get_t_proton_bethe_heitler(), dtype=float)
    asgard = solve_bethe_heitler(
        proton_energy_gev=proton_energy_gev,
        proton_density_per_gev=proton_density_per_gev,
        photon_energy_gev=photon_energy_gev,
        photon_density_per_gev=photon_density_per_gev,
        electron_energy_gev=electron_energy_gev,
    )
    asgard_t_bh_s = 1.0 / np.maximum(-np.asarray(asgard.proton_loss_rate, dtype=float), 1.0e-300)

    mask = (
        np.isfinite(am3_t_bh_s)
        & (am3_t_bh_s > 0.0)
        & np.isfinite(asgard_t_bh_s)
        & (asgard_t_bh_s > 0.0)
        & (am3_t_bh_s < 1.0e50)
    )
    if not np.any(mask):
        raise AssertionError("No overlapping positive Bethe-Heitler timescales were produced by AM3 and ASGARD.")

    ratio = asgard_t_bh_s[mask] / am3_t_bh_s[mask]
    median_ratio = float(np.median(ratio))
    rel = np.abs(asgard_t_bh_s[mask] / median_ratio - am3_t_bh_s[mask]) / am3_t_bh_s[mask]
    p99_rel = float(np.percentile(rel, 99.0))

    if not (0.99 <= median_ratio <= 1.01):
        raise AssertionError(f"ASGARD vs AM3 Bethe-Heitler median timescale ratio out of range: {median_ratio:.6e}")
    if p99_rel > 0.15:
        raise AssertionError(f"ASGARD vs AM3 Bethe-Heitler timescale shape mismatch: {p99_rel:.6e}")

    print("hadronic_am3_bethe_heitler_compare_smoke: ok")


if __name__ == "__main__":
    main()
