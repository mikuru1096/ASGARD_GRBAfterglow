from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

import am3

from asgard_core.hadronic_hummer import (
    CHARGED_PION_DECAY_S,
    MUON_DECAY_S,
    MUON_MASS_GEV,
    PI_PLUS_MASS_GEV,
    solve_hummer2010_pgamma,
)


def _configure_am3_photopion_only() -> am3.AM3:
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
    sim.set_process_bethe_heitler(0)
    sim.set_process_pp(0)
    sim.set_process_muon_decay(0)
    sim.set_process_pion_decay(0)
    sim.set_process_photopion(1)
    sim.set_process_photopion_photon_loss(0)
    sim.init_kernels()
    return sim


def _configure_am3_decay_only() -> am3.AM3:
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
    sim.set_process_bethe_heitler(0)
    sim.set_process_pp(0)
    sim.set_process_photopion(0)
    sim.set_process_muon_decay(0)
    sim.set_process_pion_decay(1)
    sim.init_kernels()
    return sim


def _check_photopion_proton_loss_timescale() -> None:
    sim = _configure_am3_photopion_only()
    hadron_energy_gev = np.asarray(sim.get_egrid_had(), dtype=float) * 1.0e-9
    photon_energy_gev = np.asarray(sim.get_egrid_photons(), dtype=float) * 1.0e-9
    neutrino_energy_gev = np.asarray(sim.get_egrid_neutrinos(), dtype=float) * 1.0e-9

    proton_density_per_gev = 1.0e-18 * (hadron_energy_gev / hadron_energy_gev[10]) ** -2.2
    photon_density_per_gev = 1.0e-10 * (photon_energy_gev / photon_energy_gev[50]) ** -1.5

    sim.clear_am3()
    sim.set_current_densities_protons(hadron_energy_gev * proton_density_per_gev)
    sim.set_current_densities_photons(photon_energy_gev * photon_density_per_gev)
    sim.evolve_step()

    am3_t_pg_s = np.asarray(sim.get_t_proton_photopion(), dtype=float)
    asgard = solve_hummer2010_pgamma(
        proton_energy_gev=hadron_energy_gev,
        proton_density_per_gev=proton_density_per_gev,
        photon_energy_gev=photon_energy_gev,
        photon_density_per_gev=photon_density_per_gev,
        gamma_energy_gev=photon_energy_gev,
        neutrino_energy_gev=neutrino_energy_gev,
        process_energy_gev=photon_energy_gev,
    )
    asgard_t_pg_s = 1.0 / np.maximum(np.asarray(asgard.proton_loss_rate, dtype=float), 1.0e-300)

    mask = np.isfinite(am3_t_pg_s) & (am3_t_pg_s > 0.0) & np.isfinite(asgard_t_pg_s) & (asgard_t_pg_s > 0.0)
    if not np.any(mask):
        raise AssertionError("No overlapping positive photopion loss timescales were produced by AM3 and ASGARD.")

    ratio = asgard_t_pg_s[mask] / am3_t_pg_s[mask]
    median_ratio = float(np.median(ratio))
    shape_rel = np.abs(asgard_t_pg_s[mask] / median_ratio - am3_t_pg_s[mask]) / am3_t_pg_s[mask]
    max_shape_rel = float(np.max(shape_rel))

    if not (0.95 <= median_ratio <= 1.05):
        raise AssertionError(f"ASGARD vs AM3 photopion proton-loss median timescale ratio out of range: {median_ratio:.6e}")
    if max_shape_rel > 1.0e-10:
        raise AssertionError(f"ASGARD vs AM3 photopion proton-loss timescale shape mismatch: {max_shape_rel:.6e}")


def _check_decay_timescales() -> None:
    sim = _configure_am3_decay_only()
    hadron_energy_gev = np.asarray(sim.get_egrid_had(), dtype=float) * 1.0e-9

    am3_t_pion_decay_s = np.asarray(sim.get_t_pion_decay(), dtype=float)
    am3_t_muon_decay_s = np.asarray(sim.get_t_muon_decay(), dtype=float)

    asgard_t_pion_decay_s = (hadron_energy_gev / PI_PLUS_MASS_GEV) * CHARGED_PION_DECAY_S
    asgard_t_muon_decay_s = (hadron_energy_gev / MUON_MASS_GEV) * MUON_DECAY_S

    pion_rel = np.abs(asgard_t_pion_decay_s - am3_t_pion_decay_s) / am3_t_pion_decay_s
    muon_rel = np.abs(asgard_t_muon_decay_s - am3_t_muon_decay_s) / am3_t_muon_decay_s

    max_pion_rel = float(np.max(pion_rel))
    max_muon_rel = float(np.max(muon_rel))
    if max_pion_rel > 3.0e-4:
        raise AssertionError(f"ASGARD vs AM3 pion decay timescale mismatch: {max_pion_rel:.6e}")
    if max_muon_rel > 3.0e-5:
        raise AssertionError(f"ASGARD vs AM3 muon decay timescale mismatch: {max_muon_rel:.6e}")


def main() -> None:
    _check_photopion_proton_loss_timescale()
    _check_decay_timescales()
    print("hadronic_am3_reference_compare_smoke: ok")


if __name__ == "__main__":
    main()
