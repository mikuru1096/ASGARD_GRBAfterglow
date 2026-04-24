from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import am3

from asgard_core.hadronic_acceleration import PROTON_MASS_GEV, estimate_max_gamma, proton_acceleration_timescale_s


LIGHT_SPEED_CGS = 2.99792458e10


def _configure_am3_proton_synch_acceleration() -> am3.AM3:
    sim = am3.AM3()
    sim.set_verbosity_level(0)
    sim.set_process_parse_sed(0)
    sim.set_process_hadronic(1)
    sim.set_process_escape(0)
    sim.set_process_expansion(0)
    sim.set_process_adiabatic_cooling(0)
    sim.set_process_electron_syn(0)
    sim.set_process_proton_syn(1)
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
    sim.set_process_photopion(0)
    sim.set_process_photopion_photon_loss(0)
    sim.set_mag_field(30.0)
    sim.set_proton_accel_factor(5.0)
    sim.init_kernels()
    return sim


def main() -> None:
    sim = _configure_am3_proton_synch_acceleration()
    proton_energy_gev = np.asarray(sim.get_egrid_had(), dtype=float) * 1.0e-9
    proton_gamma = proton_energy_gev / PROTON_MASS_GEV

    am3_t_acc_s = np.asarray(sim.get_t_proton_acceleration(), dtype=float)
    asgard_t_acc_s = proton_acceleration_timescale_s(
        gamma=proton_gamma,
        b_field_g=float(sim.get_mag_field()),
        eta_acc=float(sim.get_proton_accel_factor()),
    )
    mask = np.isfinite(am3_t_acc_s) & (am3_t_acc_s > 0.0) & np.isfinite(asgard_t_acc_s) & (asgard_t_acc_s > 0.0)
    if not np.any(mask):
        raise AssertionError("No overlapping positive acceleration timescales were produced by AM3 and ASGARD.")

    acc_ratio = asgard_t_acc_s[mask] / am3_t_acc_s[mask]
    acc_median_ratio = float(np.median(acc_ratio))
    acc_rel = np.abs(asgard_t_acc_s[mask] / acc_median_ratio - am3_t_acc_s[mask]) / am3_t_acc_s[mask]
    if not (0.999 <= acc_median_ratio <= 1.001):
        raise AssertionError(f"ASGARD vs AM3 proton acceleration median timescale ratio out of range: {acc_median_ratio:.6e}")
    if float(np.max(acc_rel)) > 1.0e-12:
        raise AssertionError(f"ASGARD vs AM3 proton acceleration timescale mismatch: {float(np.max(acc_rel)):.6e}")

    sim.set_current_densities_protons(np.full_like(proton_energy_gev, 1.0e-100))
    sim.evolve_step()
    am3_t_syn_s = np.asarray(sim.get_t_proton_syn(), dtype=float)
    valid = np.isfinite(am3_t_syn_s) & (am3_t_syn_s > 0.0)
    if not np.any(valid):
        raise AssertionError("AM3 proton synchrotron cooling timescale did not populate after one evolve_step().")

    crossing = np.flatnonzero(am3_t_acc_s >= am3_t_syn_s)
    if crossing.size == 0:
        raise AssertionError("AM3 proton acceleration and synchrotron timescales never crossed.")
    am3_emax_gev = float(proton_energy_gev[int(crossing[0])])

    radius_cm = 3.0 * float(sim.get_expansion_timescale()) * LIGHT_SPEED_CGS
    asgard = estimate_max_gamma(
        species="proton",
        b_field_g=float(sim.get_mag_field()),
        radius_cm=radius_cm,
        gamma_bulk=1.0,
        eta_acc=float(sim.get_proton_accel_factor()),
    )
    asgard_emax_gev = float(asgard.gamma_max * PROTON_MASS_GEV)
    rel_emax = abs(asgard_emax_gev - am3_emax_gev) / am3_emax_gev
    if rel_emax > 0.1:
        raise AssertionError(f"ASGARD vs AM3 proton maximum energy mismatch: {rel_emax:.6e}")

    print("hadronic_am3_acceleration_compare_smoke: ok")


if __name__ == "__main__":
    main()
