from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from asgard_core.hadronic_processes import MP_GEV, solve_pgamma
from src import constants


HUMMER_GROUPS: tuple[str, ...] = ("photopion", "pion_decay", "muon_decay")
GEV2ERG = constants.para_gev2erg
HUMMER_BACKEND = "fortran_wrapped_response"


@dataclass(frozen=True, slots=True)
class HadronicProcessOutput:
    """Process-resolved p-gamma output for a named hadronic backend."""

    l_had_pg_gamma: np.ndarray
    neutrino_frequency_hz: np.ndarray
    neutrino_luminosity: np.ndarray
    am3_process_power: np.ndarray
    process_energy_gev: np.ndarray
    process_luminosity: np.ndarray
    hadron_energy_gev: np.ndarray | None = None
    proton_reinjection_rate_per_gev: np.ndarray | None = None
    neutron_reinjection_rate_per_gev: np.ndarray | None = None
    proton_loss_rate: np.ndarray | None = None
    neutron_loss_rate: np.ndarray | None = None
    photon_loss_frequency_hz: np.ndarray | None = None
    photon_loss_rate: np.ndarray | None = None


def solve_hummer(
    radius_cm: np.ndarray,
    gam_p: np.ndarray,
    d_n_gam_p: np.ndarray,
    v_seed_hz: np.ndarray,
    seed_target_hz: np.ndarray,
    num_nu_nu: int,
    process_energy_gev: np.ndarray | None = None,
    *,
    include_pg: bool,
    include_neutrino: bool,
) -> HadronicProcessOutput:
    """Canonical backend for pgamma_scheme='hummer_2010_response'."""
    radius_arr = np.asarray(radius_cm, dtype=float)
    gam_p_arr = np.asarray(gam_p, dtype=float)
    d_n_gam_p_arr = np.asarray(d_n_gam_p, dtype=float)
    v_seed_arr = np.asarray(v_seed_hz, dtype=float)
    seed_target_arr = np.asarray(seed_target_hz, dtype=float)

    if radius_arr.ndim != 1:
        raise ValueError("radius_cm must be a 1d array.")
    if gam_p_arr.ndim != 1:
        raise ValueError("gam_p must be a 1d array.")
    if d_n_gam_p_arr.ndim != 2:
        raise ValueError("d_n_gam_p must be a 2d array.")
    if v_seed_arr.ndim != 1:
        raise ValueError("v_seed_hz must be a 1d array.")
    if seed_target_arr.ndim != 2:
        raise ValueError("seed_target_hz must be a 2d array.")
    if d_n_gam_p_arr.shape[0] != gam_p_arr.size:
        raise ValueError("d_n_gam_p first axis must match gam_p.")
    if seed_target_arr.shape[0] != v_seed_arr.size:
        raise ValueError("seed_target_hz first axis must match v_seed_hz.")
    if seed_target_arr.shape[1] != radius_arr.size:
        raise ValueError("seed_target_hz second axis must match radius_cm.")
    if int(num_nu_nu) < 2:
        raise ValueError("num_nu_nu must be >= 2.")

    num_r = radius_arr.size
    num_gam_p = gam_p_arr.size
    photon_energy_gev = constants.para_h_gev * v_seed_arr
    process_energy_arr = photon_energy_gev if process_energy_gev is None else np.asarray(process_energy_gev, dtype=float)
    if process_energy_arr.ndim != 1:
        raise ValueError("process_energy_gev must be a 1d array.")
    if np.any(process_energy_arr <= 0.0):
        raise ValueError("process_energy_gev must be positive.")
    neutrino_frequency_hz = np.logspace(
        np.log10(1.0e-3 * constants.para_gev2hz),
        np.log10(1.0e8 * constants.para_gev2hz),
        int(num_nu_nu),
    )
    neutrino_energy_gev = constants.para_h_gev * neutrino_frequency_hz

    l_had_pg_gamma = np.zeros((v_seed_arr.size, num_r), dtype=float)
    neutrino_luminosity = np.zeros((num_nu_nu, num_r), dtype=float)
    am3_process_power = np.zeros((len(HUMMER_GROUPS), num_gam_p, num_r), dtype=float)
    process_luminosity = np.zeros((len(HUMMER_GROUPS), process_energy_arr.size, num_r), dtype=float)
    (
        proton_reinj_rate, neutron_reinj_rate,
        proton_loss_rate, neutron_loss_rate,
    ) = np.zeros((4, num_gam_p, num_r), dtype=float)
    photon_loss_rate = np.zeros((v_seed_arr.size, num_r), dtype=float)

    proton_energy_gev = gam_p_arr * MP_GEV
    shell_volume_arr = _shellvolumes(radius_arr)

    for i_r in range(num_r):
        if not include_pg and not include_neutrino:
            continue
        shell_volume_cm3 = float(shell_volume_arr[i_r])
        if shell_volume_cm3 <= 0.0:
            raise ValueError("radius_cm must be positive.")

        proton_density_per_gev = d_n_gam_p_arr[:, i_r] / (shell_volume_cm3 * MP_GEV)
        photon_density_per_gev = seed_target_arr[:, i_r] / constants.para_h_gev
        backend = solve_pgamma(
            proton_gev=proton_energy_gev,
            proton_density=proton_density_per_gev,
            photon_gev=photon_energy_gev,
            photon_density=photon_density_per_gev,
            gamma_gev=photon_energy_gev,
            neutrino_gev=neutrino_energy_gev,
            process_gev=process_energy_arr,
            neutron_density=None,
        )
        proton_reinj_rate[:, i_r] = backend.proton_reinj
        neutron_reinj_rate[:, i_r] = backend.neutron_reinj
        proton_loss_rate[:, i_r] = backend.proton_loss
        neutron_loss_rate[:, i_r] = backend.neutron_loss
        photon_loss_rate[:, i_r] = backend.photon_loss

        if include_pg:
            l_had_pg_gamma[:, i_r] = _ratelnu(
                backend.gamma_gev,
                backend.gamma_rate,
                shell_volume_cm3,
            )
            for i_proc, process_rate in enumerate(backend.process_rate):
                process_luminosity[i_proc, :, i_r] = _ratelnu(
                    backend.process_gev,
                    process_rate,
                    shell_volume_cm3,
                )

        if include_neutrino:
            neutrino_luminosity[:, i_r] = _ratelnu(
                backend.neutrino_gev,
                backend.neutrino_rate,
                shell_volume_cm3,
            )

        proton_energy_weight = d_n_gam_p_arr[:, i_r] * proton_energy_gev
        total_weight = float(np.trapezoid(proton_energy_weight, proton_energy_gev))
        if total_weight <= 0.0 or not np.isfinite(total_weight):
            raise ValueError("proton energy distribution must contain positive finite energy.")
        normalized_weight = proton_energy_weight / total_weight

        photopion_total = float(np.trapezoid(
            _ratepower(backend.process_gev, backend.process_rate[0], shell_volume_cm3),
            backend.process_gev,
        ))
        pion_decay_total = float(np.trapezoid(
            _ratepower(
                backend.neutrino_gev,
                backend.pion_nu,
                shell_volume_cm3,
            ),
            backend.neutrino_gev,
        ))
        muon_decay_total = float(np.trapezoid(
            _ratepower(
                backend.neutrino_gev,
                backend.muon_nu,
                shell_volume_cm3,
            ),
            backend.neutrino_gev,
        ))
        muon_decay_total += float(np.trapezoid(
            _ratepower(
                backend.process_gev,
                backend.muon_electron,
                shell_volume_cm3,
            ),
            backend.process_gev,
        ))
        am3_process_power[0, :, i_r] = normalized_weight * photopion_total
        am3_process_power[1, :, i_r] = normalized_weight * pion_decay_total
        am3_process_power[2, :, i_r] = normalized_weight * muon_decay_total

    return HadronicProcessOutput(
        l_had_pg_gamma=l_had_pg_gamma,
        neutrino_frequency_hz=neutrino_frequency_hz,
        neutrino_luminosity=neutrino_luminosity,
        am3_process_power=am3_process_power,
        process_energy_gev=process_energy_arr,
        process_luminosity=process_luminosity,
        hadron_energy_gev=proton_energy_gev,
        proton_reinjection_rate_per_gev=proton_reinj_rate,
        neutron_reinjection_rate_per_gev=neutron_reinj_rate,
        proton_loss_rate=proton_loss_rate,
        neutron_loss_rate=neutron_loss_rate,
        photon_loss_frequency_hz=v_seed_arr,
        photon_loss_rate=photon_loss_rate,
    )


def _ratepower(
    energy_gev: np.ndarray,
    spectrum: np.ndarray,
    shell_volume_cm3: float,
) -> np.ndarray:
    energy = np.asarray(energy_gev, dtype=float)
    spec = np.asarray(spectrum, dtype=float)
    if energy.ndim != 1 or spec.ndim != 1 or energy.shape != spec.shape:
        raise ValueError("energy_gev and spectrum must be 1d arrays with matching shapes.")
    if np.any(energy <= 0.0):
        raise ValueError("energy_gev must be positive.")
    if shell_volume_cm3 <= 0.0:
        raise ValueError("shell_volume_cm3 must be positive.")
    return shell_volume_cm3 * spec * energy * GEV2ERG


def _ratelnu(
    energy_gev: np.ndarray,
    spectrum: np.ndarray,
    shell_volume_cm3: float,
) -> np.ndarray:
    return _ratepower(energy_gev, spectrum, shell_volume_cm3) * constants.para_h_gev


def _shellvolumes(radius_cm: np.ndarray) -> np.ndarray:
    radius = np.asarray(radius_cm, dtype=float)
    if radius.ndim != 1:
        raise ValueError("radius_cm must be a 1d array.")
    if np.any(radius <= 0.0) or np.any(np.diff(radius) <= 0.0):
        raise ValueError("radius_cm must be positive and strictly increasing.")
    prev_radius = np.empty_like(radius)
    prev_radius[0] = 0.0
    prev_radius[1:] = radius[:-1]
    return (4.0 / 3.0) * math.pi * (radius**3 - prev_radius**3)
