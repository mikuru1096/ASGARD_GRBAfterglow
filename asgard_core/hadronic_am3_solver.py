from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from asgard_core.hadronic_hummer import PROTON_MASS_GEV, solve_hummer2010_pgamma
from asgard_core.hadronic_pgamma import (
    kelner_aharonian_2008_secondary_spectrum,
    photon_density_hz_to_gev,
)
from src import constants


HUMMER_PROCESS_GROUP_LABELS: tuple[str, ...] = ("photopion", "pion_decay", "muon_decay")
GEV_TO_ERG = constants.para_gev2erg
HUMMER2010_RESPONSE_BACKEND = "fortran_wrapped_response"
KA2008_REFERENCE_BACKEND = "python_reference"


@dataclass(frozen=True)
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


def solve_hummer_2010_response_processes(
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
    return _solve_hummer_backend(
        radius_cm,
        gam_p,
        d_n_gam_p,
        v_seed_hz,
        seed_target_hz,
        num_nu_nu,
        process_energy_gev=process_energy_gev,
        include_pg=include_pg,
        include_neutrino=include_neutrino,
    )


def solve_ka2008_reference_processes(
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
    """Canonical KA2008 stable-secondary reference backend."""
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
    photon_energy_gev, _ = photon_density_hz_to_gev(v_seed_arr, np.ones_like(v_seed_arr))
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
    am3_process_power = np.zeros((len(HUMMER_PROCESS_GROUP_LABELS), num_gam_p, num_r), dtype=float)
    process_luminosity = np.zeros((len(HUMMER_PROCESS_GROUP_LABELS), process_energy_arr.size, num_r), dtype=float)
    proton_energy_gev = gam_p_arr * PROTON_MASS_GEV
    shell_volume_arr = _shell_volumes_from_radius(radius_arr)

    for i_r in range(num_r):
        if not include_pg and not include_neutrino:
            continue
        shell_volume_cm3 = float(shell_volume_arr[i_r])
        if shell_volume_cm3 <= 0.0:
            raise ValueError("radius_cm must be positive.")

        proton_density_per_gev = d_n_gam_p_arr[:, i_r] / (shell_volume_cm3 * PROTON_MASS_GEV)
        _, photon_density_per_gev = photon_density_hz_to_gev(v_seed_arr, seed_target_arr[:, i_r])

        gamma_spec = np.zeros_like(photon_energy_gev)
        pion_lepton_spec = np.zeros_like(photon_energy_gev)
        muon_lepton_spec = np.zeros_like(photon_energy_gev)
        neutrino_spec = np.zeros_like(neutrino_energy_gev)
        gamma_process_spec = np.zeros_like(process_energy_arr)
        pion_process_spec = np.zeros_like(process_energy_arr)
        muon_process_spec = np.zeros_like(process_energy_arr)

        if include_pg:
            gamma_spec = kelner_aharonian_2008_secondary_spectrum(
                "gamma",
                photon_energy_gev,
                proton_energy_gev,
                proton_density_per_gev,
                photon_energy_gev,
                photon_density_per_gev,
            )
            pion_lepton_spec = (
                kelner_aharonian_2008_secondary_spectrum(
                    "e_plus",
                    photon_energy_gev,
                    proton_energy_gev,
                    proton_density_per_gev,
                    photon_energy_gev,
                    photon_density_per_gev,
                )
                + kelner_aharonian_2008_secondary_spectrum(
                    "nu_mu_bar",
                    photon_energy_gev,
                    proton_energy_gev,
                    proton_density_per_gev,
                    photon_energy_gev,
                    photon_density_per_gev,
                )
                + kelner_aharonian_2008_secondary_spectrum(
                    "nu_mu",
                    photon_energy_gev,
                    proton_energy_gev,
                    proton_density_per_gev,
                    photon_energy_gev,
                    photon_density_per_gev,
                )
            )
            gamma_process_spec = kelner_aharonian_2008_secondary_spectrum(
                "gamma",
                process_energy_arr,
                proton_energy_gev,
                proton_density_per_gev,
                photon_energy_gev,
                photon_density_per_gev,
            )
            pion_process_spec = (
                kelner_aharonian_2008_secondary_spectrum(
                    "e_plus",
                    process_energy_arr,
                    proton_energy_gev,
                    proton_density_per_gev,
                    photon_energy_gev,
                    photon_density_per_gev,
                )
                + kelner_aharonian_2008_secondary_spectrum(
                    "nu_mu_bar",
                    process_energy_arr,
                    proton_energy_gev,
                    proton_density_per_gev,
                    photon_energy_gev,
                    photon_density_per_gev,
                )
                + kelner_aharonian_2008_secondary_spectrum(
                    "nu_mu",
                    process_energy_arr,
                    proton_energy_gev,
                    proton_density_per_gev,
                    photon_energy_gev,
                    photon_density_per_gev,
                )
            )
            muon_process_spec = (
                kelner_aharonian_2008_secondary_spectrum(
                    "nu_e",
                    process_energy_arr,
                    proton_energy_gev,
                    proton_density_per_gev,
                    photon_energy_gev,
                    photon_density_per_gev,
                )
                + kelner_aharonian_2008_secondary_spectrum(
                    "e_minus",
                    process_energy_arr,
                    proton_energy_gev,
                    proton_density_per_gev,
                    photon_energy_gev,
                    photon_density_per_gev,
                )
                + kelner_aharonian_2008_secondary_spectrum(
                    "nu_e_bar",
                    process_energy_arr,
                    proton_energy_gev,
                    proton_density_per_gev,
                    photon_energy_gev,
                    photon_density_per_gev,
                )
            )
            muon_lepton_spec = (
                kelner_aharonian_2008_secondary_spectrum(
                    "nu_e",
                    photon_energy_gev,
                    proton_energy_gev,
                    proton_density_per_gev,
                    photon_energy_gev,
                    photon_density_per_gev,
                )
                + kelner_aharonian_2008_secondary_spectrum(
                    "e_minus",
                    photon_energy_gev,
                    proton_energy_gev,
                    proton_density_per_gev,
                    photon_energy_gev,
                    photon_density_per_gev,
                )
                + kelner_aharonian_2008_secondary_spectrum(
                    "nu_e_bar",
                    photon_energy_gev,
                    proton_energy_gev,
                    proton_density_per_gev,
                    photon_energy_gev,
                    photon_density_per_gev,
                )
            )

        if include_neutrino:
            nu_mu = kelner_aharonian_2008_secondary_spectrum(
                "nu_mu",
                neutrino_energy_gev,
                proton_energy_gev,
                proton_density_per_gev,
                photon_energy_gev,
                photon_density_per_gev,
            )
            nu_mu_bar = kelner_aharonian_2008_secondary_spectrum(
                "nu_mu_bar",
                neutrino_energy_gev,
                proton_energy_gev,
                proton_density_per_gev,
                photon_energy_gev,
                photon_density_per_gev,
            )
            nu_e = kelner_aharonian_2008_secondary_spectrum(
                "nu_e",
                neutrino_energy_gev,
                proton_energy_gev,
                proton_density_per_gev,
                photon_energy_gev,
                photon_density_per_gev,
            )
            nu_e_bar = kelner_aharonian_2008_secondary_spectrum(
                "nu_e_bar",
                neutrino_energy_gev,
                proton_energy_gev,
                proton_density_per_gev,
                photon_energy_gev,
                photon_density_per_gev,
            )
            neutrino_spec = nu_mu + nu_mu_bar + nu_e + nu_e_bar

        l_had_pg_gamma[:, i_r] = _energy_luminosity_from_spectrum(
            photon_energy_gev,
            gamma_spec,
            shell_volume_cm3,
        )
        neutrino_luminosity[:, i_r] = _energy_luminosity_from_spectrum(
            neutrino_energy_gev,
            neutrino_spec,
            shell_volume_cm3,
        )
        process_luminosity[0, :, i_r] = _energy_luminosity_from_spectrum(process_energy_arr, gamma_process_spec, shell_volume_cm3)
        process_luminosity[1, :, i_r] = _energy_luminosity_from_spectrum(process_energy_arr, pion_process_spec, shell_volume_cm3)
        process_luminosity[2, :, i_r] = _energy_luminosity_from_spectrum(process_energy_arr, muon_process_spec, shell_volume_cm3)

        proton_energy_weight = d_n_gam_p_arr[:, i_r] * proton_energy_gev
        total_weight = float(np.trapezoid(proton_energy_weight, proton_energy_gev))
        if total_weight > 0.0 and np.isfinite(total_weight):
            normalized_weight = proton_energy_weight / total_weight
        else:
            normalized_weight = np.zeros_like(proton_energy_weight)

        if include_pg:
            photopion_total = float(np.trapezoid(l_had_pg_gamma[:, i_r], v_seed_arr))
            pion_decay_total = float(
                np.trapezoid(
                    _energy_luminosity_from_spectrum(photon_energy_gev, pion_lepton_spec, shell_volume_cm3),
                    photon_energy_gev,
                )
            )
            muon_decay_total = float(
                np.trapezoid(
                    _energy_luminosity_from_spectrum(photon_energy_gev, muon_lepton_spec, shell_volume_cm3),
                    photon_energy_gev,
                )
            )
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
    )




def _solve_hummer_backend(
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
    photon_energy_gev, _ = photon_density_hz_to_gev(v_seed_arr, np.ones_like(v_seed_arr))
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
    am3_process_power = np.zeros((len(HUMMER_PROCESS_GROUP_LABELS), num_gam_p, num_r), dtype=float)
    process_luminosity = np.zeros((len(HUMMER_PROCESS_GROUP_LABELS), process_energy_arr.size, num_r), dtype=float)
    proton_reinj_rate = np.zeros((num_gam_p, num_r), dtype=float)
    neutron_reinj_rate = np.zeros((num_gam_p, num_r), dtype=float)
    proton_loss_rate = np.zeros((num_gam_p, num_r), dtype=float)
    neutron_loss_rate = np.zeros((num_gam_p, num_r), dtype=float)
    photon_loss_rate = np.zeros((v_seed_arr.size, num_r), dtype=float)

    proton_energy_gev = gam_p_arr * PROTON_MASS_GEV
    shell_volume_arr = _shell_volumes_from_radius(radius_arr)

    for i_r in range(num_r):
        if not include_pg and not include_neutrino:
            continue
        shell_volume_cm3 = float(shell_volume_arr[i_r])
        if shell_volume_cm3 <= 0.0:
            raise ValueError("radius_cm must be positive.")

        proton_density_per_gev = d_n_gam_p_arr[:, i_r] / (shell_volume_cm3 * PROTON_MASS_GEV)
        _, photon_density_per_gev = photon_density_hz_to_gev(v_seed_arr, seed_target_arr[:, i_r])
        backend = solve_hummer2010_pgamma(
            proton_energy_gev=proton_energy_gev,
            proton_density_per_gev=proton_density_per_gev,
            photon_energy_gev=photon_energy_gev,
            photon_density_per_gev=photon_density_per_gev,
            gamma_energy_gev=photon_energy_gev,
            neutrino_energy_gev=neutrino_energy_gev,
            process_energy_gev=process_energy_arr,
            neutron_density_per_gev=None,
        )
        proton_reinj_rate[:, i_r] = backend.proton_reinjection_rate_per_gev
        neutron_reinj_rate[:, i_r] = backend.neutron_reinjection_rate_per_gev
        proton_loss_rate[:, i_r] = backend.proton_loss_rate
        neutron_loss_rate[:, i_r] = backend.neutron_loss_rate
        photon_loss_rate[:, i_r] = backend.photon_loss_rate

        if include_pg:
            l_had_pg_gamma[:, i_r] = _energy_luminosity_from_rate_spectrum(
                backend.gamma_energy_gev,
                backend.gamma_rate_per_gev,
                shell_volume_cm3,
            )
            process_luminosity[0, :, i_r] = _energy_luminosity_from_rate_spectrum(
                backend.process_energy_gev,
                backend.process_rate_per_gev[0],
                shell_volume_cm3,
            )
            process_luminosity[1, :, i_r] = _energy_luminosity_from_rate_spectrum(
                backend.process_energy_gev,
                backend.process_rate_per_gev[1],
                shell_volume_cm3,
            )
            process_luminosity[2, :, i_r] = _energy_luminosity_from_rate_spectrum(
                backend.process_energy_gev,
                backend.process_rate_per_gev[2],
                shell_volume_cm3,
            )

        if include_neutrino:
            neutrino_luminosity[:, i_r] = _energy_luminosity_from_rate_spectrum(
                backend.neutrino_energy_gev,
                backend.neutrino_rate_per_gev,
                shell_volume_cm3,
            )

        proton_energy_weight = d_n_gam_p_arr[:, i_r] * proton_energy_gev
        total_weight = float(np.trapezoid(proton_energy_weight, proton_energy_gev))
        normalized_weight = proton_energy_weight / total_weight if total_weight > 0.0 and np.isfinite(total_weight) else np.zeros_like(proton_energy_weight)

        photopion_total = float(np.trapezoid(process_luminosity[0, :, i_r], process_energy_arr))
        pion_decay_total = float(np.trapezoid(
            _energy_luminosity_from_rate_spectrum(
                backend.neutrino_energy_gev,
                backend.prompt_pion_neutrino_rate_per_gev,
                shell_volume_cm3,
            ),
            backend.neutrino_energy_gev,
        ))
        muon_decay_total = float(np.trapezoid(
            _energy_luminosity_from_rate_spectrum(
                backend.neutrino_energy_gev,
                backend.muon_neutrino_rate_per_gev,
                shell_volume_cm3,
            ),
            backend.neutrino_energy_gev,
        ))
        muon_decay_total += float(np.trapezoid(
            _energy_luminosity_from_rate_spectrum(
                backend.process_energy_gev,
                backend.muon_electron_rate_per_gev,
                shell_volume_cm3,
            ),
            backend.process_energy_gev,
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


def _energy_luminosity_from_spectrum(
    energy_gev: np.ndarray,
    spectrum: np.ndarray,
    shell_volume_cm3: float,
) -> np.ndarray:
    """Convert a KA2008 differential production spectrum to spectral luminosity."""
    energy = np.asarray(energy_gev, dtype=float)
    spec = np.asarray(spectrum, dtype=float)
    if energy.ndim != 1 or spec.ndim != 1 or energy.shape != spec.shape:
        raise ValueError("energy_gev and spectrum must be 1d arrays with matching shapes.")
    if np.any(energy <= 0.0):
        raise ValueError("energy_gev must be positive.")
    if shell_volume_cm3 <= 0.0:
        raise ValueError("shell_volume_cm3 must be positive.")

    # KA2008 returns a number-production kernel in dN / (dt dV dE)-like form after
    # convolution with the proton and photon densities.  Convert to spectral power by
    # multiplying by the shell volume and by E dE/dnu = E h, then to cgs.
    rate = constants.para_c * shell_volume_cm3 * spec
    return rate * energy * constants.para_h_gev * GEV_TO_ERG


def _energy_luminosity_from_rate_spectrum(
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
    rate = shell_volume_cm3 * spec
    return rate * energy * constants.para_h_gev * GEV_TO_ERG


def _shell_volumes_from_radius(radius_cm: np.ndarray) -> np.ndarray:
    radius = np.asarray(radius_cm, dtype=float)
    if radius.ndim != 1:
        raise ValueError("radius_cm must be a 1d array.")
    if np.any(radius <= 0.0) or np.any(np.diff(radius) <= 0.0):
        raise ValueError("radius_cm must be positive and strictly increasing.")
    prev_radius = np.empty_like(radius)
    prev_radius[0] = 0.0
    prev_radius[1:] = radius[:-1]
    return (4.0 / 3.0) * math.pi * (radius**3 - prev_radius**3)
