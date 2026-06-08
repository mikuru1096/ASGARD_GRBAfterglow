"""γγ pair/synch cascade helpers and shell-sequence time evolution."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from asgard_core.hadronic_pair_production import ELECTRON_MASS_GEV, solve_pair_production
from src import constants


@dataclass(frozen=True)
class TimeDependentPairCascadeOutput:
    photon_energy_gev: np.ndarray
    electron_energy_gev: np.ndarray
    frequency_hz: np.ndarray
    photon_loss_rate_local_s_inv: np.ndarray
    tau_pair_path: np.ndarray
    pair_density_per_gamma: np.ndarray
    pair_syn_luminosity_hz: np.ndarray
    pair_syn_seed_per_hz: np.ndarray
    cascade_photon_density_per_gev: np.ndarray
    absorbed_power_gev_per_cm3_s: np.ndarray
    injected_power_gev_per_cm3_s: np.ndarray


def shell_path_time_seconds(radius_cm: float, gamma_bulk: float) -> float:
    if radius_cm <= 0.0:
        raise ValueError("pair cascade shell radius must be positive.")
    if gamma_bulk < 1.0:
        raise ValueError("pair cascade bulk Lorentz factor must be >= 1.")
    return float(radius_cm) / (12.0 * float(gamma_bulk) * constants.para_c)


def compute_time_dependent_pair_cascade_sequence(
    photon_energy_gev: np.ndarray,
    primary_photon_density_per_gev: np.ndarray,
    electron_energy_gev: np.ndarray,
    frequency_hz: np.ndarray,
    radius_cm: np.ndarray,
    gamma_bulk: np.ndarray,
    observer_time_s: np.ndarray,
    b_field_g: np.ndarray,
    *,
    num_threads: int,
    index_syn_integr: int,
    substeps_per_shell: int,
) -> TimeDependentPairCascadeOutput:
    """Advance the gamma-gamma pair/synch cascade as a time-dependent shell sequence."""
    from src.Electron.electron_radiation import electron_radiation_kernel as electron_radiation_module

    e_ph = _as_positive_grid(photon_energy_gev, "photon_energy_gev")
    e_e = _as_positive_grid(electron_energy_gev, "electron_energy_gev")
    nu = _as_positive_grid(frequency_hz, "frequency_hz")
    primary = np.asarray(primary_photon_density_per_gev, dtype=float)
    radius = _as_positive_grid(radius_cm, "radius_cm")
    gamma = np.asarray(gamma_bulk, dtype=float)
    tobs = np.asarray(observer_time_s, dtype=float)
    b_field = np.asarray(b_field_g, dtype=float)
    if primary.shape != (e_ph.size, radius.size):
        raise ValueError("primary photon density must have shape (num_photon, num_shell).")
    if e_ph.shape != nu.shape:
        raise ValueError("photon_energy_gev and frequency_hz must have matching shapes.")
    if gamma.shape != radius.shape or tobs.shape != radius.shape or b_field.shape != radius.shape:
        raise ValueError("gamma_bulk, observer_time_s, and b_field_g must match radius_cm.")
    if np.any(primary < 0.0) or np.any(gamma < 1.0) or np.any(b_field < 0.0):
        raise ValueError("pair cascade sequence received non-physical inputs.")
    if int(substeps_per_shell) < 1:
        raise ValueError("substeps_per_shell must be positive.")

    gam_e = e_e / ELECTRON_MASS_GEV
    gam_edge = _build_gamma_edges(gam_e)
    volume = _shell_volumes_from_radius(radius)
    num_ph = e_ph.size
    num_e = e_e.size
    num_shell = radius.size
    pair_prev = np.zeros(num_e, dtype=float)

    photon_loss = np.zeros((num_ph, num_shell), dtype=float)
    tau_pair = np.zeros((num_ph, num_shell), dtype=float)
    pair_density = np.zeros((num_e, num_shell), dtype=float)
    pair_lum = np.zeros((num_ph, num_shell), dtype=float)
    pair_seed = np.zeros((num_ph, num_shell), dtype=float)
    cascade_photon_density = np.zeros((num_ph, num_shell), dtype=float)
    absorbed_power = np.zeros(num_shell, dtype=float)
    injected_power = np.zeros(num_shell, dtype=float)

    for i_shell in range(num_shell):
        dt_shell = _shell_dt(tobs, i_shell)
        t_escape = shell_path_time_seconds(float(radius[i_shell]), float(gamma[i_shell]))
        dt_sub = dt_shell / float(substeps_per_shell)
        photon_density = np.array(primary[:, i_shell], dtype=float, copy=True)
        pair_current = np.array(pair_prev, dtype=float, copy=True)

        for _ in range(int(substeps_per_shell)):
            ppair = solve_pair_production(e_ph, photon_density, e_e)
            q_pair = (
                float(volume[i_shell])
                * np.asarray(ppair.pair_injection_rate_per_gev_total, dtype=float)
                * ELECTRON_MASS_GEV
                * dt_sub
            )
            pair_current = _advance_energy_loggamma(
                gam_e,
                gam_edge,
                pair_current,
                q_pair,
                _electron_loss_rates(gam_e, float(b_field[i_shell]), _dynamical_time(radius[i_shell], gamma[i_shell])),
                dt_sub,
            )
            p_syn_i, seed_syn_i = _pair_synchrotron_state(
                electron_radiation_module,
                int(index_syn_integr),
                float(radius[i_shell]),
                float(b_field[i_shell]),
                int(num_threads),
                gam_e,
                pair_current,
                nu,
            )
            q_syn_density = (np.asarray(seed_syn_i, dtype=float) / constants.para_h_gev) / t_escape
            photon_density = _advance_photon_density(
                photon_density,
                np.asarray(primary[:, i_shell], dtype=float) / t_escape + q_syn_density,
                np.asarray(ppair.photon_loss_rate, dtype=float) + 1.0 / t_escape,
                dt_sub,
            )
            absorbed_power[i_shell] = float(ppair.absorbed_power_gev_per_cm3_s)
            injected_power[i_shell] = float(ppair.injected_power_gev_per_cm3_s)

        final_pair = solve_pair_production(e_ph, photon_density, e_e)
        p_syn_i, seed_syn_i = _pair_synchrotron_state(
            electron_radiation_module,
            int(index_syn_integr),
            float(radius[i_shell]),
            float(b_field[i_shell]),
            int(num_threads),
            gam_e,
            pair_current,
            nu,
        )
        photon_loss[:, i_shell] = np.asarray(final_pair.photon_loss_rate, dtype=float)
        tau_pair[:, i_shell] = photon_loss[:, i_shell] * t_escape
        pair_density[:, i_shell] = pair_current
        pair_lum[:, i_shell] = np.asarray(p_syn_i, dtype=float)
        pair_seed[:, i_shell] = np.asarray(seed_syn_i, dtype=float)
        cascade_photon_density[:, i_shell] = np.asarray(seed_syn_i, dtype=float) / constants.para_h_gev
        pair_prev = pair_current

    return TimeDependentPairCascadeOutput(
        photon_energy_gev=e_ph,
        electron_energy_gev=e_e,
        frequency_hz=nu,
        photon_loss_rate_local_s_inv=photon_loss,
        tau_pair_path=tau_pair,
        pair_density_per_gamma=pair_density,
        pair_syn_luminosity_hz=pair_lum,
        pair_syn_seed_per_hz=pair_seed,
        cascade_photon_density_per_gev=cascade_photon_density,
        absorbed_power_gev_per_cm3_s=absorbed_power,
        injected_power_gev_per_cm3_s=injected_power,
    )


def _as_positive_grid(values: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1D grid.")
    if arr.size < 2:
        raise ValueError(f"{name} must contain at least two points.")
    if np.any(arr <= 0.0) or np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be positive and strictly increasing.")
    return arr


def _shell_dt(observer_time_s: np.ndarray, i_shell: int) -> float:
    if i_shell <= 0:
        dt = float(observer_time_s[0])
    else:
        dt = float(observer_time_s[i_shell] - observer_time_s[i_shell - 1])
    if dt <= 0.0:
        raise ValueError("pair cascade shell times must be positive and strictly increasing.")
    return dt


def _shell_volumes_from_radius(radius_cm: np.ndarray) -> np.ndarray:
    radius = np.asarray(radius_cm, dtype=float)
    prev = np.concatenate(([0.0], radius[:-1]))
    return (4.0 / 3.0) * np.pi * (radius**3 - prev**3)


def _dynamical_time(radius_cm: float, gamma_bulk: float) -> float:
    return float(radius_cm) / (float(gamma_bulk) * constants.para_c)


def _build_gamma_edges(gamma: np.ndarray) -> np.ndarray:
    gam = np.asarray(gamma, dtype=float)
    edge = np.empty(gam.size + 1, dtype=float)
    edge[1:-1] = np.sqrt(gam[:-1] * gam[1:])
    edge[0] = gam[0] * np.sqrt(gam[0] / gam[1])
    edge[-1] = gam[-1] * np.sqrt(gam[-1] / gam[-2])
    return edge


def _electron_loss_rates(gamma: np.ndarray, b_field_g: float, t_dyn_s: float) -> np.ndarray:
    coeff_syn = constants.para_sigmat * float(b_field_g) ** 2 / (6.0 * np.pi * constants.para_m_e * constants.para_c)
    return coeff_syn * gamma * gamma + gamma / float(t_dyn_s)


def _advance_energy_loggamma(
    gamma: np.ndarray,
    gamma_edge: np.ndarray,
    density_prev: np.ndarray,
    source_content: np.ndarray,
    loss_rate: np.ndarray,
    dt_s: float,
) -> np.ndarray:
    dgamma = gamma_edge[1:] - gamma_edge[:-1]
    content = np.asarray(density_prev, dtype=float) * dgamma
    cooled_gamma = gamma - np.asarray(loss_rate, dtype=float) * float(dt_s)
    target = np.searchsorted(gamma_edge, cooled_gamma, side="right") - 1
    next_content = np.zeros_like(content)
    valid = (target >= 0) & (target < gamma.size)
    np.add.at(next_content, target[valid], content[valid])
    return next_content / dgamma + np.asarray(source_content, dtype=float)


def _advance_photon_density(
    density_prev: np.ndarray,
    source_rate: np.ndarray,
    loss_rate: np.ndarray,
    dt_s: float,
) -> np.ndarray:
    attenuation = np.exp(-np.asarray(loss_rate, dtype=float) * float(dt_s))
    active = np.asarray(loss_rate, dtype=float) > 0.0
    out = np.asarray(density_prev, dtype=float) * attenuation
    out[active] += np.asarray(source_rate, dtype=float)[active] * (1.0 - attenuation[active]) / np.asarray(loss_rate, dtype=float)[active]
    out[~active] += np.asarray(source_rate, dtype=float)[~active] * float(dt_s)
    return out


def _pair_synchrotron_state(
    electron_radiation_module,
    index_syn_integr: int,
    radius_cm: float,
    b_field_g: float,
    num_threads: int,
    gamma_e: np.ndarray,
    pair_density: np.ndarray,
    frequency_hz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if b_field_g == 0.0:
        return np.zeros_like(frequency_hz, dtype=float), np.zeros_like(frequency_hz, dtype=float)
    p_syn_i, seed_syn_i = electron_radiation_module.get_syn_selected(
        int(index_syn_integr),
        float(radius_cm),
        float(b_field_g),
        int(num_threads),
        gamma_e,
        pair_density,
        frequency_hz,
    )
    return np.asarray(p_syn_i, dtype=float), np.asarray(seed_syn_i, dtype=float)
