from __future__ import annotations

import os
from typing import Callable

import numpy as np

from asgard_core.asgard_physics_utils import compute_doppler, compute_maximum_synchrotron_frequency
from asgard_core.asgard_state import make_query_setup
from src import Structured, constants


AXISYMMETRIC_JET_KINDS = {"gaussian", "powerlaw", "twocomponent", "steppowerlaw"}
HUMMER_SCHEMES = {"hummer_2010_response"}


def solve_structured_jet_fortran(model, times_s: np.ndarray, nu_hz: np.ndarray, build_patch_config: Callable):
    from ASGARD.api_model import FluxResult

    _assert_supported_structured_fortran(model)
    times = np.asarray(times_s, dtype=float)
    frequencies = np.asarray(nu_hz, dtype=float)
    base_config = _base_patch_config(model, build_patch_config)
    setup = make_query_setup(base_config, _solve_time_grid(model, times), frequencies)
    sampled = _sample_structured_grid(model)

    outputs = Structured.structured_jet_flux_1d(
        *_structured_kernel_args(model, base_config, setup, sampled, times, frequencies)
    )
    flux = _flux_from_kernel_outputs(FluxResult, outputs)
    details = _details_from_kernel_outputs(model, base_config, sampled, outputs)
    return flux, details


def _base_patch_config(model, build_patch_config: Callable):
    return build_patch_config(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_max,
        e_iso=float(_first_active_value(model, "energy_iso")),
        gamma0=float(_first_active_value(model, "gamma0")),
        theta_center=0.0,
    )


def _structured_kernel_args(model, base_config, setup, sampled, times: np.ndarray, frequencies: np.ndarray) -> tuple:
    theta_centers, phi_centers, e_iso, gamma0, active, axisymmetric = sampled
    outer_threads, inner_threads = _structured_threads(model)
    include_reverse = bool(model.setups.rvs_shock and model.rvs_rad is not None)
    reverse_rad = model.rvs_rad

    return (
        setup.boundary,
        np.asfortranarray(e_iso),
        np.asfortranarray(gamma0),
        np.asfortranarray(active),
        np.asarray(setup.seed_frequency_hz, dtype=float),
        frequencies,
        times,
        int(base_config.num_r),
        int(base_config.num_gam_e),
        int(base_config.num_theta),
        int(base_config.num_phi),
        int(base_config.index_dyn),
        int(base_config.index_y),
        int(base_config.index_syn_integr),
        int(include_reverse),
        int(bool(model.fwd_rad.ssc)),
        int(_include_forward_hadronic(model)),
        int(bool(model.fwd_rad.proton_synch)),
        int(bool(model.fwd_rad.pg)),
        int(bool(model.fwd_rad.neutrino)),
        int(model.setups.num_gam_p),
        int(model.setups.num_nu_nu),
        float(model.setups.reverse_delta_t_s if include_reverse else 0.0),
        float(reverse_rad.eps_e if include_reverse else 0.0),
        float(reverse_rad.eps_B if include_reverse else 0.0),
        float(reverse_rad.p if include_reverse else 0.0),
        float(reverse_rad.xi_N if include_reverse else 0.0),
        float(model.setups.reverse_sigma if include_reverse else 0.0),
        float(model.fwd_rad.p),
        float(model.fwd_rad.epsilon_p),
        float(model.fwd_rad.eta_acc),
        int(axisymmetric),
        int(outer_threads),
        int(inner_threads),
        int(bool(model.setups.electron_adaptive_substeps)),
        float(model.setups.electron_substep_rtol),
        int(model.setups.electron_substep_min),
        int(model.setups.electron_substep_max),
        int(bool(model.fwd_rad.thermal_electrons)),
    )


def _flux_from_kernel_outputs(FluxResult, outputs):
    from ASGARD.api_model import FluxPair

    fwd_sync, fwd_ssc, _fwd_hadronic, rev_sync, total = (np.asarray(value, dtype=float) for value in outputs[:5])
    zero = np.zeros_like(total)
    return FluxResult(
        total=total,
        fwd=FluxPair(sync=fwd_sync, ssc=fwd_ssc),
        rev=FluxPair(sync=rev_sync, ssc=zero),
        cross_ic=None,
    )


def _details_from_kernel_outputs(model, base_config, sampled, outputs):
    from ASGARD.api_model import CharTrack, TrackBundle

    theta_centers, phi_centers, e_iso, gamma0, active, axisymmetric = sampled
    track_tobs, track_gamma, track_radius, track_mass = (np.asarray(value, dtype=float) for value in outputs[5:9])
    track_bfield, track_nu_m, track_nu_c, track_nu_a = (np.asarray(value, dtype=float) for value in outputs[9:13])
    return TrackBundle(
        fwd=CharTrack(
            t_obs=track_tobs,
            radius=track_radius,
            Gamma=track_gamma,
            N_p=track_mass / constants.para_m_p,
            Doppler=np.asarray(compute_doppler(track_gamma, model.observer.z), dtype=float),
            B_comv=track_bfield,
            nu_m=track_nu_m,
            nu_c=track_nu_c,
            nu_a=track_nu_a,
            nu_M=np.asarray(compute_maximum_synchrotron_frequency(track_gamma, track_radius, base_config), dtype=float),
        ),
        rev=None,
        patches=_patch_metadata(theta_centers, phi_centers, e_iso, gamma0, active, axisymmetric, model),
    )


def _assert_supported_structured_fortran(model) -> None:
    if str(model.setups.electron_solver).lower() != "fullhide_1d":
        raise NotImplementedError("structured_backend='fortran_1d' requires electron_solver='fullhide_1d'.")
    if bool(model.setups.rvs_shock) and model.rvs_rad is None:
        raise NotImplementedError("structured_backend='fortran_1d' requires rvs_rad when reverse shock is enabled.")
    if model.rvs_rad is not None and (bool(model.rvs_rad.ssc) or bool(model.setups.rvs_ssc)):
        raise NotImplementedError("structured_backend='fortran_1d' migrates reverse synchrotron only, not RS SSC.")
    if bool(model.setups.rvs_shock) and float(model.fwd_rad.reverse_epsilon_p) > 0.0:
        raise NotImplementedError("structured_backend='fortran_1d' does not migrate reverse-shock hadronic branches.")
    if bool(model.setups.include_cross_zone_ic):
        raise NotImplementedError("structured_backend='fortran_1d' does not support cross-zone IC.")
    if bool(model.fwd_rad.bethe_heitler or model.fwd_rad.pp or model.fwd_rad.hadronic_inverse_compton):
        raise NotImplementedError("structured_backend='fortran_1d' does not migrate BH, pp, or hadronic IC branches.")
    if bool(model.fwd_rad.pair_production) or int(model.setups.pair_cascade_iterations) > 1:
        raise NotImplementedError("structured_backend='fortran_1d' does not migrate pair-cascade branches.")
    _assert_supported_hadronic_branch(model)
    if getattr(model.jet, "spreading", False):
        raise NotImplementedError("Jet spreading is not implemented in the structured Fortran backend.")


def _assert_supported_hadronic_branch(model) -> None:
    solver = str(model.setups.hadronic_solver).lower()
    if bool(model.fwd_rad.pg or model.fwd_rad.neutrino):
        if solver != "am3_1d":
            raise NotImplementedError("structured p-gamma/neutrino output requires hadronic_solver='am3_1d'.")
        scheme = str(model.setups.pgamma_scheme if model.setups.pgamma_scheme != "disabled" else model.fwd_rad.pgamma_scheme)
        if scheme.lower() not in HUMMER_SCHEMES:
            raise NotImplementedError("structured p-gamma/neutrino output supports only the Hummer2010 response kernel.")
    elif solver not in {"legacy_1d", "am3_1d"}:
        raise NotImplementedError("structured_backend='fortran_1d' supports hadronic_solver='legacy_1d' or 'am3_1d'.")


def _sample_structured_grid(model):
    axisymmetric = str(model.jet.kind).lower() in AXISYMMETRIC_JET_KINDS
    theta_centers = _cell_centers(0.0, float(model.jet.theta_max), int(model.setups.patch_theta))
    phi_centers = np.array([0.0], dtype=float) if axisymmetric else _cell_centers(0.0, 2.0 * np.pi, int(model.setups.patch_phi))
    e_iso, gamma0 = _sample_energy_gamma(model, theta_centers, phi_centers)
    active = np.asarray((e_iso > 0.0) & (gamma0 > 1.0), dtype=np.int32)
    if not np.any(active):
        raise ValueError("No active jet elements were found for the requested structured jet.")
    return theta_centers, phi_centers, e_iso, gamma0, active, axisymmetric


def _sample_energy_gamma(model, theta_centers: np.ndarray, phi_centers: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    e_iso = np.zeros((theta_centers.size, phi_centers.size), dtype=float)
    gamma0 = np.ones_like(e_iso)
    for i_theta, theta in enumerate(theta_centers):
        for i_phi, phi in enumerate(phi_centers):
            e_iso[i_theta, i_phi] = float(model.jet.energy_iso(float(phi), float(theta)))
            gamma0[i_theta, i_phi] = float(model.jet.gamma0(float(phi), float(theta)))
    return e_iso, gamma0


def _structured_threads(model) -> tuple[int, int]:
    mode = str(getattr(model.setups, "structured_parallel_mode", "outer")).lower()
    total = int(model.setups.num_threads)
    outer = model.setups.structured_outer_threads
    inner = model.setups.structured_inner_threads
    cpu_count = os.cpu_count()
    if total < 1:
        raise ValueError("num_threads must be positive for structured jet execution.")
    if cpu_count is not None and total > int(cpu_count):
        raise ValueError("num_threads exceeds the available CPU thread count for structured jet execution.")
    if mode == "outer":
        resolved_outer = int(total if outer is None else outer)
        _validate_structured_thread_budget(resolved_outer, 1, total, cpu_count)
        return resolved_outer, 1
    if mode == "inner":
        resolved_inner = int(total if inner is None else inner)
        _validate_structured_thread_budget(1, resolved_inner, total, cpu_count)
        return 1, resolved_inner
    if mode == "nested":
        if outer is None or inner is None:
            raise ValueError("structured_parallel_mode='nested' requires structured_outer_threads and structured_inner_threads.")
        resolved_outer = int(outer)
        resolved_inner = int(inner)
        _validate_structured_thread_budget(resolved_outer, resolved_inner, total, cpu_count)
        return resolved_outer, resolved_inner
    raise ValueError("structured_parallel_mode must be 'outer', 'inner', or 'nested'.")


def _validate_structured_thread_budget(outer: int, inner: int, total: int, cpu_count: int | None) -> None:
    if outer < 1 or inner < 1:
        raise ValueError("structured outer and inner thread counts must be positive.")
    requested = outer * inner
    if requested > total:
        raise ValueError("structured outer_threads * inner_threads exceeds num_threads.")
    if cpu_count is not None and requested > int(cpu_count):
        raise ValueError("structured outer_threads * inner_threads exceeds available CPU threads.")


def _solve_time_grid(model, requested_times: np.ndarray) -> np.ndarray:
    requested = np.asarray(requested_times, dtype=float)
    base_count = max(int(model.setups.num_tobs), int(np.unique(requested).size))
    if requested.size <= 1:
        return np.logspace(
            np.log10(float(model.setups.observer_time_min_s)),
            np.log10(float(model.setups.observer_time_max_s)),
            base_count,
        )
    solve_t_min = min(float(model.setups.observer_time_min_s), float(np.min(requested)))
    solve_t_max = float(np.max(requested))
    solve_count = max(base_count, model._detail_time_count(solve_t_min, solve_t_max))
    log_t_min = np.log10(solve_t_min)
    log_t_max = np.log10(solve_t_max)
    log_step = (log_t_max - log_t_min) / float(solve_count - 2) if solve_count > 2 else 0.0
    return np.logspace(log_t_min, log_t_max + log_step, solve_count)


def _patch_metadata(theta_centers, phi_centers, e_iso, gamma0, active, axisymmetric: bool, model):
    phi_values = np.linspace(0.0, 2.0 * np.pi, int(model.setups.patch_phi), endpoint=False) if axisymmetric else phi_centers
    patches = []
    for i_theta, theta in enumerate(theta_centers):
        for i_phi, phi in enumerate(phi_values):
            source_phi_index = 0 if axisymmetric else i_phi
            if int(active[i_theta, source_phi_index]) == 0:
                continue
            patches.append(_patch_entry(model, theta, phi, e_iso[i_theta, source_phi_index], gamma0[i_theta, source_phi_index]))
    return patches


def _patch_entry(model, theta: float, phi: float, e_iso: float, gamma0: float) -> dict[str, float]:
    return {
        "phi": float(phi),
        "theta": float(theta),
        "theta_v": _angular_separation(float(theta), float(phi), float(model.observer.theta_obs), float(model.observer.phi_obs)),
        "E_iso": float(e_iso),
        "Gamma0": float(gamma0),
    }


def _cell_centers(start: float, stop: float, count: int) -> np.ndarray:
    edges = np.linspace(start, stop, count + 1)
    return 0.5 * (edges[:-1] + edges[1:])


def _angular_separation(theta1: float, phi1: float, theta2: float, phi2: float) -> float:
    cos_alpha = np.cos(theta1) * np.cos(theta2) + np.sin(theta1) * np.sin(theta2) * np.cos(phi1 - phi2)
    return float(np.arccos(np.clip(cos_alpha, -1.0, 1.0)))


def _first_active_value(model, name: str) -> float:
    for theta in _cell_centers(0.0, float(model.jet.theta_max), int(model.setups.patch_theta)):
        value = getattr(model.jet, name)(0.0, float(theta))
        if name == "energy_iso" and value > 0.0:
            return float(value)
        if name == "gamma0" and value > 1.0:
            return float(value)
    return float(getattr(model.jet, name)(0.0, 0.0))


def _include_forward_hadronic(model) -> bool:
    return bool(model.setups.hadronic_enabled and model.fwd_rad.epsilon_p > 0.0)
