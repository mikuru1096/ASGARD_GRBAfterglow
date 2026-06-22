from __future__ import annotations

import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor
from typing import Callable

import numpy as np

from asgard_core.angular_sampling import angular_separation, is_axisymmetric_jet
from asgard_core.asgard_physics_utils import compute_doppler
from asgard_core.asgard_state import make_query_cfg, make_query_setup, solve_state_from_setup
from src import Interpolation, Structured, constants


HUMMER_SCHEMES = {"hummer_2010_response"}
ELECTRON_1D_TRANSPORT_IDS = {"fullhide_1d": 1, "dg_1d": 2}
STRUCTURED_FORTRAN_MIN_GAMMA0 = 2.0
STRUCTURED_CHI_2D_MIN_GAMMA0 = 1.4


def solve_structured_jet_fortran(model, times_s: np.ndarray, nu_hz: np.ndarray, build_patch_config: Callable):
    from asgard_core.api_model import FluxPair, FluxResult

    if _uses_structured_chi_2d(model):
        return solve_structured_jet_chi_2d(model, times_s, nu_hz, build_patch_config)
    _assert_supported_structured_fortran(model)
    times = np.asarray(times_s, dtype=float)
    frequencies = np.asarray(nu_hz, dtype=float)
    sampled = _sample_structured_grid(model)
    _theta_centers, _phi_centers, e_iso, gamma0, active, _axisymmetric = sampled
    i_theta, i_phi = np.argwhere(active > 0)[0]
    solve_times = _solve_time_grid(model, times)
    base_config = build_patch_config(
        model,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_max,
        e_iso=float(e_iso[i_theta, i_phi]),
        gamma0=float(gamma0[i_theta, i_phi]),
        theta_center=0.0,
    )
    query_config = make_query_cfg(base_config, solve_times)
    query_config.num_r = max(int(query_config.num_r), int(solve_times.size))
    setup = make_query_setup(query_config, solve_times, frequencies)

    outputs = Structured.structured_jet_flux_1d(
        *_structured_kernel_args(model, query_config, setup, sampled, times, frequencies)
    )
    fwd_sync, fwd_ssc, _fwd_hadronic, rev_sync, total = (np.asarray(value, dtype=float) for value in outputs[:5])
    if base_config.nu_callback is not None:
        base_config.nu_callback(
            "structured_jet_1d",
            np.asarray(outputs[10], dtype=float),
            np.asarray(outputs[11], dtype=float),
            np.asarray(outputs[12], dtype=float),
        )
    zero = np.zeros_like(total)
    flux = FluxResult(
        total=total,
        fwd=FluxPair(sync=fwd_sync, ssc=fwd_ssc),
        rev=FluxPair(sync=rev_sync, ssc=zero),
        cross_ic=None,
    )
    details = _details_from_kernel_outputs(model, sampled, outputs)
    return flux, details


def solve_structured_jet_chi_2d(model, times_s: np.ndarray, nu_hz: np.ndarray, build_patch_config: Callable):
    from asgard_core.api_model import FluxPair, FluxResult, _make_details

    times = np.asarray(times_s, dtype=float)
    frequencies = np.asarray(nu_hz, dtype=float)
    sampled = _sample_structured_grid(
        model,
        min_gamma0=STRUCTURED_CHI_2D_MIN_GAMMA0,
        reject_transrelativistic=False,
    )
    theta_centers, _phi_centers, e_iso, gamma0, active, axisymmetric = sampled
    _assert_supported_structured_chi_2d(model, axisymmetric)
    solve_times = _solve_time_grid(model, times)
    ring_states = _solve_structured_chi_ring_states(
        model,
        build_patch_config,
        theta_centers,
        e_iso[:, 0],
        gamma0[:, 0],
        active[:, 0],
        solve_times,
        frequencies,
    )
    first_state = next((state for state in ring_states if state is not None), None)
    if first_state is None:
        raise ValueError("No active jet elements were found for the requested structured jet.")
    fwd_sync = _project_structured_chi_sync_once(model, ring_states, first_state, times, frequencies)
    zero = np.zeros_like(fwd_sync)
    flux = FluxResult(
        total=fwd_sync,
        fwd=FluxPair(sync=fwd_sync, ssc=zero),
        rev=FluxPair(sync=zero, ssc=zero),
        cross_ic=None,
    )
    details = _make_details(first_state.components, _patch_metadata(*sampled, model), state=first_state)
    return flux, details


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
        int(base_config.eats_num_theta),
        int(base_config.eats_num_phi),
        int(base_config.index_dyn),
        int(base_config.index_y),
        int(base_config.index_syn_integr),
        int(include_reverse),
        int(bool(model.fwd_rad.ssc)),
        int(bool(model.setups.hadronic_enabled and model.fwd_rad.epsilon_p > 0.0)),
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
        ELECTRON_1D_TRANSPORT_IDS[str(model.setups.electron_solver).lower()],
    )


def _solve_structured_chi_ring_states(
    model,
    build_patch_config: Callable,
    theta_centers: np.ndarray,
    e_iso: np.ndarray,
    gamma0: np.ndarray,
    active: np.ndarray,
    solve_times: np.ndarray,
    frequencies: np.ndarray,
):
    outer_threads, inner_threads = _structured_threads(model)
    ring_states = [None] * int(theta_centers.size)
    payloads = []
    theta_width = float(model.jet.theta_max) / float(theta_centers.size)
    for i_theta, theta_center in enumerate(theta_centers):
        if int(active[i_theta]) == 0:
            continue
        config = build_patch_config(
            model,
            theta_v=0.0,
            opening_angle_jet=theta_width,
            e_iso=float(e_iso[i_theta]),
            gamma0=float(gamma0[i_theta]),
            theta_center=float(theta_center),
        )
        query_config = make_query_cfg(config, solve_times)
        query_config.num_r = max(int(query_config.num_r), int(solve_times.size))
        query_config.num_threads = int(inner_threads)
        setup = make_query_setup(query_config, solve_times, frequencies)
        if int(outer_threads) == 1:
            ring_states[i_theta] = solve_state_from_setup(
                query_config,
                setup,
                requested_frequencies_hz=frequencies,
            )
        else:
            if query_config.nu_callback is not None:
                raise NotImplementedError("parallel structured chi_eats_2d ring solves do not support nu_callback.")
            payloads.append((int(i_theta), query_config, setup, np.asarray(frequencies, dtype=float)))
    if payloads:
        if "fork" not in mp.get_all_start_methods():
            raise NotImplementedError("parallel structured chi_eats_2d ring solves require a POSIX fork multiprocessing context.")
        context = mp.get_context("fork")
        worker_count = min(int(outer_threads), len(payloads))
        with ProcessPoolExecutor(max_workers=worker_count, mp_context=context) as executor:
            for i_theta, state in executor.map(_solve_structured_chi_ring_payload, payloads):
                ring_states[i_theta] = state
    return ring_states


def _solve_structured_chi_ring_payload(payload):
    i_theta, query_config, setup, frequencies = payload
    return i_theta, solve_state_from_setup(
        query_config,
        setup,
        requested_frequencies_hz=frequencies,
    )


def _project_structured_chi_sync_once(model, ring_states, first_state, times: np.ndarray, frequencies: np.ndarray) -> np.ndarray:
    order = np.argsort(frequencies)
    sorted_frequencies = frequencies[order]
    boundary = np.array(first_state.setup.boundary, dtype=float, copy=True)
    boundary[8] = float(model.jet.theta_max)
    boundary[9] = float(model.observer.theta_obs)
    arrays = _structured_chi_projection_arrays(ring_states, first_state)
    if arrays["direct_electron"]:
        flux_sorted = Interpolation.sed_interpolation_chi_structured_axisym_electron(
            boundary,
            arrays["r_tobs"],
            arrays["radius"],
            arrays["dne_chi"],
            arrays["b_chi"],
            arrays["chi_radius"],
            arrays["chi_gamma"],
            arrays["chi_weight"],
            arrays["gam_e"],
            sorted_frequencies,
            times,
            int(model.setups.structured_num_phi),
            int(model.setups.num_threads),
        )
    else:
        flux_sorted = Interpolation.sed_interpolation_chi_structured_axisym(
            boundary,
            arrays["r_tobs"],
            arrays["radius"],
            arrays["source_chi"],
            arrays["tau_chi"],
            arrays["chi_radius"],
            arrays["chi_gamma"],
            arrays["chi_weight"],
            first_state.setup.seed_frequency_hz,
            sorted_frequencies,
            times,
            int(model.setups.structured_num_phi),
            int(model.setups.num_threads),
        )
    if np.array_equal(order, np.arange(order.shape[0])):
        return flux_sorted
    flux_matrix = np.empty_like(flux_sorted)
    flux_matrix[order] = flux_sorted
    return flux_matrix


def _structured_chi_projection_arrays(ring_states, first_state) -> dict[str, np.ndarray]:
    electron = first_state.electron
    direct_electron = electron.b_chi_g is not None and not np.any(np.asarray(electron.l_syn_spec_chi, dtype=float))
    if direct_electron:
        num_gam, num_chi, num_r = np.asarray(electron.d_n_gam_e_chi, dtype=float).shape
    else:
        num_nu, num_chi, num_r = np.asarray(electron.l_syn_spec_chi, dtype=float).shape
    num_theta = len(ring_states)
    r_tobs = np.zeros((num_r, num_theta), dtype=float, order="F")
    radius = np.zeros((num_r, num_theta), dtype=float, order="F")
    source_chi = None if direct_electron else np.zeros((num_nu, num_chi, num_r, num_theta), dtype=float, order="F")
    tau_chi = None if direct_electron else np.zeros_like(source_chi, order="F")
    dne_chi = np.zeros((num_gam, num_chi, num_r, num_theta), dtype=float, order="F") if direct_electron else None
    b_chi = np.zeros((num_chi, num_r, num_theta), dtype=float, order="F") if direct_electron else None
    chi_radius = np.zeros((num_chi, num_r, num_theta), dtype=float, order="F")
    chi_gamma = np.ones_like(chi_radius, order="F")
    chi_weight = np.zeros_like(chi_radius, order="F")
    for i_theta, state in enumerate(ring_states):
        if state is None:
            continue
        r_tobs[:, i_theta] = np.asarray(state.components.fwd.characteristic_time_s, dtype=float)
        radius[:, i_theta] = np.asarray(state.components.fwd.radius_cm, dtype=float)
        if direct_electron:
            dne_chi[:, :, :, i_theta] = np.asarray(state.electron.d_n_gam_e_chi, dtype=float)
            b_chi[:, :, i_theta] = np.asarray(state.electron.b_chi_g, dtype=float)
        else:
            source = np.asarray(state.electron.l_syn_spec_chi, dtype=float) * np.asarray(state.observer.prefactor, dtype=float)[:, None, :]
            source_chi[:, :, :, i_theta] = source
            tau_chi[:, :, :, i_theta] = np.asarray(state.electron.tau_syn_chi, dtype=float)
        chi_radius[:, :, i_theta] = np.asarray(state.electron.chi_radius_cm, dtype=float)
        chi_gamma[:, :, i_theta] = np.asarray(state.electron.chi_gamma_bulk, dtype=float)
        chi_weight[:, :, i_theta] = np.asarray(state.electron.chi_dvolume_weight, dtype=float)
    return {
        "direct_electron": direct_electron,
        "r_tobs": np.asfortranarray(r_tobs),
        "radius": np.asfortranarray(radius),
        "source_chi": None if source_chi is None else np.asfortranarray(source_chi),
        "tau_chi": None if tau_chi is None else np.asfortranarray(tau_chi),
        "dne_chi": None if dne_chi is None else np.asfortranarray(dne_chi),
        "b_chi": None if b_chi is None else np.asfortranarray(b_chi),
        "gam_e": np.asfortranarray(np.asarray(electron.gam_e, dtype=float)),
        "chi_radius": np.asfortranarray(chi_radius),
        "chi_gamma": np.asfortranarray(chi_gamma),
        "chi_weight": np.asfortranarray(chi_weight),
    }


def _details_from_kernel_outputs(model, sampled, outputs):
    from asgard_core.api_model import CharTrack, TrackBundle

    theta_centers, phi_centers, e_iso, gamma0, active, axisymmetric = sampled
    track_tobs, track_gamma, track_radius, track_mass = (np.asarray(value, dtype=float) for value in outputs[5:9])
    track_bfield = np.asarray(outputs[9], dtype=float)
    return TrackBundle(
        fwd=CharTrack(
            t_obs=track_tobs,
            radius=track_radius,
            Gamma=track_gamma,
            N_p=track_mass / constants.para_m_p,
            Doppler=np.asarray(compute_doppler(track_gamma, model.observer.z), dtype=float),
            B_comv=track_bfield,
        ),
        rev=None,
        patches=_patch_metadata(theta_centers, phi_centers, e_iso, gamma0, active, axisymmetric, model),
    )


def _assert_supported_structured_fortran(model) -> None:
    if str(model.setups.electron_solver).lower() not in ELECTRON_1D_TRANSPORT_IDS:
        raise NotImplementedError("structured_backend='fortran_1d' requires electron_solver='fullhide_1d' or 'dg_1d'.")
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
    if model.jet.spreading:
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


def _uses_structured_chi_2d(model) -> bool:
    return (
        str(model.setups.geometry_kernel).lower() == "chi_eats_2d"
        or str(model.setups.electron_solver).lower().endswith("_2d")
    )


def _assert_supported_structured_chi_2d(model, axisymmetric: bool) -> None:
    if not axisymmetric:
        raise NotImplementedError("structured chi_eats_2d batch projection currently requires an axisymmetric jet profile.")
    if str(model.setups.geometry_kernel).lower() != "chi_eats_2d":
        raise NotImplementedError("structured 2d electron transport requires geometry_projection='chi_eats_2d'.")
    if str(model.setups.electron_solver).lower() != "fullhide_2d":
        raise NotImplementedError("structured chi_eats_2d batch projection currently requires electron_solver='fullhide_2d'.")
    if float(model.observer.theta_obs) != 0.0 and int(model.setups.structured_num_phi) < 2:
        raise ValueError("off-axis structured chi_eats_2d projection requires structured_num_phi >= 2.")
    if bool(model.fwd_rad.ssc):
        raise NotImplementedError("structured chi_eats_2d batch projection currently covers FS synchrotron+SSA, not SSC emission.")
    if bool(model.setups.rvs_shock):
        raise NotImplementedError("structured chi_eats_2d batch projection currently does not include reverse-shock emission.")
    if bool(model.setups.hadronic_enabled and model.fwd_rad.epsilon_p > 0.0):
        raise NotImplementedError("structured chi_eats_2d batch projection currently does not include hadronic emission.")
    if bool(model.setups.include_cross_zone_ic):
        raise NotImplementedError("structured chi_eats_2d batch projection currently does not include cross-zone IC.")
    if bool(model.fwd_rad.pair_production) or int(model.setups.pair_cascade_iterations) > 1:
        raise NotImplementedError("structured chi_eats_2d batch projection currently does not include pair-cascade emission.")
    if model.jet.spreading:
        raise NotImplementedError("Jet spreading is not implemented in the structured chi_eats_2d backend.")


def _sample_structured_grid(model, min_gamma0: float = STRUCTURED_FORTRAN_MIN_GAMMA0, reject_transrelativistic: bool = True):
    axisymmetric = is_axisymmetric_jet(model.jet)
    theta_centers = _cell_centers(0.0, float(model.jet.theta_max), int(model.setups.structured_num_theta))
    phi_centers = np.array([0.0], dtype=float) if axisymmetric else _cell_centers(0.0, 2.0 * np.pi, int(model.setups.structured_num_phi))
    e_iso = np.zeros((theta_centers.size, phi_centers.size), dtype=float)
    gamma0 = np.ones_like(e_iso)
    for i_theta, theta in enumerate(theta_centers):
        for i_phi, phi in enumerate(phi_centers):
            e_iso[i_theta, i_phi] = float(model.jet.energy_iso(float(phi), float(theta)))
            gamma0[i_theta, i_phi] = float(model.jet.gamma0(float(phi), float(theta)))
    transrelativistic = (e_iso > 0.0) & (gamma0 > 1.0) & (gamma0 < float(min_gamma0))
    if reject_transrelativistic and np.any(transrelativistic):
        gamma_min = float(np.min(gamma0[transrelativistic]))
        raise ValueError(
            "structured_backend='fortran_1d' requires Gamma0 >= "
            f"{min_gamma0:g} for every positive-energy active patch; "
            f"found Gamma0={gamma_min:.6g}. Trim the mildly/Newtonian cocoon before using "
            "the relativistic structured-jet afterglow backend."
        )
    active = np.asarray((e_iso > 0.0) & (gamma0 >= float(min_gamma0)), dtype=np.int32)
    if not np.any(active):
        raise ValueError("No active jet elements were found for the requested structured jet.")
    return theta_centers, phi_centers, e_iso, gamma0, active, axisymmetric


def _structured_threads(model) -> tuple[int, int]:
    mode = str(model.setups.structured_parallel_mode).lower()
    total = int(model.setups.num_threads)
    outer = model.setups.structured_outer_threads
    inner = model.setups.structured_inner_threads
    cpu_count = os.cpu_count()
    if total < 1:
        raise ValueError("num_threads must be positive for structured jet execution.")
    if cpu_count is not None and total > int(cpu_count):
        raise ValueError("num_threads exceeds the available CPU thread count for structured jet execution.")
    if bool(model.setups.rvs_shock and model.rvs_rad is not None):
        resolved_inner = int(total if inner is None else inner)
        if outer is not None and int(outer) != 1:
            raise ValueError("structured reverse-shock execution requires structured_outer_threads=1.")
        _validate_structured_thread_budget(1, resolved_inner, total, cpu_count)
        return 1, resolved_inner
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
    phi_values = np.linspace(0.0, 2.0 * np.pi, int(model.setups.structured_num_phi), endpoint=False) if axisymmetric else phi_centers
    theta_obs = float(model.observer.theta_obs)
    phi_obs = float(model.observer.phi_obs)
    patches = []
    for i_theta, theta in enumerate(theta_centers):
        theta_value = float(theta)
        for i_phi, phi in enumerate(phi_values):
            phi_value = float(phi)
            source_phi_index = 0 if axisymmetric else i_phi
            if int(active[i_theta, source_phi_index]) == 0:
                continue
            patches.append(
                {
                    "phi": phi_value,
                    "theta": theta_value,
                    "theta_v": float(angular_separation(theta_value, phi_value, theta_obs, phi_obs)),
                    "E_iso": float(e_iso[i_theta, source_phi_index]),
                    "Gamma0": float(gamma0[i_theta, source_phi_index]),
                }
            )
    return patches


def _cell_centers(start: float, stop: float, count: int) -> np.ndarray:
    edges = np.linspace(start, stop, count + 1)
    return 0.5 * (edges[:-1] + edges[1:])
