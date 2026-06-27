from __future__ import annotations

import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor
from typing import Callable

import numpy as np

from asgard_core.angular_sampling import angular_separation, build_patch_grid, is_axisymmetric_jet
from asgard_core.asgard_physics_utils import compute_doppler
from asgard_core.asgard_state import make_query_cfg, make_query_setup, solve_state_from_setup
from src import Interpolation, Structured, constants


from asgard_core.asgard_runtime import ELECTRON_1D_TRANSPORT_IDS

HUMMER_SCHEMES = {"hummer_2010_response"}
STRUCTURED_FORTRAN_MIN_GAMMA0 = 2.0
STRUCTURED_CHI_2D_MIN_GAMMA0 = 1.4


def _unsort_flux(flux_sorted: np.ndarray, order: np.ndarray) -> np.ndarray:
    """Undo frequency sorting applied before a Fortran call."""
    if np.array_equal(order, np.arange(order.shape[0])):
        return flux_sorted
    flux_matrix = np.empty_like(flux_sorted)
    flux_matrix[order] = flux_sorted
    return flux_matrix


def solve_structured_jet_fortran(model, times_s: np.ndarray, nu_hz: np.ndarray, build_patch_config: Callable):
    from asgard_core.api_model import FluxPair, FluxResult

    if _uses_structured_chi_2d(model):
        return solve_structured_jet_chi_2d(model, times_s, nu_hz, build_patch_config)
    _assert_supported_structured_fortran(model)
    times = np.asarray(times_s, dtype=float)
    frequencies = np.asarray(nu_hz, dtype=float)
    sampled = _sample_structured_grid(model)
    _theta_centers, _theta_edges, _phi_centers, e_iso, gamma0, active, _axisymmetric = sampled
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
        observer_time_s=times,
    )
    theta_centers, theta_edges, _phi_centers, e_iso, gamma0, active, axisymmetric = sampled
    _assert_supported_structured_chi_2d(model, axisymmetric)
    solve_times = _solve_time_grid(model, times)
    adaptive_rtol = float(getattr(model.setups, "structured_adaptive_rtol", 0.0))

    if adaptive_rtol > 0.0 and axisymmetric:
        fwd_sync, first_state = _solve_project_structured_chi_ring_flux_adaptive(
            model, build_patch_config,
            theta_centers, e_iso[:, 0], gamma0[:, 0], active[:, 0],
            theta_edges, solve_times, times, frequencies,
            adaptive_rtol=adaptive_rtol,
            adaptive_max_depth=int(getattr(model.setups, "structured_adaptive_max_depth", 4)),
        )
    else:
        fwd_sync, first_state = _solve_project_structured_chi_ring_flux(
            model, build_patch_config,
            theta_centers, e_iso[:, 0], gamma0[:, 0], active[:, 0],
            theta_edges, solve_times, times, frequencies,
        )

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
    theta_centers, _theta_edges, phi_centers, e_iso, gamma0, active, axisymmetric = sampled
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
    i_theta, query_config, setup, frequencies = payload
    return i_theta, solve_state_from_setup(
        query_config,
        setup,
        requested_frequencies_hz=frequencies,
    )


def _solve_project_structured_chi_ring_flux(
    model,
    build_patch_config: Callable,
    theta_centers: np.ndarray,
    e_iso: np.ndarray,
    gamma0: np.ndarray,
    active: np.ndarray,
    theta_edges: np.ndarray,
    solve_times: np.ndarray,
    times: np.ndarray,
    frequencies: np.ndarray,
):
    outer_threads, inner_threads = _structured_threads(model)
    active_indices = [int(i) for i, flag in enumerate(active) if int(flag) != 0]
    if not active_indices:
        raise ValueError("No active jet elements were found for the requested structured jet.")
    order = np.argsort(frequencies)
    sorted_frequencies = frequencies[order]
    flux_sorted = np.zeros((sorted_frequencies.size, times.size), dtype=float)

    first_index = active_indices[0]
    first_config, first_setup = _structured_chi_ring_config_setup(
        model,
        build_patch_config,
        float(theta_edges[first_index + 1] - theta_edges[first_index]),
        theta_centers,
        e_iso,
        gamma0,
        solve_times,
        sorted_frequencies,
        first_index,
        int(inner_threads),
    )
    first_state = solve_state_from_setup(
        first_config,
        first_setup,
        requested_frequencies_hz=sorted_frequencies,
        assemble_observer=False,
    )
    flux_sorted += _project_structured_chi_ring_state(
        first_state,
        float(theta_edges[first_index]),
        float(theta_edges[first_index + 1]),
        float(model.observer.theta_obs),
        int(model.setups.structured_num_phi),
        times,
        sorted_frequencies,
    )

    payloads = []
    for i_theta in active_indices[1:]:
        query_config, setup = _structured_chi_ring_config_setup(
            model,
            build_patch_config,
            float(theta_edges[i_theta + 1] - theta_edges[i_theta]),
            theta_centers,
            e_iso,
            gamma0,
            solve_times,
            sorted_frequencies,
            i_theta,
            int(inner_threads),
        )
        payloads.append((
            int(i_theta), query_config, setup,
            sorted_frequencies, times,
            float(theta_edges[i_theta]), float(theta_edges[i_theta + 1]),
            float(model.observer.theta_obs), int(model.setups.structured_num_phi),
        ))
    if payloads:
        if int(outer_threads) == 1:
            for _i_theta, flux_ring in map(_solve_project_structured_chi_ring_payload, payloads):
                flux_sorted += flux_ring
        else:
            if "fork" not in mp.get_all_start_methods():
                raise NotImplementedError("parallel structured chi_eats_2d ring solves require a POSIX fork multiprocessing context.")
            context = mp.get_context("fork")
            worker_count = min(int(outer_threads), len(payloads))
            with ProcessPoolExecutor(max_workers=worker_count, mp_context=context) as executor:
                for _i_theta, flux_ring in executor.map(_solve_project_structured_chi_ring_payload, payloads):
                    flux_sorted += flux_ring
    return _unsort_flux(flux_sorted, order), first_state


def _solid_angle_mid(theta_lo: float, theta_hi: float) -> float:
    """Theta such that cos(theta) bisects [cos(hi), cos(lo)]."""
    return float(np.arccos(0.5 * (np.cos(theta_lo) + np.cos(theta_hi))))


def _eval_jet_at_theta(model, theta: float) -> tuple[float, float]:
    """Evaluate E_iso and Gamma0 at theta using the jet model."""
    return float(model.jet.energy_iso(0.0, theta)), float(model.jet.gamma0(0.0, theta))


def _solve_project_structured_chi_ring_flux_adaptive(
    model,
    build_patch_config: Callable,
    theta_centers: np.ndarray,
    e_iso: np.ndarray,
    gamma0: np.ndarray,
    active: np.ndarray,
    theta_edges: np.ndarray,
    solve_times: np.ndarray,
    times: np.ndarray,
    frequencies: np.ndarray,
    adaptive_rtol: float = 0.02,
    adaptive_max_depth: int = 4,
):
    """Adaptive theta integration with solid-angle spacing and gradient-based refinement.

    1. Build a base grid with uniform solid-angle spacing (n_base = n_full // 6).
    2. Solve all base rings in parallel → {F_i} at solid-angle midpoints.
    3. For each adjacent pair (i, i+1): if the relative flux gradient
       |F_i - F_{i+1}| / |F_i + F_{i+1}| exceeds *adaptive_rtol*, add a
       refinement ring at the solid-angle midpoint between them.
    4. Recurse on refinement rings up to *adaptive_max_depth* levels.

    Midpoints use ``arccos(0.5*(cos(lo)+cos(hi)))`` for equal solid angle.
    E_iso / Gamma0 are re-evaluated at each actual midpoint via the jet model.
    """
    outer_threads, inner_threads = _structured_threads(model)
    active_indices = [int(i) for i, flag in enumerate(active) if int(flag) != 0]
    if not active_indices:
        raise ValueError("No active jet elements found.")

    theta_max = float(theta_edges[active_indices[-1] + 1])
    order = np.argsort(frequencies)
    sorted_frequencies = frequencies[order]
    flux_sorted = np.zeros((sorted_frequencies.size, times.size), dtype=float)
    theta_obs = float(model.observer.theta_obs)
    num_phi = int(model.setups.structured_num_phi)

    # ---- Build base grid (uniform solid-angle spacing) ----
    n_base = max(4, len(active_indices) // 6)
    cos_edges = np.linspace(1.0, np.cos(theta_max), n_base + 1)
    base_t_edges = np.arccos(cos_edges)

    # Each cell: (theta_lo, theta_hi, epsilon) where epsilon tracks the ring's
    # angular width relative to its neighbour for gradient comparison.
    cells = [(float(base_t_edges[i]), float(base_t_edges[i + 1]), 0)
             for i in range(n_base)]
    all_fluxes = {}   # (tlo, thi) -> flux_array
    first_state = None

    next_cells = []
    for level in range(adaptive_max_depth + 1):
        if not cells and level > 0:
            break

        # ---- Solve all cells at this level ----
        payloads = []
        cell_keys = []
        for tlo, thi, depth in cells:
            if (tlo, thi) in all_fluxes:
                continue
            tmid = _solid_angle_mid(tlo, thi)
            e_val, g_val = _eval_jet_at_theta(model, tmid)
            if g_val < 1.4:
                continue
            config = build_patch_config(
                model, theta_v=0.0, opening_angle_jet=thi - tlo,
                e_iso=e_val, gamma0=g_val, theta_center=tmid,
            )
            qc = make_query_cfg(config, solve_times)
            qc.num_r = max(int(qc.num_r), int(solve_times.size))
            qc.num_threads = int(inner_threads)
            setup = make_query_setup(qc, solve_times, sorted_frequencies)
            payloads.append((tlo, thi, qc, setup))
            cell_keys.append((tlo, thi))

        if not payloads:
            break

        # First payload: serial (provides first_state)
        tlo0, thi0, qc0, setup0 = payloads[0]
        if first_state is None:
            first_state = solve_state_from_setup(
                qc0, setup0, requested_frequencies_hz=sorted_frequencies, assemble_observer=False,
            )
        else:
            state0 = solve_state_from_setup(
                qc0, setup0, requested_frequencies_hz=sorted_frequencies, assemble_observer=False,
            )
        all_fluxes[(tlo0, thi0)] = _project_structured_chi_ring_state(
            first_state if first_state is not None else state0,
            tlo0, thi0, theta_obs, num_phi, times, sorted_frequencies,
        )

        # Remaining: serial or parallel
        remaining = payloads[1:]
        if remaining:
            rkeys = cell_keys[1:]
            if int(outer_threads) == 1 or len(remaining) <= 1:
                for (tlo, thi, qc, setup), key in zip(remaining, rkeys):
                    state = solve_state_from_setup(qc, setup, requested_frequencies_hz=sorted_frequencies,
                                                   assemble_observer=False)
                    all_fluxes[key] = _project_structured_chi_ring_state(
                        state, tlo, thi, theta_obs, num_phi, times, sorted_frequencies,
                    )
            else:
                if "fork" not in mp.get_all_start_methods():
                    raise NotImplementedError("Requires POSIX fork.")
                ctx = mp.get_context("fork")
                wc = min(int(outer_threads), len(remaining))
                proj_payloads = [(tlo, thi, qc, setup, sorted_frequencies, times, theta_obs, num_phi)
                                 for tlo, thi, qc, setup in remaining]
                with ProcessPoolExecutor(max_workers=wc, mp_context=ctx) as ex:
                    futures = [ex.submit(_adaptive_ring_solve_project, *p) for p in proj_payloads]
                    for f in futures:
                        tlo, thi, flux = f.result()
                        all_fluxes[(tlo, thi)] = flux

        # ---- Gradient-based refinement ----
        if level >= adaptive_max_depth:
            break

        # Track which cells get subdivided (their flux is replaced by children)
        subdivided = set()
        sorted_cells = sorted(cells, key=lambda x: x[0])
        for idx in range(len(sorted_cells) - 1):
            tlo_a, thi_a, d_a = sorted_cells[idx]
            tlo_b, thi_b, d_b = sorted_cells[idx + 1]
            if abs(thi_a - tlo_b) > 1e-12:
                continue
            flux_a = all_fluxes.get((tlo_a, thi_a))
            flux_b = all_fluxes.get((tlo_b, thi_b))
            if flux_a is None or flux_b is None:
                continue
            na = float(np.sum(np.abs(flux_a)))
            nb = float(np.sum(np.abs(flux_b)))
            total = na + nb
            if total < 1e-100:
                continue
            grad = abs(na - nb) / total
            if grad > adaptive_rtol:
                # Split cell A at its own solid-angle midpoint
                tmid_a = _solid_angle_mid(tlo_a, thi_a)
                next_cells.append((tlo_a, tmid_a, d_a + 1))
                next_cells.append((tmid_a, thi_a, d_a + 1))
                subdivided.add((tlo_a, thi_a))
                # Split cell B similarly
                tmid_b = _solid_angle_mid(tlo_b, thi_b)
                next_cells.append((tlo_b, tmid_b, d_b + 1))
                next_cells.append((tmid_b, thi_b, d_b + 1))
                subdivided.add((tlo_b, thi_b))

        cells = next_cells
        # Remove subdivided cells' fluxes from all_fluxes (children will replace them)
        for key in subdivided:
            all_fluxes.pop(key, None)

    if first_state is None:
        raise ValueError("No active jet elements found.")

    # ---- Sum all leaf fluxes ----
    for _, flux in all_fluxes.items():
        flux_sorted += flux

    return _unsort_flux(flux_sorted, order), first_state


def _adaptive_ring_solve_project(tlo, thi, qc, setup, freq, times, theta_obs, num_phi):
    """Picklable helper for ProcessPoolExecutor."""
    state = solve_state_from_setup(qc, setup, requested_frequencies_hz=freq, assemble_observer=False)
    flux = _project_structured_chi_ring_state(state, tlo, thi, theta_obs, num_phi, times, freq)
    return tlo, thi, flux


def _structured_chi_ring_config_setup(
    model,
    build_patch_config: Callable,
    theta_width: float,
    theta_centers: np.ndarray,
    e_iso: np.ndarray,
    gamma0: np.ndarray,
    solve_times: np.ndarray,
    frequencies: np.ndarray,
    i_theta: int,
    inner_threads: int,
):
    config = build_patch_config(
        model,
        theta_v=0.0,
        opening_angle_jet=float(theta_width),
        e_iso=float(e_iso[i_theta]),
        gamma0=float(gamma0[i_theta]),
        theta_center=float(theta_centers[i_theta]),
    )
    query_config = make_query_cfg(config, solve_times)
    query_config.num_r = max(int(query_config.num_r), int(solve_times.size))
    query_config.num_threads = int(inner_threads)
    setup = make_query_setup(query_config, solve_times, frequencies)
    return query_config, setup


def _solve_project_structured_chi_ring_payload(payload):
    (
        i_theta,
        query_config,
        setup,
        frequencies,
        times,
        theta_lo,
        theta_hi,
        theta_obs,
        num_phi,
    ) = payload
    state = solve_state_from_setup(
        query_config,
        setup,
        requested_frequencies_hz=frequencies,
        assemble_observer=False,
    )
    flux_ring = _project_structured_chi_ring_state(
        state,
        float(theta_lo), float(theta_hi), float(theta_obs), int(num_phi),
        times, frequencies,
    )
    return int(i_theta), flux_ring


def _project_structured_chi_ring_state(
    state,
    theta_lo: float,
    theta_hi: float,
    theta_obs: float,
    num_phi: int,
    times: np.ndarray,
    frequencies: np.ndarray,
    projection_seed_frequency_hz: np.ndarray | None = None,
) -> np.ndarray:
    boundary = np.array(state.setup.boundary, dtype=float, copy=True)
    boundary[9] = float(theta_obs)
    src_freq = projection_seed_frequency_hz if projection_seed_frequency_hz is not None else state.setup.seed_frequency_hz
    seed_frequency, selected = _trim_structured_ring_projection_seed(
        src_freq, frequencies, state.electron.chi_gamma_bulk,
        float(theta_lo), float(theta_hi), float(theta_obs), int(num_phi), float(state.config.z),
    )
    e = state.electron
    z = float(state.config.z)
    DL = float(state.setup.luminosity_distance_cm)
    flux_prefactor = (1.0 + z) / (4.0 * np.pi * DL * DL)
    F_ring = np.asfortranarray(
        np.asarray(e.l_syn_spec_chi, dtype=float)[selected, :, :] * flux_prefactor
    )
    Tau_ring = np.asfortranarray(
        np.asarray(e.tau_syn_chi, dtype=float)[selected, :, :]
    )
    return Interpolation.sed_interpolation_chi_structured_axisym_ring_precomputed(
        boundary,
        state.components.fwd.characteristic_time_s,
        state.components.fwd.radius_cm,
        F_ring, Tau_ring,
        e.chi_radius_cm, e.chi_gamma_bulk,
        e.chi_dvolume_weight,
        seed_frequency, frequencies, times,
        float(theta_lo), float(theta_hi), int(num_phi),
    )


def _trim_structured_ring_projection_seed(
    seed_frequency_hz: np.ndarray,
    frequencies_hz: np.ndarray,
    chi_gamma_bulk: np.ndarray,
    theta_lo: float,
    theta_hi: float,
    theta_obs: float,
    num_phi: int,
    redshift: float,
) -> tuple[np.ndarray, np.ndarray]:
    seed = seed_frequency_hz
    frequency = frequencies_hz
    beta = np.sqrt(1.0 - chi_gamma_bulk**-2)
    theta_center = 0.5 * (float(theta_lo) + float(theta_hi))
    phi = (np.arange(int(num_phi), dtype=float) + 0.5) * (2.0 * np.pi / float(num_phi))
    mu = np.cos(theta_obs) * np.cos(theta_center) + np.sin(theta_obs) * np.sin(theta_center) * np.cos(phi)
    doppler_min = np.inf
    doppler_max = 0.0
    for mu_phi in mu:
        doppler = chi_gamma_bulk * (1.0 - beta * mu_phi)
        doppler_min = min(doppler_min, float(np.min(doppler)))
        doppler_max = max(doppler_max, float(np.max(doppler)))
    selected = np.zeros(seed.shape, dtype=bool)
    redshift_factor = 1.0 + float(redshift)
    for nu_obs in frequency:
        source_min = float(nu_obs) * doppler_min * redshift_factor
        source_max = float(nu_obs) * doppler_max * redshift_factor
        lo = max(0, int(np.searchsorted(seed, source_min, side="left")) - 1)
        hi = min(seed.size - 1, int(np.searchsorted(seed, source_max, side="left")))
        selected[lo : hi + 1] = True
    return np.asfortranarray(seed[selected]), selected


def _details_from_kernel_outputs(model, sampled, outputs):
    from asgard_core.api_model import CharTrack, TrackBundle

    theta_centers, theta_edges, phi_centers, e_iso, gamma0, active, axisymmetric = sampled
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
        patches=_patch_metadata(theta_centers, theta_edges, phi_centers, e_iso, gamma0, active, axisymmetric, model),
    )


def _assert_no_unsupported_structured_features(model, backend: str) -> None:
    if not model.jet.spreading:
        return
    raise NotImplementedError(f"Jet spreading is not implemented in the structured {backend} backend.")


def _assert_no_fancy_physics(model, backend: str) -> None:
    fancy = []
    if bool(model.setups.include_cross_zone_ic):
        fancy.append("cross-zone IC")
    if bool(model.fwd_rad.pair_production) or int(model.setups.pair_cascade_iterations) > 1:
        fancy.append("pair-cascade")
    if fancy:
        raise NotImplementedError(f"structured_backend='{backend}' does not support: {', '.join(fancy)}.")


def _assert_supported_structured_fortran(model) -> None:
    if str(model.setups.electron_solver).lower() not in ELECTRON_1D_TRANSPORT_IDS:
        raise NotImplementedError("structured_backend='fortran_1d' requires electron_solver='fullhide_1d' or 'dg_1d'.")
    if bool(model.setups.rvs_shock) and model.rvs_rad is None:
        raise NotImplementedError("structured_backend='fortran_1d' requires rvs_rad when reverse shock is enabled.")
    if model.rvs_rad is not None and (bool(model.rvs_rad.ssc) or bool(model.setups.rvs_ssc)):
        raise NotImplementedError("structured_backend='fortran_1d' migrates reverse synchrotron only, not RS SSC.")
    if bool(model.setups.rvs_shock) and float(model.fwd_rad.reverse_epsilon_p) > 0.0:
        raise NotImplementedError("structured_backend='fortran_1d' does not migrate reverse-shock hadronic branches.")
    if bool(model.fwd_rad.bethe_heitler or model.fwd_rad.pp or model.fwd_rad.hadronic_inverse_compton):
        raise NotImplementedError("structured_backend='fortran_1d' does not migrate BH, pp, or hadronic IC branches.")
    _assert_supported_hadronic_branch(model)
    _assert_no_fancy_physics(model, "fortran_1d")
    _assert_no_unsupported_structured_features(model, "Fortran")


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
        raise NotImplementedError("structured chi_eats_2d ring projection currently requires an axisymmetric jet profile.")
    if str(model.setups.geometry_kernel).lower() != "chi_eats_2d":
        raise NotImplementedError("structured 2d electron transport requires geometry_projection='chi_eats_2d'.")
    if str(model.setups.electron_solver).lower() != "fullhide_2d":
        raise NotImplementedError("structured chi_eats_2d ring projection currently requires electron_solver='fullhide_2d'.")
    if float(model.observer.theta_obs) != 0.0 and int(model.setups.structured_num_phi) < 2:
        raise ValueError("off-axis structured chi_eats_2d projection requires structured_num_phi >= 2.")
    if bool(model.fwd_rad.ssc):
        raise NotImplementedError("structured chi_eats_2d ring projection currently covers FS synchrotron+SSA, not SSC emission.")
    if bool(model.setups.rvs_shock):
        raise NotImplementedError("structured chi_eats_2d ring projection currently does not include reverse-shock emission.")
    if bool(model.setups.hadronic_enabled and model.fwd_rad.epsilon_p > 0.0):
        raise NotImplementedError("structured chi_eats_2d ring projection currently does not include hadronic emission.")
    _assert_no_fancy_physics(model, "chi_eats_2d")
    _assert_no_unsupported_structured_features(model, "chi_eats_2d")


def _sample_structured_grid(
    model,
    min_gamma0: float = STRUCTURED_FORTRAN_MIN_GAMMA0,
    reject_transrelativistic: bool = True,
    observer_time_s: np.ndarray | None = None,
):
    axisymmetric = is_axisymmetric_jet(model.jet)
    grid = build_patch_grid(model, observer_time_s)
    theta_edges = np.asarray(grid.theta_edges, dtype=float)
    theta_centers = np.asarray(grid.theta_centers, dtype=float)
    phi_centers = np.array([0.0], dtype=float) if axisymmetric else np.asarray(grid.phi_centers, dtype=float)
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
    return theta_centers, theta_edges, phi_centers, e_iso, gamma0, active, axisymmetric


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
            np.log10(float(10**model.setups.t_obs_min_log10)),
            np.log10(float(10**model.setups.t_obs_max_log10)),
            base_count,
        )
    solve_t_min = min(float(10**model.setups.t_obs_min_log10), float(np.min(requested)))
    solve_t_max = float(np.max(requested))
    solve_count = max(base_count, model._detail_time_count(solve_t_min, solve_t_max))
    log_t_min = np.log10(solve_t_min)
    log_t_max = np.log10(solve_t_max)
    log_step = (log_t_max - log_t_min) / float(solve_count - 2) if solve_count > 2 else 0.0
    return np.logspace(log_t_min, log_t_max + log_step, solve_count)


def _patch_metadata(theta_centers, theta_edges, phi_centers, e_iso, gamma0, active, axisymmetric: bool, model):
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
                    "theta_lo": float(theta_edges[i_theta]),
                    "theta_hi": float(theta_edges[i_theta + 1]),
                    "theta_v": float(angular_separation(theta_value, phi_value, theta_obs, phi_obs)),
                    "patch_sampling": str(model.setups.patch_sampling).lower(),
                    "E_iso": float(e_iso[i_theta, source_phi_index]),
                    "Gamma0": float(gamma0[i_theta, source_phi_index]),
                }
            )
    return patches


