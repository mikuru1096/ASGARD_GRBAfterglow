from __future__ import annotations

import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor
from typing import Callable

import numpy as np

from asgard_core.angular_sampling import angular_separation as angsep, build_patch_grid as patchgrid, is_axisymmetric_jet as axisjet
from asgard_core.asgard_physics_utils import doppler_factor
from asgard_core.asgard_setup import build_setup
from asgard_core.asgard_state import query_cfg, solve_setup
from src import Interpolation, Structured, constants


from asgard_core.asgard_runtime import ELECTRONTRANSPORT_IDS

HUMMER_SCHEMES = {"hummer_2010_response"}
MINGAMMA1D = 2.0
MINGAMMA2D = 1.4
AMR_THETAFACTOR = 5.0e-2
AMR_PHIFACTOR = 6.0e-1
AMR_MINTHETA = 16
AMR_MINPHI = 12
AMR_STEP = 4


def _unsort(flux_sorted: np.ndarray, order: np.ndarray) -> np.ndarray:
    """Undo frequency sorting applied before a Fortran call."""
    if np.array_equal(order, np.arange(order.shape[0])):
        return flux_sorted
    flux_matrix = np.empty_like(flux_sorted)
    flux_matrix[order] = flux_sorted
    return flux_matrix


def solve_structured(model, times_s: np.ndarray, nu_hz: np.ndarray, patchbuild: Callable):
    from asgard_core.api_model import FluxPair, FluxResult

    if _usechi(model):
        return solve_chi(model, times_s, nu_hz, patchbuild)
    _checkfortran(model)
    times = np.asarray(times_s, dtype=float)
    frequencies = np.asarray(nu_hz, dtype=float)
    sampled = _samplegrid(model)
    _, _, _, e_iso, gamma0, active, _ = sampled
    i_theta, i_phi = np.argwhere(active > 0)[0]
    solve_times = _timegrid(model, times)
    base_config = patchbuild(
        model,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_max,
        e_iso=float(e_iso[i_theta, i_phi]),
        gamma0=float(gamma0[i_theta, i_phi]),
        theta_center=0.0,
    )
    query_config = query_cfg(base_config, solve_times)
    query_config.num_r = max(int(query_config.num_r), int(solve_times.size))
    setup = build_setup(query_config, frequencies, observer_time_s=solve_times)

    outputs = Structured.jet_flux_1d(
        *_kernelargs(model, query_config, setup, sampled, times, frequencies)
    )
    fwd_sync, fwd_ssc, _, rev_sync, total = (np.asarray(value, dtype=float) for value in outputs[:5])
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
    details = _details(model, sampled, outputs)
    return flux, details


def solve_chi(model, times_s: np.ndarray, nu_hz: np.ndarray, patchbuild: Callable):
    from asgard_core.api_model import FluxPair, FluxResult, _make_details as makedetails

    times = np.asarray(times_s, dtype=float)
    frequencies = np.asarray(nu_hz, dtype=float)
    adaptive_rtol = float(getattr(model.setups, "structured_adaptive_rtol", 0.0))
    theta_count, phi_count = _chicounts(model, adaptive_rtol)
    sampled = _samplegrid(
        model,
        min_gamma0=MINGAMMA2D,
        reject_transrelativistic=False,
        obstime=times,
        theta_count=theta_count,
        phi_count=phi_count,
    )
    theta_centers, theta_edges, phi_centers, e_iso, gamma0, active, axisymmetric = sampled
    _checkchi(model, axisymmetric, phi_count)
    solve_times = _timegrid(model, times)

    fwd_sync, first_state = _ringflux(
        model, patchbuild,
        theta_centers, e_iso[:, 0], gamma0[:, 0], active[:, 0],
        theta_edges, solve_times, times, frequencies,
        num_phi=phi_count,
    )

    zero = np.zeros_like(fwd_sync)
    flux = FluxResult(
        total=fwd_sync,
        fwd=FluxPair(sync=fwd_sync, ssc=zero),
        rev=FluxPair(sync=zero, ssc=zero),
        cross_ic=None,
    )
    details = makedetails(
        first_state.components,
        _patchmeta(theta_centers, theta_edges, phi_centers, e_iso, gamma0, active, axisymmetric, model, phi_count),
        state=first_state,
    )
    return flux, details


def _chicounts(model, adaptive_rtol: float) -> tuple[int, int]:
    # 中文：这里只用初始喷流参数选择全局角网格；SSA 遮挡是非局域量，不能在 Python 里做局部 flux-AMR 累加。
    # English: This only selects a global angular grid from initial jet parameters; SSA occultation is nonlocal, so Python-side local flux-AMR summation is not valid.
    thetaCount = int(model.setups.structured_num_theta)
    phiCount = int(model.setups.structured_num_phi)
    gridMode = str(getattr(model.setups, "adaptive_grid", "manual")).lower()
    if float(adaptive_rtol) <= 0.0 or gridMode in {"manual", "fixed", "none"}:
        return thetaCount, phiCount
    if gridMode not in {"initial", "auto"}:
        raise ValueError("adaptive_grid must be 'manual', 'initial', or 'auto'.")
    return _amrcounts(model, int(getattr(model.setups, "structured_adaptive_max_depth", 4)))


def _amrcounts(model, depth: int) -> tuple[int, int]:
    thetaMax = float(model.jet.theta_max)
    gammaMax = max(float(getattr(model.jet, "lf", 1.0)), MINGAMMA2D)
    depthScale = 2.0 ** int(depth)
    thetaNeed = thetaMax * gammaMax / (AMR_THETAFACTOR * depthScale)
    thetaNeed = max(thetaNeed, _thetaneed(model, thetaMax))
    thetaCount = _roundamr(thetaNeed, AMR_MINTHETA)
    thetaObs = abs(float(model.observer.theta_obs))
    if thetaObs == 0.0:
        return thetaCount, 1
    thetaRef = _phiref(thetaMax, thetaObs, gammaMax)
    phiNeed = 2.0 * np.pi * gammaMax * np.sqrt(thetaRef * thetaObs) / (AMR_PHIFACTOR * depthScale)
    phiCount = _roundamr(phiNeed, AMR_MINPHI)
    return thetaCount, phiCount


def _roundamr(value: float, minimum: int) -> int:
    return int(AMR_STEP * np.ceil(max(float(value), float(minimum)) / AMR_STEP))


def _thetaneed(model, thetaMax: float) -> float:
    scales = []
    thetaCore = float(getattr(model.jet, "theta_c", 0.0))
    if thetaCore > 0.0 and str(getattr(model.jet, "kind", "")).lower() != "tophat":
        scales.append(0.5 * thetaCore)
    thetaTable = np.asarray(getattr(model.jet, "theta_table", ()), dtype=float)
    if thetaTable.size > 1:
        scales.append(float(np.min(np.diff(thetaTable))))
    if not scales:
        return 0.0
    return float(thetaMax) / min(scales)


def _phiref(thetaMax: float, thetaObs: float, gammaMax: float) -> float:
    beamTheta = 1.0 / float(gammaMax)
    if float(thetaObs) <= float(thetaMax):
        return max(float(thetaObs), beamTheta)
    return max(float(thetaMax), beamTheta)


def _kernelargs(model, base_config, setup, sampled, times: np.ndarray, frequencies: np.ndarray) -> tuple:
    _, _, _, e_iso, gamma0, active, axisymmetric = sampled
    outer_threads, inner_threads = _threadplan(model)
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
        ELECTRONTRANSPORT_IDS[str(model.setups.electron_solver).lower()],
    )


def _ringflux(
    model,
    patchbuild: Callable,
    theta_centers: np.ndarray,
    e_iso: np.ndarray,
    gamma0: np.ndarray,
    active: np.ndarray,
    theta_edges: np.ndarray,
    solve_times: np.ndarray,
    times: np.ndarray,
    frequencies: np.ndarray,
    num_phi: int,
):
    outer_threads, inner_threads = _threadplan(model)
    active_indices = [int(i) for i, flag in enumerate(active) if int(flag) != 0]
    if not active_indices:
        raise ValueError("No active jet elements were found for the requested structured jet.")
    order = np.argsort(frequencies)
    sorted_frequencies = frequencies[order]
    state_cache = _reusecache(e_iso, gamma0, active)
    first_index = active_indices[0]
    first_config, first_setup = _ringsetup(
        model,
        patchbuild,
        float(theta_edges[first_index + 1] - theta_edges[first_index]),
        float(theta_centers[first_index]),
        float(e_iso[first_index]),
        float(gamma0[first_index]),
        solve_times,
        sorted_frequencies,
        int(inner_threads),
    )
    first_state = solve_setup(
        first_config,
        first_setup,
        requested_frequencies_hz=sorted_frequencies,
        assemble_observer=False,
    )
    ring_data = [_ringdata(
        first_state,
        float(theta_edges[first_index]),
        float(theta_edges[first_index + 1]),
    )]

    payloads = []
    for i_theta in active_indices[1:]:
        if state_cache is not None:
            ring_data.append(_ringdata(
                first_state,
                float(theta_edges[i_theta]),
                float(theta_edges[i_theta + 1]),
            ))
            continue
        query_config, setup = _ringsetup(
            model,
            patchbuild,
            float(theta_edges[i_theta + 1] - theta_edges[i_theta]),
            float(theta_centers[i_theta]),
            float(e_iso[i_theta]),
            float(gamma0[i_theta]),
            solve_times,
            sorted_frequencies,
            int(inner_threads),
        )
        payloads.append((
            int(i_theta), query_config, setup,
            sorted_frequencies,
            float(theta_edges[i_theta]), float(theta_edges[i_theta + 1]),
        ))
    if payloads:
        if int(outer_threads) == 1:
            for itheta, ring_sample in map(_ringpayload, payloads):
                ring_data.append(ring_sample)
        else:
            if "fork" not in mp.get_all_start_methods():
                raise NotImplementedError("parallel structured chi_eats_2d ring solves require a POSIX fork multiprocessing context.")
            context = mp.get_context("fork")
            worker_count = min(int(outer_threads), len(payloads))
            with ProcessPoolExecutor(max_workers=worker_count, mp_context=context) as executor:
                for itheta, ring_sample in executor.map(_ringpayload, payloads):
                    ring_data.append(ring_sample)
    flux_sorted = _ringbatch(
        first_state,
        ring_data,
        float(model.observer.theta_obs),
        int(num_phi),
        times,
        sorted_frequencies,
    )
    return _unsort(flux_sorted, order), first_state


def _reusecache(e_iso: np.ndarray, gamma0: np.ndarray, active: np.ndarray):
    active_mask = np.asarray(active, dtype=bool)
    if not np.any(active_mask):
        return None
    active_energy = np.asarray(e_iso, dtype=float)[active_mask]
    active_gamma = np.asarray(gamma0, dtype=float)[active_mask]
    if np.all(active_energy == active_energy[0]) and np.all(active_gamma == active_gamma[0]):
        return {}
    return None


def _ringsetup(
    model,
    patchbuild: Callable,
    theta_width: float,
    theta_center: float,
    e_iso: float,
    gamma0: float,
    solve_times: np.ndarray,
    frequencies: np.ndarray,
    inner_threads: int,
):
    config = patchbuild(
        model,
        theta_v=0.0,
        opening_angle_jet=float(theta_width),
        e_iso=float(e_iso),
        gamma0=float(gamma0),
        theta_center=float(theta_center),
    )
    query_config = query_cfg(config, solve_times)
    query_config.num_r = max(int(query_config.num_r), int(solve_times.size))
    query_config.num_threads = int(inner_threads)
    setup = build_setup(query_config, frequencies, observer_time_s=solve_times)
    return query_config, setup


def _ringpayload(payload):
    (
        i_theta,
        query_config,
        setup,
        frequencies,
        theta_lo,
        theta_hi,
    ) = payload
    state = solve_setup(
        query_config,
        setup,
        requested_frequencies_hz=frequencies,
        assemble_observer=False,
    )
    return int(i_theta), _ringdata(state, float(theta_lo), float(theta_hi))


def _ringdata(state, theta_lo: float, theta_hi: float):
    e = state.electron
    return (
        float(theta_lo),
        float(theta_hi),
        np.asarray(state.components.fwd.characteristic_time_s, dtype=float),
        np.asarray(state.components.fwd.radius_cm, dtype=float),
        np.asarray(e.l_syn_spec_chi, dtype=float),
        np.asarray(e.tau_syn_chi, dtype=float),
        np.asarray(e.chi_radius_cm, dtype=float),
        np.asarray(e.chi_gamma_bulk, dtype=float),
        np.asarray(e.chi_dvolume_weight, dtype=float),
    )


def _ringbatch(
    reference_state,
    ring_data,
    theta_obs: float,
    num_phi: int,
    times: np.ndarray,
    frequencies: np.ndarray,
) -> np.ndarray:
    boundary = np.array(reference_state.setup.boundary, dtype=float, copy=True)
    boundary[9] = float(theta_obs)
    theta_lo = np.asfortranarray(np.array([sample[0] for sample in ring_data], dtype=float))
    theta_hi = np.asfortranarray(np.array([sample[1] for sample in ring_data], dtype=float))
    rtobs = np.asfortranarray(np.stack([sample[2] for sample in ring_data], axis=1))
    radius = np.asfortranarray(np.stack([sample[3] for sample in ring_data], axis=1))
    flux = np.asfortranarray(np.stack([sample[4] for sample in ring_data], axis=3))
    tau = np.asfortranarray(np.stack([sample[5] for sample in ring_data], axis=3))
    chi_radius = np.asfortranarray(np.stack([sample[6] for sample in ring_data], axis=2))
    chi_gamma = np.asfortranarray(np.stack([sample[7] for sample in ring_data], axis=2))
    chi_weight = np.asfortranarray(np.stack([sample[8] for sample in ring_data], axis=2))
    return Interpolation.sed_chiring_batchlum_ray(
        boundary,
        rtobs,
        radius,
        flux,
        tau,
        chi_radius,
        chi_gamma,
        chi_weight,
        reference_state.setup.seed_frequency_hz,
        frequencies,
        times,
        theta_lo,
        theta_hi,
        int(num_phi),
    )


def _details(model, sampled, outputs):
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
            Doppler=np.asarray(doppler_factor(track_gamma, model.observer.z), dtype=float),
            B_comv=track_bfield,
        ),
        rev=None,
        patches=_patchmeta(theta_centers, theta_edges, phi_centers, e_iso, gamma0, active, axisymmetric, model),
    )


def _checkspread(model, backend: str) -> None:
    if not model.jet.spreading:
        return
    raise NotImplementedError(f"Jet spreading is not implemented in the structured {backend} backend.")


def _checkphysics(model, backend: str) -> None:
    fancy = []
    if bool(model.setups.include_cross_zone_ic):
        fancy.append("cross-zone IC")
    if bool(model.fwd_rad.pair_production) or int(model.setups.pair_cascade_iterations) > 1:
        fancy.append("pair-cascade")
    if fancy:
        raise NotImplementedError(f"structured_backend='{backend}' does not support: {', '.join(fancy)}.")


def _checkfortran(model) -> None:
    if str(model.setups.electron_solver).lower() not in ELECTRONTRANSPORT_IDS:
        raise NotImplementedError("structured_backend='fortran_1d' requires electron_solver='fullhide_1d' or 'dg_1d'.")
    if bool(model.setups.rvs_shock) and model.rvs_rad is None:
        raise NotImplementedError("structured_backend='fortran_1d' requires rvs_rad when reverse shock is enabled.")
    if model.rvs_rad is not None and (bool(model.rvs_rad.ssc) or bool(model.setups.rvs_ssc)):
        raise NotImplementedError("structured_backend='fortran_1d' migrates reverse synchrotron only, not RS SSC.")
    if bool(model.setups.rvs_shock) and float(model.fwd_rad.reverse_epsilon_p) > 0.0:
        raise NotImplementedError("structured_backend='fortran_1d' does not migrate reverse-shock hadronic branches.")
    if bool(model.fwd_rad.bethe_heitler or model.fwd_rad.pp or model.fwd_rad.hadronic_inverse_compton):
        raise NotImplementedError("structured_backend='fortran_1d' does not migrate BH, pp, or hadronic IC branches.")
    _checkhadronic(model)
    _checkphysics(model, "fortran_1d")
    _checkspread(model, "Fortran")


def _checkhadronic(model) -> None:
    solver = str(model.setups.hadronic_solver).lower()
    if bool(model.fwd_rad.pg or model.fwd_rad.neutrino):
        if solver != "am3_1d":
            raise NotImplementedError("structured p-gamma/neutrino output requires hadronic_solver='am3_1d'.")
        scheme = str(model.setups.pgamma_scheme if model.setups.pgamma_scheme != "disabled" else model.fwd_rad.pgamma_scheme)
        if scheme.lower() not in HUMMER_SCHEMES:
            raise NotImplementedError("structured p-gamma/neutrino output supports only the Hummer2010 response kernel.")
    elif solver not in {"legacy_1d", "am3_1d"}:
        raise NotImplementedError("structured_backend='fortran_1d' supports hadronic_solver='legacy_1d' or 'am3_1d'.")


def _usechi(model) -> bool:
    return (
        str(model.setups.geometry_kernel).lower() == "chi_eats_2d"
        or str(model.setups.electron_solver).lower().endswith("_2d")
    )


def _checkchi(model, axisymmetric: bool, num_phi: int) -> None:
    if not axisymmetric:
        raise NotImplementedError("structured chi_eats_2d ring projection currently requires an axisymmetric jet profile.")
    if str(model.setups.geometry_kernel).lower() != "chi_eats_2d":
        raise NotImplementedError("structured 2d electron transport requires geometry_projection='chi_eats_2d'.")
    if str(model.setups.electron_solver).lower() != "fullhide_2d":
        raise NotImplementedError("structured chi_eats_2d ring projection currently requires electron_solver='fullhide_2d'.")
    if float(model.observer.theta_obs) != 0.0 and int(num_phi) < 2:
        raise ValueError("off-axis structured chi_eats_2d projection requires structured_num_phi >= 2.")
    if bool(model.fwd_rad.ssc):
        raise NotImplementedError("structured chi_eats_2d ring projection currently covers FS synchrotron+SSA, not SSC emission.")
    if bool(model.setups.rvs_shock):
        raise NotImplementedError("structured chi_eats_2d ring projection currently does not include reverse-shock emission.")
    if bool(model.setups.hadronic_enabled and model.fwd_rad.epsilon_p > 0.0):
        raise NotImplementedError("structured chi_eats_2d ring projection currently does not include hadronic emission.")
    _checkphysics(model, "chi_eats_2d")
    _checkspread(model, "chi_eats_2d")


def _samplegrid(
    model,
    min_gamma0: float = MINGAMMA1D,
    reject_transrelativistic: bool = True,
    obstime: np.ndarray | None = None,
    theta_count: int | None = None,
    phi_count: int | None = None,
):
    axisymmetric = axisjet(model.jet)
    grid = patchgrid(model, obstime, theta_count=theta_count, phi_count=phi_count)
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


def _threadplan(model) -> tuple[int, int]:
    mode = str(model.setups.structured_parallel_mode).lower()
    total = int(model.setups.num_threads)
    outer = model.setups.structured_outer_threads
    inner = model.setups.structured_inner_threads
    cpu_count = os.cpu_count()
    if total < 1:
        raise ValueError("num_threads must be positive for structured jet execution.")
    if cpu_count is not None and total > int(cpu_count):
        raise ValueError("num_threads exceeds the available CPU thread count for structured jet execution.")
    if mode not in {"outer", "inner", "nested"}:
        raise ValueError("structured_parallel_mode must be 'outer', 'inner', or 'nested'.")
    if mode == "nested" and (outer is None or inner is None):
        raise ValueError("structured_parallel_mode='nested' requires structured_outer_threads and structured_inner_threads.")
    reverse_hadronic = bool(
        model.setups.rvs_shock
        and model.rvs_rad is not None
        and model.setups.hadronic_enabled
        and model.fwd_rad.epsilon_p > 0.0
    )
    if reverse_hadronic:
        resolved_inner = int(total if inner is None else inner)
        if outer is not None and int(outer) != 1:
            raise ValueError("structured reverse-shock plus forward-hadronic execution requires structured_outer_threads=1.")
        _checkthreads(1, resolved_inner, total, cpu_count)
        return 1, resolved_inner
    if mode == "outer":
        resolved_outer = int(total if outer is None else outer)
        _checkthreads(resolved_outer, 1, total, cpu_count)
        return resolved_outer, 1
    if mode == "inner":
        resolved_inner = int(total if inner is None else inner)
        _checkthreads(1, resolved_inner, total, cpu_count)
        return 1, resolved_inner
    resolved_outer = int(outer)
    resolved_inner = int(inner)
    _checkthreads(resolved_outer, resolved_inner, total, cpu_count)
    return resolved_outer, resolved_inner


def _checkthreads(outer: int, inner: int, total: int, cpu_count: int | None) -> None:
    if outer < 1 or inner < 1:
        raise ValueError("structured outer and inner thread counts must be positive.")
    requested = outer * inner
    if requested > total:
        raise ValueError("structured outer_threads * inner_threads exceeds num_threads.")
    if cpu_count is not None and requested > int(cpu_count):
        raise ValueError("structured outer_threads * inner_threads exceeds available CPU threads.")


def _timegrid(model, requested_times: np.ndarray) -> np.ndarray:
    requested = np.asarray(requested_times, dtype=float)
    base_count = max(int(model.setups.num_tobs), int(np.unique(requested).size))
    if requested.size <= 1:
        return np.logspace(
            np.log10(float(10**model.setups.t_obs_min_log10)),
            np.log10(float(10**model.setups.t_obs_max_log10)),
            base_count,
        )
    tmin = min(float(10**model.setups.t_obs_min_log10), float(np.min(requested)))
    tmax = float(np.max(requested))
    solve_count = max(base_count, model._detail_time_count(tmin, tmax))
    logmin = np.log10(tmin)
    logmax = np.log10(tmax)
    logstep = (logmax - logmin) / float(solve_count - 2) if solve_count > 2 else 0.0
    return np.logspace(logmin, logmax + logstep, solve_count)


def _patchmeta(theta_centers, theta_edges, phi_centers, e_iso, gamma0, active, axisymmetric: bool, model, num_phi: int | None = None):
    count = int(model.setups.structured_num_phi if num_phi is None else num_phi)
    phi_values = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False) if axisymmetric else phi_centers
    theta_obs = float(model.observer.theta_obs)
    phi_obs = float(model.observer.phi_obs)
    patches = []
    for i_theta, theta in enumerate(theta_centers):
        theta_value = float(theta)
        for i_phi, phi in enumerate(phi_values):
            phi_value = float(phi)
            sourcephi = 0 if axisymmetric else i_phi
            if int(active[i_theta, sourcephi]) == 0:
                continue
            patches.append(
                {
                    "phi": phi_value,
                    "theta": theta_value,
                    "theta_lo": float(theta_edges[i_theta]),
                    "theta_hi": float(theta_edges[i_theta + 1]),
                    "theta_v": float(angsep(theta_value, phi_value, theta_obs, phi_obs)),
                    "patch_sampling": str(model.setups.patch_sampling).lower(),
                    "E_iso": float(e_iso[i_theta, sourcephi]),
                    "Gamma0": float(gamma0[i_theta, sourcephi]),
                }
            )
    return patches
