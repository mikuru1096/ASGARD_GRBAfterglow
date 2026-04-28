from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from time import perf_counter
from typing import Optional

import numpy as np

from src.Electron.electron_radiation import electron_radiation_kernel as electron_radiation_module
from asgard_core.asgard_coupling import build_coupled_shock_geometry, build_cross_zone_seed_fields
from asgard_core.asgard_config import ExecutionPolicy, FitConfig
from asgard_core.asgard_models import PhysicalSolution, SimulationSetup
from asgard_core.hadronic_pair_production import solve_pair_production
from asgard_core.hadronic_pgamma import photon_density_hz_to_gev
from asgard_core.asgard_types import (
    BranchState,
    FluxComponents,
    SolveState,
    ObsState,
    DynamicsSolution,
    ElectronSolution,
    PhotonFieldState,
    ObserverState,
    SolverAdapterReport,
    ReverseShockEmission,
)
from asgard_core.asgard_physics_utils import ambient_density, compute_doppler, compute_magnetic_field, compute_maximum_synchrotron_frequency
from asgard_core.asgard_postprocess import interpolate_observed_flux
from asgard_core.asgard_runtime import (
    _hadronic_advance_energy_loggamma,
    _hadronic_build_gamma_edges,
    _hadronic_electron_loss_rates,
    _hadronic_shell_dt,
    solve_dynamics,
    solve_electron,
    solve_hadronic,
    solve_reverse_shock_emission,
)
from asgard_core.asgard_ssc import compute_forward_ssc_seed_adaptive, compute_ssc_auxiliary_grid
from asgard_core.asgard_setup import build_simulation_setup
from src import Radiation, constants


_ELECTRON_MASS_GEV = constants.para_m_e_gev


def make_policy(config: FitConfig) -> ExecutionPolicy:
    return ExecutionPolicy(num_threads=config.num_threads)


def _requested_frequency_bounds(
    requested_frequencies_hz: np.ndarray | None,
) -> tuple[Optional[float], Optional[float]]:
    if requested_frequencies_hz is None:
        return None, None
    requested = np.asarray(requested_frequencies_hz, dtype=float)
    requested = requested[np.isfinite(requested) & (requested > 0.0)]
    if requested.size == 0:
        return None, None
    return float(np.min(requested)), float(np.max(requested))


def _required_solve_time_count(observer_time_s: np.ndarray, default_num_tobs: int) -> int:
    observer_time_s = np.asarray(observer_time_s, dtype=float)
    if observer_time_s.size == 0:
        return int(default_num_tobs)
    return max(int(default_num_tobs), int(np.unique(observer_time_s).size))


def make_tgrid(observer_time_s: np.ndarray, default_num_tobs: int) -> np.ndarray:
    observer_time_s = np.asarray(observer_time_s, dtype=float)
    if observer_time_s.size == 0:
        raise ValueError("observer_time_s must be non-empty.")
    if np.any(observer_time_s <= 0.0):
        raise ValueError("observer_time_s must be positive.")
    if observer_time_s.size == 1:
        return observer_time_s.copy()
    num_tobs = _required_solve_time_count(observer_time_s, default_num_tobs)
    t_min = float(np.min(observer_time_s))
    t_max = float(np.max(observer_time_s))
    return np.logspace(np.log10(t_min), np.log10(t_max), num_tobs)


def make_query_cfg(
    config: FitConfig,
    observer_time_s: np.ndarray,
) -> FitConfig:
    query = deepcopy(config)
    observer_time_s = np.asarray(observer_time_s, dtype=float)
    query.num_tobs = observer_time_s.shape[0]
    query.t_obs_min_log10 = float(np.log10(observer_time_s.min()))
    query.t_obs_max_log10 = float(np.log10(observer_time_s.max()))
    return query


def make_query_setup(
    config: FitConfig,
    observer_time_s: np.ndarray,
    requested_frequencies_hz: np.ndarray | None = None,
):
    setup = build_simulation_setup(make_query_cfg(config, observer_time_s))
    setup.observer_time_s = np.asarray(observer_time_s, dtype=float)
    # Keep the internal solve frequency grid invariant to query frequencies.
    # Otherwise, merely adding an output frequency can change electron/SSC
    # evolution and produce non-physical query-dependent solutions.
    _ = requested_frequencies_hz
    return setup


def solve_state(
    config: FitConfig,
    observer_time_s: np.ndarray,
    requested_frequencies_hz: np.ndarray | None = None,
    timings: Optional[dict[str, float]] = None,
    policy: Optional[ExecutionPolicy] = None,
) -> SolveState:
    setup = make_query_setup(config, observer_time_s, requested_frequencies_hz)
    return solve_state_from_setup(
        config,
        setup,
        timings=timings,
        policy=policy,
        requested_frequencies_hz=requested_frequencies_hz,
    )


def solve_spectra(
    config: FitConfig,
    observer_time_s: np.ndarray,
    requested_frequencies_hz: np.ndarray | None = None,
) -> FluxComponents:
    return solve_state(config, observer_time_s, requested_frequencies_hz).components


def _solver_label(config: FitConfig, stage: str) -> str:
    if stage == "dynamics":
        if config.reverse or config.reverse_shock.enabled:
            return "Dynamics.dynamics_reverse"
        return "Dynamics.dynamics_forward"
    if stage == "electron":
        solver_name = config.electron_solver.lower()
        electron_label_map = {
            "fullhide": "Electron.fs_electron_fullhide_1d",
            "fullhide_1d": "Electron.fs_electron_fullhide_1d",
            "fullhide_2d": "Electron.fs_electron_fullhide_2d",
            "t2g1": "Electron.fs_electron_t2g1_1d",
            "t2g1_1d": "Electron.fs_electron_t2g1_1d",
            "slc1": "Electron.fs_electron_slc1_1d",
            "slc1_1d": "Electron.fs_electron_slc1_1d",
            "charint": "Electron.fs_electron_charint_1d",
            "charint_1d": "Electron.fs_electron_charint_1d",
            "charint_2d": "Electron.fs_electron_charint_2d",
            "weno5": "Electron.fs_electron_weno5_1d",
            "weno5_1d": "Electron.fs_electron_weno5_1d",
        }
        return electron_label_map.get(solver_name, f"Electron.{solver_name}")
    if stage == "hadronic":
        return "Hadronic.fs_hadronic_1d"
    raise ValueError(f"Unsupported solver stage: {stage}")


def _solve_dynamics_stage(
    config: FitConfig,
    setup,
    timings: Optional[dict[str, float]],
) -> tuple[DynamicsSolution, SolverAdapterReport]:
    return _timed_call(
        timings,
        _solver_label(config, "dynamics"),
        solve_dynamics,
        setup.boundary,
        config,
        return_report=True,
    )


def _solve_electron_stage(
    config: FitConfig,
    setup,
    dynamics: DynamicsSolution,
    timings: Optional[dict[str, float]],
) -> tuple[ElectronSolution, SolverAdapterReport]:
    return _timed_call(
        timings,
        _solver_label(config, "electron"),
        solve_electron,
        setup.boundary,
        dynamics,
        setup.seed_frequency_hz,
        config,
        return_report=True,
    )


def _build_photon_field_stage(
    config: FitConfig,
    setup,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    timings: Optional[dict[str, float]],
) -> PhotonFieldState:
    forward_syn_seed = np.array(electron.seed_syn, dtype=float, copy=True)
    hadronic_forward_ssc_seed = np.zeros_like(forward_syn_seed)
    hadronic_target_seed = np.array(forward_syn_seed, dtype=float, copy=True)
    if config.hadronic.enabled and config.hadronic.epsilon_p > 0.0 and config.include_forward_ssc:
        _, hadronic_forward_ssc_seed = _ssc_spectrum(
            dynamics.radius,
            electron,
            setup.seed_frequency_hz,
            forward_syn_seed,
            config,
            config.num_threads,
            timings,
            "Radiation.ssc_spec [FS-Hadronic]",
        )
        hadronic_target_seed += hadronic_forward_ssc_seed
    return PhotonFieldState(
        seed_frequency_hz=np.asarray(setup.seed_frequency_hz, dtype=float),
        forward_syn_seed=forward_syn_seed,
        hadronic_forward_ssc_seed=hadronic_forward_ssc_seed,
        hadronic_target_seed=hadronic_target_seed,
        absorption_syn_seed=np.array(forward_syn_seed, dtype=float, copy=True),
        absorption_ssc_seed=np.zeros_like(forward_syn_seed),
    )


def _solve_hadronic_stage(
    config: FitConfig,
    setup,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    photon_field: PhotonFieldState,
    timings: Optional[dict[str, float]],
) -> tuple[ElectronSolution, Optional[object], SolverAdapterReport]:
    hadronic, report = _timed_call(
        timings,
        _solver_label(config, "hadronic"),
        solve_hadronic,
        setup.boundary,
        dynamics,
        electron,
        photon_field.seed_frequency_hz,
        photon_field.hadronic_target_seed,
        config,
        return_report=True,
    )
    if hadronic is not None and hadronic.d_n_gam_e_bh is not None:
        electron = _merge_bh_into_forward_electrons(
            electron,
            hadronic,
            dynamics.radius,
            compute_magnetic_field(dynamics.r_gamma, dynamics.radius, config),
            setup.seed_frequency_hz,
            config,
        )
        hadronic.l_had_bethe_heitler = None
        hadronic.seed_had_bethe_heitler = None
    if hadronic is not None and hadronic.pg_photon_survival is not None:
        _apply_hadronic_pg_local_closure(photon_field, hadronic.pg_photon_survival)
    return electron, hadronic, report


def _apply_hadronic_pg_local_closure(
    photon_field: PhotonFieldState,
    pg_photon_survival: np.ndarray,
) -> None:
    survival = np.asarray(pg_photon_survival, dtype=float)
    photon_field.hadronic_target_seed = np.asarray(photon_field.hadronic_target_seed, dtype=float) * survival
    photon_field.absorption_syn_seed = np.asarray(photon_field.absorption_syn_seed, dtype=float) * survival
    photon_field.absorption_ssc_seed = np.asarray(photon_field.absorption_ssc_seed, dtype=float) * survival
    photon_field.hadronic_forward_ssc_seed = np.asarray(photon_field.hadronic_forward_ssc_seed, dtype=float) * survival


def _solve_reverse_stage(
    config: FitConfig,
    setup,
    dynamics: DynamicsSolution,
    timings: Optional[dict[str, float]],
) -> Optional[ReverseShockEmission]:
    if not (config.reverse or config.reverse_shock.enabled):
        return None
    return _timed_call(
        timings,
        "Radiation.seed_reverse",
        solve_reverse_shock_emission,
        setup.boundary,
        dynamics,
        setup.seed_frequency_hz,
        config,
    )


def solve_state_from_setup(
    config: FitConfig,
    setup,
    timings: Optional[dict[str, float]] = None,
    policy: Optional[ExecutionPolicy] = None,
    requested_frequencies_hz: np.ndarray | None = None,
) -> SolveState:
    execution_policy = make_policy(config) if policy is None else policy
    dynamics, dynamics_report = _solve_dynamics_stage(config, setup, timings)
    electron, electron_report = _solve_electron_stage(config, setup, dynamics, timings)
    photon_field = _build_photon_field_stage(config, setup, dynamics, electron, timings)
    electron, hadronic, hadronic_report = _solve_hadronic_stage(
        config,
        setup,
        dynamics,
        electron,
        photon_field,
        timings,
    )
    reverse_emission = _solve_reverse_stage(config, setup, dynamics, timings)
    observer = _assemble_observer_stage(
        setup,
        config,
        dynamics,
        electron,
        photon_field,
        hadronic,
        reverse_emission,
        timings=timings,
    )
    freq_min, freq_max = _requested_frequency_bounds(requested_frequencies_hz)
    return SolveState(
        config=config,
        setup=setup,
        policy=execution_policy,
        dynamics=dynamics,
        electron=electron,
        photon_field=photon_field,
        hadronic=hadronic,
        observer=observer,
        reverse_emission=reverse_emission,
        components=observer.components,
        requested_frequency_min_hz=freq_min,
        requested_frequency_max_hz=freq_max,
        timings={} if timings is None else dict(timings),
        adapter_reports={
            "dynamics": dynamics_report,
            "electron": electron_report,
            "hadronic": hadronic_report,
        },
    )


def solve_spectra_from_setup(
    config: FitConfig,
    setup,
    timings: Optional[dict[str, float]] = None,
) -> FluxComponents:
    return solve_state_from_setup(config, setup, timings=timings).components


def _build_observer_setup_from_state(
    state: SolveState,
    observer_time_s: np.ndarray,
) -> SimulationSetup:
    return SimulationSetup(
        luminosity_distance_cm=state.setup.luminosity_distance_cm,
        boundary=state.setup.boundary,
        seed_frequency_hz=state.setup.seed_frequency_hz,
        observer_time_s=np.asarray(observer_time_s, dtype=float),
    )


def covers(
    state: SolveState,
    observer_time_s: np.ndarray,
    frequencies_hz: np.ndarray | None = None,
) -> bool:
    observer_time_s = np.asarray(observer_time_s, dtype=float)
    if observer_time_s.size > 0 and float(np.min(observer_time_s)) < float(np.min(state.setup.observer_time_s)):
        return False
    if observer_time_s.size > 0 and float(np.max(observer_time_s)) > float(np.max(state.setup.observer_time_s)):
        return False
    if observer_time_s.size > 1:
        required_count = _required_solve_time_count(observer_time_s, state.config.num_tobs)
        if int(state.setup.observer_time_s.size) < required_count:
            return False
    else:
        if observer_time_s.size == 1 and int(state.setup.observer_time_s.size) == 1:
            if not np.isclose(float(state.setup.observer_time_s[0]), float(observer_time_s[0]), rtol=0.0, atol=0.0):
                return False
    freq_min, freq_max = _requested_frequency_bounds(frequencies_hz)
    if freq_min is None or freq_max is None:
        return True
    if state.requested_frequency_min_hz is None or state.requested_frequency_max_hz is None:
        return False
    return freq_min >= state.requested_frequency_min_hz and freq_max <= state.requested_frequency_max_hz


def project_flux_grid(
    state: SolveState,
    observer_time_s: np.ndarray,
    frequencies_hz: np.ndarray,
    timings: Optional[dict[str, float]] = None,
    mode: str = "full_components",
) -> ObsState:
    setup = _build_observer_setup_from_state(state, observer_time_s)
    observed = observe_components_from_setup(
        state.config,
        state.components,
        setup,
        frequencies_hz,
        timings=timings,
        mode=mode,
    )
    return ObsState(
        state=state,
        setup=setup,
        frequencies_hz=np.asarray(frequencies_hz, dtype=float),
        components=observed,
    )


def project_spec(
    state: SolveState,
    time_s: float,
    frequencies_hz: np.ndarray,
    timings: Optional[dict[str, float]] = None,
    mode: str = "full_components",
) -> ObsState:
    return project_flux_grid(
        state,
        np.array([time_s], dtype=float),
        frequencies_hz,
        timings=timings,
        mode=mode,
    )


def state_details(state: SolveState) -> dict[str, Optional[BranchState]]:
    return {"fwd": state.components.fwd, "rev": state.components.rev}


def to_physical_solution(state: SolveState) -> PhysicalSolution:
    rev = state.components.rev
    return PhysicalSolution(
        characteristic_time_s=state.components.fwd.characteristic_time_s,
        gamma=state.components.fwd.gamma,
        radius_cm=state.components.fwd.radius_cm,
        absorbed_spectral_flux=state.components.total,
        nu_m=state.components.fwd.nu_m,
        nu_c=state.components.fwd.nu_c,
        nu_a=state.components.fwd.nu_a,
        rs_nu_m=None if rev is None else rev.nu_m,
        rs_nu_c=None if rev is None else rev.nu_c,
        rs_nu_a=None if rev is None else rev.nu_a,
    )


def observe_components(
    config: FitConfig,
    components: FluxComponents,
    observer_time_s: np.ndarray,
    frequencies_hz: np.ndarray,
) -> dict[str, np.ndarray]:
    setup = make_query_setup(config, observer_time_s, frequencies_hz)
    return observe_components_from_setup(config, components, setup, frequencies_hz)


def observe_components_from_setup(
    config: FitConfig,
    components: FluxComponents,
    setup,
    frequencies_hz: np.ndarray,
    timings: Optional[dict[str, float]] = None,
    mode: str = "full_components",
) -> dict[str, np.ndarray | None]:
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    if mode not in {"full_components", "total_only"}:
        raise ValueError(f"Unsupported observe mode: {mode}")
    if mode == "total_only":
        total = _project_component(
            setup,
            components.fwd.characteristic_time_s,
            components.fwd.gamma,
            components.fwd.radius_cm,
            setup.seed_frequency_hz,
            components.total,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [total]",
        )
        return {
            "total": total,
            "fwd_sync": None,
            "fwd_ssc": None,
            "fwd_hadronic": None,
            "fwd_hadronic_bethe_heitler": None,
            "fwd_hadronic_inverse_compton": None,
            "fwd_hadronic_pair_production": None,
            "rev_sync": None,
            "rev_ssc": None,
            "cross_ic": None,
        }

    observed = {
        "fwd_sync": _project_component(
            setup,
            components.fwd.characteristic_time_s,
            components.fwd.gamma,
            components.fwd.radius_cm,
            setup.seed_frequency_hz,
            components.fwd_sync,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [fwd_sync]",
        ),
        "fwd_ssc": _project_component(
            setup,
            components.fwd.characteristic_time_s,
            components.fwd.gamma,
            components.fwd.radius_cm,
            setup.seed_frequency_hz,
            components.fwd_ssc,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [fwd_ssc]",
        ),
    }
    observed["fwd_hadronic"] = None
    if components.fwd_hadronic_gamma is not None:
        observed["fwd_hadronic"] = _project_component(
            setup,
            components.fwd.characteristic_time_s,
            components.fwd.gamma,
            components.fwd.radius_cm,
            setup.seed_frequency_hz,
            components.fwd_hadronic_gamma,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [fwd_hadronic]",
        )
    observed["fwd_hadronic_bethe_heitler"] = None
    if components.fwd_hadronic_bethe_heitler is not None:
        observed["fwd_hadronic_bethe_heitler"] = _project_component(
            setup,
            components.fwd.characteristic_time_s,
            components.fwd.gamma,
            components.fwd.radius_cm,
            setup.seed_frequency_hz,
            components.fwd_hadronic_bethe_heitler,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [fwd_hadronic_bethe_heitler]",
        )
    observed["fwd_hadronic_inverse_compton"] = None
    if components.fwd_hadronic_inverse_compton is not None:
        observed["fwd_hadronic_inverse_compton"] = _project_component(
            setup,
            components.fwd.characteristic_time_s,
            components.fwd.gamma,
            components.fwd.radius_cm,
            setup.seed_frequency_hz,
            components.fwd_hadronic_inverse_compton,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [fwd_hadronic_inverse_compton]",
        )
    observed["fwd_hadronic_pair_production"] = None
    if components.fwd_hadronic_pair_production is not None:
        observed["fwd_hadronic_pair_production"] = _project_component(
            setup,
            components.fwd.characteristic_time_s,
            components.fwd.gamma,
            components.fwd.radius_cm,
            setup.seed_frequency_hz,
            components.fwd_hadronic_pair_production,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [fwd_hadronic_pair_production]",
        )
    observed["rev_sync"] = None
    observed["rev_ssc"] = None
    observed["cross_ic"] = None
    if components.rev_sync is not None:
        observed["rev_sync"] = _project_component(
            setup,
            components.fwd.characteristic_time_s,
            components.fwd.gamma,
            components.fwd.radius_cm,
            setup.seed_frequency_hz,
            components.rev_sync,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [rev_sync]",
        )
    if components.rev_ssc is not None:
        observed["rev_ssc"] = _project_component(
            setup,
            components.fwd.characteristic_time_s,
            components.fwd.gamma,
            components.fwd.radius_cm,
            setup.seed_frequency_hz,
            components.rev_ssc,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [rev_ssc]",
        )
    if components.cross_ic is not None:
        observed["cross_ic"] = _project_component(
            setup,
            components.fwd.characteristic_time_s,
            components.fwd.gamma,
            components.fwd.radius_cm,
            setup.seed_frequency_hz,
            components.cross_ic,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [cross_ic]",
        )
    observed["total"] = np.array(observed["fwd_sync"], dtype=float, copy=True)
    observed["total"] += np.array(observed["fwd_ssc"], dtype=float, copy=False)
    if observed["fwd_hadronic"] is not None:
        observed["total"] += np.array(observed["fwd_hadronic"], dtype=float, copy=False)
    if observed["rev_sync"] is not None:
        observed["total"] += np.array(observed["rev_sync"], dtype=float, copy=False)
    if observed["rev_ssc"] is not None:
        observed["total"] += np.array(observed["rev_ssc"], dtype=float, copy=False)
    if observed["cross_ic"] is not None:
        observed["total"] += np.array(observed["cross_ic"], dtype=float, copy=False)
    if timings is not None:
        timings.setdefault("Interpolation.sed_interpolation [total]", 0.0)
    return observed


def _assemble_observer_stage(
    setup,
    config: FitConfig,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    photon_field: PhotonFieldState,
    hadronic,
    reverse_emission: Optional[ReverseShockEmission],
    timings: Optional[dict[str, float]] = None,
) -> ObserverState:
    pg_photon_survival = None if hadronic is None else hadronic.pg_photon_survival
    fwd_sync = np.asarray(electron.l_syn_spec, dtype=float)
    if pg_photon_survival is not None:
        fwd_sync = fwd_sync * np.asarray(pg_photon_survival, dtype=float)
    seed_syn_absorption = np.array(photon_field.absorption_syn_seed, dtype=float, copy=True)
    zero_component = np.zeros_like(fwd_sync)
    fwd_ssc = zero_component
    seed_ssc_total = np.array(photon_field.absorption_ssc_seed, dtype=float, copy=True)
    if config.include_forward_ssc:
        fwd_ssc, seed_ssc_total = _ssc_spectrum(
            dynamics.radius,
            electron,
            setup.seed_frequency_hz,
            seed_syn_absorption,
            config,
            config.num_threads,
            timings,
            "Radiation.ssc_spec [FS]",
        )
        if pg_photon_survival is not None:
            fwd_ssc = np.asarray(fwd_ssc, dtype=float) * np.asarray(pg_photon_survival, dtype=float)
            seed_ssc_total = np.asarray(seed_ssc_total, dtype=float) * np.asarray(pg_photon_survival, dtype=float)

    rev_sync = None
    rev_ssc = None
    cross_ic = None
    rev_details = None
    if reverse_emission is not None:
        rev_sync = reverse_emission.l_syn_spec
        seed_syn_absorption = seed_syn_absorption + reverse_emission.seed_syn
        if reverse_emission.rs_hadronic is not None:
            seed_syn_absorption = seed_syn_absorption + reverse_emission.rs_hadronic.seed_had_syn
            rev_sync = rev_sync + reverse_emission.rs_hadronic.l_had_syn_spec
        rev_details = BranchState(
            characteristic_time_s=dynamics.r_tobs,
            gamma=dynamics.r_gamma,
            radius_cm=dynamics.radius,
            swept_mass_g=dynamics.reverse_shock.swept_mass_g,
            doppler=_compute_doppler(dynamics.r_gamma, config.z),
            magnetic_field_g=reverse_emission.magnetic_field_g,
            nu_m=reverse_emission.nu_m,
            nu_c=reverse_emission.nu_c,
            nu_a=reverse_emission.nu_a,
            nu_M=reverse_emission.nu_M,
        )
        if config.reverse_shock.include_ssc:
            rev_ssc, seed_ssc_rs = _timed_call(
                timings,
                "Radiation.ssc_spec [RS-SSC]",
                Radiation.ssc_spec,
                dynamics.radius,
                dynamics.reverse_shock.gam_e,
                dynamics.reverse_shock.d_n_gam_e,
                setup.seed_frequency_hz,
                reverse_emission.seed_syn,
                config.num_threads,
            )
            seed_ssc_total = seed_ssc_total + seed_ssc_rs
        if config.reverse_shock.include_cross_zone_ic:
            coupling_geometry = build_coupled_shock_geometry(dynamics, config)
            seed_fs_to_rs, seed_rs_to_fs = build_cross_zone_seed_fields(
                electron.seed_syn,
                reverse_emission.seed_syn,
                coupling_geometry,
            )
            l_cic_fs_spec, seed_cic_fs = _ssc_spectrum(
                dynamics.radius,
                electron,
                setup.seed_frequency_hz,
                seed_rs_to_fs,
                config,
                config.num_threads,
                timings,
                "Radiation.ssc_spec [CIC-FS]",
            )
            l_cic_rs_spec, seed_cic_rs = _timed_call(
                timings,
                "Radiation.ssc_spec [CIC-RS]",
                Radiation.ssc_spec,
                dynamics.radius,
                dynamics.reverse_shock.gam_e,
                dynamics.reverse_shock.d_n_gam_e,
                setup.seed_frequency_hz,
                seed_fs_to_rs,
                config.num_threads,
            )
            cross_ic = l_cic_fs_spec + l_cic_rs_spec
            seed_ssc_total = seed_ssc_total + seed_cic_fs + seed_cic_rs

    tau_extra = np.zeros_like(fwd_sync)
    absorbed_fwd_hadronic_gamma = None
    absorbed_fwd_hadronic_bethe_heitler = None
    absorbed_fwd_hadronic_inverse_compton = None
    absorbed_fwd_hadronic_pair_production = None
    magnetic_field_g = compute_magnetic_field(dynamics.r_gamma, dynamics.radius, config)
    hadronic_ssa_transfer = np.ones_like(fwd_sync)
    if hadronic is not None:
        hadronic_ssa_transfer = _forward_synchrotron_absorption_transfer(
            electron=electron,
            radius_cm=dynamics.radius,
            magnetic_field_g=magnetic_field_g,
            seed_frequency_hz=setup.seed_frequency_hz,
            config=config,
        )
        seed_syn_absorption = seed_syn_absorption + _hadronic_absorbed_seed_density(
            hadronic=hadronic,
            radius_cm=dynamics.radius,
            seed_frequency_hz=setup.seed_frequency_hz,
            ssa_transfer=hadronic_ssa_transfer,
        )

    pair_seed_total = np.zeros_like(seed_syn_absorption)
    pair_lum_total = np.zeros_like(fwd_sync)
    tau_pair = np.zeros_like(fwd_sync)
    if bool(config.hadronic.include_pair_production):
        pair_lum_total, pair_seed_total, tau_pair = _compute_pair_production_branch(
            dynamics=dynamics,
            electron=electron,
            combined_seed_field_hz=seed_syn_absorption + seed_ssc_total,
            seed_frequency_hz=setup.seed_frequency_hz,
            magnetic_field_g=magnetic_field_g,
            config=config,
        )
        seed_syn_absorption = seed_syn_absorption + pair_seed_total
    tau_extra = tau_extra + tau_pair
    annihilation_seed_syn = seed_syn_absorption
    annihilation_seed_ssc = seed_ssc_total
    if bool(config.hadronic.include_pair_production):
        annihilation_seed_syn = np.zeros_like(seed_syn_absorption)
        annihilation_seed_ssc = np.zeros_like(seed_ssc_total)
    absorption = _timed_call(
        timings,
        "Radiation.annihilation",
        Radiation.annihilation,
        dynamics.r_gamma,
        dynamics.radius,
        setup.seed_frequency_hz,
        annihilation_seed_syn,
        annihilation_seed_ssc,
        tau_extra,
        config.num_threads,
    )
    prefactor = absorption / (4.0 * np.pi * setup.luminosity_distance_cm**2) * (1.0 + config.z)

    absorbed_fwd_sync = fwd_sync * prefactor
    absorbed_fwd_ssc = fwd_ssc * prefactor
    if hadronic is not None:
        hadronic_luminosity = _hadronic_absorbed_luminosity(hadronic, hadronic_ssa_transfer)
        hadronic_gamma_total = hadronic_luminosity["total"]
        absorbed_fwd_hadronic_bethe_heitler = _project_optional_luminosity(
            hadronic_luminosity["bethe_heitler"],
            prefactor,
        )
        absorbed_fwd_hadronic_inverse_compton = _project_optional_luminosity(
            hadronic_luminosity["inverse_compton"],
            prefactor,
        )
        absorbed_fwd_hadronic_pair_production = _project_optional_luminosity(
            hadronic_luminosity["pair_production"],
            prefactor,
        )
        if bool(config.hadronic.include_pair_production):
            hadronic_gamma_total += pair_lum_total
            absorbed_fwd_hadronic_pair_production = pair_lum_total * prefactor
        absorbed_fwd_hadronic_gamma = hadronic_gamma_total * prefactor
    absorbed_rev_sync = None if rev_sync is None else rev_sync * prefactor
    absorbed_rev_ssc = None if rev_ssc is None else rev_ssc * prefactor
    absorbed_cross_ic = None if cross_ic is None else cross_ic * prefactor

    total = absorbed_fwd_sync + absorbed_fwd_ssc
    if absorbed_rev_sync is not None:
        total = total + absorbed_rev_sync
    if absorbed_rev_ssc is not None:
        total = total + absorbed_rev_ssc
    if absorbed_cross_ic is not None:
        total = total + absorbed_cross_ic
    if absorbed_fwd_hadronic_gamma is not None:
        total = total + absorbed_fwd_hadronic_gamma

    return ObserverState(
        prefactor=np.asarray(prefactor, dtype=float),
        tau_extra=np.asarray(tau_extra, dtype=float),
        tau_pair=np.asarray(tau_pair, dtype=float),
        components=FluxComponents(
            total=total,
            fwd_sync=absorbed_fwd_sync,
            fwd_ssc=absorbed_fwd_ssc,
            fwd_hadronic_gamma=absorbed_fwd_hadronic_gamma,
            fwd_hadronic_bethe_heitler=absorbed_fwd_hadronic_bethe_heitler,
            fwd_hadronic_inverse_compton=absorbed_fwd_hadronic_inverse_compton,
            fwd_hadronic_pair_production=absorbed_fwd_hadronic_pair_production,
            rev_sync=absorbed_rev_sync,
            rev_ssc=absorbed_rev_ssc,
            cross_ic=absorbed_cross_ic,
            fwd=BranchState(
                characteristic_time_s=dynamics.r_tobs,
                gamma=dynamics.r_gamma,
                radius_cm=dynamics.radius,
                swept_mass_g=dynamics.swept_mass_g,
                doppler=_compute_doppler(dynamics.r_gamma, config.z),
                magnetic_field_g=_compute_forward_magnetic_field(dynamics.r_gamma, dynamics.radius, config),
                nu_m=electron.nu_m,
                nu_c=electron.nu_c,
                nu_a=electron.nu_a,
                nu_M=_compute_maximum_synchrotron_frequency(dynamics.r_gamma, dynamics.radius, config),
                cooling_timescale_s=electron.cooling_timescale_s,
                dynamical_timescale_s=electron.dynamical_timescale_s,
            ),
            rev=rev_details,
        ),
    )


def _compute_pair_production_branch(
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    combined_seed_field_hz: np.ndarray,
    seed_frequency_hz: np.ndarray,
    magnetic_field_g: np.ndarray,
    config: FitConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    v_seed = np.asarray(seed_frequency_hz, dtype=float)
    seed_field = np.asarray(combined_seed_field_hz, dtype=float)
    gam_e = np.asarray(electron.gam_e, dtype=float)
    e_e_gev = gam_e * _ELECTRON_MASS_GEV
    gam_edge = _hadronic_build_gamma_edges(gam_e)
    num_nu = int(v_seed.size)
    num_r = int(np.asarray(dynamics.radius, dtype=float).size)
    pair_lum = np.zeros((num_nu, num_r), dtype=float)
    pair_seed = np.zeros((num_nu, num_r), dtype=float)
    tau_pair = np.zeros((num_nu, num_r), dtype=float)
    d_n_pair_prev = np.zeros(gam_e.size, dtype=float)
    photon_energy_gev, _ = photon_density_hz_to_gev(v_seed, np.ones_like(v_seed))

    use_iterative_cascade = (
        int(getattr(config.hadronic, 'pair_cascade_iterations', 1)) > 1
    )

    for i_r in range(num_r):
        _, photon_density_per_gev = photon_density_hz_to_gev(v_seed, seed_field[:, i_r])

        if use_iterative_cascade:
            from asgard_core.hadronic_cascade import compute_iterative_pair_cascade
            c_out = compute_iterative_pair_cascade(
                photon_energy_gev=photon_energy_gev,
                photon_density_per_gev=photon_density_per_gev,
                electron_energy_gev=e_e_gev,
                radius_cm=float(dynamics.radius[i_r]),
                gamma_bulk=float(dynamics.r_gamma[i_r]),
                b_field_g=float(magnetic_field_g[i_r]),
                max_iterations=int(config.hadronic.pair_cascade_iterations),
            )
            pair_lum_shell = np.asarray(c_out.pair_syn_luminosity_hz, dtype=float) * (
                4.0 * np.pi * float(dynamics.radius[i_r])**2
            )
            pair_lum[:, i_r] = pair_lum_shell
            pair_seed[:, i_r] = pair_lum_shell / (
                max(float(dynamics.radius[i_r])**2, 1e-60)
                * 4.0 * np.pi * constants.para_c * constants.para_h
            )
            tau_pair[:, i_r] = np.asarray(c_out.tau_pair_path, dtype=float)
        else:
            ppair = solve_pair_production(
                photon_energy_gev=photon_energy_gev,
                photon_density_per_gev=photon_density_per_gev,
                electron_energy_gev=e_e_gev,
            )
            tau_pair[:, i_r] = np.maximum(np.asarray(ppair.photon_loss_rate, dtype=float), 0.0) * (
                float(dynamics.radius[i_r]) / (12.0 * max(float(dynamics.r_gamma[i_r]), 1.0) * constants.para_c)
            )
            q_pair = np.asarray(ppair.pair_injection_rate_per_gev_total, dtype=float) * (
                (4.0 / 3.0) * np.pi * (
                    float(dynamics.radius[i_r]) ** 3
                    - (0.0 if i_r == 0 else float(dynamics.radius[i_r - 1]) ** 3)
                )
            ) * _ELECTRON_MASS_GEV
            d_n_pair = _hadronic_advance_energy_loggamma(
                gam_e, gam_edge, d_n_pair_prev, q_pair,
                _hadronic_electron_loss_rates(
                    gam_e, float(magnetic_field_g[i_r]),
                    float(dynamics.radius[i_r]) / (max(float(dynamics.r_gamma[i_r]), 1.0) * constants.para_c),
                ),
                _hadronic_shell_dt(np.asarray(dynamics.r_tobs, dtype=float), i_r),
            )
            p_syn_i, seed_syn_i = electron_radiation_module.get_syn_selected(
                int(config.index_syn_integr),
                float(dynamics.radius[i_r]),
                float(max(magnetic_field_g[i_r], 1.0e-30)),
                int(config.num_threads),
                gam_e, d_n_pair, v_seed,
            )
            pair_lum[:, i_r] = np.asarray(p_syn_i, dtype=float)
            pair_seed[:, i_r] = np.asarray(seed_syn_i, dtype=float)
            d_n_pair_prev = d_n_pair
    return pair_lum, pair_seed, tau_pair


def _forward_synchrotron_absorption_transfer(
    *,
    electron: ElectronSolution,
    radius_cm: np.ndarray,
    magnetic_field_g: np.ndarray,
    seed_frequency_hz: np.ndarray,
    config: FitConfig,
) -> np.ndarray:
    radius = np.asarray(radius_cm, dtype=float)
    magnetic_field = np.asarray(magnetic_field_g, dtype=float)
    frequency = np.asarray(seed_frequency_hz, dtype=float)
    if np.any(radius <= 0.0):
        raise ValueError("forward synchrotron absorption transfer requires positive shell radii.")
    if np.any(frequency <= 0.0):
        raise ValueError("forward synchrotron absorption transfer requires positive photon frequencies.")

    transfer = np.ones((frequency.size, radius.size), dtype=float)
    gam_e = np.asarray(electron.gam_e, dtype=float)
    d_n_gam_e = np.asarray(electron.d_n_gam_e, dtype=float)
    for i_shell in range(radius.size):
        if magnetic_field[i_shell] <= 0.0:
            continue
        transfer[:, i_shell] = np.asarray(
            electron_radiation_module.get_syn_transfer(
                float(radius[i_shell]),
                float(magnetic_field[i_shell]),
                int(config.num_threads),
                gam_e,
                d_n_gam_e[:, i_shell],
                frequency,
            ),
            dtype=float,
        )
    return transfer


def _apply_local_synchrotron_absorption(luminosity: np.ndarray, ssa_transfer: np.ndarray) -> np.ndarray:
    return np.asarray(luminosity, dtype=float) * np.asarray(ssa_transfer, dtype=float)


def _project_optional_luminosity(luminosity: np.ndarray | None, prefactor: np.ndarray) -> np.ndarray | None:
    if luminosity is None:
        return None
    return np.asarray(luminosity, dtype=float) * np.asarray(prefactor, dtype=float)


def _hadronic_optional_luminosities(hadronic) -> tuple[tuple[str, np.ndarray | None], ...]:
    return (
        ("bethe_heitler", hadronic.l_had_bethe_heitler),
        ("inverse_compton", hadronic.l_had_hadronic_inverse_compton),
        ("pair_production", hadronic.l_had_pair_production),
        ("pion_synch", hadronic.l_had_pion_synch),
        ("muon_synch", hadronic.l_had_muon_synch),
        ("pion_inverse_compton", hadronic.l_had_pion_inverse_compton),
        ("muon_inverse_compton", hadronic.l_had_muon_inverse_compton),
    )


def _hadronic_absorbed_luminosity(hadronic, ssa_transfer: np.ndarray) -> dict[str, np.ndarray | None]:
    base = np.asarray(hadronic.l_had_syn_spec + hadronic.l_had_pg_gamma, dtype=float)
    out: dict[str, np.ndarray | None] = {"total": _apply_local_synchrotron_absorption(base, ssa_transfer)}
    for name, luminosity in _hadronic_optional_luminosities(hadronic):
        out[name] = None
        if luminosity is not None:
            out[name] = _apply_local_synchrotron_absorption(luminosity, ssa_transfer)
            out["total"] = np.asarray(out["total"], dtype=float) + out[name]
    return out


def _seed_density_from_luminosity(
    luminosity: np.ndarray,
    *,
    radius_cm: np.ndarray,
    seed_frequency_hz: np.ndarray,
) -> np.ndarray:
    radius = np.asarray(radius_cm, dtype=float)
    frequency = np.asarray(seed_frequency_hz, dtype=float)
    lum = np.asarray(luminosity, dtype=float)
    if lum.shape != (frequency.size, radius.size):
        raise ValueError("photon luminosity grid must have shape (num_frequency, num_shell).")
    if np.any(radius <= 0.0):
        raise ValueError("photon seed conversion requires positive shell radii.")
    if np.any(frequency <= 0.0):
        raise ValueError("photon seed conversion requires positive photon frequencies.")
    denominator = (radius[None, :] ** 2) * frequency[:, None] * (4.0 * np.pi * constants.para_c * constants.para_h)
    return lum / denominator


def _hadronic_absorbed_seed_density(
    *,
    hadronic,
    radius_cm: np.ndarray,
    seed_frequency_hz: np.ndarray,
    ssa_transfer: np.ndarray,
) -> np.ndarray:
    seed_total = np.asarray(hadronic.seed_had_syn, dtype=float) * np.asarray(ssa_transfer, dtype=float)
    if hadronic.seed_had_bethe_heitler is not None:
        seed_total = seed_total + np.asarray(hadronic.seed_had_bethe_heitler, dtype=float) * np.asarray(ssa_transfer, dtype=float)

    for luminosity in (hadronic.l_had_pg_gamma, *[item[1] for item in _hadronic_optional_luminosities(hadronic)]):
        if luminosity is None:
            continue
        escaped_luminosity = _apply_local_synchrotron_absorption(luminosity, ssa_transfer)
        seed_total = seed_total + _seed_density_from_luminosity(
            escaped_luminosity,
            radius_cm=radius_cm,
            seed_frequency_hz=seed_frequency_hz,
        )
    return seed_total


def _merge_bh_into_forward_electrons(
    electron: ElectronSolution,
    hadronic,
    radius_cm: np.ndarray,
    magnetic_field_g: np.ndarray,
    seed_frequency_hz: np.ndarray,
    config: FitConfig,
) -> ElectronSolution:
    bh_distribution = np.asarray(hadronic.d_n_gam_e_bh, dtype=float)
    total_distribution = np.asarray(electron.d_n_gam_e, dtype=float) + bh_distribution
    num_shell = total_distribution.shape[1]
    num_nu = int(np.asarray(seed_frequency_hz, dtype=float).size)
    l_syn_total = np.zeros((num_nu, num_shell), dtype=float)
    seed_syn_total = np.zeros((num_nu, num_shell), dtype=float)
    for i_shell in range(num_shell):
        p_syn_i, seed_syn_i = electron_radiation_module.get_syn_selected(
            int(config.index_syn_integr),
            float(radius_cm[i_shell]),
            float(max(magnetic_field_g[i_shell], 1.0e-30)),
            int(config.num_threads),
            np.asarray(electron.gam_e, dtype=float),
            np.asarray(total_distribution[:, i_shell], dtype=float),
            np.asarray(seed_frequency_hz, dtype=float),
        )
        l_syn_total[:, i_shell] = np.asarray(p_syn_i, dtype=float)
        seed_syn_total[:, i_shell] = np.asarray(seed_syn_i, dtype=float)
    return ElectronSolution(
        gam_e=np.asarray(electron.gam_e, dtype=float),
        d_n_gam_e=np.asarray(total_distribution, dtype=float),
        l_syn_spec=l_syn_total,
        seed_syn=seed_syn_total,
        nu_m=np.asarray(electron.nu_m, dtype=float),
        nu_c=np.asarray(electron.nu_c, dtype=float),
        nu_a=np.asarray(electron.nu_a, dtype=float),
        d_n_gam_e_bh=bh_distribution,
        d_n_gam_e_chi=None if electron.d_n_gam_e_chi is None else np.asarray(electron.d_n_gam_e_chi, dtype=float),
        chi_grid=None if electron.chi_grid is None else np.asarray(electron.chi_grid, dtype=float),
        cooling_timescale_s=None if electron.cooling_timescale_s is None else np.asarray(electron.cooling_timescale_s, dtype=float),
        dynamical_timescale_s=None if electron.dynamical_timescale_s is None else np.asarray(electron.dynamical_timescale_s, dtype=float),
        work_x_edge_log10=None if electron.work_x_edge_log10 is None else np.asarray(electron.work_x_edge_log10, dtype=float),
        work_d_n_x=None if electron.work_d_n_x is None else np.asarray(electron.work_d_n_x, dtype=float),
    )


def _ssc_spectrum(
    radius_cm: np.ndarray,
    electron: ElectronSolution,
    seed_frequency_hz: np.ndarray,
    seed_field: np.ndarray,
    config: FitConfig,
    num_threads: int,
    timings: Optional[dict[str, float]],
    label: str,
) -> tuple[np.ndarray, np.ndarray]:
    if electron.work_x_edge_log10 is not None and electron.work_d_n_x is not None:
        if config.include_forward_ssc:
            return _timed_call(
                timings,
                label,
                compute_forward_ssc_seed_adaptive,
                radius_cm,
                electron.work_x_edge_log10,
                electron.work_d_n_x,
                seed_frequency_hz,
                seed_field,
                electron.gamma,
                electron.nu_a,
                electron.nu_m,
                electron.nu_c,
                config,
            )
        return _timed_call(
            timings,
            label,
            compute_ssc_auxiliary_grid,
            radius_cm,
            electron.work_x_edge_log10,
            electron.work_d_n_x,
            seed_frequency_hz,
            seed_field,
            num_threads,
        )
    return _timed_call(
        timings,
        label,
        Radiation.ssc_spec,
        radius_cm,
        electron.gam_e,
        electron.d_n_gam_e,
        seed_frequency_hz,
        seed_field,
        num_threads,
    )


def _project_component(
    setup,
    characteristic_time_s: np.ndarray,
    gamma: np.ndarray,
    radius_cm: np.ndarray,
    seed_frequency_hz: np.ndarray,
    absorbed_spectral_flux: np.ndarray,
    frequencies_hz: np.ndarray,
    config: FitConfig,
    timings: Optional[dict[str, float]] = None,
    label: Optional[str] = None,
) -> np.ndarray:
    if not np.any(absorbed_spectral_flux):
        if timings is not None and label is not None:
            timings.setdefault(label, 0.0)
        return np.zeros((frequencies_hz.shape[0], setup.observer_time_s.shape[0]), dtype=float)
    return _timed_call(
        timings,
        label,
        interpolate_observed_flux,
        setup,
        characteristic_time_s,
        gamma,
        radius_cm,
        absorbed_spectral_flux,
        frequencies_hz,
        config,
    )


def profile_observe_setup(
    config: FitConfig,
    setup,
    frequencies_hz: np.ndarray,
) -> tuple[FluxComponents, dict[str, np.ndarray], dict[str, float]]:
    timings: dict[str, float] = {}
    components = solve_spectra_from_setup(config, setup, timings=timings)
    observed = observe_components_from_setup(config, components, setup, frequencies_hz, timings=timings)
    return components, observed, timings


def _timed_call(timings: Optional[dict[str, float]], label: Optional[str], func, *args, **kwargs):
    start = perf_counter()
    result = func(*args, **kwargs)
    elapsed = perf_counter() - start
    if timings is not None and label is not None:
        timings[label] = timings.get(label, 0.0) + elapsed
    return result


def _compute_doppler(gamma: np.ndarray, redshift: float) -> np.ndarray:
    """DEPRECATED: Use asgard_physics_utils.compute_doppler instead."""
    return compute_doppler(gamma, redshift)


def _ambient_density(radius_cm: np.ndarray, config: FitConfig) -> np.ndarray:
    """DEPRECATED: Use asgard_physics_utils.ambient_density instead."""
    return ambient_density(radius_cm, config)


def _compute_forward_magnetic_field(gamma: np.ndarray, radius_cm: np.ndarray, config: FitConfig) -> np.ndarray:
    """DEPRECATED: Use asgard_physics_utils.compute_magnetic_field instead."""
    return compute_magnetic_field(gamma, radius_cm, config)


def _compute_maximum_synchrotron_frequency(gamma: np.ndarray, radius_cm: np.ndarray, config: FitConfig) -> np.ndarray:
    """DEPRECATED: Use asgard_physics_utils.compute_maximum_synchrotron_frequency instead."""
    return compute_maximum_synchrotron_frequency(gamma, radius_cm, config)
