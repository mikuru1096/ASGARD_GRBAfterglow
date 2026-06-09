from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from time import perf_counter
from typing import Optional

import numpy as np

from src.Electron.electron_radiation import electron_radiation_kernel as electron_radiation_module
from asgard_core.asgard_coupling import build_coupled_shock_geometry, build_cross_zone_seed_fields
from asgard_core.asgard_config import ExecutionPolicy, FitConfig, SimulationSetup
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
from src import Interpolation, Radiation, constants


_ELECTRON_MASS_GEV = constants.para_m_e_gev
_PROJECTION_KINDS = {"lightcurve", "sed"}


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
    setup = build_simulation_setup(make_query_cfg(config, observer_time_s), requested_frequencies_hz)
    setup.observer_time_s = np.asarray(observer_time_s, dtype=float)
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
    projection_kind: str = "lightcurve",
) -> ObsState:
    projection_kind = _normalize_projection_kind(projection_kind)
    setup = _build_observer_setup_from_state(state, observer_time_s)
    if _uses_chi_eats_2d(state.config) and projection_kind == "lightcurve":
        observed = _observe_components_chi_eats_2d(state, setup, frequencies_hz, timings=timings, mode=mode)
    else:
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


def _normalize_projection_kind(projection_kind: str) -> str:
    kind = str(projection_kind).lower()
    if kind not in _PROJECTION_KINDS:
        raise ValueError("projection_kind must be 'lightcurve' or 'sed'.")
    return kind


def _uses_chi_eats_2d(config: FitConfig) -> bool:
    return str(config.geometry_kernel).lower() == "chi_eats_2d"


def _require_chi_eats_electron_state(state: SolveState) -> None:
    if not str(state.config.electron_solver).lower().endswith("_2d"):
        raise ValueError("geometry_kernel='chi_eats_2d' requires a 2d electron solver.")
    missing = [
        name
        for name in (
            "l_syn_spec_chi",
            "tau_syn_chi",
            "chi_radius_cm",
            "chi_gamma_bulk",
            "chi_dvolume_weight",
        )
        if getattr(state.electron, name) is None
    ]
    if missing:
        raise RuntimeError("chi_eats_2d electron state is missing: " + ", ".join(missing))


def _require_top_hat_phi_grid(config: FitConfig) -> None:
    if float(config.theta_v) != 0.0 and int(config.num_phi) == 1:
        raise ValueError("off-axis EATS projection requires num_phi >= 2; num_phi=1 is only valid for on-axis axial collapse.")


def _project_chi_fwd_sync(
    state: SolveState,
    setup,
    frequencies_hz: np.ndarray,
    timings: Optional[dict[str, float]],
) -> np.ndarray:
    _require_chi_eats_electron_state(state)
    _require_top_hat_phi_grid(state.config)
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    order = np.argsort(frequencies_hz)
    sorted_frequencies = frequencies_hz[order]
    source_chi = np.asarray(state.electron.l_syn_spec_chi, dtype=float) * np.asarray(state.observer.prefactor, dtype=float)[:, None, :]
    num_phi = 1 if float(state.config.theta_v) == 0.0 else int(state.config.num_phi)
    flux_sorted = _timed_call(
        timings,
        "Interpolation.sed_interpolation_chi [fwd_sync]",
        Interpolation.sed_interpolation_chi,
        setup.boundary,
        state.components.fwd.characteristic_time_s,
        state.components.fwd.radius_cm,
        source_chi,
        np.asarray(state.electron.tau_syn_chi, dtype=float),
        np.asarray(state.electron.chi_radius_cm, dtype=float),
        np.asarray(state.electron.chi_gamma_bulk, dtype=float),
        np.asarray(state.electron.chi_dvolume_weight, dtype=float),
        setup.seed_frequency_hz,
        sorted_frequencies,
        setup.observer_time_s,
        state.config.num_theta,
        num_phi,
        state.config.num_threads,
    )
    if np.array_equal(order, np.arange(order.shape[0])):
        return flux_sorted
    flux_matrix = np.empty_like(flux_sorted)
    flux_matrix[order] = flux_sorted
    return flux_matrix


def _observe_components_chi_eats_2d(
    state: SolveState,
    setup,
    frequencies_hz: np.ndarray,
    timings: Optional[dict[str, float]],
    mode: str,
) -> dict[str, np.ndarray | None]:
    chi_fwd_sync = _project_chi_fwd_sync(state, setup, frequencies_hz, timings)
    if mode == "total_only":
        non_chi_total = np.asarray(state.components.total, dtype=float) - np.asarray(state.components.fwd_sync, dtype=float)
        shell_total_without_fwd_sync = _project_component(
            setup,
            state.components.fwd.characteristic_time_s,
            state.components.fwd.gamma,
            state.components.fwd.radius_cm,
            setup.seed_frequency_hz,
            non_chi_total,
            frequencies_hz,
            state.config,
            timings=timings,
            label="Interpolation.sed_interpolation [total_without_fwd_sync]",
        )
        return {
            "total": shell_total_without_fwd_sync + chi_fwd_sync,
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
    observed = observe_components_from_setup(
        state.config,
        state.components,
        setup,
        frequencies_hz,
        timings=timings,
        mode="full_components",
    )
    observed["fwd_sync"] = chi_fwd_sync
    observed["total"] = sum(
        (value for key, value in observed.items() if key != "total" and value is not None),
        start=np.zeros_like(observed["fwd_sync"]),
    )
    if mode != "full_components":
        raise ValueError(f"Unsupported observe mode: {mode}")
    return observed


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

    observed = {}
    for key, attr in [
        ("fwd_sync", "fwd_sync"),
        ("fwd_ssc", "fwd_ssc"),
        ("fwd_hadronic", "fwd_hadronic_gamma"),
        ("fwd_hadronic_bethe_heitler", "fwd_hadronic_bethe_heitler"),
        ("fwd_hadronic_inverse_compton", "fwd_hadronic_inverse_compton"),
        ("fwd_hadronic_pair_production", "fwd_hadronic_pair_production"),
        ("rev_sync", "rev_sync"),
        ("rev_ssc", "rev_ssc"),
        ("cross_ic", "cross_ic"),
    ]:
        source = getattr(components, attr)
        observed[key] = None
        if source is not None:
            observed[key] = _project_component(
                setup,
                components.fwd.characteristic_time_s,
                components.fwd.gamma,
                components.fwd.radius_cm,
                setup.seed_frequency_hz,
                source,
                frequencies_hz,
                config,
                timings=timings,
                label=f"Interpolation.sed_interpolation [{key}]",
            )
    observed["total"] = sum(
        (value for value in observed.values() if value is not None),
        start=np.zeros_like(observed["fwd_sync"]),
    )
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
    s = dict(
        pg_photon_survival=pg_photon_survival,
        fwd_sync=fwd_sync,
        fwd_ssc=np.zeros_like(fwd_sync),
        rev_sync=None,
        rev_ssc=None,
        cross_ic=None,
        rev_details=None,
        seed_syn_absorption=np.array(photon_field.absorption_syn_seed, dtype=float, copy=True),
        seed_ssc_total=np.array(photon_field.absorption_ssc_seed, dtype=float, copy=True),
        tau_extra=np.zeros_like(fwd_sync),
        tau_pair=np.zeros_like(fwd_sync),
        pair_lum_total=np.zeros_like(fwd_sync),
        pair_seed_total=np.zeros_like(photon_field.absorption_syn_seed, dtype=float),
        hadronic_ssa_transfer=np.ones_like(fwd_sync),
        absorbed_fwd_hadronic_gamma=None,
        absorbed_fwd_hadronic_bethe_heitler=None,
        absorbed_fwd_hadronic_inverse_compton=None,
        absorbed_fwd_hadronic_pair_production=None,
    )
    s = _stage_forward_ssc(s, setup, config, dynamics, electron, timings)
    s = _stage_reverse_emission(s, setup, config, dynamics, electron, reverse_emission, timings)
    s = _stage_hadronic_absorption(s, setup, config, dynamics, electron, hadronic)
    s = _stage_pair_production(s, setup, config, dynamics, electron)
    s = _stage_annihilation(s, setup, config, dynamics, timings)
    return _stage_assemble_result(s, setup, config, dynamics, electron, hadronic)


def _stage_forward_ssc(
    s: dict,
    setup,
    config: FitConfig,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    timings: Optional[dict[str, float]],
) -> dict:
    if config.include_forward_ssc:
        s["fwd_ssc"], s["seed_ssc_total"] = _ssc_spectrum(
            dynamics.radius,
            electron,
            setup.seed_frequency_hz,
            s["seed_syn_absorption"],
            config,
            config.num_threads,
            timings,
            "Radiation.ssc_spec [FS]",
        )
        if s["pg_photon_survival"] is not None:
            survival = np.asarray(s["pg_photon_survival"], dtype=float)
            s["fwd_ssc"] = np.asarray(s["fwd_ssc"], dtype=float) * survival
            s["seed_ssc_total"] = np.asarray(s["seed_ssc_total"], dtype=float) * survival
    return s


def _stage_reverse_emission(
    s: dict,
    setup,
    config: FitConfig,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    reverse_emission: Optional[ReverseShockEmission],
    timings: Optional[dict[str, float]],
) -> dict:
    if reverse_emission is not None:
        s["rev_sync"] = reverse_emission.l_syn_spec
        s["seed_syn_absorption"] = s["seed_syn_absorption"] + reverse_emission.seed_syn
        if reverse_emission.rs_hadronic is not None:
            s["seed_syn_absorption"] = s["seed_syn_absorption"] + _reverse_hadronic_seed_density(
                reverse_emission.rs_hadronic,
                radius_cm=dynamics.radius,
                seed_frequency_hz=setup.seed_frequency_hz,
            )
            s["rev_sync"] = s["rev_sync"] + _reverse_hadronic_luminosity(reverse_emission.rs_hadronic)
        s["rev_details"] = BranchState(
            characteristic_time_s=dynamics.r_tobs,
            gamma=dynamics.r_gamma,
            radius_cm=dynamics.radius,
            swept_mass_g=dynamics.reverse_shock.swept_mass_g,
            doppler=compute_doppler(dynamics.r_gamma, config.z),
            magnetic_field_g=reverse_emission.magnetic_field_g,
            nu_m=reverse_emission.nu_m,
            nu_c=reverse_emission.nu_c,
            nu_a=reverse_emission.nu_a,
            nu_M=reverse_emission.nu_M,
        )
        if config.reverse_shock.include_ssc:
            s["rev_ssc"], seed_ssc_rs = _timed_call(
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
            s["seed_ssc_total"] = s["seed_ssc_total"] + seed_ssc_rs
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
            s["cross_ic"] = l_cic_fs_spec + l_cic_rs_spec
            s["seed_ssc_total"] = s["seed_ssc_total"] + seed_cic_fs + seed_cic_rs
    return s


def _stage_hadronic_absorption(
    s: dict,
    setup,
    config: FitConfig,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    hadronic,
) -> dict:
    s["magnetic_field_g"] = compute_magnetic_field(dynamics.r_gamma, dynamics.radius, config)
    if hadronic is not None:
        s["hadronic_ssa_transfer"] = _forward_synchrotron_absorption_transfer(
            electron=electron,
            radius_cm=dynamics.radius,
            magnetic_field_g=s["magnetic_field_g"],
            seed_frequency_hz=setup.seed_frequency_hz,
            config=config,
        )
        s["seed_syn_absorption"] = s["seed_syn_absorption"] + _hadronic_absorbed_seed_density(
            hadronic=hadronic,
            radius_cm=dynamics.radius,
            seed_frequency_hz=setup.seed_frequency_hz,
            ssa_transfer=s["hadronic_ssa_transfer"],
        )
    return s


def _stage_pair_production(
    s: dict,
    setup,
    config: FitConfig,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
) -> dict:
    if bool(config.hadronic.include_pair_production):
        s["pair_lum_total"], s["pair_seed_total"], s["tau_pair"] = _compute_pair_production_branch(
            dynamics=dynamics,
            electron=electron,
            combined_seed_field_hz=s["seed_syn_absorption"] + s["seed_ssc_total"],
            seed_frequency_hz=setup.seed_frequency_hz,
            magnetic_field_g=s["magnetic_field_g"],
            config=config,
        )
        s["seed_syn_absorption"] = s["seed_syn_absorption"] + s["pair_seed_total"]
    s["tau_extra"] = s["tau_extra"] + s["tau_pair"]
    return s


def _stage_annihilation(
    s: dict,
    setup,
    config: FitConfig,
    dynamics: DynamicsSolution,
    timings: Optional[dict[str, float]],
) -> dict:
    annihilation_seed_syn = s["seed_syn_absorption"]
    annihilation_seed_ssc = s["seed_ssc_total"]
    if bool(config.hadronic.include_pair_production):
        annihilation_seed_syn = np.zeros_like(s["seed_syn_absorption"])
        annihilation_seed_ssc = np.zeros_like(s["seed_ssc_total"])
    absorption = _timed_call(
        timings,
        "Radiation.annihilation",
        Radiation.annihilation,
        dynamics.r_gamma,
        dynamics.radius,
        setup.seed_frequency_hz,
        annihilation_seed_syn,
        annihilation_seed_ssc,
        s["tau_extra"],
        config.num_threads,
    )
    s["prefactor"] = absorption / (4.0 * np.pi * setup.luminosity_distance_cm**2) * (1.0 + config.z)
    return s


def _stage_assemble_result(
    s: dict,
    setup,
    config: FitConfig,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    hadronic,
) -> ObserverState:
    prefactor = s["prefactor"]
    absorbed_fwd_sync = s["fwd_sync"] * prefactor
    absorbed_fwd_ssc = s["fwd_ssc"] * prefactor
    if hadronic is not None:
        hadronic_luminosity = _hadronic_absorbed_luminosity(hadronic, s["hadronic_ssa_transfer"])
        hadronic_gamma_total = hadronic_luminosity["total"]
        s["absorbed_fwd_hadronic_bethe_heitler"] = _project_optional_luminosity(
            hadronic_luminosity["bethe_heitler"],
            prefactor,
        )
        s["absorbed_fwd_hadronic_inverse_compton"] = _project_optional_luminosity(
            hadronic_luminosity["inverse_compton"],
            prefactor,
        )
        s["absorbed_fwd_hadronic_pair_production"] = _project_optional_luminosity(
            hadronic_luminosity["pair_production"],
            prefactor,
        )
        if bool(config.hadronic.include_pair_production):
            hadronic_gamma_total += s["pair_lum_total"]
            s["absorbed_fwd_hadronic_pair_production"] = s["pair_lum_total"] * prefactor
        s["absorbed_fwd_hadronic_gamma"] = hadronic_gamma_total * prefactor
    absorbed_rev_sync = None if s["rev_sync"] is None else s["rev_sync"] * prefactor
    absorbed_rev_ssc = None if s["rev_ssc"] is None else s["rev_ssc"] * prefactor
    absorbed_cross_ic = None if s["cross_ic"] is None else s["cross_ic"] * prefactor

    total = absorbed_fwd_sync + absorbed_fwd_ssc
    if absorbed_rev_sync is not None:
        total = total + absorbed_rev_sync
    if absorbed_rev_ssc is not None:
        total = total + absorbed_rev_ssc
    if absorbed_cross_ic is not None:
        total = total + absorbed_cross_ic
    if s["absorbed_fwd_hadronic_gamma"] is not None:
        total = total + s["absorbed_fwd_hadronic_gamma"]

    return ObserverState(
        prefactor=np.asarray(prefactor, dtype=float),
        tau_extra=np.asarray(s["tau_extra"], dtype=float),
        tau_pair=np.asarray(s["tau_pair"], dtype=float),
        components=FluxComponents(
            total=total,
            fwd_sync=absorbed_fwd_sync,
            fwd_ssc=absorbed_fwd_ssc,
            fwd_hadronic_gamma=s["absorbed_fwd_hadronic_gamma"],
            fwd_hadronic_bethe_heitler=s["absorbed_fwd_hadronic_bethe_heitler"],
            fwd_hadronic_inverse_compton=s["absorbed_fwd_hadronic_inverse_compton"],
            fwd_hadronic_pair_production=s["absorbed_fwd_hadronic_pair_production"],
            rev_sync=absorbed_rev_sync,
            rev_ssc=absorbed_rev_ssc,
            cross_ic=absorbed_cross_ic,
            fwd=BranchState(
                characteristic_time_s=dynamics.r_tobs,
                gamma=dynamics.r_gamma,
                radius_cm=dynamics.radius,
                swept_mass_g=dynamics.swept_mass_g,
                doppler=compute_doppler(dynamics.r_gamma, config.z),
                magnetic_field_g=compute_magnetic_field(dynamics.r_gamma, dynamics.radius, config),
                nu_m=electron.nu_m,
                nu_c=electron.nu_c,
                nu_a=electron.nu_a,
                nu_M=compute_maximum_synchrotron_frequency(dynamics.r_gamma, dynamics.radius, config),
                cooling_timescale_s=electron.cooling_timescale_s,
                dynamical_timescale_s=electron.dynamical_timescale_s,
            ),
            rev=s["rev_details"],
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
    radius = np.asarray(dynamics.radius, dtype=float)
    gamma_bulk = np.asarray(dynamics.r_gamma, dtype=float)
    magnetic_field = np.asarray(magnetic_field_g, dtype=float)
    _validate_pair_production_branch_inputs(v_seed, seed_field, gam_e, radius, gamma_bulk, magnetic_field)
    e_e_gev = gam_e * _ELECTRON_MASS_GEV
    gam_edge = _hadronic_build_gamma_edges(gam_e)
    num_nu = int(v_seed.size)
    num_r = int(radius.size)
    pair_lum = np.zeros((num_nu, num_r), dtype=float)
    pair_seed = np.zeros((num_nu, num_r), dtype=float)
    tau_pair = np.zeros((num_nu, num_r), dtype=float)
    d_n_pair_prev = np.zeros(gam_e.size, dtype=float)
    photon_energy_gev, _ = photon_density_hz_to_gev(v_seed, np.ones_like(v_seed))
    from asgard_core.hadronic_cascade import shell_path_time_seconds

    use_iterative_cascade = (
        int(getattr(config.hadronic, 'pair_cascade_iterations', 1)) > 1
    )
    if use_iterative_cascade:
        from asgard_core.hadronic_cascade import compute_time_dependent_pair_cascade_sequence
        cascade = compute_time_dependent_pair_cascade_sequence(
            photon_energy_gev=photon_energy_gev,
            primary_photon_density_per_gev=seed_field / constants.para_h_gev,
            electron_energy_gev=e_e_gev,
            frequency_hz=v_seed,
            radius_cm=radius,
            gamma_bulk=gamma_bulk,
            observer_time_s=np.asarray(dynamics.r_tobs, dtype=float),
            b_field_g=magnetic_field,
            num_threads=int(config.num_threads),
            index_syn_integr=int(config.index_syn_integr),
            substeps_per_shell=int(config.hadronic.pair_cascade_iterations),
        )
        return (
            np.asarray(cascade.pair_syn_luminosity_hz, dtype=float),
            np.asarray(cascade.pair_syn_seed_per_hz, dtype=float),
            np.asarray(cascade.tau_pair_path, dtype=float),
        )

    for i_r in range(num_r):
        _, photon_density_per_gev = photon_density_hz_to_gev(v_seed, seed_field[:, i_r])

        ppair = solve_pair_production(
            photon_energy_gev=photon_energy_gev,
            photon_density_per_gev=photon_density_per_gev,
            electron_energy_gev=e_e_gev,
        )
        photon_loss_rate = np.asarray(ppair.photon_loss_rate, dtype=float)
        if np.any(photon_loss_rate < 0.0):
            raise RuntimeError("pair production Fortran kernel returned negative photon loss rate.")
        tau_pair[:, i_r] = photon_loss_rate * shell_path_time_seconds(float(radius[i_r]), float(gamma_bulk[i_r]))
        q_pair = np.asarray(ppair.pair_injection_rate_per_gev_total, dtype=float) * (
            (4.0 / 3.0) * np.pi * (
                float(radius[i_r]) ** 3
                - (0.0 if i_r == 0 else float(radius[i_r - 1]) ** 3)
            )
        ) * _ELECTRON_MASS_GEV
        d_n_pair = _hadronic_advance_energy_loggamma(
            gam_e, gam_edge, d_n_pair_prev, q_pair,
            _hadronic_electron_loss_rates(
                gam_e, float(magnetic_field[i_r]),
                float(radius[i_r]) / (float(gamma_bulk[i_r]) * constants.para_c),
            ),
            _hadronic_shell_dt(np.asarray(dynamics.r_tobs, dtype=float), i_r),
        )
        if magnetic_field[i_r] > 0.0:
            p_syn_i, seed_syn_i = electron_radiation_module.get_syn_selected(
                int(config.index_syn_integr),
                float(radius[i_r]),
                float(magnetic_field[i_r]),
                int(config.num_threads),
                gam_e, d_n_pair, v_seed,
            )
            pair_lum[:, i_r] = np.asarray(p_syn_i, dtype=float)
            pair_seed[:, i_r] = np.asarray(seed_syn_i, dtype=float)
        d_n_pair_prev = d_n_pair
    return pair_lum, pair_seed, tau_pair


def _validate_pair_production_branch_inputs(
    frequency_hz: np.ndarray,
    seed_field_hz: np.ndarray,
    gam_e: np.ndarray,
    radius_cm: np.ndarray,
    gamma_bulk: np.ndarray,
    magnetic_field_g: np.ndarray,
) -> None:
    if frequency_hz.ndim != 1 or gam_e.ndim != 1:
        raise ValueError("pair-production branch requires 1D frequency and electron grids.")
    if radius_cm.ndim != 1 or gamma_bulk.shape != radius_cm.shape or magnetic_field_g.shape != radius_cm.shape:
        raise ValueError("pair-production branch requires matching radius, gamma, and magnetic-field arrays.")
    if seed_field_hz.shape != (frequency_hz.size, radius_cm.size):
        raise ValueError("pair-production seed field shape must be (num_frequency, num_radius).")
    arrays = (frequency_hz, seed_field_hz, gam_e, radius_cm, gamma_bulk, magnetic_field_g)
    if not all(np.all(np.isfinite(arr)) for arr in arrays):
        raise ValueError("pair-production branch inputs must be finite.")
    if np.any(frequency_hz <= 0.0) or np.any(np.diff(frequency_hz) <= 0.0):
        raise ValueError("pair-production frequency grid must be positive and strictly increasing.")
    if np.any(gam_e < 1.0) or np.any(np.diff(gam_e) <= 0.0):
        raise ValueError("pair-production electron gamma grid must start at gamma >= 1 and be strictly increasing.")
    if np.any(seed_field_hz < 0.0):
        raise ValueError("pair-production seed field must be non-negative.")
    if np.any(radius_cm <= 0.0) or np.any(np.diff(radius_cm) <= 0.0):
        raise ValueError("pair-production shell radii must be positive and strictly increasing.")
    if np.any(gamma_bulk < 1.0):
        raise ValueError("pair-production bulk Lorentz factors must be >= 1.")
    if np.any(magnetic_field_g < 0.0):
        raise ValueError("pair-production magnetic fields must be non-negative.")


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
        ("bethe_heitler", getattr(hadronic, "l_had_bethe_heitler", None)),
        ("inverse_compton", getattr(hadronic, "l_had_hadronic_inverse_compton", None)),
        ("pair_production", getattr(hadronic, "l_had_pair_production", None)),
        ("pion_synch", getattr(hadronic, "l_had_pion_synch", None)),
        ("muon_synch", getattr(hadronic, "l_had_muon_synch", None)),
        ("pion_inverse_compton", getattr(hadronic, "l_had_pion_inverse_compton", None)),
        ("muon_inverse_compton", getattr(hadronic, "l_had_muon_inverse_compton", None)),
    )


def _reverse_hadronic_luminosity(rs_hadronic) -> np.ndarray:
    total = np.asarray(rs_hadronic.l_had_syn_spec, dtype=float)
    l_pg = getattr(rs_hadronic, "l_had_pg_gamma", None)
    if l_pg is not None:
        total = total + np.asarray(l_pg, dtype=float)
    for _name, luminosity in _hadronic_optional_luminosities(rs_hadronic):
        if luminosity is not None:
            total = total + np.asarray(luminosity, dtype=float)
    return total


def _reverse_hadronic_seed_density(
    rs_hadronic,
    *,
    radius_cm: np.ndarray,
    seed_frequency_hz: np.ndarray,
) -> np.ndarray:
    seed_total = np.asarray(rs_hadronic.seed_had_syn, dtype=float)
    seed_bh = getattr(rs_hadronic, "seed_had_bethe_heitler", None)
    if seed_bh is not None:
        seed_total = seed_total + np.asarray(seed_bh, dtype=float)
    l_pg = getattr(rs_hadronic, "l_had_pg_gamma", None)
    luminosities = (l_pg, *[item[1] for item in _hadronic_optional_luminosities(rs_hadronic)])
    for luminosity in luminosities:
        if luminosity is not None:
            seed_total = seed_total + _seed_density_from_luminosity(
                luminosity,
                radius_cm=radius_cm,
                seed_frequency_hz=seed_frequency_hz,
            )
    return seed_total


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
        if radius_cm[i_shell] <= 0.0:
            raise ValueError("BH electron merge requires positive shell radii.")
        if magnetic_field_g[i_shell] < 0.0:
            raise ValueError("BH electron merge requires non-negative magnetic fields.")
        if magnetic_field_g[i_shell] == 0.0:
            continue
        p_syn_i, seed_syn_i = electron_radiation_module.get_syn_selected(
            int(config.index_syn_integr),
            float(radius_cm[i_shell]),
            float(magnetic_field_g[i_shell]),
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
        l_syn_spec_chi=None if electron.l_syn_spec_chi is None else np.asarray(electron.l_syn_spec_chi, dtype=float),
        seed_syn_chi=None if electron.seed_syn_chi is None else np.asarray(electron.seed_syn_chi, dtype=float),
        tau_syn_chi=None if electron.tau_syn_chi is None else np.asarray(electron.tau_syn_chi, dtype=float),
        chi_radius_cm=None if electron.chi_radius_cm is None else np.asarray(electron.chi_radius_cm, dtype=float),
        chi_gamma_bulk=None if electron.chi_gamma_bulk is None else np.asarray(electron.chi_gamma_bulk, dtype=float),
        chi_dvolume_weight=None if electron.chi_dvolume_weight is None else np.asarray(electron.chi_dvolume_weight, dtype=float),
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
    _require_top_hat_phi_grid(config)
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


def _timed_call(timings: Optional[dict[str, float]], label: Optional[str], func, *args, **kwargs):
    start = perf_counter()
    result = func(*args, **kwargs)
    elapsed = perf_counter() - start
    if timings is not None and label is not None:
        timings[label] = timings.get(label, 0.0) + elapsed
    return result
