from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from time import perf_counter
from typing import Optional

import numpy as np

from asgard_coupling import build_coupled_shock_geometry, build_cross_zone_seed_fields
from asgard_models import ExecutionPolicy, FitConfig, PhysicalSolution, SimulationSetup
from asgard_postprocess import interpolate_observed_flux
from asgard_runtime import (
    DynamicsSolution,
    ElectronSolution,
    ReverseShockEmission,
    solve_dynamics,
    solve_electron,
    solve_reverse_shock_emission,
)
from asgard_ssc import compute_ssc_auxiliary_grid
from asgard_setup import build_seed_frequency_grid, build_simulation_setup
from src import Radiation, constants


@dataclass
class BranchDetails:
    characteristic_time_s: np.ndarray
    gamma: np.ndarray
    radius_cm: np.ndarray
    swept_mass_g: np.ndarray
    doppler: np.ndarray
    magnetic_field_g: np.ndarray
    nu_m: np.ndarray
    nu_c: np.ndarray
    nu_a: np.ndarray
    nu_M: np.ndarray


@dataclass
class ComponentSpectra:
    total: np.ndarray
    fwd_sync: np.ndarray
    fwd_ssc: np.ndarray
    rev_sync: Optional[np.ndarray]
    rev_ssc: Optional[np.ndarray]
    cross_ic: Optional[np.ndarray]
    fwd: BranchDetails
    rev: Optional[BranchDetails]


@dataclass
class ModelState:
    config: FitConfig
    setup: SimulationSetup
    policy: ExecutionPolicy
    dynamics: DynamicsSolution
    electron: ElectronSolution
    reverse_emission: Optional[ReverseShockEmission]
    component_spectra: ComponentSpectra
    requested_frequency_min_hz: Optional[float]
    requested_frequency_max_hz: Optional[float]
    timings: dict[str, float] = field(default_factory=dict)


@dataclass
class ObservedState:
    state: ModelState
    setup: SimulationSetup
    frequencies_hz: np.ndarray
    components: dict[str, np.ndarray | None]


def build_execution_policy(config: FitConfig) -> ExecutionPolicy:
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


def build_solve_time_grid(observer_time_s: np.ndarray, default_num_tobs: int) -> np.ndarray:
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


def build_query_config(
    config: FitConfig,
    observer_time_s: np.ndarray,
) -> FitConfig:
    query = deepcopy(config)
    observer_time_s = np.asarray(observer_time_s, dtype=float)
    query.num_tobs = observer_time_s.shape[0]
    query.t_obs_min_log10 = float(np.log10(observer_time_s.min()))
    query.t_obs_max_log10 = float(np.log10(observer_time_s.max()))
    return query


def build_query_setup(
    config: FitConfig,
    observer_time_s: np.ndarray,
    requested_frequencies_hz: np.ndarray | None = None,
):
    setup = build_simulation_setup(build_query_config(config, observer_time_s))
    setup.observer_time_s = np.asarray(observer_time_s, dtype=float)
    if requested_frequencies_hz is not None:
        setup.seed_frequency_hz = build_seed_frequency_grid(config, requested_frequencies_hz)
    return setup


def solve_model_state(
    config: FitConfig,
    observer_time_s: np.ndarray,
    requested_frequencies_hz: np.ndarray | None = None,
    timings: Optional[dict[str, float]] = None,
    policy: Optional[ExecutionPolicy] = None,
) -> ModelState:
    setup = build_query_setup(config, observer_time_s, requested_frequencies_hz)
    return solve_model_state_from_setup(
        config,
        setup,
        timings=timings,
        policy=policy,
        requested_frequencies_hz=requested_frequencies_hz,
    )


def solve_component_spectra(
    config: FitConfig,
    observer_time_s: np.ndarray,
    requested_frequencies_hz: np.ndarray | None = None,
) -> ComponentSpectra:
    return solve_model_state(config, observer_time_s, requested_frequencies_hz).component_spectra


def solve_model_state_from_setup(
    config: FitConfig,
    setup,
    timings: Optional[dict[str, float]] = None,
    policy: Optional[ExecutionPolicy] = None,
    requested_frequencies_hz: np.ndarray | None = None,
) -> ModelState:
    execution_policy = build_execution_policy(config) if policy is None else policy
    if config.reverse or config.reverse_shock.enabled:
        dynamics_label = "Dynamics.dynamics_reverse"
    else:
        dynamics_label = "Dynamics.dynamics_forward"
    dynamics = _timed_call(timings, dynamics_label, solve_dynamics, setup.boundary, config)
    solver_name = config.electron_solver.lower()
    electron_label_map = {
        "fullhide": "Electron.fs_electron_fullhide",
        "t2g1": "Electron.fs_electron_t2g1",
        "slc1": "Electron.fs_electron_slc1",
        "weno5": "Electron.fs_electron_weno5",
    }
    electron_label = electron_label_map.get(solver_name, f"Electron.{solver_name}")
    electron = _timed_call(timings, electron_label, solve_electron, setup.boundary, dynamics, setup.seed_frequency_hz, config)
    reverse_emission = None
    if config.reverse or config.reverse_shock.enabled:
        reverse_emission = _timed_call(
            timings,
            "Radiation.seed_reverse",
            solve_reverse_shock_emission,
            setup.boundary,
            dynamics,
            setup.seed_frequency_hz,
            config,
        )
    component_spectra = _assemble_component_spectra(setup, config, dynamics, electron, reverse_emission, timings=timings)
    freq_min, freq_max = _requested_frequency_bounds(requested_frequencies_hz)
    return ModelState(
        config=config,
        setup=setup,
        policy=execution_policy,
        dynamics=dynamics,
        electron=electron,
        reverse_emission=reverse_emission,
        component_spectra=component_spectra,
        requested_frequency_min_hz=freq_min,
        requested_frequency_max_hz=freq_max,
        timings={} if timings is None else dict(timings),
    )


def solve_component_spectra_from_setup(
    config: FitConfig,
    setup,
    timings: Optional[dict[str, float]] = None,
) -> ComponentSpectra:
    return solve_model_state_from_setup(config, setup, timings=timings).component_spectra


def _build_observer_setup_from_state(
    state: ModelState,
    observer_time_s: np.ndarray,
) -> SimulationSetup:
    return SimulationSetup(
        luminosity_distance_cm=state.setup.luminosity_distance_cm,
        boundary=state.setup.boundary,
        seed_frequency_hz=state.setup.seed_frequency_hz,
        observer_time_s=np.asarray(observer_time_s, dtype=float),
    )


def state_covers_request(
    state: ModelState,
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


def observe_flux_grid_from_state(
    state: ModelState,
    observer_time_s: np.ndarray,
    frequencies_hz: np.ndarray,
    timings: Optional[dict[str, float]] = None,
    mode: str = "full_components",
) -> ObservedState:
    setup = _build_observer_setup_from_state(state, observer_time_s)
    observed = observe_spectra_from_setup(
        state.config,
        state.component_spectra,
        setup,
        frequencies_hz,
        timings=timings,
        mode=mode,
    )
    return ObservedState(
        state=state,
        setup=setup,
        frequencies_hz=np.asarray(frequencies_hz, dtype=float),
        components=observed,
    )


def observe_spectrum_from_state(
    state: ModelState,
    time_s: float,
    frequencies_hz: np.ndarray,
    timings: Optional[dict[str, float]] = None,
    mode: str = "full_components",
) -> ObservedState:
    return observe_flux_grid_from_state(
        state,
        np.array([time_s], dtype=float),
        frequencies_hz,
        timings=timings,
        mode=mode,
    )


def extract_details_from_state(state: ModelState) -> dict[str, Optional[BranchDetails]]:
    return {"fwd": state.component_spectra.fwd, "rev": state.component_spectra.rev}


def extract_physical_solution_from_state(state: ModelState) -> PhysicalSolution:
    rev = state.component_spectra.rev
    return PhysicalSolution(
        characteristic_time_s=state.component_spectra.fwd.characteristic_time_s,
        gamma=state.component_spectra.fwd.gamma,
        radius_cm=state.component_spectra.fwd.radius_cm,
        absorbed_spectral_flux=state.component_spectra.total,
        nu_m=state.component_spectra.fwd.nu_m,
        nu_c=state.component_spectra.fwd.nu_c,
        nu_a=state.component_spectra.fwd.nu_a,
        rs_nu_m=None if rev is None else rev.nu_m,
        rs_nu_c=None if rev is None else rev.nu_c,
        rs_nu_a=None if rev is None else rev.nu_a,
    )


def observe_spectra(
    config: FitConfig,
    component_spectra: ComponentSpectra,
    observer_time_s: np.ndarray,
    frequencies_hz: np.ndarray,
) -> dict[str, np.ndarray]:
    setup = build_query_setup(config, observer_time_s, frequencies_hz)
    return observe_spectra_from_setup(config, component_spectra, setup, frequencies_hz)


def observe_spectra_from_setup(
    config: FitConfig,
    component_spectra: ComponentSpectra,
    setup,
    frequencies_hz: np.ndarray,
    timings: Optional[dict[str, float]] = None,
    mode: str = "full_components",
) -> dict[str, np.ndarray | None]:
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    if mode not in {"full_components", "total_only"}:
        raise ValueError(f"Unsupported observe mode: {mode}")
    if mode == "total_only":
        total = _observe_component(
            setup,
            component_spectra.fwd.characteristic_time_s,
            component_spectra.fwd.gamma,
            component_spectra.fwd.radius_cm,
            setup.seed_frequency_hz,
            component_spectra.total,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [total]",
        )
        return {
            "total": total,
            "fwd_sync": None,
            "fwd_ssc": None,
            "rev_sync": None,
            "rev_ssc": None,
            "cross_ic": None,
        }

    observed = {
        "fwd_sync": _observe_component(
            setup,
            component_spectra.fwd.characteristic_time_s,
            component_spectra.fwd.gamma,
            component_spectra.fwd.radius_cm,
            setup.seed_frequency_hz,
            component_spectra.fwd_sync,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [fwd_sync]",
        ),
        "fwd_ssc": _observe_component(
            setup,
            component_spectra.fwd.characteristic_time_s,
            component_spectra.fwd.gamma,
            component_spectra.fwd.radius_cm,
            setup.seed_frequency_hz,
            component_spectra.fwd_ssc,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [fwd_ssc]",
        ),
    }
    observed["rev_sync"] = None
    observed["rev_ssc"] = None
    observed["cross_ic"] = None
    if component_spectra.rev_sync is not None:
        observed["rev_sync"] = _observe_component(
            setup,
            component_spectra.fwd.characteristic_time_s,
            component_spectra.fwd.gamma,
            component_spectra.fwd.radius_cm,
            setup.seed_frequency_hz,
            component_spectra.rev_sync,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [rev_sync]",
        )
    if component_spectra.rev_ssc is not None:
        observed["rev_ssc"] = _observe_component(
            setup,
            component_spectra.fwd.characteristic_time_s,
            component_spectra.fwd.gamma,
            component_spectra.fwd.radius_cm,
            setup.seed_frequency_hz,
            component_spectra.rev_ssc,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [rev_ssc]",
        )
    if component_spectra.cross_ic is not None:
        observed["cross_ic"] = _observe_component(
            setup,
            component_spectra.fwd.characteristic_time_s,
            component_spectra.fwd.gamma,
            component_spectra.fwd.radius_cm,
            setup.seed_frequency_hz,
            component_spectra.cross_ic,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [cross_ic]",
        )
    observed["total"] = np.array(observed["fwd_sync"], dtype=float, copy=True)
    observed["total"] += np.array(observed["fwd_ssc"], dtype=float, copy=False)
    if observed["rev_sync"] is not None:
        observed["total"] += np.array(observed["rev_sync"], dtype=float, copy=False)
    if observed["rev_ssc"] is not None:
        observed["total"] += np.array(observed["rev_ssc"], dtype=float, copy=False)
    if observed["cross_ic"] is not None:
        observed["total"] += np.array(observed["cross_ic"], dtype=float, copy=False)
    if timings is not None:
        timings.setdefault("Interpolation.sed_interpolation [total]", 0.0)
    return observed


def _assemble_component_spectra(
    setup,
    config: FitConfig,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    reverse_emission: Optional[ReverseShockEmission],
    timings: Optional[dict[str, float]] = None,
) -> ComponentSpectra:
    fwd_sync = electron.l_syn_spec
    seed_syn_absorption = electron.seed_syn.copy()
    zero_component = np.zeros_like(fwd_sync)
    fwd_ssc = zero_component
    seed_ssc_total = zero_component
    if config.include_forward_ssc:
        fwd_ssc, seed_ssc_total = _compute_ssc_spectrum(
            dynamics.radius,
            electron,
            setup.seed_frequency_hz,
            electron.seed_syn,
            config.num_threads,
            timings,
            "Radiation.ssc_spec [FS]",
        )

    rev_sync = None
    rev_ssc = None
    cross_ic = None
    rev_details = None
    if reverse_emission is not None:
        rev_sync = reverse_emission.l_syn_spec
        seed_syn_absorption = seed_syn_absorption + reverse_emission.seed_syn
        rev_details = BranchDetails(
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
            l_cic_fs_spec, seed_cic_fs = _compute_ssc_spectrum(
                dynamics.radius,
                electron,
                setup.seed_frequency_hz,
                seed_rs_to_fs,
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

    absorption = _timed_call(
        timings,
        "Radiation.annihilation",
        Radiation.annihilation,
        dynamics.r_gamma,
        dynamics.radius,
        setup.seed_frequency_hz,
        seed_syn_absorption,
        seed_ssc_total,
        config.num_threads,
    )
    prefactor = absorption / (4.0 * np.pi * setup.luminosity_distance_cm**2) * (1.0 + config.z)

    absorbed_fwd_sync = fwd_sync * prefactor
    absorbed_fwd_ssc = fwd_ssc * prefactor
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

    return ComponentSpectra(
        total=total,
        fwd_sync=absorbed_fwd_sync,
        fwd_ssc=absorbed_fwd_ssc,
        rev_sync=absorbed_rev_sync,
        rev_ssc=absorbed_rev_ssc,
        cross_ic=absorbed_cross_ic,
        fwd=BranchDetails(
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
        ),
        rev=rev_details,
    )


def _compute_ssc_spectrum(
    radius_cm: np.ndarray,
    electron: ElectronSolution,
    seed_frequency_hz: np.ndarray,
    seed_field: np.ndarray,
    num_threads: int,
    timings: Optional[dict[str, float]],
    label: str,
) -> tuple[np.ndarray, np.ndarray]:
    if electron.work_x_edge_log10 is not None and electron.work_d_n_x is not None:
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


def _observe_component(
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


def profile_observe_spectra_from_setup(
    config: FitConfig,
    setup,
    frequencies_hz: np.ndarray,
) -> tuple[ComponentSpectra, dict[str, np.ndarray], dict[str, float]]:
    timings: dict[str, float] = {}
    component_spectra = solve_component_spectra_from_setup(config, setup, timings=timings)
    observed = observe_spectra_from_setup(config, component_spectra, setup, frequencies_hz, timings=timings)
    return component_spectra, observed, timings


def _timed_call(timings: Optional[dict[str, float]], label: Optional[str], func, *args):
    start = perf_counter()
    result = func(*args)
    elapsed = perf_counter() - start
    if timings is not None and label is not None:
        timings[label] = timings.get(label, 0.0) + elapsed
    return result


def _compute_doppler(gamma: np.ndarray, redshift: float) -> np.ndarray:
    beta = np.sqrt(1.0 - gamma**-2)
    return 1.0 / (gamma * (1.0 - beta) * (1.0 + redshift))


def _ambient_density(radius_cm: np.ndarray, config: FitConfig) -> np.ndarray:
    radius = np.asarray(radius_cm, dtype=float)
    if config.a_star > 0.0:
        d_ne_wind = config.a_star * 3.0e35 / radius**2
        density = np.where(d_ne_wind <= config.d_ne / 4.0, config.d_ne, d_ne_wind)
    else:
        density = config.d_ne * (
            1.0
            + (config.f_jump - 1.0)
            * np.exp(-(np.log10(radius) - np.log10(config.r_tr)) ** 2 / (2.0 * config.f_wide**2))
        )
    if np.any(radius < config.r0) and config.a_star > 0.0:
        density = density.copy()
        density[radius < config.r0] = config.a_star * 3.0e35 / config.r0**2
    return density


def _compute_forward_magnetic_field(gamma: np.ndarray, radius_cm: np.ndarray, config: FitConfig) -> np.ndarray:
    density = _ambient_density(radius_cm, config)
    return 0.39 * np.sqrt(config.epsilon_b * density * gamma * np.maximum(gamma - 1.0, 0.0))


def _compute_maximum_synchrotron_frequency(gamma: np.ndarray, radius_cm: np.ndarray, config: FitConfig) -> np.ndarray:
    magnetic_field = _compute_forward_magnetic_field(gamma, radius_cm, config)
    doppler = _compute_doppler(gamma, config.z)
    gam_e_max = 3.0 * constants.para_m_energy / np.sqrt(8.0 * magnetic_field * constants.para_e**3)
    return 4.2e6 * magnetic_field * gam_e_max**2 * doppler
