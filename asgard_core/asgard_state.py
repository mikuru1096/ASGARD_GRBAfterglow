from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from time import perf_counter

import numpy as np

from src.Electron.electron_radiation import electron_radiation_kernel as electron_radiation_module
import src.Hadronic.hadronic_forward_1d as hadronic_legacy_module
from asgard_core.asgard_config import ExecutionPolicy, RuntimeConfig, SimulationSetup
from asgard_core.hadronic_processes import solve_pair_production
from asgard_core.hadronic_processes import photon_density_hz_to_gev
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
from asgard_core.asgard_physics_utils import (
    ambient_density,
    compute_doppler,
    compute_magnetic_field,
    reverse_shell_baryonic_mass,
)
from asgard_core.asgard_physics_utils import density_jump_arrays
from asgard_core.asgard_postprocess import interpolate_observed_flux
from asgard_core.asgard_runtime import (
    _hadronic_pg_survival_factor,
    _hadronic_shell_comoving_dt_from_radius,
    _solver_report,
    _use_direct_chi_projection_contract,
    solve_dynamics,
    solve_electron,
    solve_electron_with_cooling_seed,
    solve_hadronic,
    solve_reverse_shock_emission,
)
from asgard_core.asgard_setup import build_simulation_setup
from src import Interpolation, Radiation, constants


_ELECTRON_MASS_GEV = constants.para_m_e_gev
_PROJECTION_KINDS = {"lightcurve", "sed"}
_COUPLING_SEPARATED = "separated"
_COUPLING_JOINT = "joint"
_JOINT_ELECTRON_PHOTON_ITERATIONS = 2


@dataclass
class _CoupledShockGeometry:
    proper_time_s: np.ndarray
    fs_width_cm: np.ndarray
    rs_width_cm: np.ndarray
    center_delay_s: np.ndarray


def build_coupled_shock_geometry(dynamics, config: RuntimeConfig) -> _CoupledShockGeometry:
    if dynamics.reverse_shock is None:
        raise ValueError("Reverse-shock dynamics are required to build coupled-shock geometry.")
    if config.reverse_shock.delta_t_s is None:
        raise ValueError("ReverseShockConfig.delta_t_s must be set when reverse shock is enabled.")

    radius_cm = dynamics.radius
    gamma = dynamics.r_gamma
    proper_time_s = _integrate_proper_time(radius_cm, gamma)

    fs_width_cm, rs_width_cm = np.zeros((2, radius_cm.size), dtype=float)

    eta_0 = config.eta_0
    shell_mass_g = reverse_shell_baryonic_mass(config)
    delta_0_cm = config.reverse_shock.delta_t_s * constants.para_c

    for i, radius_loc in enumerate(radius_cm):
        gamma_loc = gamma[i]
        n1 = ambient_density(radius_loc, config)
        n2 = 4.0 * gamma_loc * n1
        fs_width_cm[i] = dynamics.swept_mass_g[i] / (4.0 * np.pi * radius_loc**2 * n2 * constants.para_m_p)

        delta_shell_cm = max(delta_0_cm, radius_loc / eta_0**2)
        n4 = shell_mass_g / (4.0 * np.pi * constants.para_m_p * radius_loc**2 * eta_0 * delta_shell_cm)
        u2 = np.sqrt(gamma_loc * gamma_loc - 1.0)
        u4 = np.sqrt(eta_0 * eta_0 - 1.0)
        gamma34 = (gamma_loc * gamma_loc + eta_0 * eta_0 - 1.0) / (eta_0 * gamma_loc + u2 * u4)
        n3 = (4.0 * gamma34 + 3.0) * n4
        rs_width_cm[i] = (
            dynamics.reverse_shock.swept_mass_g[i] / (4.0 * np.pi * radius_loc**2 * n3 * constants.para_m_p)
        )

    center_delay_s = 0.5 * (fs_width_cm + rs_width_cm) / constants.para_c
    return _CoupledShockGeometry(
        proper_time_s=proper_time_s,
        fs_width_cm=fs_width_cm,
        rs_width_cm=rs_width_cm,
        center_delay_s=center_delay_s,
    )


def build_cross_zone_seed_fields(
    fs_seed_syn: np.ndarray,
    rs_seed_syn: np.ndarray,
    geometry: _CoupledShockGeometry,
    angular_factor: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    tau = geometry.proper_time_s
    tau_ret = tau - geometry.center_delay_s
    seed_fs_to_rs = _retarded_seed_interpolation(fs_seed_syn, tau, tau_ret, angular_factor)
    seed_rs_to_fs = _retarded_seed_interpolation(rs_seed_syn, tau, tau_ret, angular_factor)
    return seed_fs_to_rs, seed_rs_to_fs


def _integrate_proper_time(radius_cm: np.ndarray, gamma: np.ndarray) -> np.ndarray:
    proper_time_s = np.zeros_like(radius_cm)
    for i in range(1, radius_cm.shape[0]):
        gamma_mean = 0.5 * (gamma[i - 1] + gamma[i])
        beta_mean = np.sqrt(1.0 - gamma_mean**-2)
        d_radius = radius_cm[i] - radius_cm[i - 1]
        proper_time_s[i] = proper_time_s[i - 1] + d_radius / (beta_mean * gamma_mean * constants.para_c)
    return proper_time_s


def _retarded_seed_interpolation(
    seed_syn: np.ndarray,
    source_time_s: np.ndarray,
    retarded_time_s: np.ndarray,
    angular_factor: float,
) -> np.ndarray:
    shifted_seed = np.zeros_like(seed_syn)
    for i_nu in range(seed_syn.shape[0]):
        shifted_seed[i_nu] = angular_factor * np.interp(
            retarded_time_s,
            source_time_s,
            seed_syn[i_nu],
            left=0.0,
            right=seed_syn[i_nu, -1],
        )
    return shifted_seed


def make_tgrid(observer_time_s: np.ndarray, default_num_tobs: int) -> np.ndarray:
    observer_time_s = np.asarray(observer_time_s, dtype=float)
    if observer_time_s.size == 0:
        raise ValueError("observer_time_s must be non-empty.")
    if np.any(observer_time_s <= 0.0):
        raise ValueError("observer_time_s must be positive.")
    if observer_time_s.size == 1:
        return observer_time_s.copy()
    num_tobs = max(int(default_num_tobs), int(np.unique(observer_time_s).size))
    t_min = float(np.min(observer_time_s))
    t_max = float(np.max(observer_time_s))
    return np.logspace(np.log10(t_min), np.log10(t_max), num_tobs)


def make_query_cfg(
    config: RuntimeConfig,
    observer_time_s: np.ndarray,
) -> RuntimeConfig:
    query = deepcopy(config)
    observer_time_s = np.asarray(observer_time_s, dtype=float)
    query.num_tobs = observer_time_s.shape[0]
    query.t_obs_min_log10 = float(np.log10(observer_time_s.min()))
    query.t_obs_max_log10 = float(np.log10(observer_time_s.max()))
    return query


def make_query_setup(
    config: RuntimeConfig,
    observer_time_s: np.ndarray,
    requested_frequencies_hz: np.ndarray | None = None,
):
    setup = build_simulation_setup(make_query_cfg(config, observer_time_s), requested_frequencies_hz)
    setup.observer_time_s = np.asarray(observer_time_s, dtype=float)
    return setup


def solve_state(
    config: RuntimeConfig,
    observer_time_s: np.ndarray,
    requested_frequencies_hz: np.ndarray | None = None,
    timings: dict[str, float] | None = None,
    policy: ExecutionPolicy | None = None,
) -> SolveState:
    setup = make_query_setup(config, observer_time_s, requested_frequencies_hz)
    return solve_state_from_setup(
        config,
        setup,
        timings=timings,
        policy=policy,
        requested_frequencies_hz=requested_frequencies_hz,
    )


def _solver_label(config: RuntimeConfig, stage: str) -> str:
    if stage == "dynamics":
        if config.reverse or config.reverse_shock.enabled:
            return "Dynamics.dynamics_reverse"
        return "Dynamics.dynamics_forward"
    if stage == "electron":
        solver_name = config.electron_solver.lower()
        electron_label_map = {
            "fullhide_1d": "Electron.fs_electron_fullhide_1d",
            "fullhide_2d": "Electron.fs_electron_transport_2d_core",
            "t2g1_1d": "Electron.fs_electron_t2g1_1d",
            "slc1_1d": "Electron.fs_electron_slc1_1d",
            "charint_1d": "Electron.fs_electron_charint_1d",
            "charint_2d": "Electron.fs_electron_transport_2d_core",
            "weno5_1d": "Electron.fs_electron_weno5_1d",
        }
        return electron_label_map.get(solver_name, f"Electron.{solver_name}")
    if stage == "hadronic":
        return "Hadronic.fs_hadronic_1d"
    raise ValueError(f"Unsupported solver stage: {stage}")


def _electron_photon_coupling(config: RuntimeConfig) -> str:
    coupling = str(config.electron_photon_coupling).lower()
    if coupling not in {_COUPLING_SEPARATED, _COUPLING_JOINT}:
        raise ValueError("electron_photon_coupling must be 'separated' or 'joint'.")
    return coupling


def _validate_joint_electron_photon_config(config: RuntimeConfig) -> None:
    electron_solver = str(config.electron_solver).lower()
    if electron_solver != "fullhide_1d":
        raise NotImplementedError("electron_photon_coupling='joint' currently supports only electron_solver='fullhide_1d'.")
    if config.reverse or config.reverse_shock.enabled:
        raise NotImplementedError("electron_photon_coupling='joint' does not support reverse shock in this version.")
    if str(config.geometry_kernel).lower() == "chi_eats_2d" or config.downstream_num_chi is not None:
        raise NotImplementedError("electron_photon_coupling='joint' does not support chi-resolved transport.")
    if str(config.structured_backend).lower() != "fortran_1d":
        raise NotImplementedError("electron_photon_coupling='joint' does not support structured backends.")
    if not (config.hadronic.enabled and config.hadronic.epsilon_p > 0.0 and config.hadronic.include_bethe_heitler):
        raise ValueError("electron_photon_coupling='joint' requires enabled forward Bethe-Heitler hadronic physics.")
    if str(config.hadronic.solver).lower() != "am3_1d":
        raise ValueError("electron_photon_coupling='joint' requires hadronic_solver='am3_1d'.")
    if str(config.hadronic.pgamma_scheme).lower() != "hummer_2010_response":
        raise ValueError("electron_photon_coupling='joint' requires pgamma_scheme='hummer_2010_response'.")
    if int(config.index_y) != 1:
        raise ValueError("electron_photon_coupling='joint' requires ssc_cooling_mode='numeric_ic_kn' (index_y=1).")
    if bool(config.electron_adaptive_substeps):
        raise NotImplementedError("electron_photon_coupling='joint' currently requires fixed electron substeps.")


def _validate_multi_density_reverse_config(config: RuntimeConfig) -> None:
    jump_r, _, _ = density_jump_arrays(config)
    if jump_r.size < 1 or not config.reverse_shock.enabled:
        return
    if str(config.electron_solver).lower() not in {"fullhide_1d", "dg_1d"}:
        raise NotImplementedError("multi-density reverse shock v1 requires electron_solver='fullhide_1d' or 'dg_1d'.")
    if str(config.geometry_kernel).lower() != "sed_legacy" or str(config.structured_backend).lower() != "fortran_1d":
        raise NotImplementedError("multi-density reverse shock v1 supports only the direct 1D observer path.")
    if config.hadronic.enabled or config.hadronic.reverse_enabled or float(config.hadronic.reverse_epsilon_p) > 0.0:
        raise NotImplementedError("multi-density reverse shock v1 does not include hadronic processes.")
    if bool(config.include_forward_ssc) or int(config.index_y) != 0:
        raise NotImplementedError("multi-density reverse shock v1 supports electron synchrotron only.")
    if bool(config.reverse_shock.include_ssc):
        raise NotImplementedError("multi-density reverse shock v1 does not include RS SSC.")
    if bool(config.reverse_shock.include_cross_zone_ic):
        raise NotImplementedError("multi-density reverse shock v1 does not include cross-zone IC.")
    if _electron_photon_coupling(config) != _COUPLING_SEPARATED:
        raise NotImplementedError("multi-density reverse shock v1 requires separated electron-photon coupling.")


def _build_photon_field_stage(
    config: RuntimeConfig,
    setup,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    timings: dict[str, float] | None,
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
    config: RuntimeConfig,
    setup,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    photon_field: PhotonFieldState,
    timings: dict[str, float] | None,
    *,
    apply_bh_photon_sink: bool = False,
    merge_secondary_pairs: bool = True,
) -> tuple[ElectronSolution, object | None, SolverAdapterReport]:
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
    if merge_secondary_pairs and hadronic is not None and hadronic.d_n_gam_e_bh is not None:
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
        _apply_hadronic_photon_survival(photon_field, hadronic.pg_photon_survival)
    if apply_bh_photon_sink and hadronic is not None and hadronic.tau_bh is not None:
        _apply_hadronic_photon_survival(photon_field, _hadronic_pg_survival_factor(hadronic.tau_bh))
    return electron, hadronic, report


def _solve_joint_forward_stage(
    config: RuntimeConfig,
    setup,
    dynamics: DynamicsSolution,
    timings: dict[str, float] | None,
) -> tuple[ElectronSolution, PhotonFieldState, object | None, SolverAdapterReport, SolverAdapterReport]:
    _validate_joint_electron_photon_config(config)
    primary_electron, electron_report = _timed_call(
        timings,
        _solver_label(config, "electron"),
        solve_electron,
        setup.boundary,
        dynamics,
        setup.seed_frequency_hz,
        config,
        return_report=True,
    )
    electron = primary_electron
    photon_field = _build_photon_field_stage(config, setup, dynamics, electron, timings)
    hadronic = None
    hadronic_report = _solver_report("hadronic_disabled", "log-gamma-1d", "disabled", backend="none")

    for _ in range(_JOINT_ELECTRON_PHOTON_ITERATIONS):
        electron, hadronic, hadronic_report = _solve_hadronic_stage(
            config,
            setup,
            dynamics,
            primary_electron,
            photon_field,
            timings,
            apply_bh_photon_sink=True,
            merge_secondary_pairs=False,
        )
        photon_field = _build_joint_photon_field_after_hadronic(
            config,
            setup,
            dynamics,
            electron,
            hadronic,
            timings,
        )
        secondary_source_r = _apply_joint_secondary_feedback(
            config,
            setup,
            dynamics,
            electron,
            photon_field,
            hadronic,
        )
        primary_electron, electron_report = _timed_call(
            timings,
            f"{_solver_label(config, 'electron')} [joint cooling]",
            solve_electron_with_cooling_seed,
            setup.boundary,
            dynamics,
            setup.seed_frequency_hz,
            photon_field.hadronic_target_seed,
            config,
            secondary_source_r=secondary_source_r,
            return_report=True,
        )
        electron = primary_electron
        photon_field = _build_joint_photon_field_after_hadronic(
            config,
            setup,
            dynamics,
            electron,
            hadronic,
            timings,
        )
        _apply_joint_secondary_feedback(config, setup, dynamics, electron, photon_field, hadronic)

    return electron, photon_field, hadronic, electron_report, hadronic_report


def _build_joint_photon_field_after_hadronic(
    config: RuntimeConfig,
    setup,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    hadronic,
    timings: dict[str, float] | None,
) -> PhotonFieldState:
    photon_field = _build_photon_field_stage(config, setup, dynamics, electron, timings)
    if hadronic is not None and hadronic.pg_photon_survival is not None:
        _apply_hadronic_photon_survival(photon_field, hadronic.pg_photon_survival)
    if hadronic is not None and hadronic.tau_bh is not None:
        _apply_hadronic_photon_survival(photon_field, _hadronic_pg_survival_factor(hadronic.tau_bh))
    return photon_field


def _apply_joint_secondary_feedback(
    config: RuntimeConfig,
    setup,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    photon_field: PhotonFieldState,
    hadronic,
) -> np.ndarray:
    source = np.zeros_like(np.asarray(electron.d_n_gam_e, dtype=float))
    hadronic_source = None if hadronic is None else getattr(hadronic, "secondary_electron_source_r", None)
    if hadronic_source is not None:
        source = source + np.asarray(hadronic_source, dtype=float)
    if bool(config.hadronic.include_pair_production):
        magnetic_field = compute_magnetic_field(dynamics.r_gamma, dynamics.radius, config)
        pair_lum, pair_seed, tau_pair, pair_density = _compute_pair_production_branch(
            dynamics=dynamics,
            electron=electron,
            combined_seed_field_hz=photon_field.hadronic_target_seed,
            seed_frequency_hz=setup.seed_frequency_hz,
            magnetic_field_g=magnetic_field,
            config=config,
        )
        if hadronic is not None:
            hadronic.l_had_pair_production = pair_lum
        photon_field.hadronic_target_seed = np.asarray(photon_field.hadronic_target_seed, dtype=float) + pair_seed
        photon_field.absorption_syn_seed = np.asarray(photon_field.absorption_syn_seed, dtype=float) + pair_seed
        _apply_hadronic_photon_survival(photon_field, _hadronic_pg_survival_factor(tau_pair))
        source = source + _electron_density_to_source_r(np.asarray(electron.gam_e, dtype=float), pair_density, dynamics.radius)
    return source


def _apply_hadronic_photon_survival(
    photon_field: PhotonFieldState,
    photon_survival: np.ndarray,
) -> None:
    survival = np.asarray(photon_survival, dtype=float)
    photon_field.hadronic_target_seed = np.asarray(photon_field.hadronic_target_seed, dtype=float) * survival
    photon_field.absorption_syn_seed = np.asarray(photon_field.absorption_syn_seed, dtype=float) * survival
    photon_field.absorption_ssc_seed = np.asarray(photon_field.absorption_ssc_seed, dtype=float) * survival
    photon_field.hadronic_forward_ssc_seed = np.asarray(photon_field.hadronic_forward_ssc_seed, dtype=float) * survival


def solve_state_from_setup(
    config: RuntimeConfig,
    setup,
    timings: dict[str, float] | None = None,
    policy: ExecutionPolicy | None = None,
    requested_frequencies_hz: np.ndarray | None = None,
) -> SolveState:
    _validate_multi_density_reverse_config(config)
    execution_policy = ExecutionPolicy(num_threads=config.num_threads) if policy is None else policy
    dynamics, dynamics_report = _timed_call(
        timings,
        _solver_label(config, "dynamics"),
        solve_dynamics,
        setup.boundary,
        config,
        return_report=True,
    )
    if _electron_photon_coupling(config) == _COUPLING_JOINT:
        electron, photon_field, hadronic, electron_report, hadronic_report = _solve_joint_forward_stage(
            config,
            setup,
            dynamics,
            timings,
        )
    else:
        electron, electron_report = _timed_call(
            timings,
            _solver_label(config, "electron"),
            solve_electron,
            setup.boundary,
            dynamics,
            setup.seed_frequency_hz,
            config,
            return_report=True,
        )
        photon_field = _build_photon_field_stage(config, setup, dynamics, electron, timings)
        electron, hadronic, hadronic_report = _solve_hadronic_stage(
            config,
            setup,
            dynamics,
            electron,
            photon_field,
            timings,
        )
    reverse_emission = None
    if config.reverse or config.reverse_shock.enabled:
        reverse_emission = _timed_call(
            timings,
            "ReverseShock.emission",
            solve_reverse_shock_emission,
            setup.boundary,
            dynamics,
            setup.seed_frequency_hz,
            config,
        )
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
    freq_min: float | None = None
    freq_max: float | None = None
    if requested_frequencies_hz is not None:
        requested = np.asarray(requested_frequencies_hz, dtype=float)
        requested = requested[np.isfinite(requested) & (requested > 0.0)]
        if requested.size > 0:
            freq_min = float(np.min(requested))
            freq_max = float(np.max(requested))
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


def _build_observer_setup_from_state(
    state: SolveState,
    observer_time_s: np.ndarray,
) -> SimulationSetup:
    boundary = np.array(state.setup.boundary, dtype=float, copy=True)
    boundary[8] = float(state.config.opening_angle_jet)
    boundary[9] = float(state.config.theta_v)
    return SimulationSetup(
        luminosity_distance_cm=state.setup.luminosity_distance_cm,
        boundary=boundary,
        seed_frequency_hz=state.setup.seed_frequency_hz,
        observer_time_s=np.asarray(observer_time_s, dtype=float),
    )


def project_flux_grid(
    state: SolveState,
    observer_time_s: np.ndarray,
    frequencies_hz: np.ndarray,
    timings: dict[str, float] | None = None,
    mode: str = "full_components",
    projection_kind: str = "lightcurve",
) -> ObsState:
    projection_kind = _normalize_projection_kind(projection_kind)
    setup = _build_observer_setup_from_state(state, observer_time_s)
    if str(state.config.geometry_kernel).lower() == "chi_eats_2d" and projection_kind == "lightcurve":
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


def _require_chi_eats_electron_state(state: SolveState) -> None:
    if not str(state.config.electron_solver).lower().endswith("_2d"):
        raise ValueError("geometry_kernel='chi_eats_2d' requires a 2d electron solver.")
    if _use_direct_chi_projection_contract(state.config):
        missing = [
            name
            for name in (
                "d_n_gam_e_chi",
                "b_chi_g",
                "chi_radius_cm",
                "chi_gamma_bulk",
                "chi_dvolume_weight",
            )
            if getattr(state.electron, name) is None
        ]
        if missing:
            raise RuntimeError("direct chi_eats_2d electron state is missing: " + ", ".join(missing))
        return
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


def _require_top_hat_phi_grid(config: RuntimeConfig) -> None:
    if float(config.theta_v) != 0.0 and int(config.eats_num_phi) == 1:
        raise ValueError("off-axis EATS projection requires eats_num_phi >= 2; eats_num_phi=1 is only valid for on-axis axial collapse.")


_OBSERVED_COMPONENT_ATTRS = (
    ("fwd_sync", "fwd_sync"),
    ("fwd_ssc", "fwd_ssc"),
    ("fwd_hadronic", "fwd_hadronic_gamma"),
    ("fwd_hadronic_bethe_heitler", "fwd_hadronic_bethe_heitler"),
    ("fwd_hadronic_inverse_compton", "fwd_hadronic_inverse_compton"),
    ("fwd_hadronic_pair_production", "fwd_hadronic_pair_production"),
    ("rev_sync", "rev_sync"),
    ("rev_ssc", "rev_ssc"),
    ("cross_ic", "cross_ic"),
)


def _empty_observed_components(total: np.ndarray | None = None) -> dict[str, np.ndarray | None]:
    return {"total": total, **{key: None for key, _attr in _OBSERVED_COMPONENT_ATTRS}}


def _sum_observed_components(observed: dict[str, np.ndarray | None], template: np.ndarray) -> np.ndarray:
    total = np.zeros_like(template)
    for key, _attr in _OBSERVED_COMPONENT_ATTRS:
        value = observed[key]
        if value is not None:
            total = total + value
    return total


def _project_chi_fwd_sync(
    state: SolveState,
    setup,
    frequencies_hz: np.ndarray,
    timings: dict[str, float] | None,
) -> np.ndarray:
    _require_chi_eats_electron_state(state)
    _require_top_hat_phi_grid(state.config)
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    order = np.argsort(frequencies_hz)
    sorted_frequencies = frequencies_hz[order]
    num_phi = 1 if float(state.config.theta_v) == 0.0 else int(state.config.eats_num_phi)
    if _use_direct_chi_projection_contract(state.config):
        flux_sorted = _timed_call(
            timings,
            "Interpolation.sed_interpolation_chi_electron_cached [fwd_sync]",
            Interpolation.sed_interpolation_chi_electron_cached,
            setup.boundary,
            state.components.fwd.characteristic_time_s,
            state.components.fwd.radius_cm,
            np.asarray(state.electron.d_n_gam_e_chi, dtype=float),
            np.asarray(state.electron.b_chi_g, dtype=float),
            np.asarray(state.electron.chi_radius_cm, dtype=float),
            np.asarray(state.electron.chi_gamma_bulk, dtype=float),
            np.asarray(state.electron.chi_dvolume_weight, dtype=float),
            np.asarray(state.electron.gam_e, dtype=float),
            setup.seed_frequency_hz,
            sorted_frequencies,
            setup.observer_time_s,
            state.config.eats_num_theta,
            num_phi,
            state.config.num_threads,
        )
    else:
        source_chi = np.asarray(state.electron.l_syn_spec_chi, dtype=float) * np.asarray(state.observer.prefactor, dtype=float)[:, None, :]
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
            state.config.eats_num_theta,
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
    timings: dict[str, float] | None,
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
            non_chi_total,
            frequencies_hz,
            state.config,
            timings=timings,
            label="Interpolation.sed_interpolation [total_without_fwd_sync]",
        )
        return _empty_observed_components(shell_total_without_fwd_sync + chi_fwd_sync)
    observed = observe_components_from_setup(
        state.config,
        state.components,
        setup,
        frequencies_hz,
        timings=timings,
        mode="full_components",
    )
    observed["fwd_sync"] = chi_fwd_sync
    observed["total"] = _sum_observed_components(observed, observed["fwd_sync"])
    if mode != "full_components":
        raise ValueError(f"Unsupported observe mode: {mode}")
    return observed


def observe_components_from_setup(
    config: RuntimeConfig,
    components: FluxComponents,
    setup,
    frequencies_hz: np.ndarray,
    timings: dict[str, float] | None = None,
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
            components.total,
            frequencies_hz,
            config,
            timings=timings,
            label="Interpolation.sed_interpolation [total]",
        )
        return _empty_observed_components(total)

    observed = _empty_observed_components()
    for key, attr in _OBSERVED_COMPONENT_ATTRS:
        source = getattr(components, attr)
        if source is not None:
            branch = components.rev if key.startswith("rev_") and components.rev is not None else components.fwd
            observed[key] = _project_component(
                setup,
                branch.characteristic_time_s,
                branch.gamma,
                branch.radius_cm,
                source,
                frequencies_hz,
                config,
                timings=timings,
                label=f"Interpolation.sed_interpolation [{key}]",
            )
    observed["total"] = _sum_observed_components(observed, observed["fwd_sync"])
    if timings is not None:
        timings.setdefault("Interpolation.sed_interpolation [total]", 0.0)
    return observed


def _assemble_observer_stage(
    setup,
    config: RuntimeConfig,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    photon_field: PhotonFieldState,
    hadronic,
    reverse_emission: ReverseShockEmission | None,
    timings: dict[str, float] | None = None,
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
        joint_ic_seed=None,
        hadronic_ssa_transfer=np.ones_like(fwd_sync),
        absorbed_fwd_hadronic_gamma=None,
        absorbed_fwd_hadronic_bethe_heitler=None,
        absorbed_fwd_hadronic_inverse_compton=None,
        absorbed_fwd_hadronic_pair_production=None,
    )
    if _electron_photon_coupling(config) == _COUPLING_JOINT and hadronic is not None and hadronic.tau_bh is not None:
        s["tau_extra"] = s["tau_extra"] + np.asarray(hadronic.tau_bh, dtype=float)
        s["joint_ic_seed"] = np.asarray(photon_field.hadronic_target_seed, dtype=float)
    if config.include_forward_ssc:
        seed_for_ssc = s["seed_syn_absorption"] if s["joint_ic_seed"] is None else s["joint_ic_seed"]
        s["fwd_ssc"], s["seed_ssc_total"] = _ssc_spectrum(
            dynamics.radius,
            electron,
            setup.seed_frequency_hz,
            seed_for_ssc,
            config.num_threads,
            timings,
            "Radiation.ssc_spec [FS]",
        )
        if s["pg_photon_survival"] is not None and s["joint_ic_seed"] is None:
            survival = np.asarray(s["pg_photon_survival"], dtype=float)
            s["fwd_ssc"] = np.asarray(s["fwd_ssc"], dtype=float) * survival
            s["seed_ssc_total"] = np.asarray(s["seed_ssc_total"], dtype=float) * survival
    s = _stage_reverse_emission(s, setup, config, dynamics, electron, reverse_emission, timings)
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
    if bool(config.hadronic.include_pair_production):
        s["pair_lum_total"], s["pair_seed_total"], s["tau_pair"], _pair_density = _compute_pair_production_branch(
            dynamics=dynamics,
            electron=electron,
            combined_seed_field_hz=s["seed_syn_absorption"] + s["seed_ssc_total"],
            seed_frequency_hz=setup.seed_frequency_hz,
            magnetic_field_g=s["magnetic_field_g"],
            config=config,
        )
        s["seed_syn_absorption"] = s["seed_syn_absorption"] + s["pair_seed_total"]
    s["tau_extra"] = s["tau_extra"] + s["tau_pair"]
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
    prefactor = absorption / (4.0 * np.pi * setup.luminosity_distance_cm**2) * (1.0 + config.z)
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
            ),
            rev=s["rev_details"],
        ),
    )


def _stage_reverse_emission(
    s: dict,
    setup,
    config: RuntimeConfig,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    reverse_emission: ReverseShockEmission | None,
    timings: dict[str, float] | None,
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
            magnetic_field_g=dynamics.reverse_shock.magnetic_field_g,
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


def _compute_pair_production_branch(
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    combined_seed_field_hz: np.ndarray,
    seed_frequency_hz: np.ndarray,
    magnetic_field_g: np.ndarray,
    config: RuntimeConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    v_seed = np.asarray(seed_frequency_hz, dtype=float)
    seed_field = np.asarray(combined_seed_field_hz, dtype=float)
    gam_e = np.asarray(electron.gam_e, dtype=float)
    radius = np.asarray(dynamics.radius, dtype=float)
    gamma_bulk = np.asarray(dynamics.r_gamma, dtype=float)
    magnetic_field = np.asarray(magnetic_field_g, dtype=float)
    photon_energy_gev, _ = photon_density_hz_to_gev(v_seed, np.ones_like(v_seed))
    e_pair_gev = _aligned_pair_electron_grid(photon_energy_gev)
    gam_pair = e_pair_gev / _ELECTRON_MASS_GEV
    num_nu = int(v_seed.size)
    num_r = int(radius.size)
    pair_lum, pair_seed, tau_pair = np.zeros((3, num_nu, num_r), dtype=float)
    pair_density = np.zeros((gam_e.size, num_r), dtype=float)
    d_n_pair_prev_aligned = np.zeros(gam_pair.size, dtype=float)
    from asgard_core.hadronic_cascade import shell_path_time_seconds

    use_iterative_cascade = int(config.hadronic.pair_cascade_iterations) > 1
    if use_iterative_cascade:
        from asgard_core.hadronic_cascade import compute_time_dependent_pair_cascade_sequence
        cascade = compute_time_dependent_pair_cascade_sequence(
            photon_energy_gev=photon_energy_gev,
            primary_photon_density_per_gev=seed_field / constants.para_h_gev,
            electron_energy_gev=e_pair_gev,
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
            _interp_pair_density_to_electron_grid(gam_pair, np.asarray(cascade.pair_density_per_gamma, dtype=float), gam_e),
        )

    for i_r in range(num_r):
        _, photon_density_per_gev = photon_density_hz_to_gev(v_seed, seed_field[:, i_r])

        ppair = solve_pair_production(
            photon_energy_gev=photon_energy_gev,
            photon_density_per_gev=photon_density_per_gev,
            electron_energy_gev=e_pair_gev,
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
        loss_total = np.asarray(
            hadronic_legacy_module.fs_hadronic_continuous_loss_rates(
                gam_pair,
                float(magnetic_field[i_r]),
                float(radius[i_r]) / (float(gamma_bulk[i_r]) * constants.para_c),
                constants.para_m_e_gev,
                0,
            ),
            dtype=float,
        )
        d_n_pair = np.asarray(
            hadronic_legacy_module.fs_hadronic_advance_energy_loggamma(
                gam_pair,
                d_n_pair_prev_aligned,
                q_pair,
                loss_total,
                _hadronic_shell_comoving_dt_from_radius(radius, gamma_bulk, i_r),
            ),
            dtype=float,
        )
        if magnetic_field[i_r] > 0.0:
            p_syn_i, seed_syn_i = electron_radiation_module.get_syn_selected(
                int(config.index_syn_integr),
                float(radius[i_r]),
                float(magnetic_field[i_r]),
                int(config.num_threads),
                gam_pair, d_n_pair, v_seed,
            )
            pair_lum[:, i_r] = np.asarray(p_syn_i, dtype=float)
            pair_seed[:, i_r] = np.asarray(seed_syn_i, dtype=float)
        pair_density[:, i_r] = np.asarray(
            hadronic_legacy_module.fs_hadronic_positive_loglog_interp(gam_pair, d_n_pair, gam_e),
            dtype=float,
        )
        d_n_pair_prev_aligned = d_n_pair
    return pair_lum, pair_seed, tau_pair, pair_density


def _aligned_pair_electron_grid(photon_energy_gev: np.ndarray) -> np.ndarray:
    photon_energy = np.asarray(photon_energy_gev, dtype=float)
    dln = float(np.log(photon_energy[1] / photon_energy[0]))
    offset = int(np.ceil(np.log(_ELECTRON_MASS_GEV / photon_energy[0]) / dln))
    offset = max(0, offset)
    return photon_energy[0] * np.exp((offset + np.arange(photon_energy.size, dtype=float)) * dln)


def _interp_pair_density_to_electron_grid(
    gamma_pair: np.ndarray,
    pair_density: np.ndarray,
    gam_e: np.ndarray,
) -> np.ndarray:
    density = np.asarray(pair_density, dtype=float)
    out = np.zeros((np.asarray(gam_e, dtype=float).size, density.shape[1]), dtype=float)
    for i_shell in range(density.shape[1]):
        out[:, i_shell] = np.asarray(
            hadronic_legacy_module.fs_hadronic_positive_loglog_interp(gamma_pair, density[:, i_shell], gam_e),
            dtype=float,
        )
    return out


def _electron_density_to_source_r(gam_e: np.ndarray, density_per_gamma: np.ndarray, radius_cm: np.ndarray) -> np.ndarray:
    gamma = np.asarray(gam_e, dtype=float)
    density = np.asarray(density_per_gamma, dtype=float)
    radius = np.asarray(radius_cm, dtype=float)
    source = np.zeros_like(density, dtype=float)
    for i_shell in range(radius.size):
        if i_shell == 0:
            dr = radius[1] - radius[0]
        else:
            dr = radius[i_shell] - radius[i_shell - 1]
        source[:, i_shell] = density[:, i_shell] * gamma * np.log(10.0) / dr
    return source


def _forward_synchrotron_absorption_transfer(
    *,
    electron: ElectronSolution,
    radius_cm: np.ndarray,
    magnetic_field_g: np.ndarray,
    seed_frequency_hz: np.ndarray,
    config: RuntimeConfig,
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
    config: RuntimeConfig,
) -> ElectronSolution:
    bh_distribution = np.asarray(hadronic.d_n_gam_e_bh, dtype=float)
    total_distribution = np.asarray(electron.d_n_gam_e, dtype=float) + bh_distribution
    num_shell = total_distribution.shape[1]
    num_nu = int(np.asarray(seed_frequency_hz, dtype=float).size)
    l_syn_total, seed_syn_total = np.zeros((2, num_nu, num_shell), dtype=float)
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
        d_n_gam_e_bh=bh_distribution,
        d_n_gam_e_chi=None if electron.d_n_gam_e_chi is None else np.asarray(electron.d_n_gam_e_chi, dtype=float),
        chi_grid=None if electron.chi_grid is None else np.asarray(electron.chi_grid, dtype=float),
        l_syn_spec_chi=None if electron.l_syn_spec_chi is None else np.asarray(electron.l_syn_spec_chi, dtype=float),
        seed_syn_chi=None if electron.seed_syn_chi is None else np.asarray(electron.seed_syn_chi, dtype=float),
        tau_syn_chi=None if electron.tau_syn_chi is None else np.asarray(electron.tau_syn_chi, dtype=float),
        chi_radius_cm=None if electron.chi_radius_cm is None else np.asarray(electron.chi_radius_cm, dtype=float),
        chi_gamma_bulk=None if electron.chi_gamma_bulk is None else np.asarray(electron.chi_gamma_bulk, dtype=float),
        chi_dvolume_weight=None if electron.chi_dvolume_weight is None else np.asarray(electron.chi_dvolume_weight, dtype=float),
        b_chi_g=None if electron.b_chi_g is None else np.asarray(electron.b_chi_g, dtype=float),
    )


def _ssc_spectrum(
    radius_cm: np.ndarray,
    electron: ElectronSolution,
    seed_frequency_hz: np.ndarray,
    seed_field: np.ndarray,
    num_threads: int,
    timings: dict[str, float] | None,
    label: str,
) -> tuple[np.ndarray, np.ndarray]:
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
    absorbed_spectral_flux: np.ndarray,
    frequencies_hz: np.ndarray,
    config: RuntimeConfig,
    timings: dict[str, float] | None = None,
    label: str | None = None,
) -> np.ndarray:
    if label is not None and str(config.geometry_kernel).lower() == "sed_adaptive_theta":
        label = label.replace("Interpolation.sed_interpolation", "Interpolation.sed_interpolation_adaptive_theta")
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


def _timed_call(timings: dict[str, float] | None, label: str | None, func, *args, **kwargs):
    start = perf_counter()
    result = func(*args, **kwargs)
    elapsed = perf_counter() - start
    if timings is not None and label is not None:
        timings[label] = timings.get(label, 0.0) + elapsed
    return result
