from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, replace
from time import perf_counter

import numpy as np

from src.Electron.electron_radiation import electron_radiation_kernel as electron_radiation_module
from asgard_core.asgard_config import ExecutionPolicy, RuntimeConfig, SimulationSetup
from asgard_core.hadronic_cascade import solve_paircascade
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
    doppler_factor,
    magfield,
    reverse_mass,
)
from asgard_core.asgard_physics_utils import densityjumps
from asgard_core.asgard_postprocess import observe_flux_batch
from asgard_core.asgard_runtime import (
    _bhsupport,
    pgsurvival,
    _report,
    solve_dynamics,
    solve_electron,
    solve_coolingseed,
    solve_hadronic,
    solve_rsemission,
)
from asgard_core.asgard_setup import build_setup
from src import Interpolation, Radiation, constants


PROJKINDS = {"lightcurve", "sed"}
COUPLING_SEP = "separated"
COUPLING_JOINT = "joint"
JOINT_ITERS = 2
HADRONOPTIONAL = (
    ("bethe_heitler", "l_had_bethe_heitler"),
    ("inverse_compton", "l_had_hadronic_inverse_compton"),
    ("pair_production", "l_had_pair_production"),
    ("pion_synch", "l_had_pion_synch"),
    ("muon_synch", "l_had_muon_synch"),
    ("pion_inverse_compton", "l_had_pion_inverse_compton"),
    ("muon_inverse_compton", "l_had_muon_inverse_compton"),
)


@dataclass
class _CoupledShockGeometry:
    proper_time_s: np.ndarray
    fs_width_cm: np.ndarray
    rs_width_cm: np.ndarray
    center_delay_s: np.ndarray


def build_shockgeo(dynamics, config: RuntimeConfig) -> _CoupledShockGeometry:
    if dynamics.reverse_shock is None:
        raise ValueError("Reverse-shock dynamics are required to build coupled-shock geometry.")
    if config.reverse_shock.delta_t_s is None:
        raise ValueError("ReverseShockConfig.delta_t_s must be set when reverse shock is enabled.")

    radius_cm = dynamics.radius
    gamma = dynamics.r_gamma
    gamma_mean = 0.5 * (gamma[1:] + gamma[:-1])
    beta_mean = np.sqrt(1.0 - gamma_mean**-2)
    d_radius = radius_cm[1:] - radius_cm[:-1]
    proper_time_s = np.empty_like(radius_cm)
    proper_time_s[0] = 0.0
    proper_time_s[1:] = np.cumsum(d_radius / (beta_mean * gamma_mean * constants.para_c))

    fs_width_cm, rs_width_cm = np.zeros((2, radius_cm.size), dtype=float)

    eta_0 = config.eta_0
    shell_mass_g = reverse_mass(config)
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


def crossseed_fields(
    fs_seed_syn: np.ndarray,
    rs_seed_syn: np.ndarray,
    geometry: _CoupledShockGeometry,
    angular_factor: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    tau = geometry.proper_time_s
    tau_ret = tau - geometry.center_delay_s
    seed_fs_to_rs = angular_factor * np.vstack([
        np.interp(tau_ret, tau, seed_nu, left=0.0, right=seed_nu[-1])
        for seed_nu in fs_seed_syn
    ])
    seed_rs_to_fs = angular_factor * np.vstack([
        np.interp(tau_ret, tau, seed_nu, left=0.0, right=seed_nu[-1])
        for seed_nu in rs_seed_syn
    ])
    return seed_fs_to_rs, seed_rs_to_fs


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


def query_cfg(
    config: RuntimeConfig,
    observer_time_s: np.ndarray,
) -> RuntimeConfig:
    query = deepcopy(config)
    observer_time_s = np.asarray(observer_time_s, dtype=float)
    query.num_tobs = observer_time_s.shape[0]
    query.t_obs_min_log10 = float(np.log10(observer_time_s.min()))
    query.t_obs_max_log10 = float(np.log10(observer_time_s.max()))
    return query


def query_setup(
    config: RuntimeConfig,
    observer_time_s: np.ndarray,
    requested_frequencies_hz: np.ndarray | None = None,
):
    observer_time = np.asarray(observer_time_s, dtype=float)
    return replace(
        build_setup(query_cfg(config, observer_time), requested_frequencies_hz),
        observer_time_s=observer_time,
    )


def _solverlabel(config: RuntimeConfig, stage: str) -> str:
    if stage == "dynamics":
        if config.reverse or config.reverse_shock.enabled:
            return "Dynamics.dynamics_reverse"
        return "Dynamics.dynamics_forward"
    if stage == "electron":
        solver_name = config.electron_solver.lower()
        electron_label_map = {
            "fullhide_1d": "Electron.fs_fullhide_1d",
            "fullhide_2d": "Electron.fs_transport_2d",
            "dg_1d": "Electron.fs_dg_1d",
            "t2g1_1d": "Electron.fs_t2g1_1d",
            "slc1_1d": "Electron.fs_slc1_1d",
            "charint_1d": "Electron.fs_charint_1d",
            "charint_2d": "Electron.fs_transport_2d",
            "weno5_1d": "Electron.fs_weno5_1d",
        }
        return electron_label_map.get(solver_name, f"Electron.{solver_name}")
    if stage == "hadronic":
        return "Hadronic.hadronic_1d"
    raise ValueError(f"Unsupported solver stage: {stage}")


def _coupling(config: RuntimeConfig) -> str:
    coupling = str(config.electron_photon_coupling).lower()
    if coupling not in {COUPLING_SEP, COUPLING_JOINT}:
        raise ValueError("electron_photon_coupling must be 'separated' or 'joint'.")
    return coupling


def _checkjoint(config: RuntimeConfig) -> None:
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


def _formalbh(config: RuntimeConfig) -> bool:
    hadronic = config.hadronic
    return (
        bool(hadronic.enabled)
        and float(hadronic.epsilon_p) > 0.0
        and bool(hadronic.include_bethe_heitler)
        and str(hadronic.solver).lower() == "am3_1d"
        and str(hadronic.pgamma_scheme).lower() == "hummer_2010_response"
    )


def _checkforwardbh(config: RuntimeConfig) -> None:
    if _formalbh(config) and str(config.electron_solver).lower() != "fullhide_1d":
        raise NotImplementedError("Forward Bethe-Heitler transport requires electron_solver='fullhide_1d'.")


def _checkmultirs(config: RuntimeConfig) -> None:
    jump_r, _, _ = densityjumps(config)
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
    if _coupling(config) != COUPLING_SEP:
        raise NotImplementedError("multi-density reverse shock v1 requires separated electron-photon coupling.")


def _photonfield(
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
        _, hadronic_forward_ssc_seed = _timed(
            timings,
            "Radiation.ssc_spec [FS-Hadronic]",
            Radiation.ssc_spec,
            dynamics.radius,
            electron.gam_e,
            electron.d_n_gam_e,
            setup.seed_frequency_hz,
            forward_syn_seed,
            config.num_threads,
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


def _hadronstage(
    config: RuntimeConfig,
    setup,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    photon_field: PhotonFieldState,
    timings: dict[str, float] | None,
    *,
    apply_bh_photon_sink: bool = False,
    merge_secondary_pairs: bool = True,
    exact_egrid: bool = False,
) -> tuple[ElectronSolution, object | None, SolverAdapterReport]:
    hadronic, report = _timed(
        timings,
        _solverlabel(config, "hadronic"),
        solve_hadronic,
        setup.boundary,
        dynamics,
        electron,
        photon_field.seed_frequency_hz,
        photon_field.hadronic_target_seed,
        config,
        exact_egrid=exact_egrid,
        return_report=True,
    )
    if merge_secondary_pairs and hadronic is not None and hadronic.d_n_gam_e_bh is not None:
        electron = _mergebh(
            electron,
            hadronic,
            dynamics.radius,
            magfield(dynamics.r_gamma, dynamics.radius, config),
            setup.seed_frequency_hz,
            config,
        )
        hadronic.l_had_bethe_heitler = None
        hadronic.seed_had_bethe_heitler = None
    if hadronic is not None and hadronic.pg_photon_survival is not None:
        _photonsurvive(photon_field, hadronic.pg_photon_survival)
    if apply_bh_photon_sink and hadronic is not None and hadronic.tau_bh is not None:
        _photonsurvive(photon_field, pgsurvival(hadronic.tau_bh))
    return electron, hadronic, report


def _jointstage(
    config: RuntimeConfig,
    setup,
    dynamics: DynamicsSolution,
    timings: dict[str, float] | None,
    grid_top: float | None,
) -> tuple[ElectronSolution, PhotonFieldState, object | None, SolverAdapterReport, SolverAdapterReport]:
    _checkjoint(config)
    primary_electron, electron_report = _timed(
        timings,
        _solverlabel(config, "electron"),
        solve_electron,
        setup.boundary,
        dynamics,
        setup.seed_frequency_hz,
        config,
        grid_top=grid_top,
        return_report=True,
    )
    electron = primary_electron
    photon_field = _photonfield(config, setup, dynamics, electron, timings)
    primary_electron, electron_report = _timed(
        timings,
        f"{_solverlabel(config, 'electron')} [joint cooling seed]",
        solve_coolingseed,
        setup.boundary,
        dynamics,
        setup.seed_frequency_hz,
        photon_field.hadronic_target_seed,
        config,
        grid_top=grid_top,
        return_report=True,
    )
    electron = primary_electron
    photon_field = _photonfield(config, setup, dynamics, electron, timings)
    hadronic = None
    hadronic_report = _report("hadronic_disabled", "log-gamma-1d", "disabled", backend="none")

    for _ in range(JOINT_ITERS):
        electron, hadronic, hadronic_report = _hadronstage(
            config,
            setup,
            dynamics,
            primary_electron,
            photon_field,
            timings,
            apply_bh_photon_sink=True,
            merge_secondary_pairs=False,
            exact_egrid=grid_top is not None,
        )
        photon_field = _photonfield(config, setup, dynamics, electron, timings)
        if hadronic is not None and hadronic.pg_photon_survival is not None:
            _photonsurvive(photon_field, hadronic.pg_photon_survival)
        if hadronic is not None and hadronic.tau_bh is not None:
            _photonsurvive(photon_field, pgsurvival(hadronic.tau_bh))
        secondary_source_r = _jointfeedback(
            config,
            setup,
            dynamics,
            electron,
            photon_field,
            hadronic,
        )
        primary_electron, electron_report = _timed(
            timings,
            f"{_solverlabel(config, 'electron')} [joint cooling]",
            solve_coolingseed,
            setup.boundary,
            dynamics,
            setup.seed_frequency_hz,
            photon_field.hadronic_target_seed,
            config,
            secondary_source_r=secondary_source_r,
            grid_top=grid_top,
            return_report=True,
        )
        electron = primary_electron
        photon_field = _photonfield(config, setup, dynamics, electron, timings)
        if hadronic is not None and hadronic.pg_photon_survival is not None:
            _photonsurvive(photon_field, hadronic.pg_photon_survival)
        if hadronic is not None and hadronic.tau_bh is not None:
            _photonsurvive(photon_field, pgsurvival(hadronic.tau_bh))
        _jointfeedback(config, setup, dynamics, electron, photon_field, hadronic)

    return electron, photon_field, hadronic, electron_report, hadronic_report


def _jointfeedback(
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
        magnetic_field = magfield(dynamics.r_gamma, dynamics.radius, config)
        paircascade = solve_paircascade(
            seed_hz=setup.seed_frequency_hz,
            seed_field=photon_field.hadronic_target_seed,
            electron_gamma=electron.gam_e,
            radius_cm=dynamics.radius,
            gamma_bulk=dynamics.r_gamma,
            tobs_s=dynamics.r_tobs,
            bfield_g=magnetic_field,
            threads=int(config.num_threads),
            syn_index=int(config.index_syn_integr),
            substeps=int(config.hadronic.pair_cascade_iterations),
        )
        if hadronic is not None:
            hadronic.l_had_pair_production = paircascade.syn_lum
        photon_field.hadronic_target_seed = np.asarray(photon_field.hadronic_target_seed, dtype=float) + paircascade.syn_seed
        photon_field.absorption_syn_seed = np.asarray(photon_field.absorption_syn_seed, dtype=float) + paircascade.syn_seed
        _photonsurvive(photon_field, pgsurvival(paircascade.tau_pair))
        source = source + _sourcer(np.asarray(electron.gam_e, dtype=float), paircascade.density_grid, dynamics.radius)
    return source


def _photonsurvive(
    photon_field: PhotonFieldState,
    photon_survival: np.ndarray,
) -> None:
    survival = np.asarray(photon_survival, dtype=float)
    photon_field.hadronic_target_seed = np.asarray(photon_field.hadronic_target_seed, dtype=float) * survival
    photon_field.absorption_syn_seed = np.asarray(photon_field.absorption_syn_seed, dtype=float) * survival
    photon_field.absorption_ssc_seed = np.asarray(photon_field.absorption_ssc_seed, dtype=float) * survival
    photon_field.hadronic_forward_ssc_seed = np.asarray(photon_field.hadronic_forward_ssc_seed, dtype=float) * survival


def _freqrange(frequencies_hz: np.ndarray | None) -> tuple[float | None, float | None]:
    if frequencies_hz is None:
        return None, None
    requested = np.asarray(frequencies_hz, dtype=float)
    requested = requested[np.isfinite(requested) & (requested > 0.0)]
    if requested.size == 0:
        return None, None
    return float(np.min(requested)), float(np.max(requested))


def solve_setup(
    config: RuntimeConfig,
    setup,
    timings: dict[str, float] | None = None,
    policy: ExecutionPolicy | None = None,
    requested_frequencies_hz: np.ndarray | None = None,
    assemble_observer: bool = True,
) -> SolveState:
    _checkmultirs(config)
    _checkforwardbh(config)
    execution_policy = ExecutionPolicy(num_threads=config.num_threads) if policy is None else policy

    # Physical spine: dynamics -> forward electron/photon/hadronic -> reverse -> observer.
    dynamics, dynamics_report = _timed(
        timings,
        _solverlabel(config, "dynamics"),
        solve_dynamics,
        setup.boundary,
        config,
        return_report=True,
    )
    grid_top = None
    exact_egrid = _formalbh(config)
    if exact_egrid:
        _, grid_top = _bhsupport(dynamics, config)

    if _coupling(config) == COUPLING_JOINT:
        electron, photon_field, hadronic, electron_report, hadronic_report = _jointstage(
            config,
            setup,
            dynamics,
            timings,
            grid_top,
        )
    else:
        electron, electron_report = _timed(
            timings,
            _solverlabel(config, "electron"),
            solve_electron,
            setup.boundary,
            dynamics,
            setup.seed_frequency_hz,
            config,
            grid_top=grid_top,
            return_report=True,
        )
        photon_field = _photonfield(config, setup, dynamics, electron, timings)
        electron, hadronic, hadronic_report = _hadronstage(
            config,
            setup,
            dynamics,
            electron,
            photon_field,
            timings,
            exact_egrid=exact_egrid,
        )

    reverse_emission = None
    if config.reverse or config.reverse_shock.enabled:
        reverse_emission = _timed(
            timings,
            "ReverseShock.emission",
            solve_rsemission,
            setup.boundary,
            dynamics,
            setup.seed_frequency_hz,
            config,
        )

    if assemble_observer:
        observer = _observerstage(
            setup,
            config,
            dynamics,
            electron,
            photon_field,
            hadronic,
            reverse_emission,
            timings=timings,
        )
    else:
        observer = _bareobserver(setup, config, dynamics)

    freq_min, freq_max = _freqrange(requested_frequencies_hz)
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


def _bareobserver(setup, config: RuntimeConfig, dynamics) -> ObserverState:
    frequency = np.asarray(setup.seed_frequency_hz, dtype=float)
    radius = np.asarray(dynamics.radius, dtype=float)
    zeros = np.zeros((frequency.size, radius.size), dtype=float)
    fwd = BranchState(
        characteristic_time_s=np.asarray(dynamics.r_tobs, dtype=float),
        gamma=np.asarray(dynamics.r_gamma, dtype=float),
        radius_cm=radius,
        swept_mass_g=np.asarray(dynamics.swept_mass_g, dtype=float),
        doppler=doppler_factor(dynamics.r_gamma, config.z),
        magnetic_field_g=magfield(dynamics.r_gamma, dynamics.radius, config),
    )
    components = FluxComponents(
        total=zeros,
        fwd_sync=zeros.copy(),
        fwd_ssc=zeros.copy(),
        fwd_hadronic_gamma=None,
        fwd_hadronic_bethe_heitler=None,
        fwd_hadronic_inverse_compton=None,
        fwd_hadronic_pair_production=None,
        rev_sync=None,
        rev_hadronic_gamma=None,
        rev_ssc=None,
        cross_ic=None,
        fwd=fwd,
        rev=None,
    )
    return ObserverState(
        prefactor=zeros.copy(),
        tau_extra=zeros.copy(),
        tau_pair=zeros.copy(),
        components=components,
    )


def _observersetup(
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


def project_flux(
    state: SolveState,
    observer_time_s: np.ndarray,
    frequencies_hz: np.ndarray,
    timings: dict[str, float] | None = None,
    mode: str = "full_components",
    projection_kind: str = "lightcurve",
) -> ObsState:
    projection_kind = _projkind(projection_kind)
    setup = _observersetup(state, observer_time_s)
    if str(state.config.geometry_kernel).lower() == "chi_eats_2d" and projection_kind == "lightcurve":
        observed = _observechi(state, setup, frequencies_hz, timings=timings, mode=mode)
    else:
        observed = observe_setup(
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


def _projkind(projection_kind: str) -> str:
    kind = str(projection_kind).lower()
    if kind not in PROJKINDS:
        raise ValueError("projection_kind must be 'lightcurve' or 'sed'.")
    return kind


def _needchi(state: SolveState) -> None:
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


def _needphi(config: RuntimeConfig) -> None:
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
    ("rev_hadronic", "rev_hadronic_gamma"),
    ("rev_ssc", "rev_ssc"),
    ("cross_ic", "cross_ic"),
)
_FWD_COMPONENTS = tuple(item for item in _OBSERVED_COMPONENT_ATTRS if not item[0].startswith("rev_"))
_REV_COMPONENTS = tuple(item for item in _OBSERVED_COMPONENT_ATTRS if item[0].startswith("rev_"))
_TOTAL_COMPONENTS = tuple(
    item
    for item in _OBSERVED_COMPONENT_ATTRS
    if item[0] in (
        "fwd_sync",
        "fwd_ssc",
        "fwd_hadronic",
        "rev_sync",
        "rev_hadronic",
        "rev_ssc",
        "cross_ic",
    )
)
_TOTAL_FWD = tuple(item for item in _TOTAL_COMPONENTS if not item[0].startswith("rev_"))
_TOTAL_REV = tuple(item for item in _TOTAL_COMPONENTS if item[0].startswith("rev_"))


def _emptyobs(total: np.ndarray | None = None) -> dict[str, np.ndarray | None]:
    return {"total": total, **{key: None for key, _attr in _OBSERVED_COMPONENT_ATTRS}}


def _sumobs(
    observed: dict[str, np.ndarray | None],
    template: np.ndarray,
    members: tuple[tuple[str, str], ...] = _TOTAL_COMPONENTS,
) -> np.ndarray:
    total = np.zeros_like(template)
    for key, _attr in members:
        value = observed.get(key)
        if value is not None:
            total = total + value
    return total


def _batchlabel(config: RuntimeConfig, owner: str) -> str:
    kernel = (
        "sed_adaptive_theta_batch"
        if str(config.geometry_kernel).lower() == "sed_adaptive_theta"
        else "sed_interpolation_batch"
    )
    return f"Interpolation.{kernel} [{owner}]"


def _projectbatch(
    setup,
    branch: BranchState | None,
    components: FluxComponents,
    members: tuple[tuple[str, str], ...],
    frequencies_hz: np.ndarray,
    config: RuntimeConfig,
    timings: dict[str, float] | None,
    owner: str,
) -> dict[str, np.ndarray]:
    named = [
        (key, np.asarray(source, dtype=float))
        for key, attr in members
        if (source := getattr(components, attr)) is not None
    ]
    if not named:
        return {}
    shape = (frequencies_hz.size, setup.observer_time_s.size)
    active = [(key, source) for key, source in named if np.any(source)]
    activekeys = {key for key, _source in active}
    result = {key: np.zeros(shape, dtype=float) for key, _source in named if key not in activekeys}
    if not active:
        if timings is not None:
            timings.setdefault(_batchlabel(config, owner), 0.0)
        return result
    _needphi(config)
    sourcebatch = np.asfortranarray(np.stack([source for _key, source in active], axis=2))
    projected = _timed(
        timings,
        _batchlabel(config, owner),
        observe_flux_batch,
        setup,
        branch.characteristic_time_s,
        branch.gamma,
        branch.radius_cm,
        sourcebatch,
        frequencies_hz,
        config,
    )
    result.update({key: projected[:, :, index] for index, (key, _source) in enumerate(active)})
    return result


def _projectshells(
    setup,
    components: FluxComponents,
    frequencies_hz: np.ndarray,
    config: RuntimeConfig,
    timings: dict[str, float] | None,
    fwdmembers: tuple[tuple[str, str], ...],
    revmembers: tuple[tuple[str, str], ...],
    fwdowner: str = "fwd",
    revowner: str = "rev",
) -> dict[str, np.ndarray]:
    projected = _projectbatch(
        setup,
        components.fwd,
        components,
        fwdmembers,
        frequencies_hz,
        config,
        timings,
        fwdowner,
    )
    projected.update(
        _projectbatch(setup, components.rev, components, revmembers, frequencies_hz, config, timings, revowner)
    )
    return projected


def _projectchisync(
    state: SolveState,
    setup,
    frequencies_hz: np.ndarray,
    timings: dict[str, float] | None,
) -> np.ndarray:
    _needchi(state)
    _needphi(state.config)
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    order = np.argsort(frequencies_hz)
    sorted_frequencies = frequencies_hz[order]
    if int(state.config.downstream_num_chi or 0) != 1:
        # Top-hat chi-resolved SSA uses the patch-origin ray kernel, matching the structured-jet occultation contract.
        e = state.electron
        tarr = np.asarray(state.components.fwd.characteristic_time_s, dtype=float)
        rarr = np.asarray(state.components.fwd.radius_cm, dtype=float)
        lchi = np.asarray(e.l_syn_spec_chi, dtype=float)
        tau = np.asarray(e.tau_syn_chi, dtype=float)
        rchi = np.asarray(e.chi_radius_cm, dtype=float)
        gchi = np.asarray(e.chi_gamma_bulk, dtype=float)
        wchi = np.asarray(e.chi_dvolume_weight, dtype=float)
        nr = rarr.shape[0]
        ntheta = int(state.config.eats_num_theta)
        nphi = 1 if float(state.config.theta_v) == 0.0 else 2 * int(state.config.eats_num_phi)
        edges = np.linspace(0.0, float(state.config.opening_angle_jet), ntheta + 1)
        rtobs = np.asfortranarray(np.broadcast_to(tarr[:, None], (nr, ntheta)))
        radius = np.asfortranarray(np.broadcast_to(rarr[:, None], (nr, ntheta)))
        flux = np.asfortranarray(np.broadcast_to(lchi[..., None], (*lchi.shape, ntheta)))
        tauarr = np.asfortranarray(np.broadcast_to(tau[..., None], (*tau.shape, ntheta)))
        chrad = np.asfortranarray(np.broadcast_to(rchi[..., None], (*rchi.shape, ntheta)))
        chgam = np.asfortranarray(np.broadcast_to(gchi[..., None], (*gchi.shape, ntheta)))
        chwt = np.asfortranarray(np.broadcast_to(wchi[..., None], (*wchi.shape, ntheta)))
        flux_sorted = _timed(
            timings,
            "Interpolation.sed_chiring_batchlum_ray [fwd_sync]",
            Interpolation.sed_chiring_batchlum_ray,
            setup.boundary,
            rtobs,
            radius,
            flux,
            tauarr,
            chrad,
            chgam,
            chwt,
            setup.seed_frequency_hz,
            sorted_frequencies,
            setup.observer_time_s,
            np.asfortranarray(edges[:-1]),
            np.asfortranarray(edges[1:]),
            nphi,
        )
    else:
        num_phi = 1 if float(state.config.theta_v) == 0.0 else int(state.config.eats_num_phi)
        source_chi = np.asarray(state.electron.l_syn_spec_chi, dtype=float) * np.asarray(state.observer.prefactor, dtype=float)[:, None, :]
        flux_sorted = _timed(
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


def _observechi(
    state: SolveState,
    setup,
    frequencies_hz: np.ndarray,
    timings: dict[str, float] | None,
    mode: str,
) -> dict[str, np.ndarray | None]:
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    chi_fwd_sync = _projectchisync(state, setup, frequencies_hz, timings)
    if mode not in {"full_components", "total_only"}:
        raise ValueError(f"Unsupported observe mode: {mode}")
    fwdmembers = _TOTAL_FWD if mode == "total_only" else _FWD_COMPONENTS
    revmembers = _TOTAL_REV if mode == "total_only" else _REV_COMPONENTS
    fwdmembers = tuple(item for item in fwdmembers if item[0] != "fwd_sync")
    projected = _projectshells(
        setup,
        state.components,
        frequencies_hz,
        state.config,
        timings,
        fwdmembers,
        revmembers,
        fwdowner="fwd_shell",
        revowner="rev_shell",
    )
    projected["fwd_sync"] = chi_fwd_sync
    if mode == "total_only":
        return _emptyobs(_sumobs(projected, chi_fwd_sync, _TOTAL_COMPONENTS))
    observed = _emptyobs()
    observed.update(projected)
    observed["total"] = _sumobs(observed, observed["fwd_sync"])
    return observed


def observe_setup(
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
    fwdmembers = _TOTAL_FWD if mode == "total_only" else _FWD_COMPONENTS
    revmembers = _TOTAL_REV if mode == "total_only" else _REV_COMPONENTS
    projected = _projectshells(
        setup,
        components,
        frequencies_hz,
        config,
        timings,
        fwdmembers,
        revmembers,
    )
    if mode == "total_only":
        template = np.zeros((frequencies_hz.size, setup.observer_time_s.size), dtype=float)
        total = _sumobs(projected, template, _TOTAL_COMPONENTS)
        return _emptyobs(total)

    observed = _emptyobs()
    observed.update(projected)
    observed["total"] = _sumobs(observed, observed["fwd_sync"])
    return observed


def _observerstage(
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
        rev_hadronic=None,
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
    if _coupling(config) == COUPLING_JOINT and hadronic is not None and hadronic.tau_bh is not None:
        s["tau_extra"] = s["tau_extra"] + np.asarray(hadronic.tau_bh, dtype=float)
        s["joint_ic_seed"] = np.asarray(photon_field.hadronic_target_seed, dtype=float)
    if config.include_forward_ssc:
        seed_for_ssc = s["seed_syn_absorption"] if s["joint_ic_seed"] is None else s["joint_ic_seed"]
        s["fwd_ssc"], s["seed_ssc_total"] = _timed(
            timings,
            "Radiation.ssc_spec [FS]",
            Radiation.ssc_spec,
            dynamics.radius,
            electron.gam_e,
            electron.d_n_gam_e,
            setup.seed_frequency_hz,
            seed_for_ssc,
            config.num_threads,
        )
        if s["pg_photon_survival"] is not None and s["joint_ic_seed"] is None:
            survival = np.asarray(s["pg_photon_survival"], dtype=float)
            s["fwd_ssc"] = np.asarray(s["fwd_ssc"], dtype=float) * survival
            s["seed_ssc_total"] = np.asarray(s["seed_ssc_total"], dtype=float) * survival
    s = _reversestage(s, setup, config, dynamics, electron, reverse_emission, timings)
    s["magnetic_field_g"] = magfield(dynamics.r_gamma, dynamics.radius, config)
    if hadronic is not None:
        s["hadronic_ssa_transfer"] = _fsabsorb(
            electron=electron,
            radius_cm=dynamics.radius,
            magnetic_field_g=s["magnetic_field_g"],
            seed_frequency_hz=setup.seed_frequency_hz,
            config=config,
        )
        s["seed_syn_absorption"] = s["seed_syn_absorption"] + _hadronseed(
            hadronic=hadronic,
            radius_cm=dynamics.radius,
            seed_frequency_hz=setup.seed_frequency_hz,
            ssa_transfer=s["hadronic_ssa_transfer"],
        )
    if bool(config.hadronic.include_pair_production):
        paircascade = solve_paircascade(
            seed_hz=setup.seed_frequency_hz,
            seed_field=s["seed_syn_absorption"] + s["seed_ssc_total"],
            electron_gamma=electron.gam_e,
            radius_cm=dynamics.radius,
            gamma_bulk=dynamics.r_gamma,
            tobs_s=dynamics.r_tobs,
            bfield_g=s["magnetic_field_g"],
            threads=int(config.num_threads),
            syn_index=int(config.index_syn_integr),
            substeps=int(config.hadronic.pair_cascade_iterations),
        )
        s["pair_lum_total"] = paircascade.syn_lum
        s["pair_seed_total"] = paircascade.syn_seed
        s["tau_pair"] = paircascade.tau_pair
        s["seed_syn_absorption"] = s["seed_syn_absorption"] + s["pair_seed_total"]
    s["tau_extra"] = s["tau_extra"] + s["tau_pair"]
    pair_active = bool(config.hadronic.include_pair_production)
    annihilation_seed_syn = np.zeros_like(s["seed_syn_absorption"]) if pair_active else s["seed_syn_absorption"]
    annihilation_seed_ssc = np.zeros_like(s["seed_ssc_total"]) if pair_active else s["seed_ssc_total"]
    absorption = _timed(
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
        hadronic_luminosity = _hadronlum(hadronic, s["hadronic_ssa_transfer"])
        hadronic_gamma_total = hadronic_luminosity["total"]
        bethe_heitler = hadronic_luminosity["bethe_heitler"]
        inverse_compton = hadronic_luminosity["inverse_compton"]
        pair_production = hadronic_luminosity["pair_production"]
        if bethe_heitler is not None:
            s["absorbed_fwd_hadronic_bethe_heitler"] = np.asarray(bethe_heitler, dtype=float) * prefactor
        if inverse_compton is not None:
            s["absorbed_fwd_hadronic_inverse_compton"] = np.asarray(inverse_compton, dtype=float) * prefactor
        if pair_production is not None:
            s["absorbed_fwd_hadronic_pair_production"] = np.asarray(pair_production, dtype=float) * prefactor
        if bool(config.hadronic.include_pair_production):
            hadronic_gamma_total += s["pair_lum_total"]
            s["absorbed_fwd_hadronic_pair_production"] = s["pair_lum_total"] * prefactor
        s["absorbed_fwd_hadronic_gamma"] = hadronic_gamma_total * prefactor
    absorbed_rev_sync = None if s["rev_sync"] is None else s["rev_sync"] * prefactor
    absorbed_rev_hadronic = None if s["rev_hadronic"] is None else s["rev_hadronic"] * prefactor
    absorbed_rev_ssc = None if s["rev_ssc"] is None else s["rev_ssc"] * prefactor
    absorbed_cross_ic = None if s["cross_ic"] is None else s["cross_ic"] * prefactor

    total = absorbed_fwd_sync + absorbed_fwd_ssc
    if absorbed_rev_sync is not None:
        total = total + absorbed_rev_sync
    if absorbed_rev_hadronic is not None:
        total = total + absorbed_rev_hadronic
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
            rev_hadronic_gamma=absorbed_rev_hadronic,
            rev_ssc=absorbed_rev_ssc,
            cross_ic=absorbed_cross_ic,
            fwd=BranchState(
                characteristic_time_s=dynamics.r_tobs,
                gamma=dynamics.r_gamma,
                radius_cm=dynamics.radius,
                swept_mass_g=dynamics.swept_mass_g,
                doppler=doppler_factor(dynamics.r_gamma, config.z),
                magnetic_field_g=magfield(dynamics.r_gamma, dynamics.radius, config),
            ),
            rev=s["rev_details"],
        ),
    )


def _reversestage(
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
            rs_hadronic = reverse_emission.rs_hadronic
            seed_total = np.asarray(rs_hadronic.seed_had_syn, dtype=float)
            seed_bh = getattr(rs_hadronic, "seed_had_bethe_heitler", None)
            if seed_bh is not None:
                seed_total = seed_total + np.asarray(seed_bh, dtype=float)
            l_pg = getattr(rs_hadronic, "l_had_pg_gamma", None)
            for luminosity in (l_pg, *[getattr(rs_hadronic, attr, None) for _name, attr in HADRONOPTIONAL]):
                if luminosity is not None:
                    seed_total = seed_total + _seeddensity(
                        luminosity,
                        radius_cm=dynamics.radius,
                        seed_frequency_hz=setup.seed_frequency_hz,
                    )
            s["seed_syn_absorption"] = s["seed_syn_absorption"] + seed_total
            s["rev_hadronic"] = _rshadlum(reverse_emission.rs_hadronic)
        s["rev_details"] = BranchState(
            characteristic_time_s=dynamics.r_tobs,
            gamma=dynamics.r_gamma,
            radius_cm=dynamics.radius,
            swept_mass_g=dynamics.reverse_shock.swept_mass_g,
            doppler=doppler_factor(dynamics.r_gamma, config.z),
            magnetic_field_g=dynamics.reverse_shock.magnetic_field_g,
        )
        if config.reverse_shock.include_ssc:
            s["rev_ssc"], seed_ssc_rs = _timed(
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
            coupling_geometry = build_shockgeo(dynamics, config)
            seed_fs_to_rs, seed_rs_to_fs = crossseed_fields(
                electron.seed_syn,
                reverse_emission.seed_syn,
                coupling_geometry,
            )
            l_cic_fs_spec, seed_cic_fs = _timed(
                timings,
                "Radiation.ssc_spec [CIC-FS]",
                Radiation.ssc_spec,
                dynamics.radius,
                electron.gam_e,
                electron.d_n_gam_e,
                setup.seed_frequency_hz,
                seed_rs_to_fs,
                config.num_threads,
            )
            l_cic_rs_spec, seed_cic_rs = _timed(
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


def _sourcer(gam_e: np.ndarray, density_per_gamma: np.ndarray, radius_cm: np.ndarray) -> np.ndarray:
    gamma = np.asarray(gam_e, dtype=float)
    density = np.asarray(density_per_gamma, dtype=float)
    radius = np.asarray(radius_cm, dtype=float)
    dr = np.empty_like(radius)
    dr[0] = radius[1] - radius[0]
    dr[1:] = radius[1:] - radius[:-1]
    return density * gamma[:, None] * np.log(10.0) / dr[None, :]


def _shelltransfer(args):
    i, radius, bfield, gam_e, d_n_gam_e_col, frequency = args
    if bfield <= 0.0:
        return i, None
    transfer = electron_radiation_module.syn_transfer(
        float(radius), float(bfield), 1, gam_e, d_n_gam_e_col, frequency)
    return i, transfer


def _fsabsorb(
    *,
    electron: ElectronSolution,
    radius_cm: np.ndarray,
    magnetic_field_g: np.ndarray,
    seed_frequency_hz: np.ndarray,
    config: RuntimeConfig,
) -> np.ndarray:
    radius = radius_cm
    bfield = magnetic_field_g
    frequency = seed_frequency_hz
    gam_e = electron.gam_e
    d_n_gam_e = electron.d_n_gam_e
    transfer = np.ones((frequency.size, radius.size), dtype=float)
    tasks = [(i, radius[i], bfield[i], gam_e, d_n_gam_e[:, i], frequency)
             for i in range(radius.size) if bfield[i] > 0.0]
    if not tasks:
        return transfer
    workers = min(config.num_threads, len(tasks))
    if workers > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=workers) as ex:
            for i, col in ex.map(_shelltransfer, tasks):
                if col is not None:
                    transfer[:, i] = col
    else:
        for task in tasks:
            i, col = _shelltransfer(task)
            if col is not None:
                transfer[:, i] = col
    return transfer


def _rshadlum(rs_hadronic) -> np.ndarray:
    total = np.asarray(rs_hadronic.l_had_syn_spec, dtype=float)
    l_pg = getattr(rs_hadronic, "l_had_pg_gamma", None)
    if l_pg is not None:
        total = total + np.asarray(l_pg, dtype=float)
    for _name, attr in HADRONOPTIONAL:
        luminosity = getattr(rs_hadronic, attr, None)
        if luminosity is not None:
            total = total + np.asarray(luminosity, dtype=float)
    return total


def _hadronlum(hadronic, ssa_transfer: np.ndarray) -> dict[str, np.ndarray | None]:
    base = np.asarray(hadronic.l_had_syn_spec + hadronic.l_had_pg_gamma, dtype=float)
    transfer = np.asarray(ssa_transfer, dtype=float)
    out: dict[str, np.ndarray | None] = {"total": base * transfer}
    for name, attr in HADRONOPTIONAL:
        luminosity = getattr(hadronic, attr, None)
        out[name] = None
        if luminosity is not None:
            out[name] = np.asarray(luminosity, dtype=float) * transfer
            out["total"] = np.asarray(out["total"], dtype=float) + out[name]
    return out


def _seeddensity(
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


def _hadronseed(
    *,
    hadronic,
    radius_cm: np.ndarray,
    seed_frequency_hz: np.ndarray,
    ssa_transfer: np.ndarray,
) -> np.ndarray:
    transfer = np.asarray(ssa_transfer, dtype=float)
    seed_total = np.asarray(hadronic.seed_had_syn, dtype=float) * transfer
    if hadronic.seed_had_bethe_heitler is not None:
        seed_total = seed_total + np.asarray(hadronic.seed_had_bethe_heitler, dtype=float) * transfer

    for luminosity in (hadronic.l_had_pg_gamma, *[getattr(hadronic, attr, None) for _name, attr in HADRONOPTIONAL]):
        if luminosity is None:
            continue
        escaped_luminosity = np.asarray(luminosity, dtype=float) * transfer
        seed_total = seed_total + _seeddensity(
            escaped_luminosity,
            radius_cm=radius_cm,
            seed_frequency_hz=seed_frequency_hz,
        )
    return seed_total


def _mergebh(
    electron: ElectronSolution,
    hadronic,
    radius_cm: np.ndarray,
    magnetic_field_g: np.ndarray,
    seed_frequency_hz: np.ndarray,
    config: RuntimeConfig,
) -> ElectronSolution:
    bh_distribution = np.asarray(hadronic.d_n_gam_e_bh, dtype=float)
    total_distribution = electron.d_n_gam_e + bh_distribution
    gam_e_local = electron.gam_e
    nu_seed = seed_frequency_hz
    num_shell = total_distribution.shape[1]
    num_nu = nu_seed.size
    l_syn_total, seed_syn_total = np.zeros((2, num_nu, num_shell), dtype=float)
    for i_shell in range(num_shell):
        if radius_cm[i_shell] <= 0.0 or magnetic_field_g[i_shell] <= 0.0:
            continue
        _, p_syn_i, seed_syn_i, _ = electron_radiation_module.syn_state(
            config.index_syn_integr,
            float(radius_cm[i_shell]),
            float(magnetic_field_g[i_shell]),
            config.num_threads,
            gam_e_local,
            total_distribution[:, i_shell],
            nu_seed,
        )
        l_syn_total[:, i_shell] = p_syn_i
        seed_syn_total[:, i_shell] = seed_syn_i
    return ElectronSolution(
        gam_e=gam_e_local,
        d_n_gam_e=total_distribution,
        l_syn_spec=l_syn_total,
        seed_syn=seed_syn_total,
        d_n_gam_e_bh=bh_distribution,
        d_n_gam_e_chi=electron.d_n_gam_e_chi,
        chi_grid=electron.chi_grid,
        l_syn_spec_chi=electron.l_syn_spec_chi,
        seed_syn_chi=electron.seed_syn_chi,
        tau_syn_chi=electron.tau_syn_chi,
        chi_radius_cm=electron.chi_radius_cm,
        chi_gamma_bulk=electron.chi_gamma_bulk,
        chi_dvolume_weight=electron.chi_dvolume_weight,
        b_chi_g=electron.b_chi_g,
        nu_m=electron.nu_m,
        nu_c=electron.nu_c,
        nu_a=electron.nu_a,
    )


def _timed(timings: dict[str, float] | None, label: str | None, func, *args, **kwargs):
    start = perf_counter()
    result = func(*args, **kwargs)
    elapsed = perf_counter() - start
    if timings is not None and label is not None:
        timings[label] = timings.get(label, 0.0) + elapsed
    return result
