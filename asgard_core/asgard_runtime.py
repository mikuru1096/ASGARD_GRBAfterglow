from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from importlib import import_module
import time
import numpy as np

from src.Electron.electron_radiation import electron_radiation_kernel as electron_radiation_module
import src.Hadronic.FS_hadronic_1d as hadronic_legacy_module
from asgard_core.hadronic_am3_solver import (
    AM3_PROCESS_LABELS,
    HUMMER2010_RESPONSE_BACKEND,
    KA2008_REFERENCE_BACKEND,
    solve_hummer_2010_response_processes,
    solve_ka2008_reference_processes,
)
from asgard_core.hadronic_acceleration import ACCELERATION_BACKEND, InjectionConfig, estimate_max_gamma, species_injection_operator
from asgard_core.hadronic_bethe_heitler import (
    BETHE_HEITLER_BACKEND,
    ELECTRON_MASS_GEV,
    solve_bethe_heitler,
)
from asgard_core.hadronic_cascade import shell_path_time_seconds
from asgard_core.hadronic_hadronic_ic import HADRONIC_IC_BACKEND, solve_hadronic_inverse_compton
from asgard_core.hadronic_secondary_radiation import (
    SECONDARY_RADIATION_BACKEND,
    SecondarySpeciesDistribution,
    SecondaryTargetPhotonField,
    solve_secondary_radiation_spectrum,
)
from asgard_core.hadronic_species_transport import (
    ChargedMuonDistribution,
    ChargedPionDistribution,
    HadronicSpeciesSources,
    HadronicSpeciesState,
    NeutronDistribution,
    SPECIES_TRANSPORT_BACKEND,
    advance_species_transport_explicit,
)
from asgard_core.hadronic_hummer import GEV_TO_ERG, HUMMER2010_DECAY_BACKEND, HUMMER2010_OPERATOR_BACKEND, PROTON_MASS_GEV, solve_hummer2010_pgamma
from src import constants
from asgard_core.hadronic_pp import PP_DELTA_BACKEND, solve_pp_delta
from asgard_core.hadronic_pgamma import photon_density_hz_to_gev
from asgard_core.asgard_config import FitConfig
from asgard_core.asgard_numpy import trapezoid
from asgard_core.asgard_types import (
    ReverseShockParameters,
    ReverseShockDynamics,
    DynamicsSolution,
    ElectronSolution,
    HadronicSolution,
    ReverseShockEmission,
    SolverAdapterReport,
)
from asgard_core.asgard_physics_utils import ambient_density, doppler_denominator, compute_magnetic_field
from src import Dynamics, Electron, constants


_ELECTRON_SOLVER_ALIASES = {
    "fullhide": "fullhide_1d",
    "fullhide_1d": "fullhide_1d",
    "fullhide_2d": "fullhide_2d",
    "slc1": "slc1_1d",
    "slc1_1d": "slc1_1d",
    "charint": "charint_1d",
    "charint_1d": "charint_1d",
    "charint_2d": "charint_2d",
    "t2g1": "t2g1_1d",
    "t2g1_1d": "t2g1_1d",
    "weno5": "weno5_1d",
    "weno5_1d": "weno5_1d",
}

_COOLING_KERNEL_ALIASES = {
    "legacy": "legacy",
}

_ELECTRON_MODULES = {
    "charint_1d": "src.Electron.FS_electron_charint_1d",
    "charint_2d": "src.Electron.FS_electron_charint_2d",
    "fullhide_1d": "src.Electron.FS_electron_fullhide_1d",
    "fullhide_2d": "src.Electron.FS_electron_fullhide_2d",
    "slc1_1d": "src.Electron.FS_electron_slc1_1d",
    "t2g1_1d": "src.Electron.FS_electron_t2g1_1d",
    "weno5_1d": "src.Electron.FS_electron_weno5_1d",
}


@lru_cache(maxsize=None)
def _electron_module(solver: str):
    return import_module(_ELECTRON_MODULES[solver])


@lru_cache(maxsize=1)
def _electron_reverse_module():
    return import_module("src.Electron.electron_reverse_kernel")

_HADRONIC_SOLVER_ALIASES = {
    "legacy": "legacy_1d",
    "legacy_1d": "legacy_1d",
    "am3": "am3_1d",
    "am3_1d": "am3_1d",
}

_PGAMMA_SCHEME_DISABLED = "disabled"
_PGAMMA_SCHEME_HUMMER2010_RESPONSE = "hummer_2010_response"
_PGAMMA_SCHEME_KA2008_REFERENCE = "ka2008_reference"

_PGAMMA_SCHEME_ALIASES = {
    "disabled": _PGAMMA_SCHEME_DISABLED,
    "hummer_2010_response": _PGAMMA_SCHEME_HUMMER2010_RESPONSE,
    "ka2008_reference": _PGAMMA_SCHEME_KA2008_REFERENCE,
}
_PGAMMA_SCHEME_LEGACY_ALIASES = {
    "hummer_2010": _PGAMMA_SCHEME_HUMMER2010_RESPONSE,
    "hummer2010": _PGAMMA_SCHEME_HUMMER2010_RESPONSE,
    "ka2008": _PGAMMA_SCHEME_KA2008_REFERENCE,
    "kelner_aharonian_2008": _PGAMMA_SCHEME_KA2008_REFERENCE,
    "aharonian_2008": _PGAMMA_SCHEME_KA2008_REFERENCE,
    "am3_reference": _PGAMMA_SCHEME_HUMMER2010_RESPONSE,
    "am3_numeric": _PGAMMA_SCHEME_HUMMER2010_RESPONSE,
    "am3_numerical": _PGAMMA_SCHEME_HUMMER2010_RESPONSE,
    "am3": _PGAMMA_SCHEME_HUMMER2010_RESPONSE,
}


def _solver_report(
    solver: str,
    grid_semantics: str,
    status: str,
    **diagnostics,
) -> SolverAdapterReport:
    return SolverAdapterReport(
        solver=solver,
        grid_semantics=grid_semantics,
        status=status,
        diagnostics=diagnostics,
    )


def solve_dynamics(
    boundary: np.ndarray,
    config: FitConfig,
    *,
    return_report: bool = False,
) -> DynamicsSolution | tuple[DynamicsSolution, SolverAdapterReport]:
    reverse_params = _resolve_reverse_shock_parameters(config)
    if reverse_params is None:
        r_tobs, r_gamma, radius, swept_mass_g = Dynamics.dynamics_forward(boundary, config.num_r, config.index_dyn)
        solution = DynamicsSolution(r_tobs, r_gamma, radius, swept_mass_g)
        if return_report:
            return solution, _solver_report(
                "dynamics_forward",
                "shell-radius",
                "ok",
                kernel=str(config.dynamics_kernel),
                num_r=int(config.num_r),
            )
        return solution

    (
        t_cross,
        r_cross,
        e3_cross,
        gam20,
        r_tobs,
        r_gamma,
        radius,
        swept_mass_g,
        swept_reverse_mass_g,
        magnetic_field_g,
    ) = Dynamics.dynamics_reverse(
        reverse_params.delta_t_s,
        reverse_params.epsilon_e,
        reverse_params.epsilon_b,
        reverse_params.p,
        reverse_params.f_e,
        boundary,
        config.num_r,
    )
    reverse_shock = ReverseShockDynamics(
        t_cross,
        r_cross,
        e3_cross,
        gam20,
        swept_reverse_mass_g,
        magnetic_field_g,
    )
    solution = DynamicsSolution(r_tobs, r_gamma, radius, swept_mass_g, reverse_shock=reverse_shock)
    if return_report:
        return solution, _solver_report(
            "dynamics_reverse",
            "shell-radius",
            "ok",
            kernel=str(config.dynamics_kernel),
            num_r=int(config.num_r),
        )
    return solution


def solve_electron(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: FitConfig,
    *,
    return_report: bool = False,
) -> ElectronSolution | tuple[ElectronSolution, SolverAdapterReport]:
    solver_name = _resolve_electron_solver(config)
    _resolve_cooling_kernel(config)
    if solver_name == "weno5_1d":
        electron_weno5_module = _electron_module(solver_name)
        gam_e, d_n_gam_e, l_syn_spec, seed_syn = electron_weno5_module.fs_electron_weno5_1d(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.num_threads,
        )
        nu_m, nu_c, nu_a = _compute_characteristic_frequencies_weno5(
            config,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            gam_e,
            d_n_gam_e,
        )
        solution = _build_electron_solution(
            config,
            dynamics,
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu_m,
            nu_c,
            nu_a,
        )
        if return_report:
            return solution, _solver_report(
                solver_name,
                "log-gamma-1d",
                "ok",
                num_gam_e=int(config.num_gam_e),
                num_chi=1,
            )
        return solution

    if solver_name == "t2g1_1d":
        electron_t2g1_module = _electron_module(solver_name)
        gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_t2g1_module.fs_electron_t2g1_1d(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
        )
        solution = _build_electron_solution(
            config,
            dynamics,
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu_m,
            nu_c,
            nu_a,
        )
        if return_report:
            return solution, _solver_report(
                solver_name,
                "log-gamma-1d",
                "ok",
                num_gam_e=int(config.num_gam_e),
                num_chi=1,
            )
        return solution

    if solver_name == "slc1_1d":
        electron_slc1_module = _electron_module(solver_name)
        gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_slc1_module.fs_electron_slc1_1d(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
        )
        solution = _build_electron_solution(
            config,
            dynamics,
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu_m,
            nu_c,
            nu_a,
        )
        if return_report:
            return solution, _solver_report(
                solver_name,
                "log-gamma-1d",
                "ok",
                num_gam_e=int(config.num_gam_e),
                num_chi=1,
            )
        return solution

    if solver_name == "charint_1d":
        electron_charint_module = _electron_module(solver_name)
        gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_charint_module.fs_electron_charint_1d(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            config.num_gam_e,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
            1 if config.electron_adaptive_substeps else 0,
            config.electron_substep_rtol,
            config.electron_substep_min,
            config.electron_substep_max,
        )
        solution = _build_electron_solution(
            config,
            dynamics,
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu_m,
            nu_c,
            nu_a,
        )
        if return_report:
            return solution, _solver_report(
                solver_name,
                "log-gamma-1d",
                "ok",
                num_gam_e=int(config.num_gam_e),
                num_chi=1,
            )
        return solution

    if solver_name == "charint_2d":
        electron_charint_2d_module = _electron_module(solver_name)
        num_chi = _resolve_num_chi(config, solver_name)
        gam_e, d_n_gam_e_chi, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_charint_2d_module.fs_electron_charint_2d(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            config.num_gam_e,
            num_chi,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
        )
        chi_grid = _build_log_chi_grid(dynamics.r_gamma, num_chi)
        solution = _build_electron_solution(
            config,
            dynamics,
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu_m,
            nu_c,
            nu_a,
            d_n_gam_e_chi=d_n_gam_e_chi,
            chi_grid=chi_grid,
        )
        if return_report:
            return solution, _solver_report(
                solver_name,
                "log-gamma-log-chi-2d",
                "ok",
                num_gam_e=int(config.num_gam_e),
                num_chi=int(num_chi),
            )
        return solution

    if solver_name == "fullhide_2d":
        electron_fullhide_2d_module = _electron_module(solver_name)
        num_chi = _resolve_num_chi(config, solver_name)
        gam_e, d_n_gam_e_chi, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_fullhide_2d_module.fs_electron_fullhide_2d(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            config.num_gam_e,
            num_chi,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
        )
        chi_grid = _build_log_chi_grid(dynamics.r_gamma, num_chi)
        solution = _build_electron_solution(
            config,
            dynamics,
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu_m,
            nu_c,
            nu_a,
            d_n_gam_e_chi=d_n_gam_e_chi,
            chi_grid=chi_grid,
        )
        if return_report:
            return solution, _solver_report(
                solver_name,
                "log-gamma-log-chi-2d",
                "ok",
                num_gam_e=int(config.num_gam_e),
                num_chi=int(num_chi),
            )
        return solution

    electron_fullhide_1d_module = _electron_module("fullhide_1d")
    gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = electron_fullhide_1d_module.fs_electron_fullhide_1d(
        boundary,
        dynamics.r_tobs,
        dynamics.r_gamma,
        dynamics.radius,
        v_seed,
        config.num_gam_e,
        config.index_y,
        config.index_syn_integr,
        config.num_threads,
        1 if config.electron_adaptive_substeps else 0,
        config.electron_substep_rtol,
        config.electron_substep_min,
        config.electron_substep_max,
    )
    solution = _build_electron_solution(
        config,
        dynamics,
        gam_e,
        d_n_gam_e,
        l_syn_spec,
        seed_syn,
        nu_m,
        nu_c,
        nu_a,
    )
    if return_report:
        return solution, _solver_report(
            solver_name,
            "log-gamma-1d",
            "ok",
            num_gam_e=int(config.num_gam_e),
            num_chi=1,
        )
    return solution


def _resolve_electron_solver(config: FitConfig) -> str:
    if config.weno5:
        return "weno5_1d"
    solver_name = _ELECTRON_SOLVER_ALIASES.get(config.electron_solver.lower())
    if solver_name is None:
        raise ValueError(f"Unsupported electron solver: {config.electron_solver}")
    return solver_name


def _resolve_cooling_kernel(config: FitConfig) -> str:
    cooling_kernel = _COOLING_KERNEL_ALIASES.get(config.cooling_kernel.lower())
    if cooling_kernel is None:
        raise ValueError(f"Unsupported cooling kernel: {config.cooling_kernel}")
    return cooling_kernel


def _resolve_hadronic_solver(config: FitConfig) -> str:
    solver_name = _HADRONIC_SOLVER_ALIASES.get(str(config.hadronic.solver).lower())
    if solver_name is None:
        raise ValueError(f"Unsupported hadronic solver: {config.hadronic.solver}")
    return solver_name


def _resolve_pgamma_scheme(config: FitConfig) -> str:
    key = str(config.hadronic.pgamma_scheme).lower()
    scheme_name = _PGAMMA_SCHEME_ALIASES.get(key)
    if scheme_name is None:
        scheme_name = _PGAMMA_SCHEME_LEGACY_ALIASES.get(key)
    if scheme_name is None:
        raise ValueError(f"Unsupported p-gamma scheme: {config.hadronic.pgamma_scheme}")
    return scheme_name


def _resolve_num_chi(config: FitConfig, solver_name: str | None = None) -> int:
    resolved_solver = _resolve_electron_solver(config) if solver_name is None else solver_name
    user_value = config.num_chi
    if resolved_solver.endswith("_1d"):
        return 1 if user_value is None else 1
    if user_value is None:
        return 64
    if int(user_value) < 2:
        raise ValueError("num_chi must be >= 2 for 2d electron solvers.")
    return int(user_value)


def _build_log_chi_grid(r_gamma: np.ndarray, num_chi: int) -> np.ndarray:
    gamma_arr = np.asarray(r_gamma, dtype=float)
    chi_max = 1.0 + 8.0 * np.max(gamma_arr * gamma_arr)
    deta = np.log10(chi_max) / float(num_chi)
    eta_grid = (np.arange(num_chi, dtype=float) + 0.5) * deta
    return np.power(10.0, eta_grid)


def _build_electron_solution(
    config: FitConfig,
    dynamics: DynamicsSolution,
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
    l_syn_spec: np.ndarray,
    seed_syn: np.ndarray,
    nu_m: np.ndarray,
    nu_c: np.ndarray,
    nu_a: np.ndarray,
    *,
    d_n_gam_e_chi: np.ndarray | None = None,
    chi_grid: np.ndarray | None = None,
) -> ElectronSolution:
    cooling_timescale_s, dynamical_timescale_s = _compute_forward_timescales(
        dynamics.r_gamma,
        dynamics.radius,
        nu_c,
        config,
    )
    return ElectronSolution(
        gam_e=np.asarray(gam_e, dtype=float),
        d_n_gam_e=np.asarray(d_n_gam_e, dtype=float),
        l_syn_spec=np.asarray(l_syn_spec, dtype=float),
        seed_syn=np.asarray(seed_syn, dtype=float),
        nu_m=np.asarray(nu_m, dtype=float),
        nu_c=np.asarray(nu_c, dtype=float),
        nu_a=np.asarray(nu_a, dtype=float),
        d_n_gam_e_bh=None,
        d_n_gam_e_chi=None if d_n_gam_e_chi is None else np.asarray(d_n_gam_e_chi, dtype=float),
        chi_grid=None if chi_grid is None else np.asarray(chi_grid, dtype=float),
        cooling_timescale_s=cooling_timescale_s,
        dynamical_timescale_s=dynamical_timescale_s,
    )


def solve_hadronic(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    v_seed: np.ndarray,
    seed_target: np.ndarray,
    config: FitConfig,
    *,
    return_report: bool = False,
) -> HadronicSolution | None | tuple[HadronicSolution | None, SolverAdapterReport]:
    del boundary
    report_status = "ok"
    if not bool(config.hadronic.enabled):
        report = _solver_report("hadronic_disabled", "log-gamma-1d", "disabled", backend="none")
        return (None, report) if return_report else None
    if float(config.hadronic.epsilon_p) <= 0.0:
        report = _solver_report("hadronic_disabled", "log-gamma-1d", "disabled", backend="none")
        return (None, report) if return_report else None
    hadronic_solver = _resolve_hadronic_solver(config)
    pgamma_scheme = _resolve_pgamma_scheme(config)
    electron_solver = _resolve_electron_solver(config)
    if not electron_solver.endswith("_1d"):
        raise ValueError("The current hadronic kernel only supports 1d forward-shock electron solvers.")

    v_seed_arr = np.asarray(v_seed, dtype=float)
    seed_target_arr = np.asarray(seed_target, dtype=float)
    magnetic_field_g = np.asarray(compute_magnetic_field(dynamics.r_gamma, dynamics.radius, config), dtype=float)
    num_nu = int(v_seed_arr.shape[0])
    num_r = int(np.asarray(dynamics.radius, dtype=float).shape[0])
    num_gam_p = int(config.hadronic.num_gam_p)
    num_nu_nu = int(config.hadronic.num_nu_nu)

    if hadronic_solver == "legacy_1d" and (bool(config.hadronic.include_pg) or bool(config.hadronic.include_neutrino)):
        raise NotImplementedError(
            "legacy_1d only supports proton transport plus proton synchrotron. "
            "Use hadronic_solver='am3_1d' for p-gamma and neutrino channels."
        )
    if (bool(config.hadronic.include_pg) or bool(config.hadronic.include_neutrino)) and pgamma_scheme == _PGAMMA_SCHEME_DISABLED:
        raise ValueError(
            "p-gamma and neutrino channels require an explicit pgamma_scheme. "
            "Choose 'hummer_2010_response' or 'ka2008_reference' "
            "(legacy aliases: 'hummer_2010', 'am3_reference', 'aharonian_2008', etc.)."
        )
    shell_energy_inj_erg = _hadronic_shell_injection_energy(
        dynamics.radius,
        dynamics.r_gamma,
        config,
    )

    if hadronic_solver == "am3_1d" and pgamma_scheme == _PGAMMA_SCHEME_HUMMER2010_RESPONSE:
        solution = _solve_hadronic_hummer_transport_coupled(
            dynamics,
            magnetic_field_g,
            electron.gam_e,
            v_seed_arr,
            seed_target_arr,
            shell_energy_inj_erg,
            config,
        )
        if return_report:
            return solution, _solver_report(
                hadronic_solver,
                "log-gamma-1d",
                report_status,
                backend=HUMMER2010_RESPONSE_BACKEND,
                pgamma_scheme=pgamma_scheme,
                pgamma_operator_backend=HUMMER2010_OPERATOR_BACKEND,
                decay_backend=HUMMER2010_DECAY_BACKEND,
                acceleration_backend=ACCELERATION_BACKEND,
                species_transport_backend=SPECIES_TRANSPORT_BACKEND,
                secondary_radiation_backend=SECONDARY_RADIATION_BACKEND,
                bethe_heitler_backend=BETHE_HEITLER_BACKEND if bool(config.hadronic.include_bethe_heitler) else "disabled",
                pp_backend=PP_DELTA_BACKEND if bool(config.hadronic.include_pp) else "disabled",
                hadronic_ic_backend=HADRONIC_IC_BACKEND if bool(config.hadronic.include_hadronic_inverse_compton) else "disabled",
                num_gam_p=num_gam_p,
                num_nu=num_nu,
                timings=solution.timings,
            )
        return solution

    outputs = hadronic_legacy_module.fs_hadronic_1d(
        dynamics.r_tobs,
        dynamics.r_gamma,
        dynamics.radius,
        shell_energy_inj_erg,
        magnetic_field_g,
        v_seed_arr,
        seed_target_arr,
        float(config.hadronic.p_p),
        float(config.hadronic.epsilon_p),
        float(config.hadronic.eta_acc),
        1 if config.hadronic.include_proton_synch else 0,
        0,
        0,
        num_gam_p,
        num_nu_nu,
        int(config.num_threads),
    )
    gam_p, d_n_gam_p, l_had_syn_spec, seed_had_syn, l_had_pg_gamma, neutrino_frequency_hz, neutrino_luminosity = outputs

    if hadronic_solver == "am3_1d":
        if pgamma_scheme == _PGAMMA_SCHEME_KA2008_REFERENCE:
            am3_output = solve_ka2008_reference_processes(
                dynamics.radius,
                gam_p,
                d_n_gam_p,
                v_seed_arr,
                seed_target_arr,
                num_nu_nu,
                include_pg=bool(config.hadronic.include_pg),
                include_neutrino=bool(config.hadronic.include_neutrino),
            )
        else:
            am3_output = solve_hummer_2010_response_processes(
                dynamics.radius,
                gam_p,
                d_n_gam_p,
                v_seed_arr,
                seed_target_arr,
                num_nu_nu,
                include_pg=bool(config.hadronic.include_pg),
                include_neutrino=bool(config.hadronic.include_neutrino),
            )
        if bool(config.hadronic.include_pg):
            l_had_pg_gamma = am3_output.l_had_pg_gamma
        if bool(config.hadronic.include_neutrino):
            neutrino_frequency_hz = am3_output.neutrino_frequency_hz
            neutrino_luminosity = am3_output.neutrino_luminosity
        am3_process_power = am3_output.am3_process_power
    else:
        am3_process_power = np.zeros((len(AM3_PROCESS_LABELS), num_gam_p, num_r), dtype=float)

    solution = HadronicSolution(
        solver=hadronic_solver,
        gam_p=np.asarray(gam_p, dtype=float),
        d_n_gam_p=np.asarray(d_n_gam_p, dtype=float).reshape(num_gam_p, num_r),
        l_had_syn_spec=np.asarray(l_had_syn_spec, dtype=float).reshape(num_nu, num_r),
        seed_had_syn=np.asarray(seed_had_syn, dtype=float).reshape(num_nu, num_r),
        l_had_pg_gamma=np.asarray(l_had_pg_gamma, dtype=float).reshape(num_nu, num_r),
        neutrino_frequency_hz=np.asarray(neutrino_frequency_hz, dtype=float),
        neutrino_luminosity=np.asarray(neutrino_luminosity, dtype=float).reshape(num_nu_nu, num_r),
        l_had_bethe_heitler=None,
        seed_had_bethe_heitler=None,
        l_had_hadronic_inverse_compton=None,
        l_had_pair_production=None,
        tau_pg=None,
        pg_photon_survival=None,
        am3_process_power=np.asarray(am3_process_power, dtype=float).reshape(len(AM3_PROCESS_LABELS), num_gam_p, num_r),
    )
    if return_report:
        backend = "fortran_core" if hadronic_solver == "legacy_1d" else KA2008_REFERENCE_BACKEND
        return solution, _solver_report(
            hadronic_solver,
            "log-gamma-1d",
            report_status,
            backend=backend,
            pgamma_scheme=pgamma_scheme,
            num_gam_p=num_gam_p,
            num_nu=num_nu,
        )
    return solution


def _solve_hadronic_hummer_transport_coupled(
    dynamics: DynamicsSolution,
    magnetic_field_g: np.ndarray,
    electron_gamma: np.ndarray,
    v_seed_hz: np.ndarray,
    seed_target_hz: np.ndarray,
    shell_energy_inj_erg: np.ndarray,
    config: FitConfig,
) -> HadronicSolution:
    radius = np.asarray(dynamics.radius, dtype=float)
    gamma_bulk = np.asarray(dynamics.r_gamma, dtype=float)
    tobs = np.asarray(dynamics.r_tobs, dtype=float)
    b_field = np.asarray(magnetic_field_g, dtype=float)
    v_seed_arr = np.asarray(v_seed_hz, dtype=float)
    seed_target_arr = np.asarray(seed_target_hz, dtype=float)
    num_r = int(radius.size)
    num_nu = int(v_seed_arr.size)
    num_gam_p = int(config.hadronic.num_gam_p)
    num_nu_nu = int(config.hadronic.num_nu_nu)
    gam_e = np.asarray(electron_gamma, dtype=float)
    gam_e_edge = _hadronic_build_gamma_edges(gam_e)
    electron_energy_gev = gam_e * ELECTRON_MASS_GEV

    gam_p = np.logspace(
        np.log10(1.0 + 1.0e-3),
        np.log10(_hadronic_global_gamma_p_max(radius, gamma_bulk, b_field, config)),
        num_gam_p,
    )
    gam_edge = _hadronic_build_gamma_edges(gam_p)
    neutrino_frequency_hz = np.logspace(
        np.log10(1.0e-3 * constants.para_gev2hz),
        np.log10(1.0e8 * constants.para_gev2hz),
        num_nu_nu,
    )
    neutrino_energy_gev = constants.para_h_gev * neutrino_frequency_hz
    photon_energy_gev, _ = photon_density_hz_to_gev(v_seed_arr, np.ones_like(v_seed_arr))
    gam_secondary = np.array(gam_p, copy=True)
    neutron_energy_gev = gam_secondary * constants.para_m_n_gev
    pion_energy_gev = gam_secondary * constants.para_m_pi_charged_gev
    muon_energy_gev = gam_secondary * constants.para_m_mu_gev

    d_n_gam_p = np.zeros((num_gam_p, num_r), dtype=float)
    l_had_syn_spec = np.zeros((num_nu, num_r), dtype=float)
    seed_had_syn = np.zeros((num_nu, num_r), dtype=float)
    l_had_pg_gamma = np.zeros((num_nu, num_r), dtype=float)
    l_had_bh = np.zeros((num_nu, num_r), dtype=float)
    seed_had_bh = np.zeros((num_nu, num_r), dtype=float)
    d_n_gam_e_bh = np.zeros((gam_e.size, num_r), dtype=float)
    l_had_hic = np.zeros((num_nu, num_r), dtype=float)
    tau_pg = np.zeros((num_nu, num_r), dtype=float)
    pg_photon_survival = np.ones((num_nu, num_r), dtype=float)
    neutrino_luminosity = np.zeros((num_nu_nu, num_r), dtype=float)
    am3_process_power = np.zeros((len(AM3_PROCESS_LABELS), num_gam_p, num_r), dtype=float)
    d_n_gam_n = np.zeros((gam_secondary.size, num_r), dtype=float)
    d_n_gam_pi_plus = np.zeros((gam_secondary.size, num_r), dtype=float)
    d_n_gam_pi_minus = np.zeros((gam_secondary.size, num_r), dtype=float)
    d_n_gam_mu_minus_left = np.zeros((gam_secondary.size, num_r), dtype=float)
    d_n_gam_mu_minus_right = np.zeros((gam_secondary.size, num_r), dtype=float)
    d_n_gam_mu_plus_left = np.zeros((gam_secondary.size, num_r), dtype=float)
    d_n_gam_mu_plus_right = np.zeros((gam_secondary.size, num_r), dtype=float)
    l_had_pion_synch = np.zeros((num_nu, num_r), dtype=float)
    l_had_muon_synch = np.zeros((num_nu, num_r), dtype=float)
    l_had_pion_ic = np.zeros((num_nu, num_r), dtype=float)
    l_had_muon_ic = np.zeros((num_nu, num_r), dtype=float)

    d_n_prev = np.zeros(num_gam_p, dtype=float)
    d_n_bh_prev = np.zeros(gam_e.size, dtype=float)
    species_state_prev = HadronicSpeciesState(
        neutron=NeutronDistribution(gamma=gam_secondary, density_per_gamma=np.zeros_like(gam_secondary)),
        charged_pion=ChargedPionDistribution(
            gamma=gam_secondary,
            plus_density_per_gamma=np.zeros_like(gam_secondary),
            minus_density_per_gamma=np.zeros_like(gam_secondary),
        ),
        charged_muon=ChargedMuonDistribution(
            gamma=gam_secondary,
            minus_left_density_per_gamma=np.zeros_like(gam_secondary),
            minus_right_density_per_gamma=np.zeros_like(gam_secondary),
            plus_left_density_per_gamma=np.zeros_like(gam_secondary),
            plus_right_density_per_gamma=np.zeros_like(gam_secondary),
        ),
    )
    process_energy_gev = photon_energy_gev
    shell_volume_cm3 = _hadronic_shell_volumes_from_radius(radius)

    timings = {
        "pg_interaction": 0.0,
        "bethe_heitler": 0.0,
        "pp_delta": 0.0,
        "species_transport": 0.0,
        "secondary_radiation": 0.0,
        "hadronic_ic": 0.0,
        "proton_synch": 0.0,
        "bh_electron_radiation": 0.0,
        "total": 0.0,
    }
    t_total_start = time.perf_counter()

    for i_r in range(num_r):
        dt_s = _hadronic_shell_dt(tobs, i_r)
        t_dyn_s = _hadronic_dynamical_time(radius[i_r], gamma_bulk[i_r])
        gam_p_min = max(1.1, float(gamma_bulk[i_r]))
        q_inj = dt_s * species_injection_operator(
            gam_p,
            InjectionConfig(
                species="proton",
                luminosity_erg_s=float(shell_energy_inj_erg[i_r]) / dt_s,
                spectral_index=float(config.hadronic.p_p),
                gamma_min=gam_p_min,
                gamma_max=float(gam_p[-1]),
            ),
        )
        shell_volume_loc = float(shell_volume_cm3[i_r])
        d_n_trial = _hadronic_advance_energy_loggamma(
            gam_p,
            gam_edge,
            d_n_prev,
            q_inj,
            _hadronic_continuous_loss_rates(gam_p, float(b_field[i_r]), t_dyn_s,
                quantum_syn=bool(config.hadronic.quantum_syn),
                mass_gev=constants.para_m_p_gev),
            dt_s,
        )
        proton_density_trial_per_gev = d_n_trial / (shell_volume_loc * PROTON_MASS_GEV)
        neutron_density_trial_per_gev = species_state_prev.neutron.density_per_gamma / (shell_volume_loc * constants.para_m_n_gev)
        t_pg_start = time.perf_counter()
        _, photon_density_per_gev_trial = photon_density_hz_to_gev(v_seed_arr, seed_target_arr[:, i_r])
        backend_tau = solve_hummer2010_pgamma(
            proton_energy_gev=gam_p * PROTON_MASS_GEV,
            proton_density_per_gev=proton_density_trial_per_gev,
            photon_energy_gev=photon_energy_gev,
            photon_density_per_gev=photon_density_per_gev_trial,
            gamma_energy_gev=photon_energy_gev,
            neutrino_energy_gev=neutrino_energy_gev,
            process_energy_gev=process_energy_gev,
            neutron_density_per_gev=neutron_density_trial_per_gev,
        )
        tau_pg[:, i_r], pg_photon_survival[:, i_r] = _hadronic_pg_local_closure(
            radius[i_r],
            gamma_bulk[i_r],
            np.asarray(backend_tau.photon_loss_rate, dtype=float),
        )
        local_seed_target_hz = seed_target_arr[:, i_r] * pg_photon_survival[:, i_r]
        _, photon_density_per_gev = photon_density_hz_to_gev(v_seed_arr, local_seed_target_hz)
        backend = solve_hummer2010_pgamma(
            proton_energy_gev=gam_p * PROTON_MASS_GEV,
            proton_density_per_gev=proton_density_trial_per_gev,
            photon_energy_gev=photon_energy_gev,
            photon_density_per_gev=photon_density_per_gev,
            gamma_energy_gev=photon_energy_gev,
            neutrino_energy_gev=neutrino_energy_gev,
            process_energy_gev=process_energy_gev,
            neutron_density_per_gev=neutron_density_trial_per_gev,
        )
        timings["pg_interaction"] += time.perf_counter() - t_pg_start
        bh_output = None
        bh_loss = np.zeros_like(gam_p)
        pp_pair_q = np.zeros_like(gam_e)
        pp_gamma_lum = np.zeros_like(v_seed_arr)
        pp_nu_lum = np.zeros_like(neutrino_frequency_hz)
        if bool(config.hadronic.include_bethe_heitler):
            t_bh_start = time.perf_counter()
            bh_output = solve_bethe_heitler(
                proton_energy_gev=gam_p * PROTON_MASS_GEV,
                proton_density_per_gev=proton_density_trial_per_gev,
                photon_energy_gev=photon_energy_gev,
                photon_density_per_gev=photon_density_per_gev,
                electron_energy_gev=electron_energy_gev,
            )
            bh_loss = np.maximum(-np.asarray(bh_output.proton_loss_rate, dtype=float), 0.0)
            timings["bethe_heitler"] += time.perf_counter() - t_bh_start
        pp_loss = np.zeros_like(gam_p)
        if bool(config.hadronic.include_pp):
            t_pp_start = time.perf_counter()
            target_density_cm3 = float(ambient_density(np.array([radius[i_r]], dtype=float), config)[0])
            pp_output = solve_pp_delta(
                proton_energy_gev=gam_p * PROTON_MASS_GEV,
                proton_density_per_gev=proton_density_trial_per_gev,
                target_proton_density_cm3=target_density_cm3,
                gamma_energy_gev=photon_energy_gev,
                neutrino_energy_gev=neutrino_energy_gev,
                pair_energy_gev=electron_energy_gev,
            )
            pp_loss = np.maximum(-np.asarray(pp_output.proton_loss_rate, dtype=float), 0.0)
            pp_gamma_lum = _energy_luminosity_from_rate_spectrum(
                pp_output.gamma_energy_gev,
                pp_output.gamma_rate_per_gev,
                shell_volume_loc,
            )
            pp_nu_lum = _energy_luminosity_from_rate_spectrum(
                pp_output.neutrino_energy_gev,
                pp_output.neutrino_rate_per_gev,
                shell_volume_loc,
            )
            pp_pair_q = shell_volume_loc * np.asarray(pp_output.pair_rate_per_gev, dtype=float) * ELECTRON_MASS_GEV
            timings["pp_delta"] += time.perf_counter() - t_pp_start

        loss_total = _hadronic_continuous_loss_rates(
            gam_p, float(b_field[i_r]), t_dyn_s,
            quantum_syn=bool(config.hadronic.quantum_syn),
            mass_gev=constants.para_m_p_gev,
        ) + bh_loss + pp_loss
        d_n_next = _hadronic_advance_energy_loggamma(gam_p, gam_edge, d_n_prev, q_inj, loss_total, dt_s)
        pg_loss_rate = np.asarray(backend.proton_loss_rate, dtype=float)
        proton_sink = np.exp(-dt_s * pg_loss_rate)
        pg_reinj_dt = _hadronic_interaction_effective_time(pg_loss_rate, dt_s)
        proton_reinj = pg_reinj_dt * shell_volume_loc * np.asarray(backend.proton_reinjection_rate_per_gev, dtype=float) * PROTON_MASS_GEV
        d_n_next = np.maximum(d_n_next * proton_sink + proton_reinj, 0.0)
        d_n_gam_p[:, i_r] = d_n_next

        if bool(config.hadronic.include_proton_synch):
            t_syn_start = time.perf_counter()
            p_syn_i, seed_syn_i = _hadronic_proton_syn_state(
                float(radius[i_r]),
                max(float(b_field[i_r]), 1.0e-30),
                gam_p,
                d_n_next,
                v_seed_arr,
            )
            l_had_syn_spec[:, i_r] = p_syn_i
            seed_had_syn[:, i_r] = seed_syn_i
            timings["proton_synch"] += time.perf_counter() - t_syn_start

        divergence_rate_s_inv = 3.0 / t_dyn_s
        t_species_start = time.perf_counter()
        species_sources = HadronicSpeciesSources(
            neutron_per_gamma_s=_interp_source_per_gamma(
                backend.hadron_energy_gev,
                backend.neutron_reinjection_rate_per_gev,
                neutron_energy_gev,
                constants.para_m_n_gev,
                shell_volume_loc,
            ),
            charged_pion_plus_per_gamma_s=_interp_source_per_gamma(
                backend.hadron_energy_gev,
                backend.pion_plus_source_rate_per_gev,
                pion_energy_gev,
                constants.para_m_pi_charged_gev,
                shell_volume_loc,
            ),
            charged_pion_minus_per_gamma_s=_interp_source_per_gamma(
                backend.hadron_energy_gev,
                backend.pion_minus_source_rate_per_gev,
                pion_energy_gev,
                constants.para_m_pi_charged_gev,
                shell_volume_loc,
            ),
            charged_muon_minus_left_per_gamma_s=_interp_source_per_gamma(
                backend.hadron_energy_gev,
                backend.muon_minus_left_source_rate_per_gev,
                muon_energy_gev,
                constants.para_m_mu_gev,
                shell_volume_loc,
            ),
            charged_muon_minus_right_per_gamma_s=_interp_source_per_gamma(
                backend.hadron_energy_gev,
                backend.muon_minus_right_source_rate_per_gev,
                muon_energy_gev,
                constants.para_m_mu_gev,
                shell_volume_loc,
            ),
            charged_muon_plus_left_per_gamma_s=_interp_source_per_gamma(
                backend.hadron_energy_gev,
                backend.muon_plus_left_source_rate_per_gev,
                muon_energy_gev,
                constants.para_m_mu_gev,
                shell_volume_loc,
            ),
            charged_muon_plus_right_per_gamma_s=_interp_source_per_gamma(
                backend.hadron_energy_gev,
                backend.muon_plus_right_source_rate_per_gev,
                muon_energy_gev,
                constants.para_m_mu_gev,
                shell_volume_loc,
            ),
        )
        species_state_next = advance_species_transport_explicit(
            state=species_state_prev,
            sources=species_sources,
            dt_s=dt_s,
            b_field_g=float(b_field[i_r]),
            divergence_rate_s_inv=divergence_rate_s_inv,
        )
        timings["species_transport"] += time.perf_counter() - t_species_start
        neutron_sink = np.maximum(
            1.0 - dt_s * _interp_positive_loglog(backend.hadron_energy_gev, backend.neutron_loss_rate, neutron_energy_gev),
            0.0,
        )
        species_state_next = HadronicSpeciesState(
            neutron=NeutronDistribution(
                gamma=gam_secondary,
                density_per_gamma=species_state_next.neutron.density_per_gamma * neutron_sink,
            ),
            charged_pion=species_state_next.charged_pion,
            charged_muon=species_state_next.charged_muon,
        )
        d_n_gam_n[:, i_r] = species_state_next.neutron.density_per_gamma
        d_n_gam_pi_plus[:, i_r] = species_state_next.charged_pion.plus_density_per_gamma
        d_n_gam_pi_minus[:, i_r] = species_state_next.charged_pion.minus_density_per_gamma
        d_n_gam_mu_minus_left[:, i_r] = species_state_next.charged_muon.minus_left_density_per_gamma
        d_n_gam_mu_minus_right[:, i_r] = species_state_next.charged_muon.minus_right_density_per_gamma
        d_n_gam_mu_plus_left[:, i_r] = species_state_next.charged_muon.plus_left_density_per_gamma
        d_n_gam_mu_plus_right[:, i_r] = species_state_next.charged_muon.plus_right_density_per_gamma

        t_sec_start = time.perf_counter()
        photon_energy_aligned_gev = _hadronic_aligned_photon_grid(gam_p * PROTON_MASS_GEV, photon_energy_gev)
        photon_density_aligned_per_gev = _interp_positive_loglog(
            photon_energy_gev,
            photon_density_per_gev,
            photon_energy_aligned_gev,
        )
        secondary_species = SecondarySpeciesDistribution(
            pion_plus_per_gev=_interp_distribution_per_gev(pion_energy_gev, species_state_next.charged_pion.plus_density_per_gamma / (shell_volume_loc * constants.para_m_pi_charged_gev), gam_p * PROTON_MASS_GEV),
            pion_minus_per_gev=_interp_distribution_per_gev(pion_energy_gev, species_state_next.charged_pion.minus_density_per_gamma / (shell_volume_loc * constants.para_m_pi_charged_gev), gam_p * PROTON_MASS_GEV),
            muon_minus_left_per_gev=_interp_distribution_per_gev(muon_energy_gev, species_state_next.charged_muon.minus_left_density_per_gamma / (shell_volume_loc * constants.para_m_mu_gev), gam_p * PROTON_MASS_GEV),
            muon_minus_right_per_gev=_interp_distribution_per_gev(muon_energy_gev, species_state_next.charged_muon.minus_right_density_per_gamma / (shell_volume_loc * constants.para_m_mu_gev), gam_p * PROTON_MASS_GEV),
            muon_plus_left_per_gev=_interp_distribution_per_gev(muon_energy_gev, species_state_next.charged_muon.plus_left_density_per_gamma / (shell_volume_loc * constants.para_m_mu_gev), gam_p * PROTON_MASS_GEV),
            muon_plus_right_per_gev=_interp_distribution_per_gev(muon_energy_gev, species_state_next.charged_muon.plus_right_density_per_gamma / (shell_volume_loc * constants.para_m_mu_gev), gam_p * PROTON_MASS_GEV),
        )
        secondary_target = SecondaryTargetPhotonField(
            photon_energy_gev=photon_energy_aligned_gev,
            photons_on_had_grid_per_gev=photon_density_aligned_per_gev,
            ind_min_energy_pho_hadgrid=0,
        )
        secondary_radiation = solve_secondary_radiation_spectrum(
            hadron_energy_gev=gam_p * PROTON_MASS_GEV,
            species=secondary_species,
            target=secondary_target,
            magnetic_field_g=float(b_field[i_r]),
        )
        timings["secondary_radiation"] += time.perf_counter() - t_sec_start
        l_had_pion_synch[:, i_r] = _interp_positive_loglog(
            secondary_radiation.photon_energy_gev,
            _energy_luminosity_from_rate_spectrum(
                secondary_radiation.photon_energy_gev,
                secondary_radiation.pion_synch_rate_per_gev,
                shell_volume_loc,
            ),
            photon_energy_gev,
        )
        l_had_muon_synch[:, i_r] = _interp_positive_loglog(
            secondary_radiation.photon_energy_gev,
            _energy_luminosity_from_rate_spectrum(
                secondary_radiation.photon_energy_gev,
                secondary_radiation.muon_synch_rate_per_gev,
                shell_volume_loc,
            ),
            photon_energy_gev,
        )
        l_had_pion_ic[:, i_r] = _interp_positive_loglog(
            secondary_radiation.photon_energy_gev,
            _energy_luminosity_from_rate_spectrum(
                secondary_radiation.photon_energy_gev,
                secondary_radiation.pion_ic_rate_per_gev,
                shell_volume_loc,
            ),
            photon_energy_gev,
        )
        l_had_muon_ic[:, i_r] = _interp_positive_loglog(
            secondary_radiation.photon_energy_gev,
            _energy_luminosity_from_rate_spectrum(
                secondary_radiation.photon_energy_gev,
                secondary_radiation.muon_ic_rate_per_gev,
                shell_volume_loc,
            ),
            photon_energy_gev,
        )

        if bool(config.hadronic.include_pg):
            l_had_pg_gamma[:, i_r] = _energy_luminosity_from_rate_spectrum(
                backend.gamma_energy_gev,
                backend.gamma_rate_per_gev,
                shell_volume_loc,
            )
            if bool(config.hadronic.include_pp):
                l_had_pg_gamma[:, i_r] += pp_gamma_lum
            proc_lum = np.vstack([
                _energy_luminosity_from_rate_spectrum(
                    backend.process_energy_gev,
                    backend.process_rate_per_gev[i_proc],
                    shell_volume_loc,
                )
                for i_proc in range(len(AM3_PROCESS_LABELS))
            ])
            proton_energy_weight = d_n_next * (gam_p * PROTON_MASS_GEV)
            total_weight = float(trapezoid(proton_energy_weight, gam_p * PROTON_MASS_GEV))
            if total_weight > 0.0 and np.isfinite(total_weight):
                normalized_weight = proton_energy_weight / total_weight
                am3_process_power[0, :, i_r] = normalized_weight * float(trapezoid(proc_lum[0], backend.process_energy_gev))
                am3_process_power[1, :, i_r] = normalized_weight * float(trapezoid(proc_lum[1], backend.process_energy_gev))
                am3_process_power[2, :, i_r] = normalized_weight * float(trapezoid(proc_lum[2], backend.process_energy_gev))

        if bool(config.hadronic.include_neutrino):
            neutrino_luminosity[:, i_r] = _energy_luminosity_from_rate_spectrum(
                backend.neutrino_energy_gev,
                backend.neutrino_rate_per_gev,
                shell_volume_loc,
            )
            if bool(config.hadronic.include_pp):
                neutrino_luminosity[:, i_r] += pp_nu_lum

        if bool(config.hadronic.include_hadronic_inverse_compton):
            t_hic_start = time.perf_counter()
            hic_photon_energy_gev = _hadronic_aligned_photon_grid(gam_p * PROTON_MASS_GEV, photon_energy_gev)
            hic_photon_density_per_gev = _interp_positive_loglog(
                photon_energy_gev,
                photon_density_per_gev,
                hic_photon_energy_gev,
            )
            hic_output = solve_hadronic_inverse_compton(
                hadron_energy_gev=gam_p * PROTON_MASS_GEV,
                photon_energy_gev=hic_photon_energy_gev,
                photons_on_had_grid_per_gev=hic_photon_density_per_gev,
                protons_per_gev=proton_density_trial_per_gev,
                pion_plus_per_gev=np.zeros_like(gam_p, dtype=float),
                pion_minus_per_gev=np.zeros_like(gam_p, dtype=float),
                muon_minus_left_per_gev=np.zeros_like(gam_p, dtype=float),
                muon_minus_right_per_gev=np.zeros_like(gam_p, dtype=float),
                muon_plus_left_per_gev=np.zeros_like(gam_p, dtype=float),
                muon_plus_right_per_gev=np.zeros_like(gam_p, dtype=float),
            )
            hic_rate = (
                np.asarray(hic_output.epsilon_p_ic, dtype=float)
                + np.asarray(hic_output.epsilon_pi_ic, dtype=float)
                + np.asarray(hic_output.epsilon_mu_ic, dtype=float)
            ) / hic_output.photon_energy_gev
            hic_lum_aligned = _energy_luminosity_from_rate_spectrum(
                hic_output.photon_energy_gev,
                hic_rate,
                shell_volume_loc,
            )
            l_had_hic[:, i_r] = _interp_positive_loglog(
                hic_output.photon_energy_gev,
                hic_lum_aligned,
                photon_energy_gev,
            )
            timings["hadronic_ic"] += time.perf_counter() - t_hic_start

        if bh_output is not None or bool(config.hadronic.include_pp):
            t_bhe_start = time.perf_counter()
            q_bh = np.array(pp_pair_q, copy=True)
            if bh_output is not None:
                q_bh += shell_volume_loc * np.asarray(bh_output.pair_rate_per_gev, dtype=float) * ELECTRON_MASS_GEV
            loss_bh_e = _hadronic_electron_loss_rates(
                gam_e, float(b_field[i_r]), t_dyn_s,
                quantum_syn=bool(config.hadronic.quantum_syn),
            )
            d_n_bh_next = _hadronic_advance_energy_loggamma(
                gam_e,
                gam_e_edge,
                d_n_bh_prev,
                q_bh,
                loss_bh_e,
                dt_s,
            )
            p_bh_i, seed_bh_i = electron_radiation_module.get_syn_selected(
                int(config.index_syn_integr),
                float(radius[i_r]),
                max(float(b_field[i_r]), 1.0e-30),
                int(config.num_threads),
                gam_e,
                d_n_bh_next,
                v_seed_arr,
            )
            l_had_bh[:, i_r] = np.asarray(p_bh_i, dtype=float)
            seed_had_bh[:, i_r] = np.asarray(seed_bh_i, dtype=float)
            d_n_gam_e_bh[:, i_r] = d_n_bh_next
            d_n_bh_prev = d_n_bh_next
            timings["bh_electron_radiation"] += time.perf_counter() - t_bhe_start

        d_n_prev = d_n_next
        species_state_prev = species_state_next

    timings["total"] = time.perf_counter() - t_total_start
    l_had_syn_spec *= pg_photon_survival
    seed_had_syn *= pg_photon_survival
    l_had_pg_gamma *= pg_photon_survival
    if bool(config.hadronic.include_bethe_heitler):
        l_had_bh *= pg_photon_survival
        seed_had_bh *= pg_photon_survival
    if bool(config.hadronic.include_hadronic_inverse_compton):
        l_had_hic *= pg_photon_survival
    l_had_pion_synch *= pg_photon_survival
    l_had_muon_synch *= pg_photon_survival
    l_had_pion_ic *= pg_photon_survival
    l_had_muon_ic *= pg_photon_survival

    return HadronicSolution(
        solver="am3_1d",
        gam_p=gam_p,
        d_n_gam_p=d_n_gam_p,
        l_had_syn_spec=l_had_syn_spec,
        seed_had_syn=seed_had_syn,
        l_had_pg_gamma=l_had_pg_gamma,
        neutrino_frequency_hz=neutrino_frequency_hz,
        neutrino_luminosity=neutrino_luminosity,
        l_had_bethe_heitler=l_had_bh if bool(config.hadronic.include_bethe_heitler) else None,
        seed_had_bethe_heitler=seed_had_bh if bool(config.hadronic.include_bethe_heitler) else None,
        d_n_gam_e_bh=d_n_gam_e_bh if (bool(config.hadronic.include_bethe_heitler) or bool(config.hadronic.include_pp)) else None,
        l_had_hadronic_inverse_compton=l_had_hic if bool(config.hadronic.include_hadronic_inverse_compton) else None,
        gam_secondary=gam_secondary,
        d_n_gam_n=d_n_gam_n,
        d_n_gam_pi_plus=d_n_gam_pi_plus,
        d_n_gam_pi_minus=d_n_gam_pi_minus,
        d_n_gam_mu_minus_left=d_n_gam_mu_minus_left,
        d_n_gam_mu_minus_right=d_n_gam_mu_minus_right,
        d_n_gam_mu_plus_left=d_n_gam_mu_plus_left,
        d_n_gam_mu_plus_right=d_n_gam_mu_plus_right,
        l_had_pion_synch=l_had_pion_synch,
        l_had_muon_synch=l_had_muon_synch,
        l_had_pion_inverse_compton=l_had_pion_ic,
        l_had_muon_inverse_compton=l_had_muon_ic,
        tau_pg=tau_pg if (bool(config.hadronic.include_pg) or bool(config.hadronic.include_neutrino)) else None,
        pg_photon_survival=pg_photon_survival if (bool(config.hadronic.include_pg) or bool(config.hadronic.include_neutrino)) else None,
        am3_process_power=am3_process_power,
        timings=timings,
    )


def _hadronic_global_gamma_p_max(
    radius_cm: np.ndarray,
    gamma_bulk: np.ndarray,
    b_field_g: np.ndarray,
    config: FitConfig,
) -> float:
    estimates = [
        estimate_max_gamma(
            species="proton",
            b_field_g=float(bi),
            radius_cm=float(ri),
            gamma_bulk=max(float(gi), 1.0),
            eta_acc=max(float(config.hadronic.eta_acc), 1.0e-12),
        ).gamma_max
        for ri, gi, bi in zip(np.asarray(radius_cm, dtype=float), np.asarray(gamma_bulk, dtype=float), np.asarray(b_field_g, dtype=float), strict=True)
    ]
    return float(max(10.0, np.max(np.asarray(estimates, dtype=float))))


def _hadronic_shell_volumes_from_radius(radius_cm: np.ndarray) -> np.ndarray:
    radius = np.asarray(radius_cm, dtype=float)
    prev_radius = np.empty_like(radius)
    prev_radius[0] = 0.0
    prev_radius[1:] = radius[:-1]
    return (4.0 / 3.0) * np.pi * np.maximum(radius**3 - prev_radius**3, 0.0)


def _hadronic_build_gamma_edges(gam_p: np.ndarray) -> np.ndarray:
    gam = np.asarray(gam_p, dtype=float)
    if gam.size == 1:
        return np.array([max(1.0, 0.5 * gam[0]), 2.0 * gam[0]], dtype=float)
    edge = np.empty(gam.size + 1, dtype=float)
    edge[0] = gam[0] * np.sqrt(gam[0] / gam[1])
    edge[1:-1] = np.sqrt(gam[:-1] * gam[1:])
    edge[-1] = gam[-1] * np.sqrt(gam[-1] / gam[-2])
    return edge


def _hadronic_shell_dt(r_tobs: np.ndarray, i_shell: int) -> float:
    if i_shell <= 0:
        return max(float(r_tobs[0]), 1.0)
    return max(float(r_tobs[i_shell] - r_tobs[i_shell - 1]), 1.0)


def _hadronic_dynamical_time(radius_cm: float, gamma_bulk: float) -> float:
    return max(float(radius_cm) / (max(float(gamma_bulk), 1.0) * constants.para_c), 1.0)


def _hadronic_tau_pg_shell(radius_cm: float, gamma_bulk: float, alpha_gamma_s_inv: np.ndarray) -> np.ndarray:
    path_time_s = max(shell_path_time_seconds(float(radius_cm), max(float(gamma_bulk), 1.0)), 0.0)
    return np.maximum(np.asarray(alpha_gamma_s_inv, dtype=float), 0.0) * path_time_s


def _hadronic_pg_survival_factor(tau_pg: np.ndarray) -> np.ndarray:
    tau = np.asarray(tau_pg, dtype=float)
    survival = np.ones_like(tau, dtype=float)
    small = (tau > 0.0) & (tau < 1.0e-6)
    if np.any(small):
        tau_small = tau[small]
        survival[small] = 1.0 - 0.5 * tau_small + (tau_small * tau_small) / 6.0
    large = tau >= 1.0e-6
    if np.any(large):
        tau_large = tau[large]
        survival[large] = -np.expm1(-tau_large) / tau_large
    return survival


def _hadronic_pg_local_closure(
    radius_cm: float,
    gamma_bulk: float,
    alpha_gamma_s_inv: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    tau_pg = _hadronic_tau_pg_shell(radius_cm, gamma_bulk, alpha_gamma_s_inv)
    return tau_pg, _hadronic_pg_survival_factor(tau_pg)


def _hadronic_interaction_effective_time(rate_s_inv: np.ndarray, dt_s: float) -> np.ndarray:
    rate = np.asarray(rate_s_inv, dtype=float)
    out = np.full_like(rate, float(dt_s), dtype=float)
    active = rate > 0.0
    out[active] = -np.expm1(-rate[active] * float(dt_s)) / rate[active]
    return out


def _hadronic_aligned_photon_grid(hadron_energy_gev: np.ndarray, photon_energy_gev: np.ndarray) -> np.ndarray:
    hadron_energy = np.asarray(hadron_energy_gev, dtype=float)
    photon_energy = np.asarray(photon_energy_gev, dtype=float)
    dln_had = float(np.diff(np.log(hadron_energy))[0])
    log_min = float(np.log(photon_energy[0]))
    log_max = float(np.log(photon_energy[-1]))
    nbin = int(np.ceil((log_max - log_min) / dln_had)) + 1
    return np.exp(log_min + dln_had * np.arange(nbin, dtype=float))


def _hadronic_electron_loss_rates(
    gam_e: np.ndarray,
    b_field_g: float,
    t_dyn_s: float,
    quantum_syn: bool = False,
) -> np.ndarray:
    gamma_e = np.asarray(gam_e, dtype=float)
    b2 = max(float(b_field_g), 0.0) ** 2
    coeff_syn = constants.para_sigmat * b2 / (6.0 * np.pi * constants.para_m_e * constants.para_c)
    loss_syn = coeff_syn * gamma_e * gamma_e
    if quantum_syn:
        loss_syn = loss_syn * _quantum_syn_cooling_factor(
            gamma_e, b_field_g, constants.para_m_e_gev,
        )
    loss_ad = gamma_e / max(float(t_dyn_s), 1.0)
    return loss_syn + loss_ad


def _quantum_syn_cooling_factor(
    gamma: np.ndarray, b_field_g: float, mass_gev: float,
) -> np.ndarray:
    b_crit = 4.414e13
    chi = np.asarray(gamma, dtype=float) * max(float(b_field_g), 0.0) / b_crit
    chi = chi * (constants.para_m_e_gev / max(float(mass_gev), 1e-30))
    out = np.ones_like(chi, dtype=float)
    active = chi > 1e-6
    chi23 = chi[active] ** (2.0 / 3.0)
    out[active] = 1.0 / (1.0 + np.sqrt(2.0) * chi23) ** 2
    return out


def _interp_positive_loglog(
    x_old: np.ndarray,
    y_old: np.ndarray,
    x_new: np.ndarray,
) -> np.ndarray:
    xsrc = np.asarray(x_old, dtype=float)
    ysrc = np.asarray(y_old, dtype=float)
    xdst = np.asarray(x_new, dtype=float)
    out = np.zeros_like(xdst, dtype=float)
    mask = np.isfinite(xsrc) & np.isfinite(ysrc) & (xsrc > 0.0) & (ysrc > 0.0)
    if np.count_nonzero(mask) < 2:
        return out
    lx = np.log(xsrc[mask])
    ly = np.log(ysrc[mask])
    inside = (xdst >= xsrc[mask][0]) & (xdst <= xsrc[mask][-1])
    if np.any(inside):
        out[inside] = np.exp(np.interp(np.log(xdst[inside]), lx, ly))
    return out


def _interp_distribution_per_gev(
    energy_src_gev: np.ndarray,
    density_src_per_gev: np.ndarray,
    energy_dst_gev: np.ndarray,
) -> np.ndarray:
    return _interp_positive_loglog(energy_src_gev, density_src_per_gev, energy_dst_gev)


def _interp_source_per_gamma(
    energy_src_gev: np.ndarray,
    source_src_per_gev_s: np.ndarray,
    energy_dst_gev: np.ndarray,
    mass_gev: float,
    shell_volume_cm3: float,
) -> np.ndarray:
    source_per_gev_s = _interp_positive_loglog(energy_src_gev, source_src_per_gev_s, energy_dst_gev)
    return shell_volume_cm3 * source_per_gev_s * mass_gev


def _hadronic_continuous_loss_rates(
    gam_p: np.ndarray, b_field_g: float, t_dyn_s: float,
    quantum_syn: bool = False, mass_gev: float | None = None,
) -> np.ndarray:
    coeff_syn = constants.para_sigmat * b_field_g * b_field_g / (6.0 * np.pi * constants.para_m_e * constants.para_c) / (constants.para_m_p_div_m_e ** 3)
    loss_ad = gam_p / max(float(t_dyn_s), 1.0)
    loss_syn = coeff_syn * gam_p * gam_p
    if quantum_syn:
        loss_syn = loss_syn * _quantum_syn_cooling_factor(gam_p, b_field_g, mass_gev or constants.para_m_p_gev)
    return loss_ad + loss_syn


def _hadronic_advance_energy_loggamma(
    gam_p: np.ndarray,
    gam_edge: np.ndarray,
    d_n_prev: np.ndarray,
    q_inj: np.ndarray,
    loss_total: np.ndarray,
    dt_s: float,
) -> np.ndarray:
    gam = np.asarray(gam_p, dtype=float)
    dgam = gam_edge[1:] - gam_edge[:-1]
    if np.any(dgam <= 0.0):
        raise ValueError("gam_edge must be strictly increasing.")
    loss = np.asarray(loss_total, dtype=float)
    if not (
        np.all(np.isfinite(gam))
        and np.all(np.isfinite(d_n_prev))
        and np.all(np.isfinite(q_inj))
        and np.all(np.isfinite(loss))
    ):
        raise ValueError("Hadronic energy advance received non-finite inputs.")

    content = np.asarray(d_n_prev, dtype=float) * dgam
    cooled_gamma = gam - loss * float(dt_s)
    target = np.searchsorted(gam_edge, cooled_gamma, side="right") - 1
    next_content = np.zeros_like(content)
    valid = (target >= 0) & (target < gam.size)
    np.add.at(next_content, target[valid], content[valid])
    return next_content / dgam + np.asarray(q_inj, dtype=float)


def _energy_luminosity_from_rate_spectrum(
    energy_gev: np.ndarray,
    spectrum: np.ndarray,
    shell_volume_cm3: float,
) -> np.ndarray:
    energy = np.asarray(energy_gev, dtype=float)
    spec = np.asarray(spectrum, dtype=float)
    return shell_volume_cm3 * spec * energy * constants.para_h_gev * GEV_TO_ERG


def _hadronic_proton_syn_state(
    radius_cm: float,
    b_field_g: float,
    gam_p: np.ndarray,
    d_n_gam_p: np.ndarray,
    v_seed_hz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    temp_syn = np.sqrt(3.0) * constants.para_e**3 / constants.para_m_p_e
    p_had_syn = np.zeros_like(v_seed_hz, dtype=float)
    log_gam_p = np.log(gam_p)
    for i_nu, nu in enumerate(v_seed_hz):
        fx = _hadronic_syn_kernel_ultrarel(float(nu), gam_p, b_field_g)
        p_had_syn[i_nu] = temp_syn * b_field_g * float(trapezoid(d_n_gam_p * fx * gam_p, log_gam_p))
    seed_had_syn = p_had_syn / (max(radius_cm * radius_cm, 1.0e-60) * np.maximum(v_seed_hz, 1.0e-60))
    seed_had_syn /= 4.0 * np.pi * constants.para_c * constants.para_h
    return p_had_syn, seed_had_syn


def _hadronic_syn_kernel_ultrarel(
    nu_hz: float,
    gam_p: np.ndarray,
    b_field_g: float,
) -> np.ndarray:
    b_crit = 4.41e13
    mass_p_gev = constants.para_m_p_e * constants.para_erg2ev * 1.0e-9
    mass_e_gev = constants.para_m_energy * constants.para_erg2ev * 1.0e-9
    energy_photon_gev = constants.para_h_gev * nu_hz
    energy_proton_gev = gam_p * mass_p_gev
    b_dimless = b_field_g / b_crit
    mass_ratio = mass_e_gev / mass_p_gev
    xbar = energy_photon_gev * mass_p_gev / (
        3.0 * energy_proton_gev * energy_proton_gev * b_dimless * mass_ratio * mass_ratio
    )
    x_arg = 2.0 * xbar
    out = np.zeros_like(gam_p, dtype=float)

    m1 = x_arg < 1.0e-2
    out[m1] = 1.80842 * xbar[m1] ** (1.0 / 3.0) * 2.0 ** (-2.0 / 3.0)

    m2 = (x_arg >= 1.0e-2) & (x_arg < 1.0)
    if np.any(m2):
        y = np.log10(x_arg[m2])
        poly = (
            -0.35775237
            - 0.83695385 * y
            - 1.1449608 * y * y
            - 0.68137283 * y**3
            - 0.22754737 * y**4
            - 0.031967334 * y**5
        )
        out[m2] = 10.0**poly / 2.0

    m3 = (x_arg >= 1.0) & (x_arg < 10.0)
    if np.any(m3):
        y = np.log10(x_arg[m3])
        poly = (
            -0.35842494
            - 0.79652041 * y
            - 1.6113032 * y * y
            + 0.26055213 * y**3
            - 1.6979017 * y**4
            + 0.032955035 * y**5
        )
        out[m3] = 10.0**poly / 2.0

    m4 = (x_arg >= 10.0) & (x_arg < 1.0e2)
    if np.any(m4):
        out[m4] = (np.pi / 4.0) * np.exp(-x_arg[m4]) * (1.0 - 99.0 / (162.0 * x_arg[m4]))
    return out


def _hadronic_shell_injection_energy(
    radius_cm: np.ndarray,
    gamma_bulk: np.ndarray,
    config: FitConfig,
) -> np.ndarray:
    radius = np.asarray(radius_cm, dtype=float)
    gamma = np.asarray(gamma_bulk, dtype=float)
    if radius.ndim != 1 or gamma.ndim != 1 or radius.shape != gamma.shape:
        raise ValueError("radius_cm and gamma_bulk must be 1d arrays with matching shapes.")

    prev_radius = np.empty_like(radius)
    prev_radius[0] = 0.0
    prev_radius[1:] = radius[:-1]
    shell_volume_cm3 = (4.0 / 3.0) * np.pi * np.maximum(radius**3 - prev_radius**3, 0.0)
    density_cm3 = np.asarray(ambient_density(radius, config), dtype=float)
    return (
        float(config.hadronic.epsilon_p)
        * np.maximum(gamma - 1.0, 0.0)
        * shell_volume_cm3
        * density_cm3
        * constants.para_m_p
        * constants.para_c**2
    )


def _compute_forward_timescales(
    r_gamma: np.ndarray,
    radius_cm: np.ndarray,
    nu_c: np.ndarray,
    config: FitConfig,
) -> tuple[np.ndarray, np.ndarray]:
    gamma = np.asarray(r_gamma, dtype=float)
    radius = np.asarray(radius_cm, dtype=float)
    nu_c_arr = np.asarray(nu_c, dtype=float)
    magnetic_field_g = np.asarray(compute_magnetic_field(gamma, radius, config), dtype=float)
    doppler_den = np.asarray(doppler_denominator(gamma, config.z), dtype=float)
    beta = np.zeros_like(gamma)
    valid_gamma = gamma > 1.0
    beta[valid_gamma] = np.sqrt(1.0 - gamma[valid_gamma] ** (-2.0))

    gamma_c = np.zeros_like(nu_c_arr)
    valid = (magnetic_field_g > 0.0) & (doppler_den > 0.0) & (nu_c_arr > 0.0)
    gamma_c[valid] = np.sqrt(nu_c_arr[valid] * doppler_den[valid] / (4.2e6 * magnetic_field_g[valid]))

    cooling_timescale_s = np.zeros_like(nu_c_arr)
    valid_cooling = valid & (gamma > 0.0)
    cooling_timescale_s[valid_cooling] = (
        7.7e8
        * (1.0 + float(config.z))
        / (gamma[valid_cooling] * magnetic_field_g[valid_cooling] ** 2 * gamma_c[valid_cooling])
    )

    dynamical_timescale_s = np.zeros_like(radius)
    valid_dyn = beta > 0.0
    dynamical_timescale_s[valid_dyn] = radius[valid_dyn] / (gamma[valid_dyn] * beta[valid_dyn] * constants.para_c)
    return cooling_timescale_s, dynamical_timescale_s


def solve_reverse_shock_emission(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: FitConfig,
) -> ReverseShockEmission | None:
    reverse_params = _resolve_reverse_shock_parameters(config)
    if reverse_params is None or dynamics.reverse_shock is None:
        return None

    gam_e, d_n_gam_e = _solve_reverse_shock_electrons(boundary, dynamics, v_seed, config, reverse_params)
    dynamics.reverse_shock.gam_e = gam_e
    dynamics.reverse_shock.d_n_gam_e = _renormalize_reverse_shock_distribution(
        gam_e,
        d_n_gam_e,
        dynamics.reverse_shock.swept_mass_g,
        reverse_params.f_e,
    )

    (
        nu_m,
        nu_c,
        nu_a,
        magnetic_field_g,
        nu_M,
    ) = _compute_reverse_shock_characteristic_frequencies(
        config,
        reverse_params,
        dynamics,
    )
    l_syn_spec, seed_syn = _compute_reverse_shock_synchrotron_emission(dynamics, v_seed, config)

    rs_hadronic = None
    if bool(config.hadronic.reverse_enabled) and float(config.hadronic.reverse_epsilon_p) > 0.0:
        from asgard_core.hadronic_reverse import solve_rs_hadronic_core
        rs_hadronic = solve_rs_hadronic_core(
            r_tobs_s=dynamics.r_tobs,
            r_gamma=dynamics.r_gamma,
            radius_cm=dynamics.radius,
            rs_swept_mass_g=dynamics.reverse_shock.swept_mass_g,
            rs_b_field_g=dynamics.reverse_shock.magnetic_field_g,
            v_seed_hz=v_seed,
            num_gam_p=config.hadronic.num_gam_p,
            epsilon_p=float(config.hadronic.reverse_epsilon_p),
            include_proton_synch=bool(config.hadronic.include_proton_synch),
        )

    return ReverseShockEmission(
        l_syn_spec=l_syn_spec,
        seed_syn=seed_syn,
        magnetic_field_g=magnetic_field_g,
        nu_m=nu_m,
        nu_c=nu_c,
        nu_a=nu_a,
        nu_M=nu_M,
        rs_hadronic=rs_hadronic,
    )


def _solve_reverse_shock_electrons(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: FitConfig,
    reverse_params: ReverseShockParameters,
) -> tuple[np.ndarray, np.ndarray]:
    if dynamics.reverse_shock is None:
        raise ValueError("Reverse shock dynamics are required to compute reverse electrons.")
    module = _electron_reverse_module().electron_reverse_kernel
    delta_0 = reverse_params.delta_t_s * constants.para_c
    para_m_ej = config.e_iso / config.eta_0 / constants.para_c**2
    gam_e, d_n_gam_e = module.electron_reverse_evolve(
        delta_0,
        reverse_params.epsilon_e,
        reverse_params.epsilon_b,
        reverse_params.p,
        reverse_params.f_e,
        config.eta_0,
        config.epsilon_e,
        config.epsilon_b,
        config.z,
        boundary[11],
        boundary[10],
        para_m_ej,
        dynamics.reverse_shock.t_cross,
        dynamics.reverse_shock.r_cross,
        dynamics.r_tobs,
        dynamics.r_gamma,
        dynamics.radius,
        dynamics.reverse_shock.magnetic_field_g,
        v_seed,
        config.num_gam_e,
        config.index_y,
        config.index_syn_integr,
        config.num_threads,
    )
    return np.asarray(gam_e, dtype=float), np.asarray(d_n_gam_e, dtype=float)


def _compute_reverse_shock_synchrotron_emission(
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: FitConfig,
) -> tuple[np.ndarray, np.ndarray]:
    if dynamics.reverse_shock is None:
        raise ValueError("Reverse shock dynamics are required to compute reverse emission.")
    gam_e = dynamics.reverse_shock.gam_e
    d_n_gam_e = dynamics.reverse_shock.d_n_gam_e
    if gam_e is None or d_n_gam_e is None:
        raise ValueError("Reverse shock electrons are required to compute reverse emission.")
    num_nu = v_seed.shape[0]
    num_r = dynamics.radius.shape[0]
    l_syn_spec = np.zeros((num_nu, num_r), dtype=float)
    seed_syn = np.zeros((num_nu, num_r), dtype=float)
    magnetic_field_g = np.asarray(dynamics.reverse_shock.magnetic_field_g, dtype=float)

    for i in range(num_r):
        db = float(magnetic_field_g[i])
        radius_loc = float(dynamics.radius[i])
        if not np.isfinite(db) or not np.isfinite(radius_loc) or db <= 0.0 or radius_loc <= 0.0:
            continue
        p_syn_i, seed_syn_i = electron_radiation_module.get_syn_selected(
            config.index_syn_integr,
            radius_loc,
            db,
            config.num_threads,
            gam_e,
            d_n_gam_e[:, i],
            v_seed,
        )
        l_syn_spec[:, i] = np.asarray(p_syn_i, dtype=float)
        seed_syn[:, i] = np.asarray(seed_syn_i, dtype=float)
    return l_syn_spec, seed_syn


def _resolve_reverse_shock_parameters(config: FitConfig) -> ReverseShockParameters | None:
    reverse_enabled = config.reverse or config.reverse_shock.enabled
    if not reverse_enabled:
        return None

    if config.reverse_shock.delta_t_s is None:
        raise ValueError("ReverseShockConfig.delta_t_s must be set when reverse shock is enabled.")

    return ReverseShockParameters(
        delta_t_s=config.reverse_shock.delta_t_s,
        epsilon_e=config.epsilon_e if config.reverse_shock.epsilon_e is None else config.reverse_shock.epsilon_e,
        epsilon_b=config.epsilon_b if config.reverse_shock.epsilon_b is None else config.reverse_shock.epsilon_b,
        p=config.p if config.reverse_shock.p is None else config.reverse_shock.p,
        f_e=config.f_e if config.reverse_shock.f_e is None else config.reverse_shock.f_e,
    )


def _compute_characteristic_frequencies_weno5(
    config: FitConfig,
    r_tobs: np.ndarray,
    r_gamma: np.ndarray,
    radius: np.ndarray,
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    num_r = radius.shape[0]
    nu_m = np.zeros(num_r, dtype=float)
    nu_c = np.zeros(num_r, dtype=float)
    nu_a = np.zeros(num_r, dtype=float)

    for i in range(1, num_r):
        radius_loc = radius[i - 1]
        gamma_loc = 0.5 * (r_gamma[i - 1] + r_gamma[i])
        d_ne = _ambient_density(radius_loc, config)
        db = 0.39 * np.sqrt(config.epsilon_b * d_ne * (gamma_loc * (gamma_loc - 1.0)))
        gam_e_max = 3.0 * constants.para_m_energy / np.sqrt(8.0 * db * constants.para_e**3)
        gam_e_m = _minimum_electron_lorentz_factor(config, gamma_loc, gam_e_max)
        gam_e_c = 7.7e8 * (1.0 + config.z) / gamma_loc / db**2 / r_tobs[i]
        doppler_den = _doppler_denominator(gamma_loc, config.z)

        nu_m[i - 1] = _synchrotron_frequency(db, gam_e_m, doppler_den)
        nu_c[i - 1] = _synchrotron_frequency(db, gam_e_c, doppler_den)

        nu_a_comoving = electron_radiation_module.get_nu_a(radius_loc, db, gam_e, d_n_gam_e[:, i - 1])
        nu_a[i - 1] = nu_a_comoving / doppler_den

    return nu_m, nu_c, nu_a


def _compute_reverse_shock_characteristic_frequencies(
    config: FitConfig,
    reverse_params: ReverseShockParameters,
    dynamics: DynamicsSolution,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    num_r = dynamics.radius.shape[0]
    nu_m = np.zeros(num_r, dtype=float)
    nu_c = np.zeros(num_r, dtype=float)
    nu_a = np.zeros(num_r, dtype=float)
    magnetic_field_g = np.zeros(num_r, dtype=float)
    nu_M = np.zeros(num_r, dtype=float)
    gam_e = dynamics.reverse_shock.gam_e
    d_n_gam_e = dynamics.reverse_shock.d_n_gam_e
    if gam_e is None or d_n_gam_e is None:
        raise ValueError("Reverse shock electrons are required to compute reverse characteristic frequencies.")

    eta_0 = config.eta_0
    beta4 = np.sqrt(1.0 - eta_0**-2)
    gamma_floor = float(gam_e[0])
    for i in range(1, num_r):
        radius_loc = dynamics.radius[i - 1]
        gamma2 = 0.5 * (dynamics.r_gamma[i - 1] + dynamics.r_gamma[i])
        beta2 = np.sqrt(1.0 - gamma2**-2)
        gamma34 = (1.0 - beta2 * beta4) * eta_0 * gamma2

        db = float(dynamics.reverse_shock.magnetic_field_g[i - 1])
        magnetic_field_g[i - 1] = db
        gam_e_max = 3.0 * constants.para_m_energy / np.sqrt(8.0 * db * constants.para_e**3)
        gam_e_m = _minimum_reverse_shock_electron_lorentz_factor(reverse_params, gamma34, gam_e_max)
        gam_e_m = max(gam_e_m, gamma_floor)
        gam_e_c = 7.7e8 * (1.0 + config.z) / gamma2 / db**2 / dynamics.r_tobs[i]
        doppler_den = _doppler_denominator(gamma2, config.z)

        nu_m[i - 1] = _synchrotron_frequency(db, gam_e_m, doppler_den)
        nu_c[i - 1] = _synchrotron_frequency(db, gam_e_c, doppler_den)
        nu_M[i - 1] = _synchrotron_frequency(db, gam_e_max, doppler_den)

        nu_a_comoving = electron_radiation_module.get_nu_a(
            radius_loc,
            db,
            gam_e,
            d_n_gam_e[:, i - 1],
        )
        nu_a[i - 1] = nu_a_comoving / doppler_den

    return nu_m, nu_c, nu_a, magnetic_field_g, nu_M


def _compute_synchrotron_emission_from_distribution(
    radius_cm: np.ndarray,
    magnetic_field_g: np.ndarray,
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
    v_seed: np.ndarray,
    config: FitConfig,
) -> tuple[np.ndarray, np.ndarray]:
    num_nu = v_seed.shape[0]
    num_r = radius_cm.shape[0]
    l_syn_spec = np.zeros((num_nu, num_r), dtype=float)
    seed_syn = np.zeros((num_nu, num_r), dtype=float)

    for i in range(num_r):
        db = float(magnetic_field_g[i])
        radius_loc = float(radius_cm[i])
        if not np.isfinite(db) or not np.isfinite(radius_loc) or db <= 0.0 or radius_loc <= 0.0:
            continue
        p_syn_i, seed_syn_i = electron_radiation_module.get_syn_selected(
            config.index_syn_integr,
            radius_loc,
            db,
            config.num_threads,
            gam_e,
            d_n_gam_e[:, i],
            v_seed,
        )
        l_syn_spec[:, i] = np.asarray(p_syn_i, dtype=float)
        seed_syn[:, i] = np.asarray(seed_syn_i, dtype=float)

    return l_syn_spec, seed_syn




def _renormalize_reverse_shock_distribution(
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
    swept_mass_g: np.ndarray,
    f_e: float,
) -> np.ndarray:
    gam = np.asarray(gam_e, dtype=float)
    dist = np.asarray(d_n_gam_e, dtype=float).copy()
    targets = np.asarray(swept_mass_g, dtype=float) / constants.para_m_p * float(f_e)
    for i in range(dist.shape[1]):
        total = float(trapezoid(dist[:, i], gam))
        target = float(targets[i])
        if not np.isfinite(total) or total <= 0.0:
            continue
        if not np.isfinite(target) or target <= 0.0:
            continue
        dist[:, i] *= target / total
    return dist


def _ambient_density(radius_cm: float, config: FitConfig) -> float:
    """DEPRECATED: Use asgard_physics_utils.ambient_density instead."""
    return ambient_density(radius_cm, config)


def _doppler_denominator(gamma_bulk: float, redshift: float) -> float:
    """DEPRECATED: Use asgard_physics_utils.doppler_denominator instead."""
    return doppler_denominator(gamma_bulk, redshift)


def _synchrotron_frequency(magnetic_field_g: float, electron_lorentz_factor: float, doppler_den: float) -> float:
    return 4.2e6 * magnetic_field_g * electron_lorentz_factor * electron_lorentz_factor / doppler_den


def _reverse_ambient_density(radius_cm: float, config: FitConfig) -> float:
    if config.a_star >= 0.0:
        d_ne_wind = config.a_star * 3.0e35 / radius_cm**2
        if d_ne_wind <= config.d_ne / 4.0:
            return config.d_ne
        return d_ne_wind
    return config.d_ne


def _minimum_electron_lorentz_factor(config: FitConfig, gamma_bulk: float, gam_e_max: float) -> float:
    temp_gam = config.epsilon_e / config.f_e * constants.para_m_p_div_m_e * (gamma_bulk - 1.0)

    if config.p > 2.0:
        return (config.p - 2.0) / (config.p - 1.0) * temp_gam + 1.0

    if 1.0 < config.p < 2.0:
        return ((2.0 - config.p) / (config.p - 1.0) * temp_gam * gam_e_max ** (config.p - 2.0)) ** (
            1.0 / (config.p - 1.0)
        ) + 1.0

    if np.isclose(config.p, 2.0):
        gam_e_m = 1.0
        temp = temp_gam / np.log(gam_e_max / gam_e_m)
        while abs(1.0 - gam_e_m / temp) > 1.0e-5:
            temp = temp_gam / np.log(gam_e_max / gam_e_m)
            if gam_e_m > temp:
                gam_e_m = 0.5 * (gam_e_m + temp)
            else:
                gam_e_m = 0.5 * (gam_e_m + gam_e_max)
        return gam_e_m

    raise ValueError(f"Unsupported electron index p={config.p}.")


def _minimum_reverse_shock_electron_lorentz_factor(
    reverse_params: ReverseShockParameters,
    gamma34: float,
    gam_e_max: float,
) -> float:
    temp_gam = reverse_params.epsilon_e / reverse_params.f_e * constants.para_m_p_div_m_e * (gamma34 - 1.0)

    if reverse_params.p > 2.05:
        return (reverse_params.p - 2.0) / (reverse_params.p - 1.0) * temp_gam + 1.0

    if 1.0 < reverse_params.p < 2.0:
        return ((2.0 - reverse_params.p) / (reverse_params.p - 1.0) * temp_gam * gam_e_max ** (reverse_params.p - 2.0)) ** (
            1.0 / (reverse_params.p - 1.0)
        ) + 1.0

    if reverse_params.p >= 2.0:
        return 0.05 / 1.05 * temp_gam + 1.0

    raise ValueError(f"Unsupported reverse-shock electron index p={reverse_params.p}.")
