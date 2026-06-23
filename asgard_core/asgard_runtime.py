from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from importlib import import_module
import time
import numpy as np

from src.Electron.electron_radiation import electron_radiation_kernel as electron_radiation_module
import src.Hadronic.hadronic_forward_1d as hadronic_legacy_module
import src.Hadronic.hadronic_reverse_1d as _hadronic_reverse_module
from asgard_core.hadronic_am3_solver import (
    HUMMER_PROCESS_GROUP_LABELS,
    HUMMER2010_RESPONSE_BACKEND,
    solve_hummer_2010_response_processes,
)
from asgard_core.hadronic_processes import ACCELERATION_BACKEND
from asgard_core.hadronic_processes import BETHE_HEITLER_BACKEND
from asgard_core.hadronic_processes import HADRONIC_IC_BACKEND
from asgard_core.hadronic_processes import SECONDARY_RADIATION_BACKEND
from asgard_core.hadronic_processes import SPECIES_TRANSPORT_BACKEND
from asgard_core.hadronic_processes import (
    HUMMER2010_DECAY_BACKEND,
    HUMMER2010_OPERATOR_BACKEND,
)
from asgard_core.hadronic_processes import PP_DELTA_BACKEND
from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_types import (
    ReverseShockParameters,
    ReverseShockCausalityDiagnostics,
    ReverseShockDynamics,
    DynamicsSolution,
    ElectronSolution,
    HadronicSolution,
    ReverseShockEmission,
    SolverAdapterReport,
)
from asgard_core.asgard_physics_utils import (
    ambient_density,
    compute_magnetic_field,
    density_jump_arrays,
    reverse_shell_baryonic_mass,
)
from src import Dynamics, constants


_ELECTRON_MODULES = {
    "charint_1d": "src.Electron.electron_forward_charint_1d",
    "charint_2d": "src.Electron.electron_forward_transport_2d",
    "dg_1d": "src.Electron.electron_forward_dg_1d",
    "fullhide_1d": "src.Electron.electron_forward_fullhide_1d",
    "fullhide_1d_hz": "src.Electron.electron_forward_fullhide_1d_hybrid",
    "fullhide_2d": "src.Electron.electron_forward_transport_2d",
    "fullhide_2d_pic": "src.Electron.electron_forward_transport_2d_pic",
    "slc1_1d": "src.Electron.electron_forward_slc1_1d",
    "t2g1_1d": "src.Electron.electron_forward_t2g1_1d",
    "weno5_1d": "src.Electron.electron_forward_weno5_1d",
}

ELECTRON_1D_TRANSPORT_IDS = {
    "fullhide_1d": 1,
    "fullhide_1d_hz": 1,
    "dg_1d": 2,
}

_ELECTRON_1D_NU_SOLVERS = {"t2g1_1d", "slc1_1d", "charint_1d", "dg_1d"}
_ELECTRON_1D_NU_GRID = {"dg_1d": "log-four-velocity-1d-dg"}


@cache
def _electron_module(solver: str):
    return import_module(_ELECTRON_MODULES[solver])


@cache
def _electron_reverse_module():
    return import_module("src.Electron.electron_reverse_kernel")


@dataclass
class SecondaryReverseShockState:
    luminosity_syn: np.ndarray
    branch_luminosity_syn: np.ndarray
    gam_e: np.ndarray
    d_n_gam_e: np.ndarray
    event_active: np.ndarray
    start_radius_cm: np.ndarray
    shock_end_radius_cm: np.ndarray
    start_tobs_axis_s: np.ndarray
    shock_end_tobs_axis_s: np.ndarray
    gamma_contact: np.ndarray
    pressure_3: np.ndarray
    gamma_43: np.ndarray
    beta_rs: np.ndarray
    dissipated_energy_density: np.ndarray
    dissipated_energy_erg: np.ndarray
    electron_injected_energy_erg: np.ndarray
    swept_mass_g: np.ndarray
    internal_energy_erg: np.ndarray
    comoving_volume_cm3: np.ndarray
    pressure_total: np.ndarray
    enthalpy_density_total: np.ndarray
    branch_swept_mass_g: np.ndarray
    branch_internal_energy_erg: np.ndarray
    branch_comoving_volume_cm3: np.ndarray
    branch_magnetic_field_g: np.ndarray
    branch_gamma_m: np.ndarray
    branch_gamma_contact: np.ndarray
    branch_gamma_43: np.ndarray
    branch_compression: np.ndarray
    branch_beta_rs: np.ndarray
    branch_dissipated_energy_density: np.ndarray
    branch_reacceleration_seed_energy_erg: np.ndarray
    branch_reaccelerated_energy_erg: np.ndarray
    magnetic_field_g: np.ndarray


@dataclass(frozen=True)
class ReverseShockHadronicSolution:
    gam_p: np.ndarray
    d_n_gam_p: np.ndarray
    l_had_syn_spec: np.ndarray
    seed_had_syn: np.ndarray


_PGAMMA_SCHEME_DISABLED = "disabled"
_PGAMMA_SCHEME_HUMMER2010_RESPONSE = "hummer_2010_response"


def _rs_vegas_downstream_u(gamma_rel: float, sigma: float) -> float:
    sigma_cut = 1.0e-6
    gamma_eff = max(1.0, float(gamma_rel))
    ad = 4.0 / 3.0 + 1.0 / (3.0 * gamma_eff)
    gamma_sq = gamma_eff * gamma_eff
    gm1 = gamma_eff - 1.0
    gp1 = gamma_eff + 1.0
    hm1 = ad - 1.0
    hm2 = ad - 2.0
    term1 = -ad * hm2
    term2 = gamma_sq - 1.0
    if sigma <= sigma_cut:
        u2_down = gm1 * hm1 * hm1 / (term1 * gm1 + 2.0)
    else:
        a_cub = term1 * gm1 + 2.0
        b_cub = -gp1 * (-hm2 * (ad * gamma_sq + 1.0) + ad * hm1 * gamma_eff) * sigma
        b_cub = b_cub - gm1 * (term1 * (gamma_sq - 2.0) + 2.0 * gamma_eff + 3.0)
        c_cub = gp1 * (ad * (1.0 - ad / 4.0) * term2 + 1.0) * sigma * sigma
        c_cub = c_cub + term2 * (2.0 * gamma_eff + hm2 * (ad * gamma_eff - 1.0)) * sigma
        c_cub = c_cub + gp1 * gm1 * gm1 * hm1 * hm1
        d_cub = -gm1 * gp1 * gp1 * hm2 * hm2 * sigma * sigma / 4.0
        b_red = b_cub / a_cub
        c_red = c_cub / a_cub
        d_red = d_cub / a_cub
        p_cub = c_red - b_red * b_red / 3.0
        q_cub = 2.0 * b_red * b_red * b_red / 27.0 - b_red * c_red / 3.0 + d_red
        amp = np.sqrt(-p_cub / 3.0)
        arg = 3.0 * q_cub / (2.0 * p_cub * amp)
        arg = max(-1.0, min(1.0, float(arg)))
        u2_down = 2.0 * amp * np.cos((np.arccos(arg) - 2.0 * np.pi) / 3.0) - b_red / 3.0
    return float(np.sqrt(u2_down))


def _rs_shock_upstream_u(gamma_rel: float, sigma: float) -> float:
    u_down = _rs_vegas_downstream_u(gamma_rel, sigma)
    gamma_eff = max(1.0, float(gamma_rel))
    gamma_sq_minus_one = (gamma_eff - 1.0) * (gamma_eff + 1.0)
    return float(np.sqrt((1.0 + u_down * u_down) * gamma_sq_minus_one) + u_down * gamma_eff)


def _rs_fast_shock_allowed(gamma_rel: float, sigma: float) -> bool:
    if sigma <= 0.0:
        return True
    u_up = _rs_shock_upstream_u(gamma_rel, sigma)
    return bool(u_up * u_up > sigma)


def _rs_fast_wave_relative_beta(gamma4: float, sigma: float) -> float:
    if sigma <= 0.0:
        return np.inf
    beta4 = np.sqrt(1.0 - gamma4**-2)
    beta_fast = np.sqrt(sigma / (1.0 + sigma))
    return float(beta_fast / (gamma4 * gamma4 * (1.0 - beta4 * beta_fast)))


def _rs_shell_width(radius_cm: np.ndarray, gamma4: float, delta0_cm: float) -> np.ndarray:
    return np.maximum(delta0_cm, radius_cm / (gamma4 * gamma4))


def _rs_fast_wave_depth(radius_cm: np.ndarray, gamma4: float, sigma: float) -> np.ndarray:
    if sigma <= 0.0:
        return np.zeros_like(radius_cm, dtype=float)
    beta4 = np.sqrt(1.0 - gamma4**-2)
    return radius_cm * _rs_fast_wave_relative_beta(gamma4, sigma) / beta4


def _rs_shell_contact_fraction(radius_cm: np.ndarray, gamma4: float, sigma: float, shell_width_cm: np.ndarray) -> np.ndarray:
    if sigma <= 0.0:
        return np.ones_like(radius_cm, dtype=float)
    depth_cm = _rs_fast_wave_depth(radius_cm, gamma4, sigma)
    return np.where(depth_cm >= shell_width_cm, 1.0, depth_cm / shell_width_cm)


def _rs_relative_gamma_same_direction(gamma_contact: np.ndarray, gamma4: float) -> np.ndarray:
    u_contact = np.sqrt((gamma_contact - 1.0) * (gamma_contact + 1.0))
    u4 = np.sqrt((gamma4 - 1.0) * (gamma4 + 1.0))
    return (gamma_contact * gamma_contact + gamma4 * gamma4 - 1.0) / (gamma4 * gamma_contact + u_contact * u4)


def _rs_pressure_balance_sigma_cr(gamma_contact: np.ndarray, n1: np.ndarray, n4: np.ndarray) -> np.ndarray:
    return 2.0 * (4.0 * gamma_contact + 3.0) * (gamma_contact - 1.0) * n1 / (3.0 * n4)


def _rs_pressure_balance_ratio(
    gamma_contact: np.ndarray,
    n1: np.ndarray,
    n4: np.ndarray,
    sigma: float,
) -> np.ndarray:
    if sigma <= 0.0:
        return np.full_like(gamma_contact, np.inf, dtype=float)
    return _rs_pressure_balance_sigma_cr(gamma_contact, n1, n4) / sigma


def _reverse_shock_causality_diagnostics(
    config: RuntimeConfig,
    reverse_params: ReverseShockParameters,
    dynamics: DynamicsSolution,
    reverse_shock: ReverseShockDynamics,
) -> ReverseShockCausalityDiagnostics:
    sigma = float(reverse_params.sigma)
    gamma0 = float(config.eta_0)
    delta0_cm = float(reverse_params.delta_t_s) * constants.para_c
    spreading_radius_cm = gamma0 * gamma0 * delta0_cm

    if sigma <= 0.0:
        contact_radius_cm = np.inf
        initial_super_fast = True
    else:
        contact_radius_cm = delta0_cm * gamma0 * gamma0 * (np.sqrt((1.0 + sigma) / sigma) - 1.0)
        initial_super_fast = bool(gamma0 > np.sqrt(1.0 + sigma))

    if float(config.a_star) > 0.0:
        medium = "wind"
        wind_mass_per_radius = float(config.a_star) * 3.0e35 * constants.para_m_p
        deceleration_radius_cm = float(config.e_iso) / (
            4.0 * np.pi * wind_mass_per_radius * constants.para_c**2 * gamma0 * gamma0
        )
        reference_crossing_radius_cm = np.sqrt(spreading_radius_cm * deceleration_radius_cm) / (1.0 + sigma)
    else:
        medium = "ism"
        sedov_length_cm = (3.0 * float(config.e_iso) / (4.0 * np.pi * float(config.d_ne) * constants.para_m_p * constants.para_c**2)) ** (1.0 / 3.0)
        deceleration_radius_cm = sedov_length_cm / gamma0 ** (2.0 / 3.0)
        reference_crossing_radius_cm = (spreading_radius_cm * deceleration_radius_cm**3) ** 0.25 / np.sqrt(1.0 + sigma)

    radius = np.asarray(dynamics.radius, dtype=float)
    tobs = np.asarray(dynamics.r_tobs, dtype=float)
    gamma_contact = np.asarray(dynamics.r_gamma, dtype=float)
    shell_width_cm = _rs_shell_width(radius, gamma0, delta0_cm)
    n1 = np.asarray(ambient_density(radius, config), dtype=float)
    baryonic_mass_g = reverse_shell_baryonic_mass(config)
    n4 = baryonic_mass_g / (4.0 * np.pi * constants.para_m_p * radius * radius * gamma0 * shell_width_cm)
    pressure_ratio = _rs_pressure_balance_ratio(gamma_contact, n1, n4, sigma)
    pressure_allowed = pressure_ratio >= 1.0
    contact_fraction = _rs_shell_contact_fraction(radius, gamma0, sigma, shell_width_cm)
    if np.any(pressure_allowed):
        pressure_index = int(np.flatnonzero(pressure_allowed)[0])
        pressure_radius_cm = float(radius[pressure_index])
        pressure_tobs_s = float(tobs[pressure_index])
        pressure_start_ratio = float(pressure_ratio[pressure_index])
        pressure_start_contact_fraction = float(contact_fraction[pressure_index])
        pressure_seen = True
    else:
        pressure_index = -1
        pressure_radius_cm = np.inf
        pressure_tobs_s = np.inf
        pressure_start_ratio = 0.0
        pressure_start_contact_fraction = 0.0
        pressure_seen = False

    fast_wave_depth_cm = _rs_fast_wave_depth(radius, gamma0, sigma)
    fast_crossed = fast_wave_depth_cm >= shell_width_cm
    if np.any(fast_crossed):
        fast_index = int(np.flatnonzero(fast_crossed)[0])
        fast_radius_cm = float(radius[fast_index])
        fast_tobs_s = float(tobs[fast_index])
    else:
        fast_index = -1
        fast_radius_cm = np.inf
        fast_tobs_s = np.inf

    gamma34 = _rs_relative_gamma_same_direction(gamma_contact, gamma0)
    local_allowed = np.array([_rs_fast_shock_allowed(float(gamma_rel), sigma) for gamma_rel in gamma34], dtype=bool)
    if np.any(local_allowed):
        local_start_index = int(np.flatnonzero(local_allowed)[0])
        local_start_radius_cm = float(radius[local_start_index])
        local_start_tobs_s = float(tobs[local_start_index])
        local_seen = True
    else:
        local_start_index = -1
        local_start_radius_cm = np.inf
        local_start_tobs_s = np.inf
        local_seen = False
    ignition_allowed = pressure_allowed & local_allowed
    global_allowed = bool(np.any(ignition_allowed))

    swept = np.asarray(reverse_shock.swept_mass_g, dtype=float)
    growth = swept > float(swept[0]) * (1.0 + 1.0e-6)
    if np.any(growth):
        actual_start_index = int(np.flatnonzero(growth)[0])
        actual_start_radius_cm = float(radius[actual_start_index])
        actual_start_tobs_s = float(tobs[actual_start_index])
        actual_start_pressure_ratio = float(pressure_ratio[actual_start_index])
        actual_start_contact_fraction = float(contact_fraction[actual_start_index])
        actual_started = True
    else:
        actual_start_index = -1
        actual_start_radius_cm = np.inf
        actual_start_tobs_s = np.inf
        actual_start_pressure_ratio = 0.0
        actual_start_contact_fraction = 0.0
        actual_started = False

    return ReverseShockCausalityDiagnostics(
        medium=medium,
        initial_super_fast=initial_super_fast,
        global_reverse_shock_allowed=global_allowed,
        pressure_balance_condition_seen=pressure_seen,
        local_fast_condition_seen=local_seen,
        reverse_shock_started=actual_started,
        criteria_agree=bool(global_allowed == actual_started),
        contact_radius_cm=float(contact_radius_cm),
        reference_crossing_radius_cm=float(reference_crossing_radius_cm),
        pressure_balance_start_radius_cm=pressure_radius_cm,
        pressure_balance_start_tobs_s=pressure_tobs_s,
        pressure_balance_start_index=pressure_index,
        pressure_balance_start_ratio=pressure_start_ratio,
        pressure_balance_start_contact_fraction=pressure_start_contact_fraction,
        fast_wave_crossing_radius_cm=fast_radius_cm,
        fast_wave_crossing_tobs_s=fast_tobs_s,
        fast_wave_crossing_index=fast_index,
        local_start_radius_cm=local_start_radius_cm,
        local_start_tobs_s=local_start_tobs_s,
        local_start_index=local_start_index,
        actual_start_radius_cm=actual_start_radius_cm,
        actual_start_tobs_s=actual_start_tobs_s,
        actual_start_index=actual_start_index,
        actual_start_pressure_ratio=actual_start_pressure_ratio,
        actual_start_contact_fraction=actual_start_contact_fraction,
    )


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
    config: RuntimeConfig,
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
        u3_cross_erg,
        v3_cross_cm3,
        m3_cross_g,
        gamma_m_cross,
        ordered_magnetic_cross_g,
        r_tobs,
        r_gamma,
        radius,
        swept_mass_g,
        swept_reverse_mass_g,
        magnetic_field_g,
        internal_energy_erg,
        comoving_volume_cm3,
        gamma34,
        secondary_branch_swept_mass_g,
        secondary_branch_internal_energy_erg,
        secondary_branch_comoving_volume_cm3,
        secondary_branch_magnetic_field_g,
        secondary_swept_mass_g,
        secondary_internal_energy_erg,
        secondary_comoving_volume_cm3,
        secondary_magnetic_field_g,
        secondary_pressure_total,
        secondary_enthalpy_density_total,
        secondary_gamma_contact,
        secondary_pressure_3,
        secondary_gamma_43,
        secondary_beta_rs,
        secondary_dissipated_energy_density,
        secondary_dissipated_energy_erg,
        secondary_electron_injected_energy_erg,
        secondary_branch_gamma_m,
        secondary_branch_gamma_contact,
        secondary_branch_gamma_43,
        secondary_branch_compression,
        secondary_branch_beta_rs,
        secondary_branch_dissipated_energy_density,
        _secondary_nu_m,
        _secondary_nu_c,
        secondary_event_active,
        secondary_start_radius_cm,
        secondary_shock_end_radius_cm,
        secondary_start_tobs_axis_s,
        secondary_shock_end_tobs_axis_s,
    ) = Dynamics.dynamics_reverse(
        reverse_params.delta_t_s,
        reverse_params.epsilon_e,
        reverse_params.epsilon_b,
        reverse_params.p,
        reverse_params.f_e,
        reverse_params.sigma,
        boundary,
        config.num_r,
    )
    jump_r, _, _ = density_jump_arrays(config)
    n_j = jump_r.size
    secondary_branch_swept_mass_g = secondary_branch_swept_mass_g[:n_j, :]
    secondary_branch_internal_energy_erg = secondary_branch_internal_energy_erg[:n_j, :]
    secondary_branch_comoving_volume_cm3 = secondary_branch_comoving_volume_cm3[:n_j, :]
    secondary_branch_magnetic_field_g = secondary_branch_magnetic_field_g[:n_j, :]
    secondary_branch_gamma_m = secondary_branch_gamma_m[:n_j, :]
    secondary_branch_gamma_contact = secondary_branch_gamma_contact[:n_j, :]
    secondary_branch_gamma_43 = secondary_branch_gamma_43[:n_j, :]
    secondary_branch_compression = secondary_branch_compression[:n_j, :]
    secondary_branch_beta_rs = secondary_branch_beta_rs[:n_j, :]
    secondary_branch_dissipated_energy_density = secondary_branch_dissipated_energy_density[:n_j, :]
    secondary_event_active = np.asarray(secondary_event_active, dtype=bool)[:n_j]
    secondary_start_radius_cm = secondary_start_radius_cm[:n_j]
    secondary_shock_end_radius_cm = secondary_shock_end_radius_cm[:n_j]
    secondary_start_tobs_axis_s = secondary_start_tobs_axis_s[:n_j]
    secondary_shock_end_tobs_axis_s = secondary_shock_end_tobs_axis_s[:n_j]
    reverse_shock = ReverseShockDynamics(
        t_cross,
        r_cross,
        e3_cross,
        gam20,
        u3_cross_erg,
        v3_cross_cm3,
        m3_cross_g,
        gamma_m_cross,
        ordered_magnetic_cross_g,
        swept_reverse_mass_g,
        magnetic_field_g,
        internal_energy_erg,
        comoving_volume_cm3,
        gamma34,
        secondary_branch_swept_mass_g,
        secondary_branch_internal_energy_erg,
        secondary_branch_comoving_volume_cm3,
        secondary_branch_magnetic_field_g,
        secondary_swept_mass_g,
        secondary_internal_energy_erg,
        secondary_comoving_volume_cm3,
        secondary_magnetic_field_g,
        secondary_pressure_total,
        secondary_enthalpy_density_total,
        secondary_gamma_contact,
        secondary_pressure_3,
        secondary_gamma_43,
        secondary_beta_rs,
        secondary_dissipated_energy_density,
        secondary_dissipated_energy_erg,
        secondary_electron_injected_energy_erg,
        secondary_branch_gamma_m,
        secondary_branch_gamma_contact,
        secondary_branch_gamma_43,
        secondary_branch_compression,
        secondary_branch_beta_rs,
        secondary_branch_dissipated_energy_density,
        secondary_event_active,
        secondary_start_radius_cm,
        secondary_shock_end_radius_cm,
        secondary_start_tobs_axis_s,
        secondary_shock_end_tobs_axis_s,
    )
    solution = DynamicsSolution(r_tobs, r_gamma, radius, swept_mass_g, reverse_shock=reverse_shock)
    reverse_shock.causality = _reverse_shock_causality_diagnostics(config, reverse_params, solution, reverse_shock)
    if return_report:
        return solution, _solver_report(
            "dynamics_reverse",
            "shell-radius",
            "ok",
            kernel=str(config.dynamics_kernel),
            num_r=int(config.num_r),
            reverse_shock_global_allowed=bool(reverse_shock.causality.global_reverse_shock_allowed),
            reverse_shock_pressure_balance_condition_seen=bool(
                reverse_shock.causality.pressure_balance_condition_seen
            ),
            reverse_shock_local_fast_condition_seen=bool(reverse_shock.causality.local_fast_condition_seen),
            reverse_shock_started=bool(reverse_shock.causality.reverse_shock_started),
            reverse_shock_criteria_agree=bool(reverse_shock.causality.criteria_agree),
            reverse_shock_contact_radius_cm=float(reverse_shock.causality.contact_radius_cm),
            reverse_shock_reference_crossing_radius_cm=float(reverse_shock.causality.reference_crossing_radius_cm),
            reverse_shock_pressure_balance_start_radius_cm=float(
                reverse_shock.causality.pressure_balance_start_radius_cm
            ),
            reverse_shock_pressure_balance_start_ratio=float(reverse_shock.causality.pressure_balance_start_ratio),
            reverse_shock_fast_wave_crossing_radius_cm=float(reverse_shock.causality.fast_wave_crossing_radius_cm),
            reverse_shock_local_start_radius_cm=float(reverse_shock.causality.local_start_radius_cm),
            reverse_shock_actual_start_radius_cm=float(reverse_shock.causality.actual_start_radius_cm),
            reverse_shock_actual_start_pressure_ratio=float(reverse_shock.causality.actual_start_pressure_ratio),
            reverse_shock_actual_start_contact_fraction=float(reverse_shock.causality.actual_start_contact_fraction),
        )
    return solution


def _solve_electron_1d_standard(boundary, dynamics, v_seed, config, solver_name, return_report):
    module = _electron_module(solver_name)
    func = getattr(module, f"fs_electron_{solver_name}")
    args = [
        boundary, dynamics.r_tobs, dynamics.r_gamma, dynamics.radius, v_seed,
        config.num_gam_e, config.index_y, config.index_syn_integr, config.num_threads,
    ]
    if solver_name == "charint_1d":
        args.extend([
            1 if config.electron_adaptive_substeps else 0,
            config.electron_substep_rtol, config.electron_substep_min, config.electron_substep_max,
        ])
    gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = func(*args)
    return _finish_electron_solution(
        config, solver_name,
        _ELECTRON_1D_NU_GRID.get(solver_name, "log-gamma-1d"),
        gam_e, d_n_gam_e, l_syn_spec, seed_syn,
        nu=(nu_m, nu_c, nu_a), return_report=return_report,
    )


def solve_electron(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: RuntimeConfig,
    *,
    return_report: bool = False,
) -> ElectronSolution | tuple[ElectronSolution, SolverAdapterReport]:
    solver_name = _resolve_electron_solver(config)
    if config.cooling_kernel.lower() != "legacy":
        raise ValueError(f"Unsupported cooling kernel: {config.cooling_kernel}")
    if config.thermal_electrons and (solver_name not in ("fullhide_1d", "fullhide_1d_hz")):
        raise NotImplementedError("thermal_electrons currently requires electron_solver='fullhide_1d' or 'fullhide_1d_hz'.")
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
        return _finish_electron_solution(
            config,
            solver_name,
            "log-gamma-1d",
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            return_report=return_report,
        )

    if solver_name in _ELECTRON_1D_NU_SOLVERS:
        return _solve_electron_1d_standard(boundary, dynamics, v_seed, config, solver_name, return_report)

    if solver_name == "charint_2d":
        return _solve_electron_transport_2d(
            boundary,
            dynamics,
            v_seed,
            config,
            solver_name,
            use_characteristic_integrator=True,
            return_report=return_report,
        )

    if solver_name == "fullhide_2d":
        return _solve_electron_transport_2d(
            boundary,
            dynamics,
            v_seed,
            config,
            solver_name,
            use_characteristic_integrator=False,
            return_report=return_report,
        )

    if solver_name == "fullhide_2d_pic":
        electron_fullhide_2d_pic_module = _electron_module(solver_name)
        num_chi = _resolve_num_chi(config, solver_name)
        gam_e, d_n_gam_e_chi, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = (
            electron_fullhide_2d_pic_module.fs_electron_transport_2d_pic_core(
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
                False,
                "fullhide_2d_pic",
                bool(getattr(config, "electron_pic_uniform_b", False)),
                float(getattr(config, "electron_pic_eta_acc", 1.0)),
                float(getattr(config, "electron_pic_kappa_diff_scale", 1.0)),
                float(getattr(config, "electron_pic_bw_factor", 1.0)),
            )
        )
        chi_grid = _build_log_chi_grid(dynamics.r_gamma, num_chi)
        return _finish_electron_solution(
            config,
            solver_name,
            "log-gamma-log-chi-2d-pic",
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu=(nu_m, nu_c, nu_a),
            return_report=return_report,
            num_chi=num_chi,
            d_n_gam_e_chi=d_n_gam_e_chi,
            chi_grid=chi_grid,
        )

    if solver_name == "fullhide_1d_hz":
        electron_fullhide_1d_module = _electron_module("fullhide_1d_hz")
        gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = _solve_fullhide_1d(
            electron_fullhide_1d_module.fs_electron_fullhide_1d_hz,
            boundary,
            dynamics,
            v_seed,
            config,
        )
        return _finish_electron_solution(
            config,
            solver_name,
            "log-gamma-1d",
            gam_e,
            d_n_gam_e,
            l_syn_spec,
            seed_syn,
            nu=(nu_m, nu_c, nu_a),
            return_report=return_report,
        )

    electron_fullhide_1d_module = _electron_module("fullhide_1d")
    gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = _solve_fullhide_1d(
        electron_fullhide_1d_module.fs_electron_fullhide_1d,
        boundary,
        dynamics,
        v_seed,
        config,
    )
    return _finish_electron_solution(
        config,
        solver_name,
        "log-four-velocity-1d",
        gam_e,
        d_n_gam_e,
        l_syn_spec,
        seed_syn,
        nu=(nu_m, nu_c, nu_a),
        return_report=return_report,
    )


def _finish_electron_solution(
    config: RuntimeConfig,
    solver_name: str,
    grid_semantics: str,
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
    l_syn_spec: np.ndarray,
    seed_syn: np.ndarray,
    *,
    return_report: bool,
    nu: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
    num_chi: int = 1,
    **solution_kwargs,
) -> ElectronSolution | tuple[ElectronSolution, SolverAdapterReport]:
    if nu is not None:
        _emit_nu_callback(config, solver_name, *nu)
    solution = _build_electron_solution(
        gam_e,
        d_n_gam_e,
        l_syn_spec,
        seed_syn,
        **solution_kwargs,
    )
    if return_report:
        return solution, _solver_report(
            solver_name,
            grid_semantics,
            "ok",
            num_gam_e=int(config.num_gam_e),
            num_chi=int(num_chi),
        )
    return solution


def _solve_electron_transport_2d(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: RuntimeConfig,
    solver_name: str,
    *,
    use_characteristic_integrator: bool,
    return_report: bool,
) -> ElectronSolution | tuple[ElectronSolution, SolverAdapterReport]:
    electron_2d_module = _electron_module(solver_name)
    num_chi = _resolve_num_chi(config, solver_name)
    num_threads_2d = max(1, min(int(config.num_threads), int(num_chi), 16))
    emit_full_chi_spectrum = (
        str(config.geometry_kernel).lower() == "chi_eats_2d"
        or not _use_direct_chi_projection_contract(config)
    )
    (
        gam_e,
        d_n_gam_e_chi,
        d_n_gam_e,
        l_syn_spec,
        seed_syn,
        nu_m,
        nu_c,
        nu_a,
        l_syn_spec_chi,
        seed_syn_chi,
        tau_syn_chi,
        chi_radius_cm,
        chi_gamma_bulk,
        chi_dvolume_weight,
        b_chi_g,
    ) = electron_2d_module.fs_electron_transport_2d_core(
        boundary,
        dynamics.r_tobs,
        dynamics.r_gamma,
        dynamics.radius,
        v_seed,
        config.num_gam_e,
        num_chi,
        config.index_y,
        config.index_syn_integr,
        num_threads_2d,
        1 if emit_full_chi_spectrum else 0,
        config.electron_substep_max,
        use_characteristic_integrator,
        solver_name,
    )
    if not emit_full_chi_spectrum:
        l_syn_spec_chi = None
        seed_syn_chi = None
        tau_syn_chi = None
    uses_four_velocity = (
        not use_characteristic_integrator
        and str(config.fullhide2d_transport_model).lower() != "pwn_cr_v1"
    )
    grid_semantics = "log-four-velocity-q-mass-2d" if uses_four_velocity else "log-gamma-q-mass-2d"
    return _finish_electron_solution(
        config,
        solver_name,
        grid_semantics,
        gam_e,
        d_n_gam_e,
        l_syn_spec,
        seed_syn,
        nu=(nu_m, nu_c, nu_a),
        return_report=return_report,
        num_chi=num_chi,
        d_n_gam_e_chi=d_n_gam_e_chi,
        chi_grid=_build_q_mass_chi_grid(config, num_chi),
        l_syn_spec_chi=l_syn_spec_chi,
        seed_syn_chi=seed_syn_chi,
        tau_syn_chi=tau_syn_chi,
        chi_radius_cm=chi_radius_cm,
        chi_gamma_bulk=chi_gamma_bulk,
        chi_dvolume_weight=chi_dvolume_weight,
        b_chi_g=b_chi_g,
    )


def _use_direct_chi_projection_contract(config: RuntimeConfig) -> bool:
    return (
        str(config.geometry_kernel).lower() == "chi_eats_2d"
        and not bool(config.include_forward_ssc)
        and not bool(config.hadronic.enabled)
        and not bool(config.hadronic.reverse_enabled)
        and not bool(config.reverse)
        and not bool(config.reverse_shock.enabled)
    )


def _solve_fullhide_1d(
    kernel,
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: RuntimeConfig,
) -> tuple[np.ndarray, ...]:
    return kernel(
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
        1 if config.thermal_electrons else 0,
    )


def solve_electron_with_cooling_seed(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    cooling_seed: np.ndarray,
    config: RuntimeConfig,
    *,
    secondary_source_r: np.ndarray | None = None,
    return_report: bool = False,
) -> ElectronSolution | tuple[ElectronSolution, SolverAdapterReport]:
    solver_name = _resolve_electron_solver(config)
    if solver_name != "fullhide_1d":
        raise NotImplementedError("joint electron-photon coupling requires electron_solver='fullhide_1d'.")
    if int(config.index_y) != 1:
        raise NotImplementedError("joint electron-photon coupling requires ssc_cooling_mode='numeric_ic_kn' (index_y=1).")
    if bool(config.electron_adaptive_substeps):
        raise NotImplementedError("joint electron-photon coupling currently requires fixed electron substeps.")
    cooling_seed_arr = np.asfortranarray(np.asarray(cooling_seed, dtype=float))
    v_seed_arr = np.asarray(v_seed, dtype=float)
    radius_arr = np.asarray(dynamics.radius, dtype=float)
    if cooling_seed_arr.shape != (v_seed_arr.size, radius_arr.size):
        raise ValueError("joint cooling seed must have shape (num_frequency, num_radius).")
    if secondary_source_r is None:
        secondary_source_arr = np.zeros((int(config.num_gam_e), radius_arr.size), dtype=float, order="F")
    else:
        secondary_source_arr = np.asfortranarray(np.asarray(secondary_source_r, dtype=float))
        if secondary_source_arr.shape != (int(config.num_gam_e), radius_arr.size):
            raise ValueError("joint secondary electron source must have shape (num_gam_e, num_radius).")
    electron_fullhide_1d_module = _electron_module("fullhide_1d")
    gam_e, d_n_gam_e, l_syn_spec, seed_syn, nu_m, nu_c, nu_a = (
        electron_fullhide_1d_module.fs_electron_fullhide_1d_coupled(
            boundary,
            dynamics.r_tobs,
            dynamics.r_gamma,
            dynamics.radius,
            v_seed,
            cooling_seed_arr,
            secondary_source_arr,
            config.index_y,
            config.index_syn_integr,
            config.num_threads,
            0,
            config.electron_substep_rtol,
            config.electron_substep_min,
            config.electron_substep_max,
            1 if config.thermal_electrons else 0,
        )
    )
    solution = _build_electron_solution(
        gam_e,
        d_n_gam_e,
        l_syn_spec,
        seed_syn,
    )
    _emit_nu_callback(config, solver_name, nu_m, nu_c, nu_a)
    if return_report:
        return solution, _solver_report(
            solver_name,
            "log-gamma-1d-joint-cooling",
            "ok",
            num_gam_e=int(config.num_gam_e),
            num_chi=1,
        )
    return solution


def _emit_nu_callback(config: RuntimeConfig, label: str, nu_m, nu_c, nu_a) -> None:
    if config.nu_callback is not None:
        config.nu_callback(
            label,
            np.asarray(nu_m, dtype=float),
            np.asarray(nu_c, dtype=float),
            np.asarray(nu_a, dtype=float),
        )


def _resolve_electron_solver(config: RuntimeConfig) -> str:
    solver_name = config.electron_solver.lower()
    if solver_name not in _ELECTRON_MODULES:
        raise ValueError(f"Unsupported electron solver: {config.electron_solver}")
    return solver_name


def _electron_1d_transport_solver_id(config: RuntimeConfig) -> int:
    solver_name = _resolve_electron_solver(config)
    if solver_name not in ELECTRON_1D_TRANSPORT_IDS:
        raise NotImplementedError("reverse-shock electron transport supports electron_solver='fullhide_1d' or 'dg_1d'.")
    return ELECTRON_1D_TRANSPORT_IDS[solver_name]


def _resolve_pgamma_scheme(config: RuntimeConfig) -> str:
    key = str(config.hadronic.pgamma_scheme).lower()
    if key not in (_PGAMMA_SCHEME_DISABLED, _PGAMMA_SCHEME_HUMMER2010_RESPONSE):
        raise ValueError(f"Unsupported p-gamma scheme: {config.hadronic.pgamma_scheme}")
    return key


def _resolve_num_chi(config: RuntimeConfig, solver_name: str | None = None) -> int:
    resolved_solver = _resolve_electron_solver(config) if solver_name is None else solver_name
    user_value = config.downstream_num_chi
    if resolved_solver.endswith("_1d"):
        return 1
    if user_value is None:
        return 64
    if int(user_value) < 2:
        raise ValueError("downstream_num_chi must be >= 2 for 2d electron solvers.")
    return int(user_value)


def _build_log_chi_grid(r_gamma: np.ndarray, num_chi: int) -> np.ndarray:
    gamma_arr = np.asarray(r_gamma, dtype=float)
    chi_max = 1.0 + 8.0 * np.max(gamma_arr * gamma_arr)
    deta = np.log10(chi_max) / float(num_chi)
    eta_grid = (np.arange(num_chi, dtype=float) + 0.5) * deta
    return np.power(10.0, eta_grid)


def _build_q_mass_chi_grid(config: RuntimeConfig, num_chi: int) -> np.ndarray:
    sigma = 4.0
    q_active = 1.0 - (1.0 - 1.0 / sigma) ** sigma
    q_grid = (np.arange(num_chi, dtype=float) + 0.5) * q_active / float(num_chi)
    k_medium = 2 if float(config.a_star) > 0.0 else 0
    alpha = float(4 - k_medium) / float(3 - k_medium)
    return np.power(1.0 - q_grid, -alpha)


def _build_electron_solution(
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
    l_syn_spec: np.ndarray,
    seed_syn: np.ndarray,
    *,
    d_n_gam_e_chi: np.ndarray | None = None,
    chi_grid: np.ndarray | None = None,
    l_syn_spec_chi: np.ndarray | None = None,
    seed_syn_chi: np.ndarray | None = None,
    tau_syn_chi: np.ndarray | None = None,
    chi_radius_cm: np.ndarray | None = None,
    chi_gamma_bulk: np.ndarray | None = None,
    chi_dvolume_weight: np.ndarray | None = None,
    b_chi_g: np.ndarray | None = None,
) -> ElectronSolution:
    return ElectronSolution(
        gam_e=np.asarray(gam_e, dtype=float),
        d_n_gam_e=np.asarray(d_n_gam_e, dtype=float),
        l_syn_spec=np.asarray(l_syn_spec, dtype=float),
        seed_syn=np.asarray(seed_syn, dtype=float),
        d_n_gam_e_bh=None,
        d_n_gam_e_chi=None if d_n_gam_e_chi is None else np.asarray(d_n_gam_e_chi, dtype=float),
        chi_grid=None if chi_grid is None else np.asarray(chi_grid, dtype=float),
        l_syn_spec_chi=None if l_syn_spec_chi is None else np.asarray(l_syn_spec_chi, dtype=float),
        seed_syn_chi=None if seed_syn_chi is None else np.asarray(seed_syn_chi, dtype=float),
        tau_syn_chi=None if tau_syn_chi is None else np.asarray(tau_syn_chi, dtype=float),
        chi_radius_cm=None if chi_radius_cm is None else np.asarray(chi_radius_cm, dtype=float),
        chi_gamma_bulk=None if chi_gamma_bulk is None else np.asarray(chi_gamma_bulk, dtype=float),
        chi_dvolume_weight=None if chi_dvolume_weight is None else np.asarray(chi_dvolume_weight, dtype=float),
        b_chi_g=None if b_chi_g is None else np.asarray(b_chi_g, dtype=float),
    )


def solve_hadronic(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    electron: ElectronSolution,
    v_seed: np.ndarray,
    seed_target: np.ndarray,
    config: RuntimeConfig,
    *,
    return_report: bool = False,
) -> HadronicSolution | None | tuple[HadronicSolution | None, SolverAdapterReport]:
    del boundary
    if not bool(config.hadronic.enabled) or float(config.hadronic.epsilon_p) <= 0.0:
        report = _solver_report("hadronic_disabled", "log-gamma-1d", "disabled", backend="none")
        return (None, report) if return_report else None
    hadronic_solver_key = str(config.hadronic.solver).lower()
    if hadronic_solver_key not in {"legacy_1d", "am3_1d"}:
        raise ValueError(f"Unsupported hadronic solver: {config.hadronic.solver}")
    hadronic_solver = hadronic_solver_key
    pgamma_scheme = _resolve_pgamma_scheme(config)
    electron_solver = _resolve_electron_solver(config)
    if not electron_solver.endswith("_1d"):
        raise NotImplementedError("The current hadronic kernel only supports 1d forward-shock electron solvers.")

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
            "Choose 'hummer_2010_response'."
        )
    prev_radius = np.empty_like(dynamics.radius)
    prev_radius[0] = 0.0
    prev_radius[1:] = dynamics.radius[:-1]
    shell_volume_cm3 = (4.0 / 3.0) * np.pi * (dynamics.radius**3 - prev_radius**3)
    shell_energy_inj_erg = (
        float(config.hadronic.epsilon_p)
        * (dynamics.r_gamma - 1.0)
        * shell_volume_cm3
        * np.asarray(ambient_density(dynamics.radius, config), dtype=float)
        * constants.para_m_p
        * constants.para_c**2
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
                "ok",
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
        am3_process_power = np.zeros((len(HUMMER_PROCESS_GROUP_LABELS), num_gam_p, num_r), dtype=float)

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
        d_n_gam_e_bh=None,
        secondary_electron_source_r=None,
        tau_bh=None,
        bh_photon_loss_rate=None,
        l_had_hadronic_inverse_compton=None,
        l_had_pair_production=None,
        tau_pg=None,
        pg_photon_survival=None,
        am3_process_power=np.asarray(am3_process_power, dtype=float).reshape(len(HUMMER_PROCESS_GROUP_LABELS), num_gam_p, num_r),
    )
    if return_report:
        backend = "fortran_core" if hadronic_solver == "legacy_1d" else HUMMER2010_RESPONSE_BACKEND
        return solution, _solver_report(
            hadronic_solver,
            "log-gamma-1d",
            "ok",
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
    config: RuntimeConfig,
    pp_target_density_cm3: np.ndarray | None = None,
) -> HadronicSolution:
    radius = dynamics.radius
    gamma_bulk = dynamics.r_gamma
    b_field = magnetic_field_g
    v_seed_arr = v_seed_hz
    seed_target_arr = seed_target_hz
    shell_energy_inj = shell_energy_inj_erg
    pp_target_density_arr = pp_target_density_cm3
    if pp_target_density_arr is not None and pp_target_density_arr.shape != radius.shape:
        raise ValueError("pp_target_density_cm3 must match the shell radius grid.")
    num_r = radius.size
    num_nu = v_seed_arr.size
    num_gam_p = int(config.hadronic.num_gam_p)
    num_nu_nu = int(config.hadronic.num_nu_nu)
    gam_e = electron_gamma
    if pp_target_density_arr is None:
        pp_target_density_arr = ambient_density(radius, config)
    t_total_start = time.perf_counter()
    (
        gam_p, gam_secondary, d_n_gam_p, l_had_syn_spec, seed_had_syn,
        l_had_pg_gamma, neutrino_frequency_hz, neutrino_luminosity, l_had_bh,
        seed_had_bh, d_n_gam_e_bh, secondary_electron_source_r, tau_bh,
        bh_photon_loss_rate, l_had_hic, d_n_gam_n, d_n_gam_pi_plus,
        d_n_gam_pi_minus, d_n_gam_mu_minus_left, d_n_gam_mu_minus_right,
        d_n_gam_mu_plus_left, d_n_gam_mu_plus_right, l_had_pion_synch,
        l_had_muon_synch, l_had_pion_ic, l_had_muon_ic, tau_pg,
        pg_photon_survival, am3_process_power,
    ) = hadronic_legacy_module.fs_hadronic_formal_transport_1d(
        dynamics.r_tobs,
        gamma_bulk,
        radius,
        b_field,
        v_seed_arr,
        seed_target_arr,
        gam_e,
        shell_energy_inj,
        pp_target_density_arr,
        float(config.hadronic.p_p),
        float(config.hadronic.eta_acc),
        int(config.index_syn_integr),
        1 if bool(config.hadronic.include_proton_synch) else 0,
        1 if bool(config.hadronic.include_pg) else 0,
        1 if bool(config.hadronic.include_neutrino) else 0,
        1 if bool(config.hadronic.include_bethe_heitler) else 0,
        1 if bool(config.hadronic.include_hadronic_inverse_compton) else 0,
        1 if bool(config.hadronic.include_pp) else 0,
        1 if bool(config.hadronic.quantum_syn) else 0,
        int(config.num_threads),
        num_gam_p,
        num_nu_nu,
    )
    timings = {"formal_transport_fortran": time.perf_counter() - t_total_start}
    sed_components = {
        "proton_synchrotron": l_had_syn_spec,
        "pgamma_pi0_decay": l_had_pg_gamma,
        "pion_synchrotron": l_had_pion_synch,
        "muon_synchrotron": l_had_muon_synch,
        "pion_inverse_compton": l_had_pion_ic,
        "muon_inverse_compton": l_had_muon_ic,
    }
    if bool(config.hadronic.include_bethe_heitler):
        sed_components["bethe_heitler"] = l_had_bh
    if bool(config.hadronic.include_hadronic_inverse_compton):
        sed_components["hadronic_inverse_compton"] = l_had_hic
    return HadronicSolution(
        solver="am3_1d",
        gam_p=np.asarray(gam_p, dtype=float),
        d_n_gam_p=np.asarray(d_n_gam_p, dtype=float).reshape(num_gam_p, num_r),
        l_had_syn_spec=np.asarray(l_had_syn_spec, dtype=float).reshape(num_nu, num_r),
        seed_had_syn=np.asarray(seed_had_syn, dtype=float).reshape(num_nu, num_r),
        l_had_pg_gamma=np.asarray(l_had_pg_gamma, dtype=float).reshape(num_nu, num_r),
        neutrino_frequency_hz=np.asarray(neutrino_frequency_hz, dtype=float),
        neutrino_luminosity=np.asarray(neutrino_luminosity, dtype=float).reshape(num_nu_nu, num_r),
        l_had_bethe_heitler=l_had_bh if bool(config.hadronic.include_bethe_heitler) else None,
        seed_had_bethe_heitler=seed_had_bh if bool(config.hadronic.include_bethe_heitler) else None,
        d_n_gam_e_bh=d_n_gam_e_bh if (bool(config.hadronic.include_bethe_heitler) or bool(config.hadronic.include_pp)) else None,
        secondary_electron_source_r=secondary_electron_source_r if (bool(config.hadronic.include_bethe_heitler) or bool(config.hadronic.include_pp)) else None,
        tau_bh=tau_bh if bool(config.hadronic.include_bethe_heitler) else None,
        bh_photon_loss_rate=bh_photon_loss_rate if bool(config.hadronic.include_bethe_heitler) else None,
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
        sed_components=sed_components,
    )



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


def solve_reverse_shock_emission(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: RuntimeConfig,
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

    l_syn_spec, seed_syn = _compute_reverse_shock_synchrotron_emission(dynamics, v_seed, config)
    secondary_rs = _compute_secondary_reverse_shock_synchrotron(dynamics, v_seed, config, reverse_params)
    if secondary_rs is not None:
        l_syn_spec = l_syn_spec + secondary_rs.luminosity_syn

    rs_hadronic = None
    if bool(config.hadronic.reverse_enabled) and float(config.hadronic.reverse_epsilon_p) > 0.0:
        requires_full_chain = any((
            bool(config.hadronic.include_pg),
            bool(config.hadronic.include_neutrino),
            bool(config.hadronic.include_bethe_heitler),
            bool(config.hadronic.include_hadronic_inverse_compton),
            bool(config.hadronic.include_pp),
        ))
        if requires_full_chain:
            if (bool(config.hadronic.include_pg) or bool(config.hadronic.include_neutrino)) and (
                _resolve_pgamma_scheme(config) != _PGAMMA_SCHEME_HUMMER2010_RESPONSE
            ):
                raise ValueError("Reverse-shock p-gamma currently requires pgamma_scheme='hummer_2010_response'.")
            rs_swept = np.asarray(dynamics.reverse_shock.swept_mass_g, dtype=float)
            rs_shell_mass = np.empty_like(rs_swept)
            rs_shell_mass[0] = rs_swept[0]
            rs_shell_mass[1:] = rs_swept[1:] - rs_swept[:-1]
            prev_radius = np.empty_like(dynamics.radius)
            prev_radius[0] = 0.0
            prev_radius[1:] = dynamics.radius[:-1]
            rs_shell_volume = (4.0 / 3.0) * np.pi * (dynamics.radius**3 - prev_radius**3)
            rs_shell_energy = (
                float(config.hadronic.reverse_epsilon_p)
                * rs_shell_mass
                * (np.asarray(dynamics.r_gamma, dtype=float) - 1.0)
                * constants.para_c
                * constants.para_c
            )
            rs_hadronic = _solve_hadronic_hummer_transport_coupled(
                dynamics,
                np.asarray(dynamics.reverse_shock.magnetic_field_g, dtype=float),
                gam_e,
                v_seed,
                seed_syn,
                rs_shell_energy,
                config,
                pp_target_density_cm3=rs_shell_mass / (constants.para_m_p * rs_shell_volume),
            )
        else:
            num_gam_p = int(config.hadronic.num_gam_p)
            num_nu = int(np.asarray(v_seed, dtype=float).size)
            num_r = int(np.asarray(dynamics.radius, dtype=float).size)
            shell_energy = (
                float(config.hadronic.reverse_epsilon_p)
                * np.asarray(dynamics.reverse_shock.swept_mass_g, dtype=float)
                * np.maximum(np.asarray(dynamics.r_gamma, dtype=float) - 1.0, 0.0)
                * constants.para_c
                * constants.para_c
            )
            gam_p, d_n_gam_p, l_had_syn_spec, seed_had_syn = _hadronic_reverse_module.fs_hadronic_reverse_1d(
                np.asarray(dynamics.r_tobs, dtype=float),
                np.asarray(dynamics.r_gamma, dtype=float),
                np.asarray(dynamics.radius, dtype=float),
                shell_energy,
                np.asarray(dynamics.reverse_shock.magnetic_field_g, dtype=float),
                np.asarray(v_seed, dtype=float),
                1 if bool(config.hadronic.include_proton_synch) else 0,
                num_gam_p,
            )
            rs_hadronic = ReverseShockHadronicSolution(
                gam_p=np.asarray(gam_p, dtype=float),
                d_n_gam_p=np.asarray(d_n_gam_p, dtype=float).reshape(num_gam_p, num_r),
                l_had_syn_spec=np.asarray(l_had_syn_spec, dtype=float).reshape(num_nu, num_r),
                seed_had_syn=np.asarray(seed_had_syn, dtype=float).reshape(num_nu, num_r),
            )

    return ReverseShockEmission(
        l_syn_spec=l_syn_spec,
        seed_syn=seed_syn,
        rs_hadronic=rs_hadronic,
        secondary_rs=secondary_rs,
    )


def _solve_reverse_shock_electrons(
    boundary: np.ndarray,
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: RuntimeConfig,
    reverse_params: ReverseShockParameters,
) -> tuple[np.ndarray, np.ndarray]:
    if dynamics.reverse_shock is None:
        raise ValueError("Reverse shock dynamics are required to compute reverse electrons.")
    module = _electron_reverse_module().electron_reverse_kernel
    delta_0 = reverse_params.delta_t_s * constants.para_c
    para_m_ej = reverse_shell_baryonic_mass(config)
    solver_id = _electron_1d_transport_solver_id(config)
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
        boundary[20],
        boundary[21],
        boundary[22],
        boundary[-1],
        dynamics.reverse_shock.t_cross,
        dynamics.reverse_shock.r_cross,
        dynamics.reverse_shock.u3_cross_erg,
        dynamics.reverse_shock.v3_cross_cm3,
        dynamics.reverse_shock.m3_cross_g,
        dynamics.r_tobs,
        dynamics.r_gamma,
        dynamics.radius,
        dynamics.reverse_shock.magnetic_field_g,
        dynamics.reverse_shock.swept_mass_g,
        dynamics.reverse_shock.internal_energy_erg,
        dynamics.reverse_shock.comoving_volume_cm3,
        v_seed,
        config.num_gam_e,
        config.index_y,
        config.index_syn_integr,
        config.num_threads,
        solver_id=solver_id,
    )
    return np.asarray(gam_e, dtype=float), np.asarray(d_n_gam_e, dtype=float)


def _solve_shell_syn(args):
    i, radius, bfield, index_syn, gam_e, d_n_gam_e_col, v_seed = args
    if bfield <= 0.0:
        return i, None, None
    p_syn, seed_syn = electron_radiation_module.get_syn_selected(
        index_syn, float(radius), float(bfield), 1, gam_e, d_n_gam_e_col, v_seed)
    return i, p_syn, seed_syn


def _compute_reverse_shock_synchrotron_emission(
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: RuntimeConfig,
) -> tuple[np.ndarray, np.ndarray]:
    if dynamics.reverse_shock is None:
        raise ValueError("Reverse shock dynamics are required to compute reverse emission.")
    gam_e = dynamics.reverse_shock.gam_e
    d_n_gam_e = dynamics.reverse_shock.d_n_gam_e
    if gam_e is None or d_n_gam_e is None:
        raise ValueError("Reverse shock electrons are required to compute reverse emission.")
    num_nu, num_r = v_seed.shape[0], dynamics.radius.shape[0]
    l_syn_spec, seed_syn = np.zeros((2, num_nu, num_r), dtype=float)
    bfield = dynamics.reverse_shock.magnetic_field_g
    radius = dynamics.radius

    tasks = [(i, radius[i], bfield[i], config.index_syn_integr,
              gam_e, d_n_gam_e[:, i], v_seed) for i in range(num_r) if bfield[i] > 0.0]
    if not tasks:
        return l_syn_spec, seed_syn
    workers = min(config.num_threads, len(tasks))
    if workers > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=workers) as ex:
            for i, p_syn_i, seed_syn_i in ex.map(_solve_shell_syn, tasks):
                if p_syn_i is not None:
                    l_syn_spec[:, i] = p_syn_i
                    seed_syn[:, i] = seed_syn_i
    else:
        for task in tasks:
            i, p_syn_i, seed_syn_i = _solve_shell_syn(task)
            if p_syn_i is not None:
                l_syn_spec[:, i] = p_syn_i
                seed_syn[:, i] = seed_syn_i
    return l_syn_spec, seed_syn


def _compute_secondary_reverse_shock_synchrotron(
    dynamics: DynamicsSolution,
    v_seed: np.ndarray,
    config: RuntimeConfig,
    reverse_params: ReverseShockParameters,
) -> SecondaryReverseShockState | None:
    jump_r, _, _ = density_jump_arrays(config)
    if jump_r.size == 0 or dynamics.reverse_shock is None:
        return None
    radius = np.asarray(dynamics.radius, dtype=float)
    if reverse_params.p <= 2.0:
        raise ValueError("secondary reverse shock v1 requires p > 2.")
    event_active = np.asarray(dynamics.reverse_shock.secondary_event_active, dtype=bool)
    start_radius = np.asarray(dynamics.reverse_shock.secondary_start_radius_cm, dtype=float)
    end_radius = np.asarray(dynamics.reverse_shock.secondary_shock_end_radius_cm, dtype=float)
    start_tobs = np.asarray(dynamics.reverse_shock.secondary_start_tobs_axis_s, dtype=float)
    end_tobs = np.asarray(dynamics.reverse_shock.secondary_shock_end_tobs_axis_s, dtype=float)
    if not np.any(event_active):
        return None
    rs = dynamics.reverse_shock
    if not np.any(rs.secondary_swept_mass_g > 0.0):
        return None
    parent_branch = np.zeros(rs.secondary_branch_gamma_43.shape[0], dtype=np.int32)
    consecutive = event_active[:-1] & event_active[1:]
    parent_branch[1:] = np.where(consecutive, np.arange(1, parent_branch.size), 0)
    (
        gam_e_sec, dist, branch_luminosity, luminosity,
        reacceleration_seed_energy, reaccelerated_energy,
    ) = _electron_reverse_module().electron_reverse_kernel.electron_secondary_reverse_branch_reaccelerated(
        reverse_params.epsilon_e,
        reverse_params.epsilon_b,
        reverse_params.p,
        reverse_params.f_e,
        config.z,
        dynamics.r_tobs,
        dynamics.r_gamma,
        radius,
        rs.secondary_branch_magnetic_field_g,
        rs.secondary_branch_swept_mass_g,
        rs.secondary_branch_internal_energy_erg,
        rs.secondary_branch_comoving_volume_cm3,
        rs.secondary_branch_gamma_m,
        rs.secondary_branch_gamma_43,
        rs.secondary_branch_compression,
        parent_branch,
        v_seed,
        config.num_gam_e,
        config.index_syn_integr,
        config.num_threads,
        solver_id=_electron_1d_transport_solver_id(config),
    )
    return SecondaryReverseShockState(
        luminosity_syn=luminosity,
        branch_luminosity_syn=branch_luminosity,
        gam_e=gam_e_sec,
        d_n_gam_e=dist,
        event_active=event_active,
        start_radius_cm=start_radius,
        shock_end_radius_cm=end_radius,
        start_tobs_axis_s=start_tobs,
        shock_end_tobs_axis_s=end_tobs,
        gamma_contact=rs.secondary_gamma_contact,
        pressure_3=rs.secondary_pressure_3,
        gamma_43=rs.secondary_gamma_43,
        beta_rs=rs.secondary_beta_rs,
        dissipated_energy_density=rs.secondary_dissipated_energy_density,
        dissipated_energy_erg=rs.secondary_dissipated_energy_erg,
        electron_injected_energy_erg=rs.secondary_electron_injected_energy_erg,
        swept_mass_g=rs.secondary_swept_mass_g,
        internal_energy_erg=rs.secondary_internal_energy_erg,
        comoving_volume_cm3=rs.secondary_comoving_volume_cm3,
        pressure_total=rs.secondary_pressure_total,
        enthalpy_density_total=rs.secondary_enthalpy_density_total,
        branch_swept_mass_g=rs.secondary_branch_swept_mass_g,
        branch_internal_energy_erg=rs.secondary_branch_internal_energy_erg,
        branch_comoving_volume_cm3=rs.secondary_branch_comoving_volume_cm3,
        branch_magnetic_field_g=rs.secondary_branch_magnetic_field_g,
        branch_gamma_m=rs.secondary_branch_gamma_m,
        branch_gamma_contact=rs.secondary_branch_gamma_contact,
        branch_gamma_43=rs.secondary_branch_gamma_43,
        branch_compression=rs.secondary_branch_compression,
        branch_beta_rs=rs.secondary_branch_beta_rs,
        branch_dissipated_energy_density=rs.secondary_branch_dissipated_energy_density,
        branch_reacceleration_seed_energy_erg=reacceleration_seed_energy,
        branch_reaccelerated_energy_erg=reaccelerated_energy,
        magnetic_field_g=rs.secondary_magnetic_field_g,
    )


def _resolve_reverse_shock_parameters(config: RuntimeConfig) -> ReverseShockParameters | None:
    reverse_enabled = config.reverse or config.reverse_shock.enabled
    if not reverse_enabled:
        return None

    if config.reverse_shock.delta_t_s is None:
        raise ValueError("ReverseShockConfig.delta_t_s must be set when reverse shock is enabled.")
    if config.reverse_shock.sigma < 0.0:
        raise ValueError("ReverseShockConfig.sigma must be non-negative.")

    return ReverseShockParameters(
        delta_t_s=config.reverse_shock.delta_t_s,
        sigma=config.reverse_shock.sigma,
        epsilon_e=config.epsilon_e if config.reverse_shock.epsilon_e is None else config.reverse_shock.epsilon_e,
        epsilon_b=config.epsilon_b if config.reverse_shock.epsilon_b is None else config.reverse_shock.epsilon_b,
        p=config.p if config.reverse_shock.p is None else config.reverse_shock.p,
        f_e=config.f_e if config.reverse_shock.f_e is None else config.reverse_shock.f_e,
    )


def _renormalize_reverse_shock_distribution(
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
    swept_mass_g: np.ndarray,
    f_e: float,
) -> np.ndarray:
    gam = np.asarray(gam_e, dtype=float)
    dist = np.where(np.isfinite(d_n_gam_e) & (d_n_gam_e > 0.0), d_n_gam_e, 0.0)
    targets = np.asarray(swept_mass_g, dtype=float) / constants.para_m_p * float(f_e)
    totals = np.trapezoid(dist, gam, axis=0)
    valid = np.isfinite(totals) & (totals > 0.0) & np.isfinite(targets) & (targets > 0.0)
    dist[:, valid] *= (targets[valid] / totals[valid])[None, :]
    return dist
