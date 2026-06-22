"""
Centralized type definitions for ASGARD.

This module contains all dataclass definitions used throughout the codebase,
providing a single source of truth for data structures.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


# ============================================================================
# State and Solution Types (from asgard_state.py)
# ============================================================================

@dataclass
class BranchState:
    """State information for a single shock branch (forward or reverse)."""
    characteristic_time_s: np.ndarray
    gamma: np.ndarray
    radius_cm: np.ndarray
    swept_mass_g: np.ndarray
    doppler: np.ndarray
    magnetic_field_g: np.ndarray


@dataclass
class FluxComponents:
    """Flux components from different emission processes."""
    total: np.ndarray
    fwd_sync: np.ndarray
    fwd_ssc: np.ndarray
    fwd_hadronic_gamma: np.ndarray | None
    fwd_hadronic_bethe_heitler: np.ndarray | None
    fwd_hadronic_inverse_compton: np.ndarray | None
    fwd_hadronic_pair_production: np.ndarray | None
    rev_sync: np.ndarray | None
    rev_ssc: np.ndarray | None
    cross_ic: np.ndarray | None
    fwd: BranchState
    rev: BranchState | None


@dataclass(frozen=True)
class SolverAdapterReport:
    """Thin solver-adapter report exposed to the orchestration layer."""
    solver: str
    grid_semantics: str
    status: str
    diagnostics: dict[str, Any] = field(default_factory=dict)


@dataclass
class PhotonFieldState:
    """Photon-field contract passed between electron, hadronic, and observer stages."""
    seed_frequency_hz: np.ndarray
    forward_syn_seed: np.ndarray
    hadronic_forward_ssc_seed: np.ndarray
    hadronic_target_seed: np.ndarray
    absorption_syn_seed: np.ndarray
    absorption_ssc_seed: np.ndarray


@dataclass
class ObserverState:
    """Observer-side assembly state prior to interpolation onto query grids."""
    prefactor: np.ndarray
    tau_extra: np.ndarray
    tau_pair: np.ndarray
    components: FluxComponents


@dataclass
class SolveState:
    """Complete state from a simulation solve."""
    config: Any  # RuntimeConfig - avoid circular import
    setup: Any  # SimulationSetup
    policy: Any  # ExecutionPolicy
    dynamics: Any  # DynamicsSolution
    electron: Any  # ElectronSolution
    photon_field: PhotonFieldState
    hadronic: Any | None  # HadronicSolution
    observer: ObserverState
    reverse_emission: Any | None  # ReverseShockEmission
    components: FluxComponents
    requested_frequency_min_hz: float | None
    requested_frequency_max_hz: float | None
    timings: dict[str, float] = field(default_factory=dict)
    adapter_reports: dict[str, SolverAdapterReport] = field(default_factory=dict)


@dataclass
class ObsState:
    """Observed state with frequency-dependent components."""
    state: SolveState
    setup: Any  # SimulationSetup
    frequencies_hz: np.ndarray
    components: dict[str, np.ndarray | None]


# ============================================================================
# Solver Types (from asgard_runtime.py)
# ============================================================================

@dataclass
class ReverseShockParameters:
    """Parameters for reverse shock physics."""
    delta_t_s: float
    sigma: float
    epsilon_e: float
    epsilon_b: float
    p: float
    f_e: float


@dataclass(frozen=True)
class ReverseShockCausalityDiagnostics:
    """Global and local reverse-shock causality diagnostics."""
    medium: str
    initial_super_fast: bool
    global_reverse_shock_allowed: bool
    pressure_balance_condition_seen: bool
    local_fast_condition_seen: bool
    reverse_shock_started: bool
    criteria_agree: bool
    contact_radius_cm: float
    reference_crossing_radius_cm: float
    pressure_balance_start_radius_cm: float
    pressure_balance_start_tobs_s: float
    pressure_balance_start_index: int
    pressure_balance_start_ratio: float
    pressure_balance_start_contact_fraction: float
    fast_wave_crossing_radius_cm: float
    fast_wave_crossing_tobs_s: float
    fast_wave_crossing_index: int
    local_start_radius_cm: float
    local_start_tobs_s: float
    local_start_index: int
    actual_start_radius_cm: float
    actual_start_tobs_s: float
    actual_start_index: int
    actual_start_pressure_ratio: float
    actual_start_contact_fraction: float


@dataclass
class ReverseShockDynamics:
    """Dynamics solution for reverse shock."""
    t_cross: float
    r_cross: float
    e3_cross: float
    gam20: float
    u3_cross_erg: float
    v3_cross_cm3: float
    m3_cross_g: float
    gamma_m_cross: float
    ordered_magnetic_cross_g: float
    swept_mass_g: np.ndarray
    magnetic_field_g: np.ndarray
    internal_energy_erg: np.ndarray
    comoving_volume_cm3: np.ndarray
    gamma34: np.ndarray
    secondary_branch_swept_mass_g: np.ndarray | None = None
    secondary_branch_internal_energy_erg: np.ndarray | None = None
    secondary_branch_comoving_volume_cm3: np.ndarray | None = None
    secondary_branch_magnetic_field_g: np.ndarray | None = None
    secondary_swept_mass_g: np.ndarray | None = None
    secondary_internal_energy_erg: np.ndarray | None = None
    secondary_comoving_volume_cm3: np.ndarray | None = None
    secondary_magnetic_field_g: np.ndarray | None = None
    secondary_pressure_total: np.ndarray | None = None
    secondary_enthalpy_density_total: np.ndarray | None = None
    secondary_gamma_contact: np.ndarray | None = None
    secondary_pressure_3: np.ndarray | None = None
    secondary_gamma_43: np.ndarray | None = None
    secondary_beta_rs: np.ndarray | None = None
    secondary_dissipated_energy_density: np.ndarray | None = None
    secondary_dissipated_energy_erg: np.ndarray | None = None
    secondary_electron_injected_energy_erg: np.ndarray | None = None
    secondary_branch_gamma_m: np.ndarray | None = None
    secondary_branch_gamma_contact: np.ndarray | None = None
    secondary_branch_gamma_43: np.ndarray | None = None
    secondary_branch_compression: np.ndarray | None = None
    secondary_branch_beta_rs: np.ndarray | None = None
    secondary_branch_dissipated_energy_density: np.ndarray | None = None
    secondary_event_active: np.ndarray | None = None
    secondary_start_radius_cm: np.ndarray | None = None
    secondary_shock_end_radius_cm: np.ndarray | None = None
    secondary_start_tobs_axis_s: np.ndarray | None = None
    secondary_shock_end_tobs_axis_s: np.ndarray | None = None
    causality: ReverseShockCausalityDiagnostics | None = None
    gam_e: np.ndarray | None = None
    d_n_gam_e: np.ndarray | None = None


@dataclass
class DynamicsSolution:
    """Solution from dynamics solver."""
    r_tobs: np.ndarray
    r_gamma: np.ndarray
    radius: np.ndarray
    swept_mass_g: np.ndarray
    reverse_shock: ReverseShockDynamics | None = None


@dataclass
class ElectronSolution:
    """Solution from electron solver."""
    gam_e: np.ndarray
    d_n_gam_e: np.ndarray
    l_syn_spec: np.ndarray
    seed_syn: np.ndarray
    d_n_gam_e_bh: np.ndarray | None = None
    d_n_gam_e_chi: np.ndarray | None = None
    chi_grid: np.ndarray | None = None
    l_syn_spec_chi: np.ndarray | None = None
    seed_syn_chi: np.ndarray | None = None
    tau_syn_chi: np.ndarray | None = None
    chi_radius_cm: np.ndarray | None = None
    chi_gamma_bulk: np.ndarray | None = None
    chi_dvolume_weight: np.ndarray | None = None
    b_chi_g: np.ndarray | None = None


@dataclass
class ReverseShockEmission:
    """Emission from reverse shock."""
    l_syn_spec: np.ndarray
    seed_syn: np.ndarray
    rs_hadronic: Any | None = None
    secondary_rs: Any | None = None


@dataclass
class HadronicSolution:
    """Solution from the 1d hadronic solver."""
    solver: str
    gam_p: np.ndarray
    d_n_gam_p: np.ndarray
    l_had_syn_spec: np.ndarray
    seed_had_syn: np.ndarray
    l_had_pg_gamma: np.ndarray
    neutrino_frequency_hz: np.ndarray
    neutrino_luminosity: np.ndarray
    l_had_bethe_heitler: np.ndarray | None = None
    seed_had_bethe_heitler: np.ndarray | None = None
    d_n_gam_e_bh: np.ndarray | None = None
    secondary_electron_source_r: np.ndarray | None = None
    tau_bh: np.ndarray | None = None
    bh_photon_loss_rate: np.ndarray | None = None
    l_had_hadronic_inverse_compton: np.ndarray | None = None
    l_had_pair_production: np.ndarray | None = None
    gam_secondary: np.ndarray | None = None
    d_n_gam_n: np.ndarray | None = None
    d_n_gam_pi_plus: np.ndarray | None = None
    d_n_gam_pi_minus: np.ndarray | None = None
    d_n_gam_mu_minus_left: np.ndarray | None = None
    d_n_gam_mu_minus_right: np.ndarray | None = None
    d_n_gam_mu_plus_left: np.ndarray | None = None
    d_n_gam_mu_plus_right: np.ndarray | None = None
    l_had_pion_synch: np.ndarray | None = None
    l_had_muon_synch: np.ndarray | None = None
    l_had_pion_inverse_compton: np.ndarray | None = None
    l_had_muon_inverse_compton: np.ndarray | None = None
    tau_pg: np.ndarray | None = None
    pg_photon_survival: np.ndarray | None = None
    am3_process_power: np.ndarray | None = None
    timings: dict[str, float] = field(default_factory=dict)
    sed_components: dict[str, np.ndarray] = field(default_factory=dict)


DynamicsState = DynamicsSolution
ElectronState = ElectronSolution
HadronicState = HadronicSolution


__all__ = [
    # State types
    "BranchState",
    "FluxComponents",
    "SolverAdapterReport",
    "PhotonFieldState",
    "ObserverState",
    "SolveState",
    "ObsState",
    # Solver types
    "ReverseShockParameters",
    "ReverseShockDynamics",
    "DynamicsSolution",
    "DynamicsState",
    "ElectronSolution",
    "ElectronState",
    "HadronicSolution",
    "HadronicState",
    "ReverseShockEmission",
]
