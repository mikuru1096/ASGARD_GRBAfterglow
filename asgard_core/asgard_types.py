"""
Centralized type definitions for ASGARD.

This module contains all dataclass definitions used throughout the codebase,
providing a single source of truth for data structures.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, ClassVar, Optional

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
    nu_m: np.ndarray
    nu_c: np.ndarray
    nu_a: np.ndarray
    nu_M: np.ndarray
    cooling_timescale_s: np.ndarray | None = None
    dynamical_timescale_s: np.ndarray | None = None


@dataclass
class FluxComponents:
    """Flux components from different emission processes."""
    total: np.ndarray
    fwd_sync: np.ndarray
    fwd_ssc: np.ndarray
    fwd_hadronic_gamma: Optional[np.ndarray]
    fwd_hadronic_bethe_heitler: Optional[np.ndarray]
    fwd_hadronic_inverse_compton: Optional[np.ndarray]
    fwd_hadronic_pair_production: Optional[np.ndarray]
    rev_sync: Optional[np.ndarray]
    rev_ssc: Optional[np.ndarray]
    cross_ic: Optional[np.ndarray]
    fwd: BranchState
    rev: Optional[BranchState]


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
    producer: ClassVar[str] = "photon_field_stage"
    consumers: ClassVar[tuple[str, ...]] = ("hadronic_stage", "observer_stage")
    mutable_fields: ClassVar[tuple[str, ...]] = ()

    seed_frequency_hz: np.ndarray
    forward_syn_seed: np.ndarray
    hadronic_forward_ssc_seed: np.ndarray
    hadronic_target_seed: np.ndarray
    absorption_syn_seed: np.ndarray
    absorption_ssc_seed: np.ndarray


@dataclass
class ObserverState:
    """Observer-side assembly state prior to interpolation onto query grids."""
    producer: ClassVar[str] = "observer_stage"
    consumers: ClassVar[tuple[str, ...]] = ("projection_stage", "api")
    mutable_fields: ClassVar[tuple[str, ...]] = ()

    prefactor: np.ndarray
    tau_extra: np.ndarray
    tau_pair: np.ndarray
    components: FluxComponents


@dataclass
class SolveState:
    """Complete state from a simulation solve."""
    config: Any  # FitConfig - avoid circular import
    setup: Any  # SimulationSetup
    policy: Any  # ExecutionPolicy
    dynamics: Any  # DynamicsSolution
    electron: Any  # ElectronSolution
    photon_field: PhotonFieldState
    hadronic: Optional[Any]  # HadronicSolution
    observer: ObserverState
    reverse_emission: Optional[Any]  # ReverseShockEmission
    components: FluxComponents
    requested_frequency_min_hz: Optional[float]
    requested_frequency_max_hz: Optional[float]
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
    gam_e: np.ndarray | None = None
    d_n_gam_e: np.ndarray | None = None


@dataclass
class DynamicsSolution:
    """Solution from dynamics solver."""
    producer: ClassVar[str] = "dynamics_stage"
    consumers: ClassVar[tuple[str, ...]] = ("electron_stage", "observer_stage")
    mutable_fields: ClassVar[tuple[str, ...]] = ()

    r_tobs: np.ndarray
    r_gamma: np.ndarray
    radius: np.ndarray
    swept_mass_g: np.ndarray
    reverse_shock: ReverseShockDynamics | None = None


@dataclass
class ElectronSolution:
    """Solution from electron solver."""
    producer: ClassVar[str] = "electron_stage"
    consumers: ClassVar[tuple[str, ...]] = ("photon_field_stage", "hadronic_stage", "observer_stage")
    mutable_fields: ClassVar[tuple[str, ...]] = ("d_n_gam_e_bh", "l_syn_spec", "seed_syn")

    gam_e: np.ndarray
    d_n_gam_e: np.ndarray
    l_syn_spec: np.ndarray
    seed_syn: np.ndarray
    nu_m: np.ndarray
    nu_c: np.ndarray
    nu_a: np.ndarray
    d_n_gam_e_bh: np.ndarray | None = None
    d_n_gam_e_chi: np.ndarray | None = None
    chi_grid: np.ndarray | None = None
    cooling_timescale_s: np.ndarray | None = None
    dynamical_timescale_s: np.ndarray | None = None
    work_x_edge_log10: np.ndarray | None = None
    work_d_n_x: np.ndarray | None = None


@dataclass
class ReverseShockEmission:
    """Emission from reverse shock."""
    l_syn_spec: np.ndarray
    seed_syn: np.ndarray
    magnetic_field_g: np.ndarray
    nu_m: np.ndarray
    nu_c: np.ndarray
    nu_a: np.ndarray
    nu_M: np.ndarray
    rs_hadronic: Optional[Any] = None


@dataclass
class HadronicSolution:
    """Solution from the 1d hadronic solver."""
    producer: ClassVar[str] = "hadronic_stage"
    consumers: ClassVar[tuple[str, ...]] = ("observer_stage", "api")
    mutable_fields: ClassVar[tuple[str, ...]] = ()

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
