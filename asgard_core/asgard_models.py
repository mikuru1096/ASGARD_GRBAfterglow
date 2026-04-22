"""
DEPRECATED: This module is deprecated. Use asgard_config instead.

This module is kept for backward compatibility during migration.
All new code should import from asgard_config.
"""

from __future__ import annotations

# Re-export everything from asgard_config for backward compatibility
from asgard_core.asgard_config import (
    ExecutionPolicy,
    FitConfig,
    NumericalConfig,
    OutputConfig,
    PhysicsConfig,
    ReverseShockConfig,
    SimulationConfig,
    SpectrumOutputConfig,
    default_num_threads,
)

import os
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from src import constants


# Keep only the result/setup dataclasses that aren't in asgard_config
@dataclass
class PhysicalSolution:
    characteristic_time_s: np.ndarray
    gamma: np.ndarray
    radius_cm: np.ndarray
    absorbed_spectral_flux: np.ndarray
    nu_m: np.ndarray
    nu_c: np.ndarray
    nu_a: np.ndarray
    rs_nu_m: Optional[np.ndarray] = None
    rs_nu_c: Optional[np.ndarray] = None
    rs_nu_a: Optional[np.ndarray] = None


@dataclass
class SimulationSetup:
    luminosity_distance_cm: float
    boundary: np.ndarray
    seed_frequency_hz: np.ndarray
    observer_time_s: np.ndarray


@dataclass
class FitResult:
    t_obs_s: np.ndarray
    characteristic_time_s: np.ndarray
    bands: tuple[str, ...]
    bands_flux: np.ndarray
    redchi: float
    nu_m: np.ndarray
    nu_c: np.ndarray
    nu_a: np.ndarray
    rs_nu_m: Optional[np.ndarray] = None
    rs_nu_c: Optional[np.ndarray] = None
    rs_nu_a: Optional[np.ndarray] = None
    spectrum_time_s: Optional[float] = None
    spectrum_freq_hz: Optional[np.ndarray] = None
    spectrum_fnu: Optional[np.ndarray] = None
    spectrum_redchi: Optional[float] = None


__all__ = [
    "default_num_threads",
    "SpectrumOutputConfig",
    "ReverseShockConfig",
    "ExecutionPolicy",
    "PhysicsConfig",
    "NumericalConfig",
    "OutputConfig",
    "SimulationConfig",
    "FitConfig",
    "PhysicalSolution",
    "SimulationSetup",
    "FitResult",
]
