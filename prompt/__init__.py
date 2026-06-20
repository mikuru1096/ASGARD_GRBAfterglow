"""Prompt-emission experimental branch for magnetized GRB internal shocks."""

from .internal_shock import InternalShockNumerics, InternalShockShell, fast_shock_allowed, simulate_internal_shock
from .radiation import InternalShockMicrophysics

__all__ = [
    "InternalShockMicrophysics",
    "InternalShockNumerics",
    "InternalShockShell",
    "fast_shock_allowed",
    "simulate_internal_shock",
]
