from __future__ import annotations

from importlib import import_module
from typing import Any


_BINDINGS = {
    "sed_interpolation": "SED_interpolation",
    "sed_adaptive_theta": "SED_interpolation",
    "sed_interpolation_chi": "SED_interpolation",
    "sed_chi_electron": "SED_interpolation",
    "sed_chi_ring": "SED_interpolation",
    "sed_chiring_batchlum": "SED_interpolation",
    "sed_chiring_batchlum_ray": "SED_interpolation",
}


def __getattr__(name: str) -> Any:
    module_name = _BINDINGS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f"{__name__}.{module_name}"), name)
    globals()[name] = value
    return value


__all__ = sorted(_BINDINGS)
__doc__ = "Observer-frame SED interpolation bindings."
