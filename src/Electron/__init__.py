from __future__ import annotations

from importlib import import_module
from typing import Any


_SOLVERS = {
    "fs_charint_1d": "electron_forward_charint_1d",
    "fs_dg_1d": "electron_forward_dg_1d",
    "fs_fullhide_1d": "electron_forward_fullhide_1d",
    "fs_fullhide_hz": "electron_forward_fullhide_1d_hybrid",
    "fs_transport_2d": "electron_forward_transport_2d",
    "fs_slc1_1d": "electron_forward_slc1_1d",
    "fs_t2g1_1d": "electron_forward_t2g1_1d",
    "fs_weno5_1d": "electron_forward_weno5_1d",
    "electron_reverse_kernel": "electron_reverse_kernel",
}


def __getattr__(name: str) -> Any:
    module_name = _SOLVERS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(f"{__name__}.{module_name}")
    value = getattr(module, name)
    globals()[name] = value
    return value


__all__ = sorted(_SOLVERS)
__doc__ = "Electron solver bindings."
