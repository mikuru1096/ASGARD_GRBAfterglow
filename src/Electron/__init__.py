from __future__ import annotations

from importlib import import_module
from typing import Any


_SOLVERS = {
    "fs_electron_charint_1d": "FS_electron_charint_1d",
    "fs_electron_charint_2d": "FS_electron_charint_2d",
    "fs_electron_fullhide_1d": "FS_electron_fullhide_1d",
    "fs_electron_fullhide_2d": "FS_electron_fullhide_2d",
    "fs_electron_slc1_1d": "FS_electron_slc1_1d",
    "fs_electron_t2g1_1d": "FS_electron_t2g1_1d",
    "fs_electron_weno5_1d": "FS_electron_weno5_1d",
}
_ALIASES = {
    "fs_electron_charint": "fs_electron_charint_1d",
    "fs_electron_fullhide": "fs_electron_fullhide_1d",
    "fs_electron_slc1": "fs_electron_slc1_1d",
    "fs_electron_t2g1": "fs_electron_t2g1_1d",
    "fs_electron_weno5": "fs_electron_weno5_1d",
}


def __getattr__(name: str) -> Any:
    target = _ALIASES.get(name, name)
    module_name = _SOLVERS.get(target)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(f"{__name__}.{module_name}")
    value = getattr(module, target)
    globals()[name] = value
    return value


__all__ = sorted((*_SOLVERS, *_ALIASES))
__doc__ = "Electron solver bindings."
