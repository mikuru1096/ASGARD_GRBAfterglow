from __future__ import annotations

from importlib import import_module
from typing import Any


_BINDINGS = {
    "fs_hadronic_1d": "hadronic_forward_1d",
    "fs_hadronic_reverse_1d": "hadronic_reverse_1d",
}


def __getattr__(name: str) -> Any:
    module_name = _BINDINGS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f"{__name__}.{module_name}"), name)
    globals()[name] = value
    return value


__all__ = sorted(_BINDINGS)
__doc__ = "Hadronic solver bindings."
