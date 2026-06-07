from __future__ import annotations

from importlib import import_module
from typing import Any


_BINDINGS = {
    "structured_jet_flux_1d": "structured_jet_1d",
}


def __getattr__(name: str) -> Any:
    module_name = _BINDINGS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f"{__name__}.{module_name}"), name)
    globals()[name] = value
    return value


__all__ = sorted(_BINDINGS)
__doc__ = "Structured jet Fortran aggregate bindings."
