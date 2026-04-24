from __future__ import annotations

from importlib import import_module
from typing import Any


_BINDINGS = {
    "annihilation": "Annihilation",
    "cal_ebl": "Cal_ebl",
    "seed_reverse": "Seed_reverse",
    "ssc_spec": "SSC_spec",
    "ssc_spec_nonuniform": "SSC_spec",
}


def __getattr__(name: str) -> Any:
    module_name = _BINDINGS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(f"{__name__}.{module_name}"), name)
    globals()[name] = value
    return value


__all__ = sorted(_BINDINGS)
__doc__ = "Radiation process bindings."
