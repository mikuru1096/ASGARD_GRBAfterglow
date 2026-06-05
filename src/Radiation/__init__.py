from __future__ import annotations

from importlib import import_module
from typing import Any


_BINDINGS = {
    "annihilation": "radiation_gamma_gamma_absorption",
    "cal_ebl": "Cal_ebl",
    "seed_reverse": "radiation_reverse_seed",
    "ssc_spec": "radiation_ssc_spectrum",
    "ssc_spec_nonuniform": "radiation_ssc_spectrum",
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
