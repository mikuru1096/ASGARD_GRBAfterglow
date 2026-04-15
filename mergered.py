from __future__ import annotations

import warnings

from asgard_legacy import legacy_kwargs_to_config
from asgard_observables import (
    FITTING_BANDS,
    FITTING_FREQUENCIES_HZ,
    OUTPUT_BANDS,
    ZEROPOINT_FLUX,
)
from asgard_models import (
    FitConfig,
    FitResult,
    ReverseShockConfig,
    SpectrumOutputConfig,
)
from asgard_plot import plot_characteristic_frequencies, plot_light_curve, plot_spectrum
from asgard_solver import run_fit


def fit(*args, **kwargs) -> float:
    if args:
        raise TypeError("fit() no longer accepts positional arguments. Use run_fit(FitConfig).")
    warnings.warn("fit(**kwargs) is deprecated. Use run_fit(FitConfig) instead.", DeprecationWarning, stacklevel=2)
    config = legacy_kwargs_to_config(kwargs)
    return run_fit(config).redchi


__all__ = [
    "FITTING_BANDS",
    "OUTPUT_BANDS",
    "FITTING_FREQUENCIES_HZ",
    "ZEROPOINT_FLUX",
    "SpectrumOutputConfig",
    "ReverseShockConfig",
    "FitConfig",
    "FitResult",
    "run_fit",
    "plot_light_curve",
    "plot_spectrum",
    "plot_characteristic_frequencies",
    "fit",
]
