from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from typing import Any

import numpy as np

from asgard_core.api_fit import Param
from asgard_core.api_model import Model
from asgard_core.api_observe import _as_model, _param_path, _total_matrix
from asgard_core.asgard_state import project_flux_grid, solve_state_from_setup
from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_observables import build_multiband_observer_frequencies, combine_multiband_flux
from asgard_core.asgard_postprocess import compute_light_curve_redchi
from asgard_core.asgard_setup import build_simulation_setup


@dataclass(frozen=True)
class ParamBinding:
    name: str
    target: Any
    attr_name: str
    original_value: Any
    scale: str
    fixed_value: float | None


@dataclass(frozen=True)
class FluxData:
    pair_mode: bool
    time_indices: np.ndarray
    freq_indices: np.ndarray
    flux: np.ndarray
    flux_err: np.ndarray


@dataclass(frozen=True)
class SpecData:
    time_index: int
    freq_indices: np.ndarray
    flux: np.ndarray
    flux_err: np.ndarray


@dataclass(frozen=True)
class BandData:
    time_index: int
    freq_indices: np.ndarray
    sample_frequencies_hz: np.ndarray
    flux: float
    flux_err: float


@dataclass(frozen=True)
class ObsBlock:
    observer_time_s: np.ndarray
    requested_frequencies_hz: np.ndarray
    flux_density: tuple[FluxData, ...]
    spectra: tuple[SpecData, ...]
    band_fluxes: tuple[BandData, ...]


@dataclass(frozen=True)
class ObsPlan:
    blocks: tuple[ObsBlock, ...]


@dataclass
class InferenceProblem:
    model: Model
    observations: ObsPlan
    param_bindings: tuple[ParamBinding, ...]


@dataclass(frozen=True)
class FitProblem:
    observer_time_s: np.ndarray
    requested_frequencies_hz: np.ndarray
    num_xrt: int


def compile_problem(
    data_or_config,
    model_or_config: Any | None = None,
    *,
    params: list[Param] | None = None,
):
    if isinstance(data_or_config, RuntimeConfig) and model_or_config is None:
        num_xrt, requested_frequencies_hz = build_multiband_observer_frequencies()
        setup = build_simulation_setup(data_or_config)
        return FitProblem(
            observer_time_s=np.asarray(setup.observer_time_s, dtype=float),
            requested_frequencies_hz=np.asarray(requested_frequencies_hz, dtype=float),
            num_xrt=num_xrt,
        )

    data = data_or_config if isinstance(data_or_config, dict) else None
    if data is None:
        raise TypeError("compile_problem expects either (config) or (dict, model_or_config).")
    model = model_or_config if isinstance(model_or_config, Model) else _as_model(model_or_config)
    param_bindings = []
    for param in [] if params is None else params:
        path = _param_path(model, param)
        target = model
        parts = path.split(".")
        for name in parts[:-1]:
            target = getattr(target, name)
        attr_name = parts[-1]
        param_bindings.append(
            ParamBinding(
                name=param.name,
                target=target,
                attr_name=attr_name,
                original_value=getattr(target, attr_name),
                scale=str(param.scale.value),
                fixed_value=None if not param.is_fixed() else float(param.lower),
            )
        )
    return InferenceProblem(
        model=model,
        observations=_compile_obs(data),
        param_bindings=tuple(param_bindings),
    )


def eval_loglike(
    problem,
    params_or_config,
    *,
    timings: dict[str, float] | None = None,
) -> float:
    if isinstance(problem, FitProblem):
        if not isinstance(params_or_config, RuntimeConfig):
            raise TypeError("FitProblem expects a RuntimeConfig input.")
        return _eval_cfg(problem, params_or_config, timings=timings)
    if isinstance(problem, InferenceProblem):
        if not isinstance(params_or_config, dict):
            raise TypeError("InferenceProblem expects a parameter dictionary.")
        return _eval_model(problem, params_or_config, timings=timings)
    raise TypeError(f"Unsupported compiled inference problem: {type(problem)!r}")


def _compile_obs(
    data: dict,
) -> ObsPlan:
    blocks: list[ObsBlock] = []
    for dataset in data["flux_density"]:
        times_s = np.asarray(dataset["times_s"], dtype=float)
        frequencies_hz = np.asarray(dataset["frequencies_hz"], dtype=float)
        observer_time_s = np.unique(times_s.reshape(-1))
        requested_frequencies_hz = np.unique(frequencies_hz.reshape(-1))
        blocks.append(
            ObsBlock(
                observer_time_s=observer_time_s,
                requested_frequencies_hz=requested_frequencies_hz,
                flux_density=(
                    FluxData(
                        pair_mode=times_s.shape == frequencies_hz.shape,
                        time_indices=_idx(observer_time_s, times_s.reshape(-1)),
                        freq_indices=_idx(requested_frequencies_hz, frequencies_hz.reshape(-1)),
                        flux=np.asarray(dataset["flux"], dtype=float),
                        flux_err=np.asarray(dataset["flux_err"], dtype=float),
                    ),
                ),
                spectra=(),
                band_fluxes=(),
            )
        )
    for dataset in data["spectrum"]:
        time_s = float(dataset["time_s"])
        frequencies_hz = np.asarray(dataset["frequencies_hz"], dtype=float)
        requested_frequencies_hz = np.unique(frequencies_hz)
        blocks.append(
            ObsBlock(
                observer_time_s=np.array([time_s], dtype=float),
                requested_frequencies_hz=requested_frequencies_hz,
                flux_density=(),
                spectra=(
                    SpecData(
                        time_index=0,
                        freq_indices=_idx(requested_frequencies_hz, frequencies_hz),
                        flux=np.asarray(dataset["flux"], dtype=float),
                        flux_err=np.asarray(dataset["flux_err"], dtype=float),
                    ),
                ),
                band_fluxes=(),
            )
        )
    for dataset in data["flux"]:
        time_s = float(dataset["time_s"])
        sample_frequencies_hz = np.logspace(
            np.log10(float(dataset["nu_min_hz"])),
            np.log10(float(dataset["nu_max_hz"])),
            int(dataset["num_points"]),
        )
        requested_frequencies_hz = np.unique(sample_frequencies_hz)
        blocks.append(
            ObsBlock(
                observer_time_s=np.array([time_s], dtype=float),
                requested_frequencies_hz=requested_frequencies_hz,
                flux_density=(),
                spectra=(),
                band_fluxes=(
                    BandData(
                        time_index=0,
                        freq_indices=_idx(requested_frequencies_hz, sample_frequencies_hz),
                        sample_frequencies_hz=sample_frequencies_hz,
                        flux=float(dataset["flux"]),
                        flux_err=float(dataset["flux_err"]),
                    ),
                ),
            )
        )
    if not blocks:
        blocks.append(
            ObsBlock(
                observer_time_s=np.array([1.0e2], dtype=float),
                requested_frequencies_hz=np.array([1.0e14], dtype=float),
                flux_density=(),
                spectra=(),
                band_fluxes=(),
            )
        )
    return ObsPlan(blocks=tuple(blocks))


def _idx(reference: np.ndarray, values: np.ndarray) -> np.ndarray:
    indices = np.searchsorted(reference, np.asarray(values, dtype=float))
    if np.any(indices >= reference.shape[0]):
        raise ValueError("Compiled inference request could not map values onto the unified grid.")
    if not np.allclose(reference[indices], values, rtol=0.0, atol=0.0):
        raise ValueError("Compiled inference request could not map values onto the unified grid.")
    return indices.astype(int, copy=False)


def _eval_cfg(
    problem: FitProblem,
    config: RuntimeConfig,
    *,
    timings: dict[str, float] | None = None,
) -> float:
    setup = build_simulation_setup(config)
    setup.observer_time_s = np.array(problem.observer_time_s, dtype=float, copy=True)
    state = solve_state_from_setup(
        config,
        setup,
        timings=timings,
        requested_frequencies_hz=problem.requested_frequencies_hz,
    )
    observed = project_flux_grid(
        state,
        problem.observer_time_s,
        problem.requested_frequencies_hz,
        timings=timings,
        mode="total_only",
        projection_kind="lightcurve",
    )
    band_flux_matrix = np.asarray(observed.components["total"], dtype=float)
    bands_flux = combine_multiband_flux(band_flux_matrix, problem.requested_frequencies_hz, problem.num_xrt)
    redchi = compute_light_curve_redchi(bands_flux, problem.observer_time_s, config)
    if np.isnan(redchi):
        redchi = np.inf
    return -0.5 * redchi


def _eval_model(
    problem: InferenceProblem,
    values: dict[str, float],
    *,
    timings: dict[str, float] | None = None,
) -> float:
    original_values: list[tuple[Any, str, Any]] = []
    try:
        for binding in problem.param_bindings:
            original_values.append((binding.target, binding.attr_name, getattr(binding.target, binding.attr_name)))
            raw_value = binding.fixed_value if binding.fixed_value is not None else float(values[binding.name])
            transformed = 10.0**raw_value if binding.scale in {"log", "log10"} else float(raw_value)
            setattr(binding.target, binding.attr_name, transformed)
        total_matrices = [
            _total_matrix(
                problem.model,
                block.observer_time_s,
                block.requested_frequencies_hz,
                timings=timings,
                projection_kind="lightcurve",
            )
            for block in problem.observations.blocks
        ]
    finally:
        for target, attr_name, original_value in reversed(original_values):
            setattr(target, attr_name, original_value)

    residual_start = perf_counter()
    loglike = 0.0
    for block, total_matrix in zip(problem.observations.blocks, total_matrices):
        for dataset in block.flux_density:
            if dataset.pair_mode:
                pred = total_matrix[dataset.freq_indices, dataset.time_indices]
            else:
                pred = total_matrix[np.ix_(dataset.freq_indices, dataset.time_indices)]
            resid = (pred - dataset.flux) / dataset.flux_err
            loglike -= 0.5 * float(np.sum(resid * resid))
        for dataset in block.spectra:
            pred = total_matrix[dataset.freq_indices, dataset.time_index]
            resid = (pred - dataset.flux) / dataset.flux_err
            loglike -= 0.5 * float(np.sum(resid * resid))
        for dataset in block.band_fluxes:
            pred = float(np.trapezoid(total_matrix[dataset.freq_indices, dataset.time_index], dataset.sample_frequencies_hz))
            resid = (pred - dataset.flux) / dataset.flux_err
            loglike -= 0.5 * float(resid * resid)
    if timings is not None:
        timings["Residual aggregation"] = timings.get("Residual aggregation", 0.0) + (perf_counter() - residual_start)
    return loglike
