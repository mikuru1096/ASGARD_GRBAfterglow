from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from typing import Any, Optional

import numpy as np

from ASGARD.api_fit import ParamDef
from ASGARD.api_model import Model
from ASGARD.api_observe import _as_model, _param_path, _total_matrix
from asgard_core.asgard_numpy import trapezoid
from asgard_core.asgard_state import project_flux_grid, solve_state_from_setup
from asgard_core.asgard_config import FitConfig
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
    fixed_value: Optional[float]


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
    params: Optional[list[ParamDef]] = None,
):
    if isinstance(data_or_config, FitConfig) and model_or_config is None:
        return _compile_cfg(data_or_config)

    data = data_or_config if isinstance(data_or_config, dict) else None
    if data is None:
        raise TypeError("compile_problem expects either (config) or (dict, model_or_config).")
    model = model_or_config if isinstance(model_or_config, Model) else _as_model(model_or_config)
    return _compile_model(data, model, params=params)


def eval_loglike(
    problem,
    params_or_config,
    *,
    timings: Optional[dict[str, float]] = None,
) -> float:
    if isinstance(problem, FitProblem):
        if not isinstance(params_or_config, FitConfig):
            raise TypeError("FitProblem expects a FitConfig input.")
        return _eval_cfg(problem, params_or_config, timings=timings)
    if isinstance(problem, InferenceProblem):
        if not isinstance(params_or_config, dict):
            raise TypeError("InferenceProblem expects a parameter dictionary.")
        return _eval_model(problem, params_or_config, timings=timings)
    raise TypeError(f"Unsupported compiled inference problem: {type(problem)!r}")


def _compile_cfg(config: FitConfig) -> FitProblem:
    num_xrt, requested_frequencies_hz = build_multiband_observer_frequencies()
    setup = build_simulation_setup(config)
    return FitProblem(
        observer_time_s=np.asarray(setup.observer_time_s, dtype=float),
        requested_frequencies_hz=np.asarray(requested_frequencies_hz, dtype=float),
        num_xrt=num_xrt,
    )


def _compile_model(
    data: dict,
    model: Model,
    *,
    params: Optional[list[ParamDef]] = None,
) -> InferenceProblem:
    observations = _compile_obs(data)
    bindings = _compile_params(model, [] if params is None else params)
    return InferenceProblem(
        model=model,
        observations=observations,
        param_bindings=bindings,
    )


def _compile_obs(
    data: dict,
) -> ObsPlan:
    blocks: list[ObsBlock] = []
    for dataset in data["flux_density"]:
        times_s = np.asarray(dataset["times_s"], dtype=float)
        frequencies_hz = np.asarray(dataset["frequencies_hz"], dtype=float)
        blocks.append(
            _make_flux_block(
                np.unique(times_s.reshape(-1)),
                np.unique(frequencies_hz.reshape(-1)),
                times_s,
                frequencies_hz,
                np.asarray(dataset["flux"], dtype=float),
                np.asarray(dataset["flux_err"], dtype=float),
            )
        )
    for dataset in data["spectrum"]:
        time_s = float(dataset["time_s"])
        frequencies_hz = np.asarray(dataset["frequencies_hz"], dtype=float)
        blocks.append(
            _make_spec_block(
                np.array([time_s], dtype=float),
                np.unique(frequencies_hz),
                frequencies_hz,
                np.asarray(dataset["flux"], dtype=float),
                np.asarray(dataset["flux_err"], dtype=float),
            )
        )
    for dataset in data["flux"]:
        time_s = float(dataset["time_s"])
        sample_frequencies_hz = np.logspace(
            np.log10(float(dataset["nu_min_hz"])),
            np.log10(float(dataset["nu_max_hz"])),
            int(dataset["num_points"]),
        )
        blocks.append(
            _make_band_block(
                np.array([time_s], dtype=float),
                np.unique(sample_frequencies_hz),
                sample_frequencies_hz,
                float(dataset["flux"]),
                float(dataset["flux_err"]),
            )
        )
    if not blocks:
        blocks.append(_make_default_block())
    return ObsPlan(blocks=tuple(blocks))


def _make_flux_block(
    observer_time_s: np.ndarray,
    requested_frequencies_hz: np.ndarray,
    times_s: np.ndarray,
    frequencies_hz: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray,
) -> ObsBlock:
    return ObsBlock(
        observer_time_s=observer_time_s,
        requested_frequencies_hz=requested_frequencies_hz,
        flux_density=(
            FluxData(
                pair_mode=times_s.shape == frequencies_hz.shape,
                time_indices=_idx(observer_time_s, times_s.reshape(-1)),
                freq_indices=_idx(requested_frequencies_hz, frequencies_hz.reshape(-1)),
                flux=flux,
                flux_err=flux_err,
            ),
        ),
        spectra=(),
        band_fluxes=(),
    )


def _make_spec_block(
    observer_time_s: np.ndarray,
    requested_frequencies_hz: np.ndarray,
    frequencies_hz: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray,
) -> ObsBlock:
    return ObsBlock(
        observer_time_s=observer_time_s,
        requested_frequencies_hz=requested_frequencies_hz,
        flux_density=(),
        spectra=(
            SpecData(
                time_index=0,
                freq_indices=_idx(requested_frequencies_hz, frequencies_hz),
                flux=flux,
                flux_err=flux_err,
            ),
        ),
        band_fluxes=(),
    )


def _make_band_block(
    observer_time_s: np.ndarray,
    requested_frequencies_hz: np.ndarray,
    sample_frequencies_hz: np.ndarray,
    flux: float,
    flux_err: float,
) -> ObsBlock:
    return ObsBlock(
        observer_time_s=observer_time_s,
        requested_frequencies_hz=requested_frequencies_hz,
        flux_density=(),
        spectra=(),
        band_fluxes=(
            BandData(
                time_index=0,
                freq_indices=_idx(requested_frequencies_hz, sample_frequencies_hz),
                sample_frequencies_hz=sample_frequencies_hz,
                flux=flux,
                flux_err=flux_err,
            ),
        ),
    )


def _make_default_block() -> ObsBlock:
    return ObsBlock(
        observer_time_s=np.array([1.0e2], dtype=float),
        requested_frequencies_hz=np.array([1.0e14], dtype=float),
        flux_density=(),
        spectra=(),
        band_fluxes=(),
    )


def _compile_params(model: Model, params: list[ParamDef]) -> tuple[ParamBinding, ...]:
    bindings: list[ParamBinding] = []
    for param in params:
        path = _param_path(model, param)
        target, attr_name = _find_attr(model, path)
        bindings.append(
            ParamBinding(
                name=param.name,
                target=target,
                attr_name=attr_name,
                original_value=getattr(target, attr_name),
                scale=str(param.scale.value),
                fixed_value=None if not param.is_fixed() else float(param.lower),
            )
        )
    return tuple(bindings)


def _find_attr(root: Any, path: str) -> tuple[Any, str]:
    target = root
    parts = path.split(".")
    for name in parts[:-1]:
        target = getattr(target, name)
    return target, parts[-1]


def _idx(reference: np.ndarray, values: np.ndarray) -> np.ndarray:
    indices = np.searchsorted(reference, np.asarray(values, dtype=float))
    if np.any(indices >= reference.shape[0]):
        raise ValueError("Compiled inference request could not map values onto the unified grid.")
    if not np.allclose(reference[indices], values, rtol=0.0, atol=0.0):
        raise ValueError("Compiled inference request could not map values onto the unified grid.")
    return indices.astype(int, copy=False)


def _eval_cfg(
    problem: FitProblem,
    config: FitConfig,
    *,
    timings: Optional[dict[str, float]] = None,
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
    timings: Optional[dict[str, float]] = None,
) -> float:
    original_values: list[tuple[Any, str, Any]] = []
    try:
        for binding in problem.param_bindings:
            original_values.append((binding.target, binding.attr_name, getattr(binding.target, binding.attr_name)))
            transformed = _bind_value(binding, values)
            setattr(binding.target, binding.attr_name, transformed)
        total_matrices = [
            _total_matrix(
                problem.model,
                block.observer_time_s,
                block.requested_frequencies_hz,
                timings=timings,
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
            pred = float(trapezoid(total_matrix[dataset.freq_indices, dataset.time_index], dataset.sample_frequencies_hz))
            resid = (pred - dataset.flux) / dataset.flux_err
            loglike -= 0.5 * float(resid * resid)
    if timings is not None:
        timings["Residual aggregation"] = timings.get("Residual aggregation", 0.0) + (perf_counter() - residual_start)
    return loglike


def _bind_value(binding: ParamBinding, values: dict[str, float]) -> float:
    raw_value = binding.fixed_value if binding.fixed_value is not None else float(values[binding.name])
    if binding.scale in {"log", "log10"}:
        return 10.0**raw_value
    return float(raw_value)
