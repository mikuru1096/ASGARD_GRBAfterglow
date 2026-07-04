from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from functools import cached_property
from time import perf_counter
from typing import Any, Optional

import numpy as np

from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_postprocess import (
    band_freqs,
    combine_flux,
    light_chi,
)
from asgard_core.asgard_setup import build_simulation_setup
from asgard_core.asgard_state import project_flux_grid, solve_state_from_setup
from .api_model import (
    Model,
    Scale,
    make_empty_obs,
    make_flux_density_entry,
    make_spectrum_entry,
    make_flux_entry,
)

class Param:
    def __init__(self, name: str, *args, scale: Scale = Scale.LINEAR) -> None:
        self.name = name
        self.path: Optional[str]
        if len(args) >= 3 and isinstance(args[0], str):
            self.path = args[0]
            self.lower = float(args[1])
            self.upper = float(args[2])
            if len(args) >= 4:
                scale = args[3]
        elif len(args) >= 2:
            self.path = None
            self.lower = float(args[0])
            self.upper = float(args[1])
            if len(args) >= 3:
                scale = args[2]
        else:
            raise TypeError("Param requires either (name, path, lower, upper, scale) or (name, lower, upper, scale).")
        self.scale = Scale(scale)

    def transform(self, value: float) -> float:
        if self.scale == Scale.LOG10:
            return 10.0**value
        if self.scale == Scale.FIXED:
            return float(self.lower)
        return value

    def is_fixed(self) -> bool:
        return self.scale == Scale.FIXED or np.isclose(self.lower, self.upper)


@dataclass
class FitResult:
    best_params: dict[str, float]
    best_loglike: float
    samples: Optional[np.ndarray] = None
    log_prob: Optional[np.ndarray] = None
    logz: Optional[float] = None
    logzerr: Optional[float] = None
    labels: Optional[list[str]] = None
    top_k_params: Optional[list[dict[str, float]]] = None
class Fitter:
    def __init__(
        self,
        model: Model,
        data: Optional[dict] = None,
        params: Optional[list[Param]] = None,
        num_workers: int = 1,
    ) -> None:
        self.model = model
        self.data = make_empty_obs() if data is None else data
        self.params = [] if params is None else list(params)
        self.num_workers = int(num_workers)

    @cached_property
    def _compiled_problem(self):
        return compile_problem(self.data, self.model, params=self.params)

    def build_model(self, values: dict[str, float]) -> Model:
        model = deepcopy(self.model)
        from .api_observe import _parampath

        for param in self.params:
            if param.name in values or param.is_fixed():
                raw_value = param.lower if param.is_fixed() else values[param.name]
                parts = _parampath(model, param).split(".")
                target = model
                for name in parts[:-1]:
                    target = getattr(target, name)
                setattr(target, parts[-1], param.transform(raw_value))
        return model

    def loglike(self, values: dict[str, float]) -> float:
        return eval_loglike(self._compiled_problem, values)

    def add_flux_density(self, times_s, frequencies_hz, flux, flux_err) -> None:
        entry = make_flux_density_entry(times_s, frequencies_hz, flux, flux_err)
        self.data["flux_density"].append(entry)
        self.__dict__.pop("_compiled_problem", None)

    def add_spectrum(self, time_s, frequencies_hz, flux, flux_err) -> None:
        entry = make_spectrum_entry(time_s, frequencies_hz, flux, flux_err)
        self.data["spectrum"].append(entry)
        self.__dict__.pop("_compiled_problem", None)

    def add_flux(self, nu_min_hz, nu_max_hz, time_s, flux, flux_err, num_points=64) -> None:
        entry = make_flux_entry(nu_min_hz, nu_max_hz, time_s, flux, flux_err, int(num_points))
        self.data["flux"].append(entry)
        self.__dict__.pop("_compiled_problem", None)

    def _run_multinest(self, outputfiles_basename: str, verbose: bool = True) -> FitResult:
        from pymultinest.solve import solve

        n_dims = len(self.params)

        def _prior_transform(cube):
            params = np.asarray(cube, dtype=float).copy()
            for i, param in enumerate(self.params):
                params[i] = param.lower + params[i] * (param.upper - param.lower)
            return params

        def _loglike(theta):
            return self.loglike({param.name: theta[i] for i, param in enumerate(self.params)})

        result = solve(
            LogLikelihood=_loglike,
            Prior=_prior_transform,
            n_dims=n_dims,
            outputfiles_basename=outputfiles_basename,
            verbose=verbose,
        )
        samples = np.asarray(result["samples"], dtype=float)
        best_idx = int(np.argmax(result["log_likelihood"]))
        best = {param.name: float(samples[best_idx, i]) for i, param in enumerate(self.params)}
        return FitResult(
            best_params=best,
            best_loglike=float(result["log_likelihood"][best_idx]),
            samples=samples,
            log_prob=np.asarray(result["log_likelihood"], dtype=float),
            logz=float(result["logZ"]),
            logzerr=float(result["logZerr"]),
            labels=[param.name for param in self.params],
            top_k_params=[best],
        )

    def fit(
        self,
        param_defs: Optional[list[Param]] = None,
        *,
        sampler: str = "emcee",
        total_steps: int = 128,
        burn_frac: float = 0.5,
        thin: int = 1,
        nwalkers: Optional[int] = None,
        outputfiles_basename: str = "chains/asgard-",
        verbose: bool = False,
    ) -> FitResult:
        if param_defs is not None:
            self.params = list(param_defs)
            self.__dict__.pop("_compiled_problem", None)
        if not self.params:
            raise ValueError("No parameter definitions were provided to Fitter.fit().")

        active_params = [param for param in self.params if not param.is_fixed()]
        if sampler.lower() == "pymultinest":
            if active_params != self.params:
                raise NotImplementedError("sampler='pymultinest' currently requires all fitted parameters to be active.")
            return self._run_multinest(outputfiles_basename=outputfiles_basename, verbose=verbose)

        if sampler.lower() != "emcee":
            raise ValueError("sampler must be 'emcee' or 'pymultinest'.")
        if not active_params:
            best = {param.name: param.lower for param in self.params}
            best_loglike = self.loglike(best)
            return FitResult(best_params=best, best_loglike=best_loglike, labels=[param.name for param in self.params], top_k_params=[best])

        import emcee

        ndim = len(active_params)
        walker_num = max(2 * ndim + 2, 8) if nwalkers is None else int(nwalkers)
        initial = np.zeros((walker_num, ndim), dtype=float)
        for i, param in enumerate(active_params):
            span = param.upper - param.lower
            if span <= 0.0:
                initial[:, i] = param.lower
            else:
                initial[:, i] = param.lower + span * np.random.default_rng(1234 + i).uniform(0.2, 0.8, size=walker_num)

        def _lnprob(theta):
            trial = {param.name: theta[i] for i, param in enumerate(active_params)}
            for param in active_params:
                if trial[param.name] < param.lower or trial[param.name] > param.upper:
                    return -np.inf
            return self.loglike(trial)

        sampler_obj = emcee.EnsembleSampler(walker_num, ndim, _lnprob)
        sampler_obj.run_mcmc(initial, int(total_steps), progress=False)
        chain = sampler_obj.get_chain(flat=True)
        log_prob = sampler_obj.get_log_prob(flat=True)
        burn = int(burn_frac * chain.shape[0])
        flat = chain[burn::thin]
        flat_log_prob = log_prob[burn::thin]
        best_idx = int(np.argmax(flat_log_prob))
        best = {param.name: float(flat[best_idx, i]) for i, param in enumerate(active_params)}
        for param in self.params:
            if param.is_fixed():
                best[param.name] = param.lower
        return FitResult(
            best_params=best,
            best_loglike=float(flat_log_prob[best_idx]),
            samples=flat,
            log_prob=flat_log_prob,
            labels=[param.name for param in active_params],
            top_k_params=[best],
        )


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
        num_xrt, requested_frequencies_hz = band_freqs()
        setup = build_simulation_setup(data_or_config)
        return FitProblem(
            observer_time_s=np.asarray(setup.observer_time_s, dtype=float),
            requested_frequencies_hz=np.asarray(requested_frequencies_hz, dtype=float),
            num_xrt=num_xrt,
        )

    data = data_or_config if isinstance(data_or_config, dict) else None
    if data is None:
        raise TypeError("compile_problem expects either (config) or (dict, model_or_config).")
    from .api_observe import _asmodel, _parampath

    model = model_or_config if isinstance(model_or_config, Model) else _asmodel(model_or_config)
    param_bindings = []
    for param in [] if params is None else params:
        path = _parampath(model, param)
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


def _compile_obs(data: dict) -> ObsPlan:
    blocks: list[ObsBlock] = []
    for dataset in data["flux_density"]:
        times_s = _positive("times_s", dataset["times_s"])
        frequencies_hz = _positive("frequencies_hz", dataset["frequencies_hz"])
        flux = _finite("flux", dataset["flux"])
        fluxerr = _positive("flux_err", dataset["flux_err"])
        pair_mode = times_s.shape == frequencies_hz.shape
        if pair_mode:
            flux = flux.reshape(-1)
            fluxerr = fluxerr.reshape(-1)
        expected = times_s.reshape(-1).shape if pair_mode else (frequencies_hz.size, times_s.size)
        if flux.shape != expected or fluxerr.shape != expected:
            raise ValueError("flux_density flux and flux_err must match the implied observation shape.")
        observer_time_s = np.unique(times_s.reshape(-1))
        requested_frequencies_hz = np.unique(frequencies_hz.reshape(-1))
        blocks.append(
            ObsBlock(
                observer_time_s=observer_time_s,
                requested_frequencies_hz=requested_frequencies_hz,
                flux_density=(
                    FluxData(
                        pair_mode=pair_mode,
                        time_indices=_idx(observer_time_s, times_s.reshape(-1)),
                        freq_indices=_idx(requested_frequencies_hz, frequencies_hz.reshape(-1)),
                        flux=flux,
                        flux_err=fluxerr,
                    ),
                ),
                spectra=(),
                band_fluxes=(),
            )
        )
    for dataset in data["spectrum"]:
        time_s = float(_positive("time_s", [dataset["time_s"]])[0])
        frequencies_hz = _positive("frequencies_hz", dataset["frequencies_hz"])
        flux = _finite("flux", dataset["flux"]).reshape(-1)
        fluxerr = _positive("flux_err", dataset["flux_err"]).reshape(-1)
        if flux.shape != frequencies_hz.reshape(-1).shape or fluxerr.shape != frequencies_hz.reshape(-1).shape:
            raise ValueError("spectrum flux and flux_err must match frequencies_hz.")
        frequencies_hz = frequencies_hz.reshape(-1)
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
                        flux=flux,
                        flux_err=fluxerr,
                    ),
                ),
                band_fluxes=(),
            )
        )
    for dataset in data["flux"]:
        time_s = float(_positive("time_s", [dataset["time_s"]])[0])
        nu_min = float(_positive("nu_min_hz", [dataset["nu_min_hz"]])[0])
        nu_max = float(_positive("nu_max_hz", [dataset["nu_max_hz"]])[0])
        points = int(dataset["num_points"])
        if nu_max <= nu_min:
            raise ValueError("nu_max_hz must exceed nu_min_hz.")
        if points < 2:
            raise ValueError("num_points must be at least 2 for band flux integration.")
        flux = float(_finite("flux", [dataset["flux"]])[0])
        fluxerr = float(_positive("flux_err", [dataset["flux_err"]])[0])
        sample_frequencies_hz = np.logspace(
            np.log10(nu_min),
            np.log10(nu_max),
            points,
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
                        flux=flux,
                        flux_err=fluxerr,
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


def _finite(name: str, values) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if np.any(~np.isfinite(array)):
        raise ValueError(f"{name} must contain finite values.")
    return array


def _positive(name: str, values) -> np.ndarray:
    array = _finite(name, values)
    if np.any(array <= 0.0):
        raise ValueError(f"{name} must contain positive values.")
    return array


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
    bands_flux = combine_flux(band_flux_matrix, problem.requested_frequencies_hz, problem.num_xrt)
    chi2 = light_chi(bands_flux, problem.observer_time_s, config)
    return -0.5 * chi2


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
            problem.model._total_matrix(
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
