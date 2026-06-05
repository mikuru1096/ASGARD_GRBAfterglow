from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from functools import cached_property
from typing import Any, Optional

import numpy as np

from .api_model import (
    Model,
    FluxResult,
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
        if self.scale in (Scale.LOG, Scale.LOG10):
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

ParamDef = Param


class Fitter:
    def __init__(
        self,
        model: Optional[Model] = None,
        data: Optional[dict] = None,
        params: Optional[list[ParamDef]] = None,
        cfg: Optional[Any] = None,
        num_workers: int = 1,
        **config_kwargs,
    ) -> None:
        if isinstance(model, dict) and data is None:
            data = model
            model = None
        if cfg is None and config_kwargs:
            cfg = config_kwargs
        from .api_observe import _as_model

        self.model = model if model is not None else _as_model(cfg)
        self.data = make_empty_obs() if data is None else data
        self.params = [] if params is None else list(params)
        self.num_workers = int(num_workers)

    @cached_property
    def _compiled_problem(self):
        from asgard_core.asgard_fit import compile_problem

        return compile_problem(self.data, self.model, params=self.params)

    def build_model(self, values: dict[str, float]) -> Model:
        model = deepcopy(self.model)
        from .api_observe import _param_path, _set_dotted_attr

        for param in self.params:
            if param.name in values or param.is_fixed():
                raw_value = param.lower if param.is_fixed() else values[param.name]
                path = _param_path(model, param)
                _set_dotted_attr(model, path, param.transform(raw_value))
        return model

    def loglike(self, values: dict[str, float]) -> float:
        from asgard_core.asgard_fit import eval_loglike

        return eval_loglike(self._compiled_problem, values)

    def flux_density_grid(self, values: dict[str, float], times_s: np.ndarray, nu_hz: np.ndarray) -> FluxResult:
        return self.build_model(values).flux_density_grid(times_s, nu_hz)

    def add_flux_density(self, times_s=None, frequencies_hz=None, flux=None, flux_err=None, **kwargs) -> None:
        entry = make_flux_density_entry(
            times_s if times_s is not None else kwargs.get("t"),
            frequencies_hz if frequencies_hz is not None else kwargs.get("nu"),
            flux if flux is not None else kwargs.get("f_nu"),
            flux_err if flux_err is not None else kwargs.get("err"),
        )
        self.data["flux_density"].append(entry)
        self.__dict__.pop("_compiled_problem", None)

    def add_spectrum(self, time_s=None, frequencies_hz=None, flux=None, flux_err=None, **kwargs) -> None:
        entry = make_spectrum_entry(
            time_s if time_s is not None else kwargs.get("t"),
            frequencies_hz if frequencies_hz is not None else kwargs.get("nu"),
            flux if flux is not None else kwargs.get("f_nu"),
            flux_err if flux_err is not None else kwargs.get("err"),
        )
        self.data["spectrum"].append(entry)
        self.__dict__.pop("_compiled_problem", None)

    def add_flux(self, nu_min_hz=None, nu_max_hz=None, time_s=None, flux=None, flux_err=None, num_points=64, **kwargs) -> None:
        entry = make_flux_entry(
            nu_min_hz if nu_min_hz is not None else kwargs.get("nu_min"),
            nu_max_hz if nu_max_hz is not None else kwargs.get("nu_max"),
            time_s if time_s is not None else kwargs.get("t"),
            flux if flux is not None else kwargs.get("value"),
            flux_err if flux_err is not None else kwargs.get("err"),
            int(num_points),
        )
        self.data["flux"].append(entry)
        self.__dict__.pop("_compiled_problem", None)

    def run_emcee(self, initial: np.ndarray, nwalkers: int, nsteps: int) -> FitResult:
        import emcee

        ndim = len(self.params)

        def _lnprob(theta):
            trial = {param.name: theta[i] for i, param in enumerate(self.params)}
            return self.loglike(trial)

        sampler = emcee.EnsembleSampler(nwalkers, ndim, _lnprob)
        sampler.run_mcmc(initial, nsteps, progress=False)
        flat = sampler.get_chain(flat=True)
        log_prob = sampler.get_log_prob(flat=True)
        best_idx = int(np.argmax(log_prob))
        best = {param.name: flat[best_idx, i] for i, param in enumerate(self.params)}
        return FitResult(
            best_params=best,
            best_loglike=float(log_prob[best_idx]),
            samples=flat,
            log_prob=log_prob,
            labels=[param.name for param in self.params],
            top_k_params=[best],
        )

    def run_multinest(self, outputfiles_basename: str, verbose: bool = True) -> FitResult:
        from pymultinest.solve import solve

        n_dims = len(self.params)

        def _prior_transform(cube):
            params = np.asarray(cube, dtype=float).copy()
            for i, param in enumerate(self.params):
                params[i] = param.lower + params[i] * (param.upper - param.lower)
            return params

        def _loglike(theta):
            trial = {param.name: theta[i] for i, param in enumerate(self.params)}
            return self.loglike(trial)

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
        param_defs: Optional[list[ParamDef]] = None,
        *,
        sampler: str = "emcee",
        total_steps: int = 128,
        burn_frac: float = 0.5,
        thin: int = 1,
        nwalkers: Optional[int] = None,
        outputfiles_basename: str = "chains/vegas-",
        verbose: bool = False,
        nsteps: Optional[int] = None,
        nburn: Optional[int] = None,
        npool: Optional[int] = None,
        resolution: Optional[tuple[float, float, int]] = None,
    ) -> FitResult:
        if param_defs is not None:
            self.params = list(param_defs)
            self._compiled_problem = None
        if not self.params:
            raise ValueError("No parameter definitions were provided to Fitter.fit().")
        if nsteps is not None:
            total_steps = int(nsteps)
        if nburn is not None:
            burn_frac = float(nburn) / float(total_steps)
        if npool is not None:
            self.num_workers = int(npool)
        if resolution is not None:
            self.model._apply_resolutions(resolution)
            self._compiled_problem = None

        active_params = [param for param in self.params if not param.is_fixed()]
        if sampler.lower() in ("multinest", "pymultinest"):
            if active_params != self.params:
                raise NotImplementedError("run_multinest currently requires all fitted parameters to be active.")
            return self.run_multinest(outputfiles_basename=outputfiles_basename, verbose=verbose)

        if sampler.lower() not in ("emcee", "mcmc"):
            raise ValueError(f"Unsupported sampler: {sampler}")
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
                if param.name not in trial:
                    continue
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
