from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from functools import lru_cache
from typing import Optional

import numpy as np

from asgard_core.asgard_models import FitConfig
from asgard_core.asgard_state import SolveState, project_flux_grid

from .api_model import FluxPair, Model, FluxResult

def _array_signature(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values, dtype=float))
    digest = hashlib.blake2b(digest_size=16)
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
    digest.update(array.view(np.uint8))
    return digest.hexdigest()


def _remember_cache_entry(cache: dict, key, value, max_items: int = 8) -> None:
    cache[key] = value
    if len(cache) > max_items:
        cache.pop(next(iter(cache)))


def _fit_config_signature(config: FitConfig) -> str:
    payload = asdict(config)
    payload.pop("num_tobs", None)
    payload.pop("t_obs_min_log10", None)
    payload.pop("t_obs_max_log10", None)
    digest = hashlib.blake2b(digest_size=16)
    digest.update(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8"))
    return digest.hexdigest()


def _model_signature(model: Model) -> str:
    payload = {
        "medium": {
            "type": type(model.medium).__name__,
            "label": getattr(model.medium, "label", None),
            "data": {k: v for k, v in vars(model.medium).items() if k != "rho"},
        },
        "jet": {"type": type(model.jet).__name__, "data": vars(model.jet)},
        "observer": vars(model.observer),
        "fwd_rad": vars(model.fwd_rad),
        "rvs_rad": None if model.rvs_rad is None else vars(model.rvs_rad),
        "setups": asdict(model.setups),
    }
    digest = hashlib.blake2b(digest_size=16)
    digest.update(json.dumps(payload, sort_keys=True, default=float, separators=(",", ":")).encode("utf-8"))
    return digest.hexdigest()


def _pack_flux(observed: dict[str, np.ndarray | None]) -> FluxResult:
    total = np.asarray(observed["total"], dtype=float)
    fwd_sync = np.zeros_like(total) if observed["fwd_sync"] is None else np.asarray(observed["fwd_sync"], dtype=float)
    fwd_ssc = np.zeros_like(total) if observed["fwd_ssc"] is None else np.asarray(observed["fwd_ssc"], dtype=float)
    rev_sync = np.zeros_like(total) if observed["rev_sync"] is None else np.asarray(observed["rev_sync"], dtype=float)
    rev_ssc = np.zeros_like(total) if observed["rev_ssc"] is None else np.asarray(observed["rev_ssc"], dtype=float)
    cross_ic = None if observed["cross_ic"] is None else np.asarray(observed["cross_ic"], dtype=float)
    return FluxResult(
        total=total,
        fwd=FluxPair(sync=fwd_sync, ssc=fwd_ssc),
        rev=FluxPair(sync=rev_sync, ssc=rev_ssc),
        cross_ic=cross_ic,
    )


def _observe_parts(
    state: SolveState,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    mode: str = "full_components",
) -> FluxResult:
    observed_state = project_flux_grid(state, times_s, nu_hz, mode=mode)
    return _pack_flux(observed_state.components)


def _observe_total(
    state: SolveState,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    timings: Optional[dict[str, float]] = None,
) -> np.ndarray:
    observed_state = project_flux_grid(state, times_s, nu_hz, timings=timings, mode="total_only")
    return np.asarray(observed_state.components["total"], dtype=float)


def _log_time_midpoint(t_left: float, t_right: float) -> float:
    if t_left <= 0.0 or t_right <= 0.0:
        return 0.5 * (t_left + t_right)
    return float(np.sqrt(t_left * t_right))


def _interpolate_segment_value(
    t_left: float,
    f_left: float,
    t_right: float,
    f_right: float,
    t_mid: float,
) -> float:
    if f_left > 0.0 and f_right > 0.0 and t_left > 0.0 and t_right > 0.0 and t_mid > 0.0:
        slope = (np.log(f_right) - np.log(f_left)) / (np.log(t_right) - np.log(t_left))
        return float(np.exp(np.log(f_left) + slope * (np.log(t_mid) - np.log(t_left))))
    weight = (t_mid - t_left) / (t_right - t_left)
    return float(f_left + weight * (f_right - f_left))


def _relative_segment_error(actual: float, estimate: float) -> float:
    scale = max(abs(actual), 1.0e-99)
    return abs(actual - estimate) / scale


def _should_use_adaptive_flux_grid(times_s: np.ndarray, nu_hz: np.ndarray) -> bool:
    times_s = np.asarray(times_s, dtype=float)
    nu_hz = np.asarray(nu_hz, dtype=float)
    if times_s.ndim != 1 or nu_hz.ndim != 1:
        return False
    if times_s.size > 96:
        return False
    if times_s.size < 48 or nu_hz.size == 0 or nu_hz.size > 8:
        return False
    # Ultra-high-energy bands are sensitive to cutoff-region curvature.
    # Use direct evaluation to avoid interpolation-induced hard truncation.
    if float(np.max(nu_hz)) > 1.0e28:
        return False
    if np.any(~np.isfinite(times_s)) or np.any(times_s <= 0.0):
        return False
    return float(np.max(times_s)) / float(np.min(times_s)) >= 2.0


def _should_use_adaptive_spectrum_grid(times_s: np.ndarray, nu_hz: np.ndarray) -> bool:
    times_s = np.asarray(times_s, dtype=float)
    nu_hz = np.asarray(nu_hz, dtype=float)
    if times_s.ndim != 1 or nu_hz.ndim != 1:
        return False
    if nu_hz.size > 256:
        return False
    if times_s.size != 1 or nu_hz.size < 64:
        return False
    if np.any(~np.isfinite(nu_hz)) or np.any(nu_hz <= 0.0):
        return False
    return float(np.max(nu_hz)) / float(np.min(nu_hz)) >= 10.0


def _adaptive_segment_error(
    t_left: float,
    t_mid: float,
    t_right: float,
    left_values: dict[str, Optional[np.ndarray]],
    mid_values: dict[str, Optional[np.ndarray]],
    right_values: dict[str, Optional[np.ndarray]],
) -> float:
    max_error = 0.0
    for key in ("total", "fwd_sync", "fwd_ssc", "rev_sync", "rev_ssc", "cross_ic"):
        left = left_values[key]
        mid = mid_values[key]
        right = right_values[key]
        if left is None or mid is None or right is None:
            continue
        estimate = np.empty_like(mid)
        for i_freq in range(mid.shape[0]):
            estimate[i_freq] = _interpolate_segment_value(
                t_left,
                float(left[i_freq]),
                t_right,
                float(right[i_freq]),
                t_mid,
            )
        scale = np.maximum(np.abs(mid), 1.0e-99)
        error = float(np.max(np.abs(mid - estimate) / scale))
        if error > max_error:
            max_error = error
    return max_error


def _adaptive_frequency_error(
    nu_left: float,
    nu_mid: float,
    nu_right: float,
    left_values: dict[str, Optional[float]],
    mid_values: dict[str, Optional[float]],
    right_values: dict[str, Optional[float]],
) -> float:
    max_error = 0.0
    for key in ("total", "fwd_sync", "fwd_ssc", "rev_sync", "rev_ssc", "cross_ic"):
        left = left_values[key]
        mid = mid_values[key]
        right = right_values[key]
        if left is None or mid is None or right is None:
            continue
        estimate = _interpolate_segment_value(nu_left, float(left), nu_right, float(right), nu_mid)
        error = _relative_segment_error(float(mid), estimate)
        if error > max_error:
            max_error = error
    return max_error


def _interpolate_time_series(
    source_times_s: np.ndarray,
    source_values: np.ndarray,
    target_times_s: np.ndarray,
) -> np.ndarray:
    source_times_s = np.asarray(source_times_s, dtype=float)
    target_times_s = np.asarray(target_times_s, dtype=float)
    result = np.empty((source_values.shape[0], target_times_s.size), dtype=float)
    for i_freq in range(source_values.shape[0]):
        y = np.asarray(source_values[i_freq], dtype=float)
        if np.all(y >= 0.0):
            positive = y > 0.0
            if np.count_nonzero(positive) >= 2:
                y_pos = y[positive]
                floor = max(float(np.min(y_pos)) * 1.0e-30, 1.0e-300)
                y_safe = np.maximum(y, floor)
                result[i_freq] = np.exp(
                    np.interp(
                        np.log(target_times_s),
                        np.log(source_times_s),
                        np.log(y_safe),
                    )
                )
                continue
            if np.count_nonzero(positive) == 1:
                result[i_freq] = np.full(target_times_s.shape, float(y[positive][0]), dtype=float)
                continue
        if np.all(y > 0.0):
            result[i_freq] = np.exp(
                np.interp(
                    np.log(target_times_s),
                    np.log(source_times_s),
                    np.log(y),
                )
            )
        else:
            result[i_freq] = np.interp(target_times_s, source_times_s, y)
    return result


def _interpolate_frequency_series(
    source_freqs_hz: np.ndarray,
    source_values: np.ndarray,
    target_freqs_hz: np.ndarray,
) -> np.ndarray:
    source_freqs_hz = np.asarray(source_freqs_hz, dtype=float)
    source_values = np.asarray(source_values, dtype=float)
    target_freqs_hz = np.asarray(target_freqs_hz, dtype=float)
    if np.all(source_values >= 0.0):
        positive = source_values > 0.0
        if np.count_nonzero(positive) >= 2:
            src_pos = source_values[positive]
            floor = max(float(np.min(src_pos)) * 1.0e-30, 1.0e-300)
            values_safe = np.maximum(source_values, floor)
            return np.exp(
                np.interp(
                    np.log(target_freqs_hz),
                    np.log(source_freqs_hz),
                    np.log(values_safe),
                )
            )
        if np.count_nonzero(positive) == 1:
            return np.full(target_freqs_hz.shape, float(source_values[positive][0]), dtype=float)
    if np.all(source_values > 0.0):
        return np.exp(
            np.interp(
                np.log(target_freqs_hz),
                np.log(source_freqs_hz),
                np.log(source_values),
            )
        )
    return np.interp(target_freqs_hz, source_freqs_hz, source_values)


def _adaptive_flux_density_grid(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    *,
    tolerance: float = 1.0e-2,
    max_depth: int = 6,
    min_ratio: float = 1.02,
    coarse_nodes: int = 16,
) -> FluxResult:
    requested_times = np.asarray(times_s, dtype=float)
    order = np.argsort(requested_times)
    sorted_times = requested_times[order]
    unique_times = np.unique(sorted_times)
    if unique_times.size <= coarse_nodes:
        return model._compute_raw(requested_times, nu_hz)

    node_times = np.logspace(
        np.log10(float(unique_times[0])),
        np.log10(float(unique_times[-1])),
        min(int(coarse_nodes), unique_times.size),
    )
    node_result = model._compute_raw(node_times, nu_hz, solve_reference_times_s=requested_times)
    component_cache = {
        float(t): {
            "total": node_result.total[:, idx],
            "fwd_sync": node_result.fwd.sync[:, idx],
            "fwd_ssc": node_result.fwd.ssc[:, idx],
            "rev_sync": node_result.rev.sync[:, idx],
            "rev_ssc": node_result.rev.ssc[:, idx],
            "cross_ic": None if node_result.cross_ic is None else node_result.cross_ic[:, idx],
        }
        for idx, t in enumerate(node_times)
    }
    segments = [(float(node_times[i]), float(node_times[i + 1])) for i in range(node_times.size - 1)]

    for _ in range(max_depth):
        midpoint_times: list[float] = []
        midpoint_meta: list[tuple[float, float, float]] = []
        next_segments: list[tuple[float, float]] = []
        for t_left, t_right in segments:
            if t_right / t_left < min_ratio:
                next_segments.append((t_left, t_right))
                continue
            t_mid = _log_time_midpoint(t_left, t_right)
            if t_mid in component_cache:
                next_segments.append((t_left, t_right))
                continue
            midpoint_times.append(t_mid)
            midpoint_meta.append((t_left, t_mid, t_right))
        if not midpoint_times:
            break
        midpoint_result = model._compute_raw(
            np.array(midpoint_times, dtype=float),
            nu_hz,
            solve_reference_times_s=requested_times,
        )
        for idx, t_mid in enumerate(midpoint_times):
            component_cache[float(t_mid)] = {
                "total": midpoint_result.total[:, idx],
                "fwd_sync": midpoint_result.fwd.sync[:, idx],
                "fwd_ssc": midpoint_result.fwd.ssc[:, idx],
                "rev_sync": midpoint_result.rev.sync[:, idx],
                "rev_ssc": midpoint_result.rev.ssc[:, idx],
                "cross_ic": None if midpoint_result.cross_ic is None else midpoint_result.cross_ic[:, idx],
            }
        refined = False
        for t_left, t_mid, t_right in midpoint_meta:
            error = _adaptive_segment_error(
                t_left,
                t_mid,
                t_right,
                component_cache[t_left],
                component_cache[t_mid],
                component_cache[t_right],
            )
            if error > tolerance and t_right / t_left >= min_ratio:
                next_segments.append((t_left, t_mid))
                next_segments.append((t_mid, t_right))
                refined = True
            else:
                next_segments.append((t_left, t_right))
        segments = next_segments
        if not refined:
            break

    source_times = np.array(sorted(component_cache.keys()), dtype=float)
    total_grid = np.column_stack([component_cache[float(t)]["total"] for t in source_times])
    fwd_sync_grid = np.column_stack([component_cache[float(t)]["fwd_sync"] for t in source_times])
    fwd_ssc_grid = np.column_stack([component_cache[float(t)]["fwd_ssc"] for t in source_times])
    rev_sync_grid = np.column_stack([component_cache[float(t)]["rev_sync"] for t in source_times])
    rev_ssc_grid = np.column_stack([component_cache[float(t)]["rev_ssc"] for t in source_times])
    has_cross_ic = any(component_cache[float(t)]["cross_ic"] is not None for t in source_times)
    cross_ic_grid = None
    if has_cross_ic:
        cross_ic_grid = np.column_stack([
            np.zeros_like(component_cache[float(source_times[0])]["total"])
            if component_cache[float(t)]["cross_ic"] is None
            else component_cache[float(t)]["cross_ic"]
            for t in source_times
        ])

    total_sorted = _interpolate_time_series(source_times, total_grid, sorted_times)
    fwd_sync_sorted = _interpolate_time_series(source_times, fwd_sync_grid, sorted_times)
    fwd_ssc_sorted = _interpolate_time_series(source_times, fwd_ssc_grid, sorted_times)
    rev_sync_sorted = _interpolate_time_series(source_times, rev_sync_grid, sorted_times)
    rev_ssc_sorted = _interpolate_time_series(source_times, rev_ssc_grid, sorted_times)
    cross_ic_sorted = None if cross_ic_grid is None else _interpolate_time_series(source_times, cross_ic_grid, sorted_times)

    inverse_order = np.argsort(order)
    total = total_sorted[:, inverse_order]
    fwd_sync = fwd_sync_sorted[:, inverse_order]
    fwd_ssc = fwd_ssc_sorted[:, inverse_order]
    rev_sync = rev_sync_sorted[:, inverse_order]
    rev_ssc = rev_ssc_sorted[:, inverse_order]
    cross_ic = None if cross_ic_sorted is None else cross_ic_sorted[:, inverse_order]
    return FluxResult(
        total=total,
        fwd=FluxPair(sync=fwd_sync, ssc=fwd_ssc),
        rev=FluxPair(sync=rev_sync, ssc=rev_ssc),
        cross_ic=cross_ic,
    )


def _adaptive_spectrum_grid(
    model: Model,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    *,
    tolerance: float = 1.0e-2,
    max_depth: int = 6,
    min_ratio: float = 1.02,
    coarse_nodes: int = 16,
) -> FluxResult:
    requested_freqs = np.asarray(nu_hz, dtype=float)
    order = np.argsort(requested_freqs)
    sorted_freqs = requested_freqs[order]
    unique_freqs = np.unique(sorted_freqs)
    if unique_freqs.size <= coarse_nodes:
        return model._compute_raw(times_s, nu_hz)

    node_freqs = np.logspace(
        np.log10(float(unique_freqs[0])),
        np.log10(float(unique_freqs[-1])),
        min(int(coarse_nodes), unique_freqs.size),
    )
    node_result = model._compute_raw(times_s, node_freqs)
    component_cache = {
        float(freq): {
            "total": float(node_result.total[idx, 0]),
            "fwd_sync": float(node_result.fwd.sync[idx, 0]),
            "fwd_ssc": float(node_result.fwd.ssc[idx, 0]),
            "rev_sync": float(node_result.rev.sync[idx, 0]),
            "rev_ssc": float(node_result.rev.ssc[idx, 0]),
            "cross_ic": None if node_result.cross_ic is None else float(node_result.cross_ic[idx, 0]),
        }
        for idx, freq in enumerate(node_freqs)
    }
    segments = [(float(node_freqs[i]), float(node_freqs[i + 1])) for i in range(node_freqs.size - 1)]

    for _ in range(max_depth):
        midpoint_freqs: list[float] = []
        midpoint_meta: list[tuple[float, float, float]] = []
        next_segments: list[tuple[float, float]] = []
        for nu_left, nu_right in segments:
            if nu_right / nu_left < min_ratio:
                next_segments.append((nu_left, nu_right))
                continue
            nu_mid = _log_time_midpoint(nu_left, nu_right)
            if nu_mid in component_cache:
                next_segments.append((nu_left, nu_right))
                continue
            midpoint_freqs.append(nu_mid)
            midpoint_meta.append((nu_left, nu_mid, nu_right))
        if not midpoint_freqs:
            break
        midpoint_result = model._compute_raw(times_s, np.array(midpoint_freqs, dtype=float))
        for idx, nu_mid in enumerate(midpoint_freqs):
            component_cache[float(nu_mid)] = {
                "total": float(midpoint_result.total[idx, 0]),
                "fwd_sync": float(midpoint_result.fwd.sync[idx, 0]),
                "fwd_ssc": float(midpoint_result.fwd.ssc[idx, 0]),
                "rev_sync": float(midpoint_result.rev.sync[idx, 0]),
                "rev_ssc": float(midpoint_result.rev.ssc[idx, 0]),
                "cross_ic": None if midpoint_result.cross_ic is None else float(midpoint_result.cross_ic[idx, 0]),
            }
        refined = False
        for nu_left, nu_mid, nu_right in midpoint_meta:
            error = _adaptive_frequency_error(
                nu_left,
                nu_mid,
                nu_right,
                component_cache[nu_left],
                component_cache[nu_mid],
                component_cache[nu_right],
            )
            if error > tolerance and nu_right / nu_left >= min_ratio:
                next_segments.append((nu_left, nu_mid))
                next_segments.append((nu_mid, nu_right))
                refined = True
            else:
                next_segments.append((nu_left, nu_right))
        segments = next_segments
        if not refined:
            break

    source_freqs = np.array(sorted(component_cache.keys()), dtype=float)
    total_source = np.array([component_cache[float(freq)]["total"] for freq in source_freqs], dtype=float)
    fwd_sync_source = np.array([component_cache[float(freq)]["fwd_sync"] for freq in source_freqs], dtype=float)
    fwd_ssc_source = np.array([component_cache[float(freq)]["fwd_ssc"] for freq in source_freqs], dtype=float)
    rev_sync_source = np.array([component_cache[float(freq)]["rev_sync"] for freq in source_freqs], dtype=float)
    rev_ssc_source = np.array([component_cache[float(freq)]["rev_ssc"] for freq in source_freqs], dtype=float)
    has_cross_ic = any(component_cache[float(freq)]["cross_ic"] is not None for freq in source_freqs)
    cross_ic_source = None
    if has_cross_ic:
        cross_ic_source = np.array(
            [
                0.0 if component_cache[float(freq)]["cross_ic"] is None else float(component_cache[float(freq)]["cross_ic"])
                for freq in source_freqs
            ],
            dtype=float,
        )

    total_sorted = _interpolate_frequency_series(source_freqs, total_source, sorted_freqs)
    fwd_sync_sorted = _interpolate_frequency_series(source_freqs, fwd_sync_source, sorted_freqs)
    fwd_ssc_sorted = _interpolate_frequency_series(source_freqs, fwd_ssc_source, sorted_freqs)
    rev_sync_sorted = _interpolate_frequency_series(source_freqs, rev_sync_source, sorted_freqs)
    rev_ssc_sorted = _interpolate_frequency_series(source_freqs, rev_ssc_source, sorted_freqs)
    cross_ic_sorted = None if cross_ic_source is None else _interpolate_frequency_series(source_freqs, cross_ic_source, sorted_freqs)

    inverse_order = np.argsort(order)
    total = total_sorted[inverse_order][:, None]
    fwd_sync = fwd_sync_sorted[inverse_order][:, None]
    fwd_ssc = fwd_ssc_sorted[inverse_order][:, None]
    rev_sync = rev_sync_sorted[inverse_order][:, None]
    rev_ssc = rev_ssc_sorted[inverse_order][:, None]
    cross_ic = None if cross_ic_sorted is None else cross_ic_sorted[inverse_order][:, None]
    return FluxResult(
        total=total,
        fwd=FluxPair(sync=fwd_sync, ssc=fwd_ssc),
        rev=FluxPair(sync=rev_sync, ssc=rev_ssc),
        cross_ic=cross_ic,
    )


def _batch_fetch_pair_result(
    model: Model,
    cache: dict[tuple[float, float], tuple[float, float, float, float, float, Optional[float]]],
    query_pairs: list[tuple[float, float]],
) -> None:
    missing: list[tuple[float, float]] = []
    seen: set[tuple[float, float]] = set()
    for pair in query_pairs:
        if pair not in cache and pair not in seen:
            missing.append(pair)
            seen.add(pair)
    if not missing:
        return
    times_s = np.array([pair[0] for pair in missing], dtype=float)
    frequencies_hz = np.array([pair[1] for pair in missing], dtype=float)
    result = model.flux_density(times_s, frequencies_hz)
    for idx, pair in enumerate(missing):
        cross_ic = None if result.cross_ic is None else float(result.cross_ic[idx])
        cache[pair] = (
            float(result.total[idx]),
            float(result.fwd.sync[idx]),
            float(result.fwd.ssc[idx]),
            float(result.rev.sync[idx]),
            float(result.rev.ssc[idx]),
            cross_ic,
        )


@lru_cache(maxsize=None)
def _cached_leggauss(num_subsamples: int) -> tuple[np.ndarray, np.ndarray]:
    nodes_1d, weights_1d = np.polynomial.legendre.leggauss(max(int(num_subsamples), 1))
    return np.asarray(nodes_1d, dtype=float), np.asarray(weights_1d, dtype=float)


def _adaptive_exposure_average(
    model: Model,
    times_s: np.ndarray,
    frequencies_hz: np.ndarray,
    exposures_s: np.ndarray,
    num_subsamples: int,
) -> FluxResult:
    pair_cache: dict[tuple[float, float], tuple[float, float, float, float, float, Optional[float]]] = {}
    exposure_nodes: list[np.ndarray] = []
    exposure_weights: list[np.ndarray] = []
    initial_pairs: list[tuple[float, float]] = []

    for time_s, freq_hz, exposure_s in zip(times_s, frequencies_hz, exposures_s):
        t_start = max(float(time_s) - 0.5 * float(exposure_s), 1.0e-30)
        t_stop = float(time_s) + 0.5 * float(exposure_s)
        if np.isclose(t_start, t_stop):
            nodes = np.array([float(time_s)], dtype=float)
            weights = np.array([1.0], dtype=float)
        else:
            nodes_1d, weights_1d = _cached_leggauss(int(num_subsamples))
            half_width = 0.5 * (t_stop - t_start)
            center = 0.5 * (t_stop + t_start)
            nodes = half_width * nodes_1d + center
            weights = half_width * weights_1d
        exposure_nodes.append(nodes)
        exposure_weights.append(weights)
        initial_pairs.extend((float(node), float(freq_hz)) for node in nodes)

    _batch_fetch_pair_result(model, pair_cache, initial_pairs)

    total = np.zeros(times_s.shape[0], dtype=float)
    fwd_sync = np.zeros_like(total)
    fwd_ssc = np.zeros_like(total)
    rev_sync = np.zeros_like(total)
    rev_ssc = np.zeros_like(total)
    cross_ic = np.zeros_like(total)
    has_cross_ic = False

    for idx, nodes in enumerate(exposure_nodes):
        freq_hz = float(frequencies_hz[idx])
        node_array = np.asarray(nodes, dtype=float)
        weight_array = np.asarray(exposure_weights[idx], dtype=float)
        duration = float(np.sum(weight_array))
        values = np.array([pair_cache[(float(node), freq_hz)] for node in node_array], dtype=object)
        total_values = np.array([entry[0] for entry in values], dtype=float)
        fwd_sync_values = np.array([entry[1] for entry in values], dtype=float)
        fwd_ssc_values = np.array([entry[2] for entry in values], dtype=float)
        rev_sync_values = np.array([entry[3] for entry in values], dtype=float)
        rev_ssc_values = np.array([entry[4] for entry in values], dtype=float)
        cross_values = np.array([0.0 if entry[5] is None else entry[5] for entry in values], dtype=float)
        has_cross_ic = has_cross_ic or np.any(cross_values != 0.0)

        if duration == 0.0 or node_array.size == 1:
            total[idx] = total_values[0]
            fwd_sync[idx] = fwd_sync_values[0]
            fwd_ssc[idx] = fwd_ssc_values[0]
            rev_sync[idx] = rev_sync_values[0]
            rev_ssc[idx] = rev_ssc_values[0]
            cross_ic[idx] = cross_values[0]
            continue

        inv_duration = 1.0 / duration
        total[idx] = float(np.dot(total_values, weight_array) * inv_duration)
        fwd_sync[idx] = float(np.dot(fwd_sync_values, weight_array) * inv_duration)
        fwd_ssc[idx] = float(np.dot(fwd_ssc_values, weight_array) * inv_duration)
        rev_sync[idx] = float(np.dot(rev_sync_values, weight_array) * inv_duration)
        rev_ssc[idx] = float(np.dot(rev_ssc_values, weight_array) * inv_duration)
        cross_ic[idx] = float(np.dot(cross_values, weight_array) * inv_duration)

    return FluxResult(
        total=total,
        fwd=FluxPair(sync=fwd_sync, ssc=fwd_ssc),
        rev=FluxPair(sync=rev_sync, ssc=rev_ssc),
        cross_ic=None if not has_cross_ic else cross_ic,
    )
