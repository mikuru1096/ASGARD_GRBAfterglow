from __future__ import annotations

import hashlib
from functools import lru_cache
from typing import Optional

import numpy as np

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
    projection_kind: str = "lightcurve",
) -> FluxResult:
    observed_state = project_flux_grid(state, times_s, nu_hz, mode=mode, projection_kind=projection_kind)
    return _pack_flux(observed_state.components)


def _observe_total(
    state: SolveState,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    timings: Optional[dict[str, float]] = None,
    projection_kind: str = "lightcurve",
) -> np.ndarray:
    observed_state = project_flux_grid(
        state,
        times_s,
        nu_hz,
        timings=timings,
        mode="total_only",
        projection_kind=projection_kind,
    )
    return np.asarray(observed_state.components["total"], dtype=float)


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
