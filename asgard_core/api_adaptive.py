from __future__ import annotations

import hashlib
from functools import cache

import numpy as np

from asgard_core.asgard_state import SolveState, project_flux_grid
from src import constants

from .api_model import CharTrack, FluxPair, Model, FluxResult

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


def _observe_parts(
    state: SolveState,
    times_s: np.ndarray,
    nu_hz: np.ndarray,
    mode: str = "full_components",
    projection_kind: str = "lightcurve",
) -> FluxResult:
    observed_state = project_flux_grid(state, times_s, nu_hz, mode=mode, projection_kind=projection_kind)
    observed = observed_state.components
    total = np.asarray(observed["total"], dtype=float)
    rev_sync = np.zeros_like(total) if observed["rev_sync"] is None else np.asarray(observed["rev_sync"], dtype=float)
    if observed.get("rev_hadronic") is not None:
        rev_sync = rev_sync + np.asarray(observed["rev_hadronic"], dtype=float)
    return FluxResult(
        total=total,
        fwd=FluxPair(
            sync=np.zeros_like(total) if observed["fwd_sync"] is None else np.asarray(observed["fwd_sync"], dtype=float),
            ssc=np.zeros_like(total) if observed["fwd_ssc"] is None else np.asarray(observed["fwd_ssc"], dtype=float),
        ),
        rev=FluxPair(
            sync=rev_sync,
            ssc=np.zeros_like(total) if observed["rev_ssc"] is None else np.asarray(observed["rev_ssc"], dtype=float),
        ),
        cross_ic=None if observed["cross_ic"] is None else np.asarray(observed["cross_ic"], dtype=float),
    )


def _adaptive_observer_time_grid(model: Model, times_s: np.ndarray) -> np.ndarray:
    # 用完整半径发射历史生成 EATS 解析时间网格；默认 API 的返回时间轴不变。
    user_times = _positive_unique_times(times_s)
    details = model.details(float(user_times[0]), float(user_times[-1]))
    base_log = np.logspace(np.log10(user_times[0]), np.log10(user_times[-1]), int(model.setups.num_tobs))
    knots = [user_times, base_log]
    knots.extend(_arrival_time_knots(model, track) for track in (details.fwd, details.rev) if track is not None)
    merged = _positive_unique_times(np.concatenate([item for item in knots if item.size > 0]))
    merged = merged[(merged >= user_times[0]) & (merged <= user_times[-1])]
    midpoints = np.sqrt(merged[:-1] * merged[1:]) if merged.size >= 2 else np.empty(0, dtype=float)
    return _positive_unique_times(np.concatenate((merged, midpoints)))


def _positive_unique_times(times_s: np.ndarray) -> np.ndarray:
    values = np.asarray(times_s, dtype=float).reshape(-1)
    if values.size == 0:
        raise ValueError("times_s must be non-empty.")
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("times_s must contain finite positive values.")
    return np.unique(values)


def _arrival_time_knots(model: Model, track: CharTrack) -> np.ndarray:
    radius = np.asarray(track.radius, dtype=float)
    t_axis = np.asarray(track.t_obs, dtype=float)
    mask = _emitting_radius_mask(track)
    if not np.any(mask):
        return np.empty(0, dtype=float)
    mu_values = _eats_mu_values(model)
    shell_times = [t_axis[mask]]
    shell_times.extend(t_axis[mask] + radius[mask] * (1.0 - mu) * (1.0 + model.observer.z) / constants.para_c for mu in mu_values)
    return np.concatenate(shell_times)


def _emitting_radius_mask(track: CharTrack) -> np.ndarray:
    radius = np.asarray(track.radius, dtype=float)
    mask = np.isfinite(radius) & (radius > 0.0)
    fields = [track.B_comv, track.secondary_rs_B, track.secondary_rs_u_diss]
    active = np.zeros(radius.shape, dtype=bool)
    for field in fields:
        if field is None:
            continue
        values = np.asarray(field, dtype=float)
        if values.shape == radius.shape:
            active |= np.isfinite(values) & (values > 0.0)
    if np.any(active):
        return mask & active
    return mask


def _eats_mu_values(model: Model) -> np.ndarray:
    theta_edges = np.linspace(0.0, model.jet.theta_j, int(model.setups.eats_num_theta) + 1)
    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    if float(model.observer.theta_obs) == 0.0 or int(model.setups.eats_num_phi) == 1:
        return np.cos(theta_centers)
    phi_edges = np.linspace(0.0, np.pi, int(model.setups.eats_num_phi) + 1)
    phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    theta = theta_centers[:, None]
    phi = phi_centers[None, :]
    theta_obs = float(model.observer.theta_obs)
    mu = np.cos(theta_obs) * np.cos(theta) + np.sin(theta_obs) * np.sin(theta) * np.cos(phi)
    return np.unique(mu.reshape(-1))


def _batch_fetch_pair_result(
    model: Model,
    cache: dict[tuple[float, float], tuple[float, float, float, float, float, float | None]],
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


@cache
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
    if np.any(~np.isfinite(times_s)) or np.any(times_s <= 0.0):
        raise ValueError("times_s must contain finite positive values.")
    if np.any(~np.isfinite(frequencies_hz)) or np.any(frequencies_hz <= 0.0):
        raise ValueError("frequencies_hz must contain finite positive values.")
    if np.any(~np.isfinite(exposures_s)) or np.any(exposures_s < 0.0):
        raise ValueError("exposures_s must contain finite non-negative values.")
    if np.any(times_s - 0.5 * exposures_s <= 0.0):
        raise ValueError("exposure windows must stay at positive observer times.")

    pair_cache: dict[tuple[float, float], tuple[float, float, float, float, float, float | None]] = {}
    exposure_nodes: list[np.ndarray] = []
    exposure_weights: list[np.ndarray] = []
    initial_pairs: list[tuple[float, float]] = []

    for time_s, freq_hz, exposure_s in zip(times_s, frequencies_hz, exposures_s):
        t_start = float(time_s) - 0.5 * float(exposure_s)
        t_stop = float(time_s) + 0.5 * float(exposure_s)
        if float(exposure_s) == 0.0:
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
    fwd_sync, fwd_ssc, rev_sync, rev_ssc, cross_ic = np.zeros((5, times_s.shape[0]), dtype=float)
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
