from __future__ import annotations

import numpy as np

from asgard_models import FitConfig
from src import Radiation, constants


def _ambient_density(radius_cm: np.ndarray, config: FitConfig) -> np.ndarray:
    radius_cm = np.asarray(radius_cm, dtype=float)
    if config.a_star > 0.0:
        d_ne_wind = config.a_star * 3.0e35 / radius_cm**2
        d_ne = np.where(d_ne_wind <= config.d_ne / 4.0, config.d_ne, d_ne_wind)
    else:
        d_ne = config.d_ne * (
            1.0
            + (config.f_jump - 1.0)
            * np.exp(-(np.log10(radius_cm) - np.log10(config.r_tr)) ** 2 / (2.0 * config.f_wide**2))
        )

    if config.a_star > 0.0:
        d_ne = np.where(radius_cm < config.r0, config.a_star * 3.0e35 / config.r0**2, d_ne)
    return d_ne


def _compute_forward_nu_M(gamma: np.ndarray, radius_cm: np.ndarray, config: FitConfig) -> np.ndarray:
    magnetic_field = 0.39 * np.sqrt(config.epsilon_b * _ambient_density(radius_cm, config) * gamma * np.maximum(gamma - 1.0, 0.0))
    doppler = gamma * (1.0 - np.sqrt(1.0 - gamma**-2)) * (1.0 + config.z)
    gam_e_max = 3.0 * constants.para_m_energy / np.sqrt(8.0 * magnetic_field * constants.para_e**3)
    return 4.2e6 * magnetic_field * gam_e_max**2 / doppler


def adaptive_log_grid_with_breaks(
    grid_min: float,
    grid_max: float,
    breaks: np.ndarray,
    break_weights: np.ndarray,
    pts_per_decade: float = 3.0,
    max_refined_breaks: int = 4,
    refine_radius_decades: float = 0.5,
    refine_factor: float = 3.0,
    max_points: int = 128,
) -> np.ndarray:
    lg_min = np.log10(grid_min)
    lg_max = np.log10(grid_max)
    coarse_step = 1.0 / max(float(pts_per_decade), 1.0e-6)
    fine_step = coarse_step / float(refine_factor)
    merge_eps = 0.5 * fine_step

    valid = []
    for i, value in enumerate(np.asarray(breaks, dtype=float)):
        if not np.isfinite(value) or value <= 0.0:
            continue
        lg_value = np.log10(value)
        if lg_min < lg_value < lg_max:
            weight = 1.0
            if i < len(break_weights):
                candidate = float(break_weights[i])
                if np.isfinite(candidate) and candidate >= 0.0:
                    weight = candidate
            valid.append((lg_value, weight))

    valid.sort(key=lambda item: item[1], reverse=True)
    n_refine = min(max_refined_breaks, len(valid))

    span = max(lg_max - lg_min, 0.0)
    n_coarse = max(1, int(np.ceil(span / coarse_step)))
    points = [(lg_min + span * i / n_coarse, False) for i in range(n_coarse + 1)]
    for lg_value, _ in valid:
        points.append((lg_value, True))
    for i in range(n_refine):
        lo = max(lg_min, valid[i][0] - refine_radius_decades)
        hi = min(lg_max, valid[i][0] + refine_radius_decades)
        n_local = max(1, int((hi - lo) / fine_step))
        for k in range(n_local + 1):
            points.append((lo + (hi - lo) * k / n_local, False))

    points.sort(key=lambda item: item[0])
    merged: list[tuple[float, bool]] = []
    for point in points:
        if merged and point[0] - merged[-1][0] < merge_eps:
            if point[1] and not merged[-1][1]:
                merged[-1] = point
            else:
                merged[-1] = (merged[-1][0], merged[-1][1] or point[1])
        else:
            merged.append(point)

    if not merged:
        merged = [(lg_min, True), (lg_max, True)]
    else:
        merged[0] = (lg_min, True)
        if len(merged) == 1:
            merged.append((lg_max, True))
        else:
            merged[-1] = (lg_max, True)

    if len(merged) > max_points:
        keep = [merged[i * (len(merged) - 1) // (max_points - 1)] for i in range(max_points)]
        merged = keep

    return np.power(10.0, np.array([item[0] for item in merged], dtype=float))


def _interp_positive_loglog(source_grid: np.ndarray, values: np.ndarray, target_grid: np.ndarray) -> np.ndarray:
    source_grid = np.asarray(source_grid, dtype=float)
    target_grid = np.asarray(target_grid, dtype=float)
    values = np.asarray(values, dtype=float)
    out = np.zeros((target_grid.shape[0], values.shape[1]), dtype=float)
    log_source = np.log(source_grid)
    log_target = np.log(target_grid)

    for i in range(values.shape[1]):
        column = values[:, i]
        positive = column > 0.0
        if not np.any(positive):
            continue
        x = log_source[positive]
        y = np.log(column[positive])
        inside = (log_target >= x[0]) & (log_target <= x[-1])
        if np.any(inside):
            out[inside, i] = np.exp(np.interp(log_target[inside], x, y))
    return out


def _build_forward_ssc_grid(
    full_grid_hz: np.ndarray,
    seed_syn: np.ndarray,
    gamma_bulk: np.ndarray,
    radius_cm: np.ndarray,
    nu_a: np.ndarray,
    nu_m: np.ndarray,
    nu_c: np.ndarray,
    config: FitConfig,
) -> np.ndarray:
    nu_M = _compute_forward_nu_M(gamma_bulk, radius_cm, config)
    break_arrays = [nu_a, nu_m, nu_c, nu_M]
    break_list: list[float] = []
    weight_list: list[float] = []

    for arr in break_arrays:
        valid = np.asarray(arr, dtype=float)
        valid = valid[np.isfinite(valid) & (valid > 0.0)]
        if valid.size == 0:
            continue
        for q in (0.1, 0.5, 0.9):
            value = float(np.quantile(valid, q))
            idx = int(np.clip(np.searchsorted(full_grid_hz, value), 0, full_grid_hz.shape[0] - 1))
            weight = float(np.sum(seed_syn[idx])) / max(full_grid_hz[idx] * full_grid_hz[idx], 1.0)
            break_list.append(value)
            weight_list.append(weight)

    return adaptive_log_grid_with_breaks(
        float(full_grid_hz[0]),
        float(full_grid_hz[-1]),
        np.array(break_list, dtype=float),
        np.array(weight_list, dtype=float),
        pts_per_decade=3.0,
        max_refined_breaks=4,
        refine_radius_decades=0.5,
        refine_factor=3.0,
        max_points=min(128, full_grid_hz.shape[0]),
    )


def compute_forward_ssc_adaptive(
    radius_cm: np.ndarray,
    gam_e: np.ndarray,
    d_n_gam_e: np.ndarray,
    full_grid_hz: np.ndarray,
    seed_syn: np.ndarray,
    gamma_bulk: np.ndarray,
    nu_a: np.ndarray,
    nu_m: np.ndarray,
    nu_c: np.ndarray,
    config: FitConfig,
) -> tuple[np.ndarray, np.ndarray]:
    reduced_grid_hz = _build_forward_ssc_grid(full_grid_hz, seed_syn, gamma_bulk, radius_cm, nu_a, nu_m, nu_c, config)
    if reduced_grid_hz.shape[0] >= full_grid_hz.shape[0]:
        return Radiation.ssc_spec(radius_cm, gam_e, d_n_gam_e, full_grid_hz, seed_syn, config.num_threads)

    reduced_seed = _interp_positive_loglog(full_grid_hz, seed_syn, reduced_grid_hz)
    p_ssc_reduced, seed_ssc_reduced = Radiation.ssc_spec(
        radius_cm,
        gam_e,
        d_n_gam_e,
        reduced_grid_hz,
        reduced_seed,
        config.num_threads,
    )
    p_ssc_full = _interp_positive_loglog(reduced_grid_hz, p_ssc_reduced, full_grid_hz)
    seed_ssc_full = _interp_positive_loglog(reduced_grid_hz, seed_ssc_reduced, full_grid_hz)
    return p_ssc_full, seed_ssc_full
