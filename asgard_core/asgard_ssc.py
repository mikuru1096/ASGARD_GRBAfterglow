from __future__ import annotations

import numpy as np

from asgard_core.asgard_config import FitConfig
from asgard_core.asgard_physics_utils import ambient_density, compute_magnetic_field, compute_doppler
from src import Radiation, constants


DEFAULT_AUXILIARY_GAMMA_COUNT = 64


def _ambient_density(radius_cm: np.ndarray, config: FitConfig) -> np.ndarray:
    """DEPRECATED: Use asgard_physics_utils.ambient_density instead."""
    return ambient_density(radius_cm, config)


def _compute_forward_nu_M(gamma: np.ndarray, radius_cm: np.ndarray, config: FitConfig) -> np.ndarray:
    magnetic_field = compute_magnetic_field(gamma, radius_cm, config)
    doppler = compute_doppler(gamma, config.z)
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
    if not np.isfinite(grid_min) or not np.isfinite(grid_max) or grid_min <= 0.0 or grid_max <= grid_min:
        raise ValueError("adaptive SSC grid bounds must be finite, positive, and strictly ordered.")
    if not np.isfinite(pts_per_decade) or pts_per_decade <= 0.0:
        raise ValueError("pts_per_decade must be positive.")
    if int(max_refined_breaks) < 0:
        raise ValueError("max_refined_breaks must be non-negative.")
    if not np.isfinite(refine_radius_decades) or refine_radius_decades < 0.0:
        raise ValueError("refine_radius_decades must be non-negative.")
    if not np.isfinite(refine_factor) or refine_factor <= 0.0:
        raise ValueError("refine_factor must be positive.")
    if int(max_points) < 2:
        raise ValueError("max_points must be at least 2.")
    lg_min = np.log10(grid_min)
    lg_max = np.log10(grid_max)
    coarse_step = 1.0 / float(pts_per_decade)
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

    span = lg_max - lg_min
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


def _shell_support_bounds(work_x_edge_log10: np.ndarray, work_d_n_x: np.ndarray) -> tuple[float, float]:
    x_edge = np.asarray(work_x_edge_log10, dtype=float)
    q_work = np.asarray(work_d_n_x, dtype=float)
    if x_edge.ndim != 2 or q_work.ndim != 2:
        raise ValueError("SSC auxiliary support bounds require 2D edge and work arrays.")
    if x_edge.shape[1] != q_work.shape[1] or x_edge.shape[0] != q_work.shape[0] + 1:
        raise ValueError("SSC auxiliary support bounds require edge shape (num_cell+1, num_shell).")
    if not np.all(np.isfinite(x_edge)) or not np.all(np.isfinite(q_work)):
        raise ValueError("SSC auxiliary support bounds require finite inputs.")
    if np.any(q_work < 0.0):
        raise ValueError("SSC auxiliary support bounds require non-negative work distributions.")
    x_lo = np.inf
    x_hi = -np.inf
    for i_shell in range(q_work.shape[1]):
        shell_edge = x_edge[:, i_shell]
        if np.any(np.diff(shell_edge) <= 0.0):
            raise ValueError("SSC auxiliary support bounds require strictly increasing edge grids in every shell.")
        q_shell = q_work[:, i_shell]
        active = np.where(q_shell > 0.0)[0]
        if active.size == 0:
            continue
        x_lo = min(x_lo, float(shell_edge[active[0]]))
        x_hi = max(x_hi, float(shell_edge[active[-1] + 1]))
    if not np.isfinite(x_lo) or not np.isfinite(x_hi) or x_hi <= x_lo:
        raise ValueError("SSC auxiliary support bounds are undefined because the work distribution is identically zero.")
    return x_lo, x_hi


def _build_auxiliary_gamma_edges(
    work_x_edge_log10: np.ndarray,
    work_d_n_x: np.ndarray,
    num_auxiliary_gamma: int,
) -> np.ndarray:
    x_lo, x_hi = _shell_support_bounds(work_x_edge_log10, work_d_n_x)
    margin = 0.08
    return np.linspace(x_lo - margin, x_hi + margin, num_auxiliary_gamma + 1, dtype=float)


def _validate_nonuniform_log_grid(x_edge: np.ndarray, q: np.ndarray) -> None:
    edge = np.asarray(x_edge, dtype=float)
    values = np.asarray(q, dtype=float)
    if edge.ndim != 1 or values.ndim != 1:
        raise ValueError("nonuniform SSC helper expects 1D edge and value arrays.")
    if edge.size != values.size + 1:
        raise ValueError("nonuniform SSC edge grid must have one more element than cell values.")
    if not np.all(np.isfinite(edge)) or not np.all(np.isfinite(values)):
        raise ValueError("nonuniform SSC helper received non-finite inputs.")
    if np.any(np.diff(edge) <= 0.0):
        raise ValueError("nonuniform SSC edge grid must be strictly increasing.")
    if np.any(values < 0.0):
        raise ValueError("nonuniform SSC cell values must be non-negative.")


def _minmod(a: float, b: float) -> float:
    if a * b <= 0.0:
        return 0.0
    return a if abs(a) <= abs(b) else b


def _minmod3(a: float, b: float, c: float) -> float:
    return _minmod(a, _minmod(b, c))


def _loglinear_cell_int(qlog_center: float, slope: float, x_center: float, xa: float, xb: float) -> float:
    if xb <= xa:
        return 0.0
    if abs(slope) <= 1.0e-14:
        return (10.0**qlog_center - 1.0) * (xb - xa)
    alpha = slope
    beta = qlog_center - alpha * x_center
    ln10 = np.log(10.0)
    u_a = (beta + alpha * xa) * ln10
    u_b = (beta + alpha * xb) * ln10
    u_max = max(u_a, u_b)
    exp_span = np.exp(u_b - u_max) - np.exp(u_a - u_max)
    return np.exp(u_max) * exp_span / (alpha * ln10) - (xb - xa)


def _build_log_prefix_nonuniform(x_edge_old: np.ndarray, q_old: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    _validate_nonuniform_log_grid(x_edge_old, q_old)
    num_old = q_old.shape[0]
    x_center = 0.5 * (x_edge_old[:-1] + x_edge_old[1:])
    qlog = np.log10(1.0 + q_old)
    slope = np.zeros(num_old, dtype=float)
    for i_cell in range(1, num_old - 1):
        dl = (qlog[i_cell] - qlog[i_cell - 1]) / (x_center[i_cell] - x_center[i_cell - 1])
        dr = (qlog[i_cell + 1] - qlog[i_cell]) / (x_center[i_cell + 1] - x_center[i_cell])
        dc = (qlog[i_cell + 1] - qlog[i_cell - 1]) / (x_center[i_cell + 1] - x_center[i_cell - 1])
        slope[i_cell] = _minmod3(2.0 * dl, dc, 2.0 * dr)
    prefix = np.zeros(num_old + 1, dtype=float)
    for i_cell in range(num_old):
        prefix[i_cell + 1] = prefix[i_cell] + _loglinear_cell_int(
            qlog[i_cell],
            slope[i_cell],
            x_center[i_cell],
            x_edge_old[i_cell],
            x_edge_old[i_cell + 1],
        )
    return prefix, qlog, slope, x_center


def _log_prefix_eval_nonuniform(
    x_edge_old: np.ndarray,
    prefix: np.ndarray,
    qlog: np.ndarray,
    slope: np.ndarray,
    x_center: np.ndarray,
    x_eval: float,
) -> float:
    xa = max(float(x_edge_old[0]), min(float(x_eval), float(x_edge_old[-1])))
    if xa <= x_edge_old[0]:
        return 0.0
    if xa >= x_edge_old[-1]:
        return float(prefix[-1])
    i_cell = int(np.searchsorted(x_edge_old[1:], xa, side="left"))
    return float(prefix[i_cell]) + _loglinear_cell_int(qlog[i_cell], slope[i_cell], x_center[i_cell], x_edge_old[i_cell], xa)


def _conservative_remap_log_nonuniform(x_edge_old: np.ndarray, x_edge_new: np.ndarray, q_old: np.ndarray) -> np.ndarray:
    prefix, qlog, slope, x_center = _build_log_prefix_nonuniform(x_edge_old, q_old)
    if np.any(np.diff(np.asarray(x_edge_new, dtype=float)) <= 0.0):
        raise ValueError("target nonuniform SSC edge grid must be strictly increasing.")
    num_new = x_edge_new.shape[0] - 1
    q_new = np.zeros(num_new, dtype=float)
    for i_cell in range(num_new):
        dx_cell = float(x_edge_new[i_cell + 1] - x_edge_new[i_cell])
        right = _log_prefix_eval_nonuniform(x_edge_old, prefix, qlog, slope, x_center, float(x_edge_new[i_cell + 1]))
        left = _log_prefix_eval_nonuniform(x_edge_old, prefix, qlog, slope, x_center, float(x_edge_new[i_cell]))
        q_new[i_cell] = (right - left) / dx_cell
        if q_new[i_cell] < 0.0:
            raise ValueError("conservative nonuniform SSC remap produced a negative cell average.")
    return q_new


def _ppm_interfaces_nonuniform(x_edge: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    _validate_nonuniform_log_grid(x_edge, q)
    x_center = 0.5 * (x_edge[:-1] + x_edge[1:])
    dx_cell = x_edge[1:] - x_edge[:-1]
    slope = np.zeros_like(q, dtype=float)
    for i_cell in range(1, q.shape[0] - 1):
        dl = (q[i_cell] - q[i_cell - 1]) / (x_center[i_cell] - x_center[i_cell - 1])
        dr = (q[i_cell + 1] - q[i_cell]) / (x_center[i_cell + 1] - x_center[i_cell])
        dc = (q[i_cell + 1] - q[i_cell - 1]) / (x_center[i_cell + 1] - x_center[i_cell - 1])
        slope[i_cell] = _minmod3(2.0 * dl, dc, 2.0 * dr)

    q_left = np.array(q, copy=True)
    q_right = np.array(q, copy=True)
    for i_cell in range(1, q.shape[0] - 1):
        q_left[i_cell] = q[i_cell] - 0.5 * slope[i_cell] * dx_cell[i_cell]
        q_right[i_cell] = q[i_cell] + 0.5 * slope[i_cell] * dx_cell[i_cell]
        q_min = min(q[i_cell - 1], q[i_cell], q[i_cell + 1])
        q_max = max(q[i_cell - 1], q[i_cell], q[i_cell + 1])
        q_left[i_cell] = max(q_min, min(q_max, q_left[i_cell]))
        q_right[i_cell] = max(q_min, min(q_max, q_right[i_cell]))
        if (q_right[i_cell] - q[i_cell]) * (q[i_cell] - q_left[i_cell]) <= 0.0:
            q_left[i_cell] = q[i_cell]
            q_right[i_cell] = q[i_cell]
        else:
            q_bar = 0.5 * (q_left[i_cell] + q_right[i_cell])
            dq = q_right[i_cell] - q_left[i_cell]
            if dq * (q[i_cell] - q_bar) > dq * dq / 6.0:
                q_left[i_cell] = 3.0 * q[i_cell] - 2.0 * q_right[i_cell]
            elif dq * (q[i_cell] - q_bar) < -dq * dq / 6.0:
                q_right[i_cell] = 3.0 * q[i_cell] - 2.0 * q_left[i_cell]
    return q_left, q_right


def _ppm_point_values_nonuniform(x_src_edge: np.ndarray, q_src: np.ndarray, x_tgt: np.ndarray) -> np.ndarray:
    q_left, q_right = _ppm_interfaces_nonuniform(x_src_edge, q_src)
    q_tgt = np.zeros_like(x_tgt, dtype=float)
    for i_tgt, x_val in enumerate(x_tgt):
        if x_val < x_src_edge[0] or x_val > x_src_edge[-1]:
            continue
        i_src = int(np.searchsorted(x_src_edge[1:], x_val, side="left"))
        dx = x_src_edge[i_src + 1] - x_src_edge[i_src]
        xi = (x_val - x_src_edge[i_src]) / dx
        coeff_c = q_src[i_src] - 0.5 * (q_left[i_src] + q_right[i_src])
        q_tgt[i_tgt] = q_left[i_src] + xi * (q_right[i_src] - q_left[i_src] + 6.0 * coeff_c) + 6.0 * coeff_c * xi * (1.0 - xi)
        if q_tgt[i_tgt] < 0.0:
            raise ValueError("nonuniform SSC PPM reconstruction produced a negative point value.")
    return q_tgt


def project_work_grid_to_auxiliary_gamma(
    work_x_edge_log10: np.ndarray,
    work_d_n_x: np.ndarray,
    num_auxiliary_gamma: int = DEFAULT_AUXILIARY_GAMMA_COUNT,
) -> tuple[np.ndarray, np.ndarray]:
    x_aux_edge = _build_auxiliary_gamma_edges(work_x_edge_log10, work_d_n_x, num_auxiliary_gamma)
    x_aux = 0.5 * (x_aux_edge[:-1] + x_aux_edge[1:])
    dx_aux = x_aux_edge[1:] - x_aux_edge[:-1]
    q_aux = np.zeros((num_auxiliary_gamma, work_d_n_x.shape[1]), dtype=float)
    for i_shell in range(work_d_n_x.shape[1]):
        x_src_edge = np.asarray(work_x_edge_log10[:, i_shell], dtype=float)
        q_src = np.asarray(work_d_n_x[:, i_shell], dtype=float)
        q_aux_avg = _conservative_remap_log_nonuniform(x_src_edge, x_aux_edge, q_src)
        q_aux_shell = _ppm_point_values_nonuniform(x_aux_edge, q_aux_avg, x_aux)
        total_src = float(np.sum(q_src * np.diff(x_src_edge)))
        total_aux = float(np.sum(q_aux_shell * dx_aux))
        if total_src > 0.0 and total_aux > 0.0:
            q_aux_shell *= total_src / total_aux
        q_aux[:, i_shell] = q_aux_shell

    gam_aux = np.power(10.0, x_aux)
    d_n_gam_aux = q_aux / (gam_aux[:, None] * np.log(10.0))
    return np.asfortranarray(gam_aux), np.asfortranarray(d_n_gam_aux)


def compute_ssc_auxiliary_grid(
    radius_cm: np.ndarray,
    work_x_edge_log10: np.ndarray,
    work_d_n_x: np.ndarray,
    seed_frequency_hz: np.ndarray,
    seed_field: np.ndarray,
    num_threads: int,
    num_auxiliary_gamma: int = DEFAULT_AUXILIARY_GAMMA_COUNT,
) -> tuple[np.ndarray, np.ndarray]:
    radius_cm = np.asfortranarray(np.asarray(radius_cm, dtype=float))
    work_x_edge_log10 = np.asfortranarray(np.asarray(work_x_edge_log10, dtype=float))
    work_d_n_x = np.asfortranarray(np.asarray(work_d_n_x, dtype=float))
    seed_frequency_hz = np.asfortranarray(np.asarray(seed_frequency_hz, dtype=float))
    seed_field = np.asfortranarray(np.asarray(seed_field, dtype=float))

    if np.all(work_d_n_x == 0.0):
        zeros = np.zeros((seed_frequency_hz.shape[0], radius_cm.shape[0]), dtype=float, order="F")
        return zeros, zeros.copy(order="F")

    if work_x_edge_log10.ndim == 2 and work_d_n_x.ndim == 2:
        return Radiation.ssc_spec_nonuniform(
            radius_cm,
            work_x_edge_log10,
            work_d_n_x,
            seed_frequency_hz,
            seed_field,
            num_threads,
        )

    gam_aux, d_n_gam_aux = project_work_grid_to_auxiliary_gamma(
        work_x_edge_log10,
        work_d_n_x,
        num_auxiliary_gamma=num_auxiliary_gamma,
    )
    return Radiation.ssc_spec(
        radius_cm,
        gam_aux,
        d_n_gam_aux,
        seed_frequency_hz,
        seed_field,
        num_threads,
    )


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
            weight = float(np.sum(seed_syn[idx])) / (full_grid_hz[idx] * full_grid_hz[idx])
            break_list.append(value)
            weight_list.append(weight)

    return adaptive_log_grid_with_breaks(
        float(full_grid_hz[0]),
        float(full_grid_hz[-1]),
        np.array(break_list, dtype=float),
        np.array(weight_list, dtype=float),
        pts_per_decade=5.0,
        max_refined_breaks=6,
        refine_radius_decades=0.35,
        refine_factor=4.0,
        max_points=max(full_grid_hz.shape[0], int(np.ceil(1.5 * full_grid_hz.shape[0]))),
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


def _densify_log_grid(grid_hz: np.ndarray, target_points: int) -> np.ndarray:
    grid = np.asarray(grid_hz, dtype=float)
    grid = grid[np.isfinite(grid) & (grid > 0.0)]
    if grid.size == 0:
        raise ValueError("grid_hz must contain positive finite values.")
    grid = np.unique(np.sort(grid))
    if grid.size >= target_points:
        return grid

    while grid.size < target_points:
        if grid.size < 2:
            grid = np.array([grid[0], grid[0] * 10.0], dtype=float)
            continue
        gaps = np.diff(np.log(grid))
        idx = int(np.argmax(gaps))
        midpoint = float(np.exp(0.5 * (np.log(grid[idx]) + np.log(grid[idx + 1]))))
        if midpoint <= grid[idx] or midpoint >= grid[idx + 1]:
            break
        grid = np.insert(grid, idx + 1, midpoint)
    return grid


def compute_forward_ssc_seed_adaptive(
    radius_cm: np.ndarray,
    work_x_edge_log10: np.ndarray,
    work_d_n_x: np.ndarray,
    full_grid_hz: np.ndarray,
    seed_syn: np.ndarray,
    gamma_bulk: np.ndarray,
    nu_a: np.ndarray,
    nu_m: np.ndarray,
    nu_c: np.ndarray,
    config: FitConfig,
) -> tuple[np.ndarray, np.ndarray]:
    reduced_grid_hz = _build_forward_ssc_grid(full_grid_hz, seed_syn, gamma_bulk, radius_cm, nu_a, nu_m, nu_c, config)
    target_points = max(full_grid_hz.shape[0], int(np.ceil(1.5 * full_grid_hz.shape[0])))
    if reduced_grid_hz.shape[0] < target_points:
        reduced_grid_hz = _densify_log_grid(reduced_grid_hz, target_points)
    reduced_seed = _interp_positive_loglog(full_grid_hz, seed_syn, reduced_grid_hz)
    p_ssc_reduced, seed_ssc_reduced = compute_ssc_auxiliary_grid(
        radius_cm,
        work_x_edge_log10,
        work_d_n_x,
        reduced_grid_hz,
        reduced_seed,
        config.num_threads,
    )
    p_ssc_full = _interp_positive_loglog(reduced_grid_hz, p_ssc_reduced, full_grid_hz)
    seed_ssc_full = _interp_positive_loglog(reduced_grid_hz, seed_ssc_reduced, full_grid_hz)
    return p_ssc_full, seed_ssc_full
