from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src import Interpolation


@dataclass(frozen=True)
class EATSGeometry:
    redshift: float
    opening_angle_rad: float
    viewing_angle_rad: float = 0.0


@dataclass(frozen=True)
class EATSNumerics:
    num_theta: int = 64
    num_phi: int = 1
    num_threads: int = 1
    adaptive_rtol: float = 0.0
    adaptive_max_depth: int = 0


def project_branch_flux(
    *,
    characteristic_time_s: np.ndarray,
    gamma: np.ndarray,
    radius_cm: np.ndarray,
    source_flux: np.ndarray,
    seed_frequency_hz: np.ndarray,
    observer_frequency_hz: np.ndarray,
    observer_time_s: np.ndarray,
    geometry: EATSGeometry,
    numerics: EATSNumerics,
) -> np.ndarray:
    if geometry.viewing_angle_rad != 0.0 and numerics.num_phi == 1:
        raise ValueError(
            "off-axis EATS projection requires num_phi >= 2; "
            "num_phi=1 is only valid for on-axis axial collapse."
        )
    boundary = np.zeros(10, dtype=float)
    boundary[0] = float(gamma[0])
    boundary[3] = float(radius_cm[0])
    boundary[7] = geometry.redshift
    boundary[8] = geometry.opening_angle_rad
    boundary[9] = geometry.viewing_angle_rad
    order = np.argsort(observer_frequency_hz)
    sorted_observer_frequency = np.asarray(observer_frequency_hz, dtype=float)[order]
    if numerics.adaptive_max_depth > 0:
        projected = Interpolation.sed_interpolation_adaptive_theta(
            boundary,
            np.asarray(characteristic_time_s, dtype=float),
            np.asarray(gamma, dtype=float),
            np.asarray(radius_cm, dtype=float),
            np.asarray(source_flux, dtype=float),
            np.asarray(seed_frequency_hz, dtype=float),
            sorted_observer_frequency,
            np.asarray(observer_time_s, dtype=float),
            numerics.adaptive_rtol,
            numerics.adaptive_max_depth,
            numerics.num_theta,
            numerics.num_phi,
            numerics.num_threads,
        )
    else:
        projected = Interpolation.sed_interpolation(
            boundary,
            np.asarray(characteristic_time_s, dtype=float),
            np.asarray(gamma, dtype=float),
            np.asarray(radius_cm, dtype=float),
            np.asarray(source_flux, dtype=float),
            np.asarray(seed_frequency_hz, dtype=float),
            sorted_observer_frequency,
            np.asarray(observer_time_s, dtype=float),
            numerics.num_theta,
            numerics.num_phi,
            numerics.num_threads,
        )
    if np.array_equal(order, np.arange(order.shape[0])):
        return projected
    restored = np.empty_like(projected)
    restored[order] = projected
    return restored
