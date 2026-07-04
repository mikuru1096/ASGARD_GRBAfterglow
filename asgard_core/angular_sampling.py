from __future__ import annotations

from dataclasses import dataclass

import numpy as np


SUPPORTED_PATCH_SAMPLING = ("uniform",)
AXISYMMETRIC_JET_KINDS = frozenset(("tophat", "gaussian", "powerlaw", "twocomponent", "steppowerlaw", "tabulated"))


def is_axisymmetric_jet(jet) -> bool:
    return str(getattr(jet, "kind", "")).lower() in AXISYMMETRIC_JET_KINDS


@dataclass(frozen=True)
class PatchGrid:
    theta_centers: np.ndarray
    theta_edges: np.ndarray
    phi_centers: np.ndarray
    phi_edges: np.ndarray
    domega: np.ndarray
    patch_half_angle: np.ndarray


def build_patch_grid(
    model,
    observer_time_s: np.ndarray | None = None,
    theta_count: int | None = None,
    phi_count: int | None = None,
) -> PatchGrid:
    return _build_uniform_grid(model, theta_count=theta_count, phi_count=phi_count)


def _patch_counts(model, theta_count: int | None = None, phi_count: int | None = None) -> tuple[int, int]:
    n_theta = int(model.setups.structured_num_theta if theta_count is None else theta_count)
    n_phi = int(model.setups.structured_num_phi if phi_count is None else phi_count)
    if n_theta <= 0 or n_phi <= 0:
        raise ValueError("structured_num_theta and structured_num_phi must be positive for patch sampling.")
    return n_theta, n_phi


def _build_uniform_grid(model, theta_count: int | None = None, phi_count: int | None = None) -> PatchGrid:
    n_theta, n_phi = _patch_counts(model, theta_count=theta_count, phi_count=phi_count)
    theta_edges = np.linspace(0.0, float(model.jet.theta_max), n_theta + 1)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, n_phi + 1)
    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    domega = _solid_angle_cells(theta_edges, phi_edges)
    return PatchGrid(
        theta_centers=theta_centers,
        theta_edges=theta_edges,
        phi_centers=phi_centers,
        phi_edges=phi_edges,
        domega=domega,
        patch_half_angle=np.sqrt(np.maximum(domega, 1.0e-12) / np.pi),
    )


def angular_separation(
    theta: float | np.ndarray,
    phi: float | np.ndarray,
    theta_obs: float,
    phi_obs: float,
) -> np.ndarray:
    cos_sep = (
        np.cos(theta) * np.cos(theta_obs)
        + np.sin(theta) * np.sin(theta_obs) * np.cos(phi - phi_obs)
    )
    return np.arccos(np.clip(cos_sep, -1.0, 1.0))


def _solid_angle_cells(theta_edges: np.ndarray, phi_edges: np.ndarray) -> np.ndarray:
    dcos = np.cos(theta_edges[:-1]) - np.cos(theta_edges[1:])
    dphi = np.diff(phi_edges)
    return dcos[:, None] * dphi[None, :]
