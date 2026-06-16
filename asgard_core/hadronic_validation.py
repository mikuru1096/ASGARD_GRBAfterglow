import numpy as np


def as_strictly_increasing_grid(values: np.ndarray, name: str, *, require_finite: bool = False) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be a 1d array.")
    if arr.size < 2:
        raise ValueError(f"{name} must contain at least two points.")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive.")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing.")
    if require_finite and not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    return arr


def as_matching(values: np.ndarray, reference: np.ndarray, name: str, *, require_finite: bool = False) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != reference.shape:
        raise ValueError(f"{name} must match grid shape.")
    if require_finite and not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite values.")
    return arr


def as_matching_nonnegative(values: np.ndarray, reference: np.ndarray, name: str) -> np.ndarray:
    arr = as_matching(values, reference, name)
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must be non-negative.")
    return arr
