from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.hadronic_pgamma import (
    PGammaMicrophysics,
    sigma_kappa_stecker_box,
    strict_isotropic_pgamma_loss_rate,
)


LOSS_REFERENCE = {
    "gamma_grid": np.array(
        [
            1.0000000000000000e02,
            4.6415888336127773e02,
            2.1544346900318824e03,
            1.0000000000000000e04,
            4.6415888336127726e04,
            2.1544346900318822e05,
            1.0000000000000000e06,
            4.6415888336127726e06,
            2.1544346900318824e07,
            1.0000000000000000e08,
        ],
        dtype=float,
    ),
    "loss_rate": np.array(
        [
            2.0840992611629027e-40,
            4.4902723459201695e-40,
            9.6741356973808052e-40,
            2.0842028519572350e-39,
            4.4903079243473938e-39,
            9.6741290071385000e-39,
            2.0842040916616518e-38,
            4.4903311800264170e-38,
            9.6740427112578652e-38,
            2.0842142615877767e-37,
        ],
        dtype=float,
    ),
}

LOSS_TOLERANCE = 5.0e-11


def _relative_error(actual: np.ndarray, reference: np.ndarray) -> np.ndarray:
    denom = np.maximum(np.abs(reference), 1.0e-300)
    return np.abs(actual - reference) / denom


def evaluate_loss_reference() -> dict[str, object]:
    gamma_grid = LOSS_REFERENCE["gamma_grid"]
    nu_hz = np.logspace(12.0, 24.0, 512)
    photon_density_hz = 1.0e-30 * (nu_hz / nu_hz[0]) ** (-1.5)

    rate = strict_isotropic_pgamma_loss_rate(
        gamma_grid,
        nu_hz,
        photon_density_hz,
        PGammaMicrophysics("stecker_box", sigma_kappa_stecker_box),
    )
    rel = _relative_error(rate, LOSS_REFERENCE["loss_rate"])
    log_rate = np.log10(np.maximum(rate, 1.0e-300))
    log_gamma = np.log10(gamma_grid)
    slope = np.diff(log_rate) / np.diff(log_gamma)

    metrics = {
        "max_rel_error": float(np.max(rel)),
        "mean_rel_error": float(np.mean(rel)),
        "min_rate": float(np.min(rate)),
        "max_rate": float(np.max(rate)),
        "slope_std": float(np.std(slope)),
    }
    return {
        "reference": LOSS_REFERENCE,
        "actual": {"gamma_grid": gamma_grid, "loss_rate": rate},
        "metrics": metrics,
    }


def assert_loss_reference(payload: dict[str, object]) -> None:
    actual = payload["actual"]["loss_rate"]
    metrics = payload["metrics"]
    if not np.all(np.isfinite(actual)):
        raise AssertionError("p-gamma loss benchmark returned non-finite values.")
    if not np.all(actual > 0.0):
        raise AssertionError("p-gamma loss benchmark expects strictly positive loss rates.")
    if not np.all(np.diff(actual) > 0.0):
        raise AssertionError("p-gamma loss benchmark expects monotonic increasing loss rate with gamma.")
    if float(metrics["max_rel_error"]) > LOSS_TOLERANCE:
        raise AssertionError(
            f"KA2008 loss regression max_rel_error={metrics['max_rel_error']:.3e} exceeds {LOSS_TOLERANCE:.3e}"
        )


def main() -> None:
    payload = evaluate_loss_reference()
    metrics = payload["metrics"]
    print(f"loss_max_rel_error={metrics['max_rel_error']:.3e}")
    print(f"loss_mean_rel_error={metrics['mean_rel_error']:.3e}")
    print(f"loss_slope_std={metrics['slope_std']:.3e}")
    assert_loss_reference(payload)


if __name__ == "__main__":
    main()
