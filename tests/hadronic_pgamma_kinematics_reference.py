from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_core.hadronic_pgamma import (
    ETA0,
    _x_pm_gamma,
    _x_two_pion_bounds,
    kelner_aharonian_2008_gamma_phi,
)


KINEMATICS_REFERENCE = {
    "eta_gamma_ratio": np.array([1.1, 1.6, 2.0, 4.0, 10.0, 30.0], dtype=float),
    "x_pm_gamma_minus": np.array(
        [
            8.443219597711255e-02,
            4.717072237800911e-02,
            3.617644083742078e-02,
            1.724602718265588e-02,
            6.816727979421467e-03,
            2.267812306298499e-03,
        ],
        dtype=float,
    ),
    "x_pm_gamma_plus": np.array(
        [
            1.876740473829344e-01,
            3.008763364872557e-01,
            3.620932145710073e-01,
            5.483450802153742e-01,
            7.563451205767892e-01,
            9.035998297781394e-01,
        ],
        dtype=float,
    ),
    "eta_two_pion_ratio": np.array([3.0, 4.0, 6.0, 10.0, 30.0], dtype=float),
    "x_two_pion_minus": np.array(
        [
            3.698682612884475e-02,
            2.346027985599794e-02,
            1.376466741636708e-02,
            7.584547953156292e-03,
            2.346242976458196e-03,
        ],
        dtype=float,
    ),
    "x_two_pion_plus": np.array(
        [
            2.969680589917278e-01,
            4.030972442322421e-01,
            5.375558806042156e-01,
            6.797766956419868e-01,
            8.733941175322963e-01,
        ],
        dtype=float,
    ),
    "eta_gamma_phi_ratio": np.array([1.1, 2.0, 10.0], dtype=float),
    "x_gamma_phi": np.array([1.0e-4, 1.0e-1, 3.0e-1], dtype=float),
    "gamma_phi": np.array(
        [
            1.128135599792708e-19,
            1.6289759016147586e-17,
            3.6353018962795205e-18,
        ],
        dtype=float,
    ),
}

KINEMATICS_TOLERANCE = 5.0e-12


def _relative_error(actual: np.ndarray, reference: np.ndarray) -> np.ndarray:
    denom = np.maximum(np.abs(reference), 1.0e-300)
    return np.abs(actual - reference) / denom


def evaluate_kinematics_reference() -> dict[str, object]:
    eta_gamma = KINEMATICS_REFERENCE["eta_gamma_ratio"] * ETA0
    eta_two = KINEMATICS_REFERENCE["eta_two_pion_ratio"] * ETA0

    x_minus_gamma, x_plus_gamma = _x_pm_gamma(eta_gamma)
    x_minus_two, x_plus_two = _x_two_pion_bounds(eta_two)

    eta_phi = KINEMATICS_REFERENCE["eta_gamma_phi_ratio"] * ETA0
    x_phi = KINEMATICS_REFERENCE["x_gamma_phi"]
    gamma_phi = kelner_aharonian_2008_gamma_phi(eta_phi, x_phi)

    rel_xminus_gamma = _relative_error(x_minus_gamma, KINEMATICS_REFERENCE["x_pm_gamma_minus"])
    rel_xplus_gamma = _relative_error(x_plus_gamma, KINEMATICS_REFERENCE["x_pm_gamma_plus"])
    rel_xminus_two = _relative_error(x_minus_two, KINEMATICS_REFERENCE["x_two_pion_minus"])
    rel_xplus_two = _relative_error(x_plus_two, KINEMATICS_REFERENCE["x_two_pion_plus"])
    rel_gamma_phi = _relative_error(gamma_phi, KINEMATICS_REFERENCE["gamma_phi"])

    rel_stack = np.concatenate([rel_xminus_gamma, rel_xplus_gamma, rel_xminus_two, rel_xplus_two, rel_gamma_phi])
    metrics = {
        "max_rel_error": float(np.max(rel_stack)),
        "mean_rel_error": float(np.mean(rel_stack)),
        "gamma_xminus_max_rel": float(np.max(rel_xminus_gamma)),
        "gamma_xplus_max_rel": float(np.max(rel_xplus_gamma)),
        "two_pion_xminus_max_rel": float(np.max(rel_xminus_two)),
        "two_pion_xplus_max_rel": float(np.max(rel_xplus_two)),
        "gamma_phi_max_rel": float(np.max(rel_gamma_phi)),
    }
    return {
        "reference": KINEMATICS_REFERENCE,
        "actual": {
            "x_pm_gamma_minus": x_minus_gamma,
            "x_pm_gamma_plus": x_plus_gamma,
            "x_two_pion_minus": x_minus_two,
            "x_two_pion_plus": x_plus_two,
            "gamma_phi": gamma_phi,
        },
        "metrics": metrics,
    }


def assert_kinematics_reference(payload: dict[str, object]) -> None:
    metrics = payload["metrics"]
    if float(metrics["max_rel_error"]) > KINEMATICS_TOLERANCE:
        raise AssertionError(
            f"KA2008 kinematics regression max_rel_error={metrics['max_rel_error']:.3e} exceeds {KINEMATICS_TOLERANCE:.3e}"
        )


def main() -> None:
    payload = evaluate_kinematics_reference()
    metrics = payload["metrics"]
    print(f"kinematics_max_rel_error={metrics['max_rel_error']:.3e}")
    print(f"kinematics_mean_rel_error={metrics['mean_rel_error']:.3e}")
    assert_kinematics_reference(payload)


if __name__ == "__main__":
    main()
