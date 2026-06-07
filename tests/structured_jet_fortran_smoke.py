from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


TIMES = np.logspace(4.0, 5.0, 4)
FREQS = np.array([1.0e10, 1.0e14])


def _run_worker(args: list[str]) -> dict:
    command = [sys.executable, __file__, "--worker", *args]
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    if result.returncode != 0:
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise SystemExit(result.returncode)
    return json.loads(result.stdout)


def _worker_model(backend: str, patch_phi: int, jet_kind: str):
    from ASGARD import Ejecta, GaussianJet, ISM, Model, Observer, Radiation, Setups

    if jet_kind == "gaussian":
        jet = GaussianJet(E_iso=1.0e52, Gamma0=120.0, theta_c=0.08, theta_max=0.24)
    elif jet_kind == "ejecta":
        jet = Ejecta(
            E_iso=lambda phi, theta: 1.0e52 * (1.0 + 0.2 * np.cos(phi)) * np.exp(-0.5 * (theta / 0.08) ** 2),
            Gamma0=lambda phi, theta: 1.0
            + 119.0 * (1.0 + 0.1 * np.sin(phi)) * np.exp(-0.5 * (theta / 0.08) ** 2),
            theta_max=0.24,
        )
    else:
        raise ValueError(f"unknown jet kind: {jet_kind}")
    return Model(
        jet,
        ISM(1.0),
        Observer(1.0e26, 0.1, 0.06),
        Radiation(0.1, 1.0e-3, 2.3, ssc=False),
        setups=Setups(
            structured_backend=backend,
            num_gam_e=21,
            num_nu=17,
            num_r=16,
            num_theta=12,
            num_phi=12,
            num_tobs=12,
            patch_theta=2,
            patch_phi=patch_phi,
            num_threads=1,
            electron_solver="fullhide_1d",
            electron_adaptive_substeps=False,
        ),
    )


def _worker(args: argparse.Namespace) -> None:
    model = _worker_model(args.backend, args.patch_phi, args.jet)
    start = time.perf_counter()
    flux = np.asarray(model.flux_density_grid(TIMES, FREQS).total, dtype=float)
    elapsed = time.perf_counter() - start
    print(json.dumps({"flux": flux.tolist(), "elapsed_s": elapsed}))


def _relative_max(lhs: np.ndarray, rhs: np.ndarray) -> float:
    mask = np.isfinite(lhs) & np.isfinite(rhs) & (lhs > 0.0) & (rhs > 0.0)
    if not np.any(mask):
        raise AssertionError("no positive finite flux samples for comparison")
    return float(np.max(np.abs(lhs[mask] - rhs[mask]) / np.maximum(np.abs(rhs[mask]), 1.0e-300)))


def _assert_smooth_positive(flux: np.ndarray) -> None:
    if not np.all(np.isfinite(flux)):
        raise AssertionError("structured flux contains non-finite values")
    if np.any(flux <= 0.0):
        raise AssertionError("structured flux contains non-positive values")
    log_flux = np.log(flux)
    adjacent_jump = np.max(np.abs(np.diff(log_flux, axis=1)))
    if adjacent_jump > 4.0:
        raise AssertionError(f"structured light curve has a large adjacent log jump: {adjacent_jump:.3g}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--backend", default="fortran_1d")
    parser.add_argument("--patch-phi", type=int, default=4)
    parser.add_argument("--jet", choices=("gaussian", "ejecta"), default="gaussian")
    args = parser.parse_args()
    if args.worker:
        _worker(args)
        return

    fortran_gaussian = np.asarray(
        _run_worker(["--backend", "fortran_1d", "--patch-phi", "4", "--jet", "gaussian"])["flux"],
        dtype=float,
    )
    python_gaussian = np.asarray(
        _run_worker(["--backend", "python_patch", "--patch-phi", "4", "--jet", "gaussian"])["flux"],
        dtype=float,
    )
    rel = _relative_max(fortran_gaussian, python_gaussian)
    if rel > 0.30:
        raise AssertionError(f"structured Fortran/Python Gaussian flux mismatch: relmax={rel:.3g}")
    _assert_smooth_positive(fortran_gaussian)

    ejecta = np.asarray(
        _run_worker(["--backend", "fortran_1d", "--patch-phi", "4", "--jet", "ejecta"])["flux"],
        dtype=float,
    )
    _assert_smooth_positive(ejecta)

    fast = _run_worker(["--backend", "fortran_1d", "--patch-phi", "4", "--jet", "gaussian"])["elapsed_s"]
    wide_phi = _run_worker(["--backend", "fortran_1d", "--patch-phi", "12", "--jet", "gaussian"])["elapsed_s"]
    if wide_phi > 2.5 * fast:
        raise AssertionError(f"axisymmetric structured Fortran runtime scales too strongly with patch_phi: {wide_phi/fast:.2f}x")


if __name__ == "__main__":
    main()
