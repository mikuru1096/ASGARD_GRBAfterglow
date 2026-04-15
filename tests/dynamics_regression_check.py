from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from asgard_presets import build_baseline_config, build_reverse_demo_config
from asgard_setup import build_simulation_setup
from src import Dynamics


def build_forward_config():
    return build_baseline_config(num_r=128)


def build_reverse_config():
    return build_reverse_demo_config(num_r=128, num_gam_e=96)


def write_baseline(path: Path) -> None:
    forward_config = build_forward_config()
    forward_setup = build_simulation_setup(forward_config)
    forward = Dynamics.dynamics_forward(forward_setup.boundary, forward_config.num_r, forward_config.index_dyn)

    reverse_config = build_reverse_config()
    reverse_setup = build_simulation_setup(reverse_config)
    reverse = Dynamics.dynamics_reverse(
        reverse_config.reverse_shock.delta_t_s,
        reverse_config.reverse_shock.epsilon_e,
        reverse_config.reverse_shock.epsilon_b,
        reverse_config.reverse_shock.p,
        reverse_config.reverse_shock.f_e,
        reverse_setup.boundary,
        reverse_config.num_r,
        reverse_config.num_gam_e,
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        forward_r_tobs=forward[0],
        forward_r_gamma=forward[1],
        forward_radius=forward[2],
        forward_mass=forward[3],
        reverse_t_cross=np.array(reverse[0]),
        reverse_r_cross=np.array(reverse[1]),
        reverse_e3_cross=np.array(reverse[2]),
        reverse_gam20=np.array(reverse[3]),
        reverse_r_tobs=reverse[4],
        reverse_r_gamma=reverse[5],
        reverse_radius=reverse[6],
        reverse_m2=reverse[7],
        reverse_m3=reverse[8],
    )


def check_against_baseline(path: Path, rtol: float, atol: float) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Baseline file not found: {path}")

    baseline = np.load(path, allow_pickle=False)

    forward_config = build_forward_config()
    forward_setup = build_simulation_setup(forward_config)
    forward = Dynamics.dynamics_forward(forward_setup.boundary, forward_config.num_r, forward_config.index_dyn)

    reverse_config = build_reverse_config()
    reverse_setup = build_simulation_setup(reverse_config)
    reverse = Dynamics.dynamics_reverse(
        reverse_config.reverse_shock.delta_t_s,
        reverse_config.reverse_shock.epsilon_e,
        reverse_config.reverse_shock.epsilon_b,
        reverse_config.reverse_shock.p,
        reverse_config.reverse_shock.f_e,
        reverse_setup.boundary,
        reverse_config.num_r,
        reverse_config.num_gam_e,
    )

    _assert_close("forward_r_tobs", forward[0], baseline["forward_r_tobs"], rtol, atol)
    _assert_close("forward_r_gamma", forward[1], baseline["forward_r_gamma"], rtol, atol)
    _assert_close("forward_radius", forward[2], baseline["forward_radius"], rtol, atol)
    _assert_close("forward_mass", forward[3], baseline["forward_mass"], rtol, atol)

    _assert_scalar_close("reverse_t_cross", reverse[0], float(baseline["reverse_t_cross"]), rtol, atol)
    _assert_scalar_close("reverse_r_cross", reverse[1], float(baseline["reverse_r_cross"]), rtol, atol)
    _assert_scalar_close("reverse_e3_cross", reverse[2], float(baseline["reverse_e3_cross"]), rtol, atol)
    _assert_scalar_close("reverse_gam20", reverse[3], float(baseline["reverse_gam20"]), rtol, atol)
    _assert_close("reverse_r_tobs", reverse[4], baseline["reverse_r_tobs"], rtol, atol)
    _assert_close("reverse_r_gamma", reverse[5], baseline["reverse_r_gamma"], rtol, atol)
    _assert_close("reverse_radius", reverse[6], baseline["reverse_radius"], rtol, atol)
    _assert_close("reverse_m2", reverse[7], baseline["reverse_m2"], rtol, atol)
    _assert_close("reverse_m3", reverse[8], baseline["reverse_m3"], rtol, atol)
    print("PASS: dynamics regression check succeeded.")


def _assert_close(name: str, current: np.ndarray, baseline: np.ndarray, rtol: float, atol: float) -> None:
    if current.shape != baseline.shape:
        raise ValueError(f"{name} shape mismatch: {current.shape} != {baseline.shape}")
    if not np.allclose(current, baseline, rtol=rtol, atol=atol):
        diff = float(np.max(np.abs(current - baseline)))
        raise AssertionError(f"{name} mismatch: max abs diff = {diff:.6e}")


def _assert_scalar_close(name: str, current: float, baseline: float, rtol: float, atol: float) -> None:
    if not np.isclose(current, baseline, rtol=rtol, atol=atol):
        diff = abs(current - baseline)
        raise AssertionError(f"{name} mismatch: abs diff = {diff:.6e}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Regression check for forward/reverse dynamics outputs.")
    parser.add_argument("--baseline", type=Path, default=Path("tests") / "baseline_dynamics.npz")
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    parser.add_argument("--atol", type=float, default=1.0e-12)
    parser.add_argument("--write-baseline", action="store_true")
    args = parser.parse_args()

    if args.write_baseline:
        write_baseline(args.baseline)
        print(f"Baseline written: {args.baseline}")
        return

    check_against_baseline(args.baseline, args.rtol, args.atol)


if __name__ == "__main__":
    main()
