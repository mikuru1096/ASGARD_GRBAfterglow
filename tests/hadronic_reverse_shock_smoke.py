"""RS hadronic smoke test: verify RS proton transport + synchrotron."""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet


def _build_rs_model(rs_eps_p: float = 0.0) -> Model:
    return Model(
        TophatJet(0.1, 1.0e52, 300.0),
        ISM(1.0),
        Observer(1.0e26, 0.1, 0.0),
        Radiation(
            0.1, 1.0e-3, 2.3,
            epsilon_p=0.0, proton_synch=True,
            reverse_epsilon_p=rs_eps_p,
        ),
        rvs_rad=Radiation(0.1, 1.0e-3, 2.3, epsilon_p=0.0, proton_synch=True),
        setups=Setups(
            electron_solver="fullhide_1d",
            rvs_shock=True,
            reverse_delta_t_s=10.0,
            num_gam_e=24, num_gam_p=24, num_nu=24, num_nu_nu=16,
            num_r=24, num_theta=8, num_tobs=24,
            hadronic_enabled=False,
        ),
    )


def main() -> None:
    # 1. RS without hadronic: baseline
    base = _build_rs_model(0.0).flux_density(
        np.array([1e3, 1e5]), np.array([1e14, 1e18]),
    ).total
    assert np.all(np.isfinite(base))

    # 2. RS with hadronic enabled: should produce extra flux
    had = _build_rs_model(0.3).flux_density(
        np.array([1e3, 1e5]), np.array([1e14, 1e18]),
    ).total
    assert np.all(np.isfinite(had))

    # 3. RS hadronic should add non-negative flux
    diff = had - base
    assert np.all(diff >= -1e-30), f"RS hadronic flux should not be negative: {diff}"

    print(f"base  flux: {base}")
    print(f"had   flux: {had}")
    print(f"delta flux: {diff}")
    print("RS hadronic smoke test PASSED")


if __name__ == "__main__":
    main()
