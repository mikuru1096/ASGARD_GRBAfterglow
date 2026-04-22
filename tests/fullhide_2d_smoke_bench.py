from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet
from asgard_core.asgard_config import FitConfig
from asgard_core.asgard_state import solve_state
from tests._bench_common import run_case


NUM_GAM_E = 8
NUM_CHI = 4
NUM_NU = 8
NUM_R = 8
NUM_THETA = 8
NUM_TOBS = 4


def _build_model(solver: str) -> Model:
    return Model(
        TophatJet(0.1, 1.0e52, 300.0),
        ISM(1.0),
        Observer(1.0e26, 0.1, 0.0),
        Radiation(0.1, 1.0e-3, 2.3),
        setups=Setups(
            electron_solver=solver,
            num_gam_e=NUM_GAM_E,
            num_chi=NUM_CHI,
            num_nu=NUM_NU,
            num_r=NUM_R,
            num_theta=NUM_THETA,
            num_tobs=NUM_TOBS,
            electron_adaptive_substeps=False,
        ),
    )


def case_basic_smoke():
    model = _build_model("charint_2d")
    flux = model.flux_density(np.array([1.0e4]), np.array([1.0e14])).total
    assert flux.shape == (1,)
    assert np.all(np.isfinite(flux))
    return {"flux_shape": list(flux.shape), "flux": float(flux[0])}


def case_electron_grid():
    state = solve_state(
        FitConfig(
            electron_solver="charint_2d",
            num_gam_e=NUM_GAM_E,
            num_chi=NUM_CHI,
            num_nu=NUM_NU,
            num_r=NUM_R,
            num_theta=NUM_THETA,
            num_tobs=NUM_TOBS,
        ),
        np.array([1.0e2, 1.1e2]),
    )
    electron = state.electron
    chi_grid = np.asarray(electron.chi_grid, dtype=float)
    d_n_gam_e_chi = np.asarray(electron.d_n_gam_e_chi, dtype=float)
    d_n_gam_e = np.asarray(electron.d_n_gam_e, dtype=float)

    assert chi_grid.ndim == 1
    assert chi_grid.shape == (NUM_CHI,)
    assert np.all(np.isfinite(chi_grid))

    assert d_n_gam_e_chi.ndim == 3
    assert d_n_gam_e_chi.shape[0] == NUM_GAM_E
    assert d_n_gam_e_chi.shape[1] == NUM_CHI
    assert d_n_gam_e_chi.shape[2] == d_n_gam_e.shape[1]
    assert np.all(np.isfinite(d_n_gam_e_chi))

    return {
        "chi_shape": list(chi_grid.shape),
        "d_n_gam_e_chi_shape": list(d_n_gam_e_chi.shape),
        "d_n_gam_e_shape": list(d_n_gam_e.shape),
    }


def main() -> None:
    results = [
        run_case("[1/2] charint_2d:basic_smoke", case_basic_smoke),
        run_case("[2/2] charint_2d:electron_grid", case_electron_grid),
    ]

    failed = [item for item in results if item["status"] == "FAIL"]
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
