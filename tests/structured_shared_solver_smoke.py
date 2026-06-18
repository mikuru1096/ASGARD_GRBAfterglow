from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from tests.public_api_builders import numerics, radiation, reverse_shock, solver_options, top_hat_model


def _model(solver: str):
    return top_hat_model(
        fwd_rad=radiation(include_ssc=False, epsilon_B=1.0e-3),
        rvs_rad=radiation(include_ssc=False, epsilon_B=1.0e-2, p=2.4, accelerated_electron_fraction=1.0),
        numerics=numerics(
            num_radius=24,
            num_theta=6,
            num_phi=1,
            num_observer_time=6,
            num_electron_gamma=121,
            num_photon_frequency=21,
            num_threads=1,
        ),
        solver_options=solver_options(
            electron_solver=solver,
            structured_backend="fortran_1d",
            ssc_cooling_mode="none",
            synchrotron_integration="fixed_grid",
        ),
        reverse_shock=reverse_shock(enabled=True, shell_duration_s=10.0),
    )


def main() -> None:
    times = np.array([1.0e2, 1.0e4], dtype=float)
    nu = np.array([1.0e9, 1.0e14], dtype=float)
    for solver in ("fullhide_1d", "dg_1d"):
        flux = _model(solver).flux_density_grid(times, nu)
        assert np.all(np.isfinite(flux.total))
        assert np.all(np.isfinite(flux.rev.sync))
        assert np.max(flux.total) > 0.0
        assert np.max(flux.rev.sync) > 0.0
    print("structured-shared-solver-smoke-ok")


if __name__ == "__main__":
    main()
