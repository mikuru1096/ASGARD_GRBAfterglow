from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core import tabulated_angular_jet
from asgard_core.structured_jet_kernel import _structured_threads
from tests.public_api_builders import numerics, radiation, reverse_shock, solver_options, top_hat_model


def _model(solver: str):
    return top_hat_model(
        fwd_rad=radiation(include_ssc=False, epsilon_B=1.0e-3),
        rvs_rad=radiation(include_ssc=False, epsilon_B=1.0e-2, p=2.4, accelerated_electron_fraction=1.0),
        numerics=numerics(
            num_radius=24,
            eats_num_theta=6,
            eats_num_phi=1,
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


def _tabulated_model():
    return top_hat_model(
        jet=tabulated_angular_jet(
            theta_rad=np.array([0.0, 0.05, 0.12], dtype=float),
            energy_iso_erg=np.array([1.0e52, 3.0e51, 1.0e50], dtype=float),
            lorentz_factor=np.array([150.0, 80.0, 5.0], dtype=float),
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        fwd_rad=radiation(include_ssc=False, epsilon_B=1.0e-3),
        rvs_rad=None,
        numerics=numerics(
            num_radius=24,
            eats_num_theta=6,
            eats_num_phi=1,
            num_observer_time=6,
            num_electron_gamma=61,
            num_photon_frequency=21,
            num_threads=1,
        ),
        solver_options=solver_options(
            structured_backend="fortran_1d",
            ssc_cooling_mode="none",
            synchrotron_integration="fixed_grid",
        ),
        reverse_shock=reverse_shock(enabled=False),
    )


def _transrelativistic_tabulated_model():
    return top_hat_model(
        jet=tabulated_angular_jet(
            theta_rad=np.array([0.0, 0.05, 0.12], dtype=float),
            energy_iso_erg=np.array([1.0e50, 8.0e49, 5.0e49], dtype=float),
            lorentz_factor=np.array([1.5, 1.4, 1.3], dtype=float),
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        fwd_rad=radiation(include_ssc=False, epsilon_B=1.0e-3),
        rvs_rad=None,
        numerics=numerics(
            num_radius=24,
            eats_num_theta=6,
            eats_num_phi=1,
            num_observer_time=6,
            num_electron_gamma=61,
            num_photon_frequency=21,
            num_threads=1,
        ),
        solver_options=solver_options(
            structured_backend="fortran_1d",
            ssc_cooling_mode="none",
            synchrotron_integration="fixed_grid",
        ),
        reverse_shock=reverse_shock(enabled=False),
    )


def main() -> None:
    times = np.array([1.0e2, 1.0e4], dtype=float)
    nu = np.array([1.0e9, 1.0e14], dtype=float)
    threaded_reverse = _model("fullhide_1d")
    threaded_reverse.setups.num_threads = 4
    assert _structured_threads(threaded_reverse) == (1, 4)

    for solver in ("fullhide_1d", "dg_1d"):
        flux = _model(solver).flux_density_grid(times, nu)
        assert np.all(np.isfinite(flux.total))
        assert np.all(np.isfinite(flux.rev.sync))
        assert np.max(flux.total) > 0.0
        assert np.max(flux.rev.sync) > 0.0
    tabulated_flux = _tabulated_model().flux_density_grid(times, nu)
    assert np.all(np.isfinite(tabulated_flux.total))
    assert np.max(tabulated_flux.total) > 0.0

    dense_times = np.logspace(2.0, 4.0, 9)
    dense_flux = _tabulated_model().flux_density_grid(dense_times, nu[:1])
    assert np.all(np.isfinite(dense_flux.total))
    assert dense_flux.total.shape == (1, dense_times.size)

    try:
        _transrelativistic_tabulated_model().flux_density_grid(times, nu[:1])
    except ValueError as exc:
        assert "Gamma0 >= 2" in str(exc)
    else:
        raise AssertionError("structured fortran backend accepted transrelativistic active patches")
    print("structured-shared-solver-smoke-ok")


if __name__ == "__main__":
    main()
