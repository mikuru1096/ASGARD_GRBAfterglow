from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core import Model
from asgard_core.api_model import _build_fit_config_for_patch, _solve_patch_state
from tests.public_api_builders import hadronic, numerics, radiation, solver_options, top_hat_model


def _build_model(enabled: bool, epsilon_p: float, include_pg: bool = False, include_neutrino: bool = False) -> Model:
    hadronic_solver = "am3_1d" if (include_pg or include_neutrino) else "legacy_1d"
    pgamma_scheme = "hummer_2010_response" if (include_pg or include_neutrino) else "disabled"
    return top_hat_model(
        fwd_rad=radiation(
            proton_energy_fraction=epsilon_p,
            include_ssc=True,
            proton_synch=True,
            include_pgamma=include_pg,
            neutrino=include_neutrino,
            pgamma_scheme=pgamma_scheme,
        ),
        numerics=numerics(
            num_electron_gamma=24,
            num_photon_frequency=32,
            num_radius=24,
            eats_num_theta=16,
            num_observer_time=24,
        ),
        solver_options=solver_options(electron_solver="fullhide_1d"),
        hadronic=hadronic(
            enabled=enabled,
            solver=hadronic_solver,
            num_proton_gamma=32,
            num_neutrino_frequency=24,
            pgamma_scheme=pgamma_scheme,
        ),
    )


def main() -> None:
    times = np.array([1.0e4, 1.0e5], dtype=float)
    freqs = np.array([1.0e14, 1.0e18], dtype=float)

    base = _build_model(False, 0.3).flux_density(times, freqs).total
    disabled = _build_model(False, 0.3).flux_density(times, freqs).total
    rel = np.max(np.abs(disabled - base) / np.maximum(np.abs(base), 1.0e-99))
    assert float(rel) < 1.0e-10

    zero_eps = _build_model(True, 0.0).flux_density(times, freqs).total
    rel_zero = np.max(np.abs(zero_eps - base) / np.maximum(np.abs(base), 1.0e-99))
    assert float(rel_zero) < 1.0e-10

    had = _build_model(True, 0.3).flux_density(times, freqs).total
    assert np.all(np.isfinite(had))
    assert np.max(np.abs(had - base)) > 0.0

    details = _build_model(True, 0.3).details(1.0e3, 1.0e6)
    assert details.fwd.gamma_p is not None
    assert details.fwd.dN_dgamma_p is not None
    assert details.fwd.hadronic_gamma is not None
    assert np.all(np.isfinite(details.fwd.dN_dgamma_p))
    assert np.all(np.isfinite(details.fwd.hadronic_gamma))
    assert details.fwd.neutrino_luminosity is None

    pg_details = _build_model(True, 0.3, include_pg=True, include_neutrino=True).details(1.0e3, 1.0e6)
    assert pg_details.fwd.neutrino_luminosity is not None
    assert np.all(np.isfinite(pg_details.fwd.neutrino_luminosity))

    pg_model = _build_model(True, 0.3, include_pg=True, include_neutrino=False)
    pg_config = _build_fit_config_for_patch(
        pg_model,
        theta_v=pg_model.observer.theta_obs,
        opening_angle_jet=pg_model.jet.theta_j,
        e_iso=pg_model.jet.E_iso,
        gamma0=pg_model.jet.lf,
        theta_center=0.0,
    )
    pg_state = _solve_patch_state(pg_model, pg_config, times, freqs)
    assert pg_state.hadronic is not None
    assert pg_state.hadronic.secondary_electron_source_r is None


if __name__ == "__main__":
    main()
