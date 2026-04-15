from __future__ import annotations

from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import VegasAfterglow.api as va_api
from VegasAfterglow import GaussianJet, ISM, Model, Observer, Radiation, Setups, TophatJet


def _build_direct_model() -> Model:
    return Model(
        medium=ISM(0.1),
        jet=TophatJet(E_iso=1.0e53, Gamma0=100.0, theta_c=0.1),
        observer=Observer(lumi_dist=1.0e28, z=0.4, theta_obs=0.0, phi_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.5, xi_N=0.1, ssc=True, kn=False),
        setups=Setups(
            num_threads=8,
            num_gam_e=121,
            num_nu=121,
            num_r=140,
            num_theta=120,
            num_phi=1,
            num_tobs=64,
            observer_time_min_s=1.0e2,
            observer_time_max_s=1.0e6,
        ),
    )


def _build_structured_model() -> Model:
    return Model(
        medium=ISM(0.1),
        jet=GaussianJet(E_iso=1.0e53, Gamma0=100.0, theta_c=0.05, theta_max=0.12),
        observer=Observer(lumi_dist=1.0e28, z=0.4, theta_obs=0.08, phi_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.5, xi_N=0.1, ssc=False, kn=False),
        setups=Setups(
            num_threads=8,
            num_gam_e=101,
            num_nu=101,
            num_r=100,
            num_theta=64,
            num_phi=1,
            num_tobs=48,
            patch_theta=2,
            patch_phi=3,
            observer_time_min_s=1.0e2,
            observer_time_max_s=1.0e5,
        ),
    )


def main() -> None:
    times_s = np.logspace(2.0, 5.0, 12)
    frequencies_hz = np.array([9.0e9, 4.84e14], dtype=float)

    direct_model = _build_direct_model()
    with patch.object(va_api, "solve_model_state_from_setup", wraps=va_api.solve_model_state_from_setup) as solve_mock:
        direct_model.flux_density_grid(times_s, frequencies_hz)
        solve_after_flux = solve_mock.call_count
        direct_model.spectrum(1.0e3, frequencies_hz)
        direct_model.details(1.0e2, 1.0e5)
        solve_after_reuse = solve_mock.call_count
    assert solve_after_flux == 1
    assert solve_after_reuse == solve_after_flux
    direct_states = next(iter(direct_model._state_cache.values()))
    assert len(direct_states) == 1
    assert direct_states[0].setup.observer_time_s.size == direct_model.setups.num_tobs

    structured_model = _build_structured_model()
    with patch.object(va_api, "solve_model_state_from_setup", wraps=va_api.solve_model_state_from_setup) as solve_mock:
        structured_model.flux_density_grid(times_s, frequencies_hz)
        solve_after_flux = solve_mock.call_count
        structured_model.details(1.0e2, 1.0e5)
        solve_after_reuse = solve_mock.call_count
    assert solve_after_flux > 0
    assert solve_after_reuse == solve_after_flux
    structured_states = next(iter(structured_model._state_cache.values()))
    assert all(state.setup.observer_time_s.size == structured_model.setups.num_tobs for state in structured_states)

    print("PASS: state reuse check succeeded.")


if __name__ == "__main__":
    main()
