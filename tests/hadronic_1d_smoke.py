from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet
from ASGARD.api_model import _build_fit_config_for_patch, _solve_patch_state


def _build_model(enabled: bool, epsilon_p: float, include_pg: bool = False, include_neutrino: bool = False) -> Model:
    hadronic_solver = "am3_1d" if (include_pg or include_neutrino) else "legacy_1d"
    pgamma_scheme = "hummer_2010_response" if (include_pg or include_neutrino) else "disabled"
    return Model(
        TophatJet(0.1, 1.0e52, 300.0),
        ISM(1.0),
        Observer(1.0e26, 0.1, 0.0),
        Radiation(
            0.1,
            1.0e-3,
            2.3,
            epsilon_p=epsilon_p,
            ssc=True,
            proton_synch=True,
            pg=include_pg,
            neutrino=include_neutrino,
            pgamma_scheme=pgamma_scheme,
        ),
        setups=Setups(
            electron_solver="fullhide_1d",
            num_gam_e=24,
            num_gam_p=32,
            num_nu=32,
            num_nu_nu=24,
            num_r=24,
            num_theta=16,
            num_tobs=24,
            hadronic_enabled=enabled,
            hadronic_solver=hadronic_solver,
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
        phi_center=0.0,
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
