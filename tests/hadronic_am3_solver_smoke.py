from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core import Model
from tests.public_api_builders import hadronic, numerics, radiation, solver_options, top_hat_model


def _build_model(*, include_pg: bool = False, include_neutrino: bool = False) -> Model:
    return top_hat_model(
        fwd_rad=radiation(
            proton_energy_fraction=0.2,
            include_ssc=True,
            proton_synch=True,
            include_pgamma=include_pg,
            neutrino=include_neutrino,
        ),
        numerics=numerics(
            num_electron_gamma=24,
            num_photon_frequency=40,
            num_radius=24,
            num_theta=16,
            num_observer_time=24,
        ),
        solver_options=solver_options(electron_solver="fullhide_1d"),
        hadronic=hadronic(
            enabled=True,
            solver="am3_1d",
            num_proton_gamma=40,
            num_neutrino_frequency=24,
        ),
    )


def main() -> None:
    details = _build_model().details(1.0e3, 1.0e6).fwd
    assert details.gamma_p is not None
    assert details.dN_dgamma_p is not None
    assert details.hadronic_gamma is not None
    assert np.all(np.isfinite(details.hadronic_gamma))
    assert details.neutrino_luminosity is None


if __name__ == "__main__":
    main()
