"""RS hadronic smoke test: verify RS proton transport + synchrotron."""
from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core import Model
from tests.public_api_builders import hadronic, numerics, radiation, reverse_shock, solver_options, top_hat_model


def _build_rs_model(rs_eps_p: float = 0.0, full_chain: bool = False) -> Model:
    return top_hat_model(
        fwd_rad=radiation(
            proton_energy_fraction=0.0,
            proton_synch=True,
            reverse_proton_energy_fraction=rs_eps_p,
            include_pgamma=full_chain,
            bethe_heitler=full_chain,
            pp=full_chain,
            neutrino=full_chain,
            pgamma_scheme="hummer_2010_response" if full_chain else "disabled",
        ),
        rvs_rad=radiation(proton_energy_fraction=0.0, proton_synch=True),
        numerics=numerics(
            num_electron_gamma=24,
            num_photon_frequency=24,
            num_radius=24,
            num_theta=8,
            num_observer_time=24,
        ),
        solver_options=solver_options(electron_solver="fullhide_1d"),
        reverse_shock=reverse_shock(enabled=True, shell_duration_s=10.0),
        hadronic=hadronic(
            enabled=False,
            num_proton_gamma=24,
            num_neutrino_frequency=16,
            pgamma_scheme="hummer_2010_response" if full_chain else "disabled",
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

    # 4. RS full hadronic chain: p-gamma/BH/pp/secondary paths should run through the formal 1D kernels.
    full = _build_rs_model(0.2, full_chain=True).flux_density(
        np.array([1e3, 1e5]), np.array([1e14, 1e18]),
    ).total
    assert np.all(np.isfinite(full))
    assert np.any(full > base)

    print(f"base  flux: {base}")
    print(f"had   flux: {had}")
    print(f"delta flux: {diff}")
    print(f"full  flux: {full}")
    print("RS hadronic smoke test PASSED")


if __name__ == "__main__":
    main()
