from __future__ import annotations

from _repo_path import ensure_repo_root_on_path


ensure_repo_root_on_path()

from tests.public_api_builders import hadronic, numerics, radiation, solver_options, top_hat_model


def main() -> None:
    model = top_hat_model(
        fwd_rad=radiation(proton_energy_fraction=0.2, proton_synch=True, include_pgamma=True, neutrino=True),
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
            pgamma_scheme="ka2008_reference",
            num_proton_gamma=40,
            num_neutrino_frequency=24,
        ),
    )
    try:
        model.details(1.0e3, 1.0e6)
    except ValueError as exc:
        message = str(exc)
        assert "tabulated" in message
        assert "eta/eta0" in message
        return
    raise AssertionError("Expected KA2008 domain guard for pgamma_scheme='ka2008_reference'.")


if __name__ == "__main__":
    main()
