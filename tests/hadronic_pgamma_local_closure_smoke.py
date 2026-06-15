from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import Model, Observer, UniformMedium, top_hat_jet
from ASGARD.api_model import _build_fit_config_for_patch
from asgard_core.asgard_physics_utils import compute_magnetic_field
from asgard_core.asgard_state import _forward_synchrotron_absorption_transfer, solve_state
from tests.public_api_builders import hadronic, numerics, radiation, observer_grid, solver_options, top_hat_model


def _build_model() -> Model:
    return top_hat_model(
        jet=top_hat_jet(
            energy_iso_erg=3.0e53,
            initial_lorentz_factor=400.0,
            opening_angle_rad=0.12,
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        medium=UniformMedium(number_density_cm3=30.0),
        observer=Observer(z=0.2, viewing_angle_rad=0.0, viewing_azimuth_rad=0.0, luminosity_distance_cm=1.0e27),
        fwd_rad=radiation(
            epsilon_e=0.15,
            epsilon_B=0.03,
            p=2.2,
            proton_energy_fraction=0.6,
            include_ssc=True,
            proton_synch=True,
            include_pgamma=True,
            neutrino=True,
            pgamma_scheme="hummer_2010_response",
        ),
        numerics=numerics(
            num_electron_gamma=28,
            num_photon_frequency=40,
            num_radius=26,
            num_theta=16,
            num_observer_time=26,
        ),
        observer_grid=observer_grid(time_min_s=1.0e2, time_max_s=1.0e5),
        solver_options=solver_options(electron_solver="fullhide_1d"),
        hadronic=hadronic(
            enabled=True,
            solver="am3_1d",
            num_proton_gamma=36,
            num_neutrino_frequency=24,
            pgamma_scheme="hummer_2010_response",
        ),
    )


def main() -> None:
    model = _build_model()
    config = _build_fit_config_for_patch(
        model,
        phi_center=0.0,
        theta_v=model.observer.theta_obs,
        opening_angle_jet=model.jet.theta_j,
        e_iso=model.jet.E_iso,
        gamma0=model.jet.lf,
        theta_center=0.0,
    )
    times_s = np.logspace(2.0, 5.0, 26)
    state = solve_state(config, times_s)

    if state.hadronic is None:
        raise AssertionError("Expected hadronic state for local photopion-closure smoke test.")

    tau_pg = state.hadronic.tau_pg
    survival = state.hadronic.pg_photon_survival
    if tau_pg is None or survival is None:
        raise AssertionError("Expected photopion local-closure diagnostics in hadronic state.")

    tau_pg = np.asarray(tau_pg, dtype=float)
    survival = np.asarray(survival, dtype=float)
    if not np.all(np.isfinite(tau_pg)):
        raise AssertionError("Photopion local tau contains non-finite values.")
    if not np.all(np.isfinite(survival)):
        raise AssertionError("Photopion local survival contains non-finite values.")
    if not np.all((survival > 0.0) & (survival <= 1.0)):
        raise AssertionError("Photopion local survival must stay within (0, 1].")
    if not np.any(tau_pg > 0.0):
        raise AssertionError("Photopion local tau never became positive.")
    if not np.allclose(state.observer.tau_extra, 0.0, rtol=0.0, atol=1.0e-14):
        raise AssertionError("Observer tau_extra still contains photopion tau after local closure migration.")

    hadronic_gamma = np.asarray(state.hadronic.l_had_pg_gamma, dtype=float)
    if not np.all(np.isfinite(hadronic_gamma)):
        raise AssertionError("Photopion gamma luminosity contains non-finite values after local closure.")
    if float(np.max(hadronic_gamma)) <= 0.0:
        raise AssertionError("Photopion gamma luminosity vanished in the local-closure smoke test.")

    transfer = _forward_synchrotron_absorption_transfer(
        electron=state.electron,
        radius_cm=state.dynamics.radius,
        magnetic_field_g=compute_magnetic_field(state.dynamics.r_gamma, state.dynamics.radius, config),
        seed_frequency_hz=state.setup.seed_frequency_hz,
        config=config,
    )
    if not np.all(np.isfinite(transfer)):
        raise AssertionError("Hadronic observer SSA transfer contains non-finite values.")
    if float(np.min(transfer)) < 0.0 or float(np.max(transfer)) > 1.0 + 1.0e-12:
        raise AssertionError("Hadronic observer SSA transfer must stay within [0, 1].")

    photon_frequency = np.asarray(state.setup.seed_frequency_hz, dtype=float)
    opt_to_gev = (photon_frequency >= 1.0e14) & (photon_frequency <= 2.5e23)
    if np.any(opt_to_gev) and float(np.min(transfer[opt_to_gev, :])) < 0.9:
        raise AssertionError("Hadronic observer SSA transfer collapsed the optical-GeV band.")

    observer_hadronic = state.observer.components.fwd_hadronic_gamma
    if observer_hadronic is None or not np.all(np.isfinite(observer_hadronic)):
        raise AssertionError("Observed hadronic gamma component contains non-finite values.")

    print("hadronic_pgamma_local_closure_smoke: ok")


if __name__ == "__main__":
    main()
