from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ASGARD import Model, Observer, top_hat_jet
from tests.public_api_builders import numerics, radiation, solver_options, top_hat_model


def _model(coupling: str | None) -> Model:
    return top_hat_model(
        jet=top_hat_jet(
            energy_iso_erg=1.0e52,
            initial_lorentz_factor=100.0,
            opening_angle_rad=0.1,
            shell_duration_s=None,
            magnetar=None,
            spreading=False,
        ),
        observer=Observer(z=0.05, viewing_angle_rad=0.0, viewing_azimuth_rad=0.0, luminosity_distance_cm=1.0e26),
        fwd_rad=radiation(include_ssc=True),
        numerics=numerics(
            num_electron_gamma=18,
            num_photon_frequency=18,
            num_radius=14,
            num_theta=8,
            num_observer_time=14,
        ),
        solver_options=solver_options(
            electron_solver="fullhide_1d",
            electron_photon_coupling="separated" if coupling is None else coupling,
        ),
    )


def main() -> None:
    times = np.array([1.0e4, 1.0e5], dtype=float)
    freqs = np.array([1.0e14, 1.0e18], dtype=float)
    default_flux = _model(None).flux_density(times, freqs).total
    explicit_flux = _model("separated").flux_density(times, freqs).total
    np.testing.assert_allclose(explicit_flux, default_flux, rtol=1.0e-14, atol=0.0)
    print("electron_photon_separated_regression_smoke: ok")


if __name__ == "__main__":
    main()
