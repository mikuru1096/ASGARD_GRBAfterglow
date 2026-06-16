from __future__ import annotations

from _repo_path import ensure_repo_root_on_path

import numpy as np


ensure_repo_root_on_path()

from asgard_core.asgard_config import RuntimeConfig
from asgard_core.asgard_state import _compute_pair_production_branch
from asgard_core.asgard_types import DynamicsSolution, ElectronSolution
from asgard_core.hadronic_pair_production import ELECTRON_MASS_GEV
from asgard_core.api_model import _build_fit_config_for_patch
from tests.public_api_builders import hadronic, radiation, top_hat_model
from src import constants


def _state_inputs() -> tuple[DynamicsSolution, ElectronSolution, np.ndarray, np.ndarray, np.ndarray]:
    dln = 0.1
    n_ph = 84
    n_e = 56
    ind_min = 28
    x_ph = np.exp((-ind_min + np.arange(n_ph, dtype=float)) * dln)
    gam_e = np.exp(np.arange(n_e, dtype=float) * dln)
    frequency_hz = x_ph * ELECTRON_MASS_GEV / constants.para_h_gev
    photon_density_per_gev = 3.0e9 * np.power(x_ph, -2.1)
    seed_field_hz = photon_density_per_gev[:, None] * constants.para_h_gev * np.array([[1.0, 0.8]])
    radius = np.array([1.0e16, 1.2e16], dtype=float)
    dynamics = DynamicsSolution(
        r_tobs=np.array([100.0, 130.0], dtype=float),
        r_gamma=np.array([100.0, 95.0], dtype=float),
        radius=radius,
        swept_mass_g=np.ones_like(radius),
    )
    electron = ElectronSolution(
        gam_e=gam_e,
        d_n_gam_e=np.zeros((n_e, radius.size), dtype=float),
        l_syn_spec=np.zeros((n_ph, radius.size), dtype=float),
        seed_syn=np.zeros((n_ph, radius.size), dtype=float),
        nu_m=np.zeros(radius.size, dtype=float),
        nu_c=np.zeros(radius.size, dtype=float),
        nu_a=np.zeros(radius.size, dtype=float),
    )
    magnetic_field_g = np.full(radius.size, 1.0e10, dtype=float)
    return dynamics, electron, seed_field_hz, frequency_hz, magnetic_field_g


def _config(iterations: int) -> RuntimeConfig:
    config = RuntimeConfig()
    config.hadronic.include_pair_production = True
    config.hadronic.pair_cascade_iterations = iterations
    return config


def test_pair_branch_single_step_tau() -> None:
    dynamics, electron, seed_field_hz, frequency_hz, magnetic_field_g = _state_inputs()

    pair_lum, pair_seed, tau_pair, _pair_density = _compute_pair_production_branch(
        dynamics, electron, seed_field_hz, frequency_hz, magnetic_field_g, _config(1),
    )

    assert np.all(np.isfinite(pair_lum))
    assert np.all(np.isfinite(pair_seed))
    assert np.all(tau_pair >= 0.0)
    assert float(np.max(tau_pair)) > 0.0


def test_pair_branch_iterative_tau_is_not_placeholder() -> None:
    dynamics, electron, seed_field_hz, frequency_hz, magnetic_field_g = _state_inputs()

    pair_lum, pair_seed, tau_pair, _pair_density = _compute_pair_production_branch(
        dynamics, electron, seed_field_hz, frequency_hz, magnetic_field_g, _config(2),
    )

    assert np.all(np.isfinite(pair_lum))
    assert np.all(np.isfinite(pair_seed))
    assert np.all(tau_pair >= 0.0)
    assert float(np.max(pair_lum)) > 0.0
    assert float(np.max(tau_pair)) > 0.0


def test_public_pair_cascade_config_mapping() -> None:
    model = top_hat_model(
        fwd_rad=radiation(proton_energy_fraction=0.1, pair_production=True),
        hadronic=hadronic(enabled=True, pair_cascade_iterations=3),
    )
    config = _build_fit_config_for_patch(
        model,
        theta_v=0.0,
        opening_angle_jet=0.1,
        e_iso=1.0e52,
        gamma0=300.0,
    )
    assert config.hadronic.include_pair_production
    assert config.hadronic.pair_cascade_iterations == 3


def main() -> None:
    test_pair_branch_single_step_tau()
    test_pair_branch_iterative_tau_is_not_placeholder()
    test_public_pair_cascade_config_mapping()


if __name__ == "__main__":
    main()
