from __future__ import annotations

from typing import Optional

import numpy as np

from asgard_component_backend import extract_physical_solution_from_state, solve_model_state_from_setup
from asgard_observables import OUTPUT_BANDS
from asgard_postprocess import compute_band_fluxes, compute_light_curve_redchi, compute_spectrum_flux
from asgard_models import FitConfig, FitResult, PhysicalSolution, SimulationSetup
from asgard_setup import build_simulation_setup


def _build_fit_result(
    setup: SimulationSetup,
    physical: PhysicalSolution,
    config: FitConfig,
) -> FitResult:
    bands_flux = compute_band_fluxes(setup, physical, config)
    redchi = compute_light_curve_redchi(bands_flux, setup.observer_time_s, config)
    spectrum_freq_hz, spectrum_fnu = compute_spectrum_flux(setup, physical, config)
    return FitResult(
        t_obs_s=setup.observer_time_s,
        characteristic_time_s=physical.characteristic_time_s,
        bands=OUTPUT_BANDS,
        bands_flux=bands_flux,
        redchi=redchi,
        nu_m=physical.nu_m,
        nu_c=physical.nu_c,
        nu_a=physical.nu_a,
        rs_nu_m=physical.rs_nu_m,
        rs_nu_c=physical.rs_nu_c,
        rs_nu_a=physical.rs_nu_a,
        spectrum_freq_hz=spectrum_freq_hz,
        spectrum_fnu=spectrum_fnu,
    )


def _solve_physical_pipeline(
    setup: SimulationSetup,
    config: FitConfig,
) -> PhysicalSolution:
    state = solve_model_state_from_setup(config, setup)
    return extract_physical_solution_from_state(state)


def run_fit(config: FitConfig) -> FitResult:
    setup = build_simulation_setup(config)
    physical = _solve_physical_pipeline(setup, config)
    result = _build_fit_result(setup, physical, config)

    if config.plot_lc:
        from asgard_plot import plot_light_curve

        plot_light_curve(result, show=config.show_plots)

    return result
