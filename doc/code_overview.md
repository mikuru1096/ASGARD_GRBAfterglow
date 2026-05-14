# ASGARD Code Overview

本文档按当前工作树整理代码结构、运行主链和关键边界。算法细节见 `doc/electron_solver_algorithms.md`。

## 1. Public API

- `ASGARD/api_model.py`: `Model`, `Medium`, `ISM`, `Wind`, `TophatJet`, `GaussianJet`, `PowerLawJet`, `TwoComponentJet`, `StepPowerLawJet`, `Ejecta`, `Observer`, `Radiation`, `Setups`
- `Model.flux_density_grid(times_s, nu_hz)`, `flux_density(times_s, nu_hz)`, `spectrum(time_s, nu_hz)`, `flux(time_s, nu_min, nu_max)`, `sky_image(t_obs, nu_obs, fov)`, `details()`
  - `Model.polarization(times_s, nu_hz, magnetic_geometry=..., local_emissivity=...)`
- Hadronic public switches: `Radiation.pair_production`, `Radiation.pg`, `Radiation.bethe_heitler`, `Radiation.pp`, `Radiation.neutrino`, `Radiation.reverse_epsilon_p`; cascade substeps use `Setups.pair_cascade_iterations`
- Reverse-shock magnetization switch: `ReverseShockConfig.sigma` / `Setups.reverse_sigma`
- `ASGARD/api_observe.py`: `observe(model, config)`, `run_fit(config)`
- `ASGARD/api_fit.py`: `Fitter`, `Param`, `FitResult`
- Electron solver names: `fullhide_1d`, `slc1_1d`, `charint_1d`, `charint_2d`, `t2g1_1d`, `weno5_1d`, `fullhide_2d`
- Aliases: `fullhide`, `slc1`, `charint`, `t2g1`, `weno5` → `*_1d`

## 2. Runtime Main Chain

```text
Model / observe / run_fit
  → FitConfig → SimulationSetup
  → solve_state_from_setup
  → solve_dynamics → solve_electron → photon_field_stage
  → solve_hadronic → solve_reverse_shock_emission
  → observer assembly → Radiation.annihilation
  → Interpolation.sed_interpolation → API result
```

Core state objects (`asgard_core/asgard_types.py`):
- `DynamicsSolution`: `r_tobs`, `r_gamma`, `radius`, `swept_mass_g`
- `ReverseShockDynamics`: `M3`, total `B3`, ordered crossing field `B3_ordered_cross`, `U3/V3`, `gamma34` and crossing thermal records; `B3` is the sum of turbulent `sqrt(8 pi epsilon_B,r U3/V3)` and the optional ordered upstream field amplified by the magnetized compression ratio
- `ElectronSolution`: `gam_e`, `d_n_gam_e`, `l_syn_spec`, `seed_syn`, `nu_m`, `nu_c`, `nu_a`; 2D adds `d_n_gam_e_chi`, `chi_grid`; BH adds `d_n_gam_e_bh`
- `PhotonFieldState`: forward synch seed, hadronic target field, absorption seed field
- `HadronicSolution`: 1D hadronic proton/secondary/radiation results
- `ObserverState`: absorption factors, `tau_pair`, flux components
- `FluxComponents`: `total`, FS synch/SSC, hadronic, RS synch/SSC, cross-zone IC

State machine (`asgard_core/asgard_state.py`):
- `solve_state_from_setup`: dynamics → electron → photon field → hadronic → reverse shock → observer
- `_build_photon_field_stage`: copy electron `seed_syn`; hadronic SSC seed → target field (when `epsilon_p > 0`)
- `_solve_hadronic_stage`: call `solve_hadronic`; BH secondary e± → merge into forward electron → recompute `l_syn_spec/seed_syn`; pγ photon survival applied to photon field
- `_assemble_observer_stage`: assemble FS synch/SSC, RS synch/SSC, cross-zone IC, hadronic components; hadronic photons use electron Fortran kernel's SSA transfer

Fit shortest path:
```text
Fitter.loglike → compile_problem → eval_loglike → solve_state_from_setup
  → project_flux_grid → combine_multiband_flux → compute_light_curve_redchi
```

## 3. Python Orchestration (`asgard_core/`)

- `asgard_setup.py`: `FitConfig → SimulationSetup`
- `asgard_runtime.py`: backend selection, Fortran extension dispatch, array wrapping
- `asgard_state.py`: main state machine, cross-stage coupling
- `asgard_ssc.py`: forward SSC auxiliary grid + seed
- `asgard_coupling.py`: FS/RS cross-zone IC geometry + seed field coupling
- `asgard_postprocess.py`: observer projection, band aggregation, fit postprocessing
- `asgard_fit.py`: fit problem compilation, likelihood path
- `asgard_types.py`: runtime dataclass contracts
- Polarization timing diagnosis against Lan 2023 is recorded in `doc/polarization_timing_diagnostic.md`; current evidence points to dynamics / jet-evolution benchmarking before any projection-layer change.

Hadronic Python modules (orchestration/wrapping/benchmark only):
- Fortran wrappers: `hadronic_hummer.py`, `hadronic_bethe_heitler.py`, `hadronic_hadronic_ic.py`, `hadronic_pp.py`, `hadronic_pair_production.py`, `hadronic_species_transport.py`, `hadronic_secondary_radiation.py`, `hadronic_acceleration.py`
- Reverse shock wrapper: `hadronic_reverse.py`; full RS hadronic chain is dispatched from `asgard_runtime.py` through the formal 1D hadronic kernels when RS full-chain flags are enabled. RS seed photons, RS `B3`, shell energy and baryon target density come from the ASGARD RS dynamics/thermal state, not from a Vegas-style averaged thermal closure.
- Reference/benchmark: `hadronic_pgamma.py`, `hadronic_am3_solver.py`, `hadronic_am3_benchmark.py`, `hadronic_cascade.py`; IC-mediated electromagnetic cascade remains outside the current pair-cascade contract, see `doc/pair_cascade_extension_boundary.md`

Final AM3-derived microphysics lives in `src/Hadronic/*.f90`.

## 4. Fortran Kernels

### Dynamics
- `src/Dynamics/Dynamics_forward.f90`: forward shock dynamics, ISM/wind, density jumps, energy injection
- `src/Dynamics/Dynamics_reverse.f90`: reverse shock dynamics with explicit region-3 `U3/V3/gamma34`; injection uses shock-front `gamma34`, turbulent field and post-crossing thermal evolution use `U3/V3`, and optional upstream `sigma_r` adds MHD-jump compression plus an ordered magnetic component; f2py returns `B3_ordered_cross` after `gamma_m_cross`
- `src/Dynamics/dynamics_common.f90`: shared dynamics auxiliaries

### Electron
- Main entries: `FS_electron_fullhide_1d.f90`, `FS_electron_fullhide_2d.f90`, `FS_electron_charint_1d.f90`, `FS_electron_charint_2d.f90`, `FS_electron_slc1_1d.f90`, `FS_electron_t2g1_1d.f90`, `FS_electron_weno5_1d.f90`
- Shared kernels: `electron_common.f90` (log-γ grid, injection, spectral init, characteristic helpers), `electron_cooling_kernel.f90` (synch/IC/SSA cooling assembly), `electron_radiation_kernel.f90` (synchrotron spectrum and seed), `electron_seed_history_kernel.f90` (historical photon fields), `electron_transport_2d_kernel.f90` (log-χ geometry, η/log-χ transport, 2D energy advance), `electron_injection_profiles.f90` (injection profiles), `electron_transport_common.f90` (implicit coefficients, PPM, characteristic transport, semi-Lagrangian), `electron_reverse_kernel.f90` (reverse shock electron evolution), `adaptive_resampling_mod.f90` (log-space adaptive resampling)

### Radiation / Interpolation
- `src/Radiation/SSC_spec.f90`: SSC spectrum and seed (uniform + nonuniform)
- `src/Radiation/Seed_reverse.f90`: reverse shock synchrotron radiation + seed
- `src/Radiation/Annihilation.f90`: γ-γ absorption
- `src/Radiation/radiation_common.f90`: Simpson weights, power-law interp, pair cross-section, synchrotron seed core, transfer factor
- `src/Radiation/synchrotron_polarization_kernel.f90`: frequency-dependent synchrotron polarization emissivity
- `src/Radiation/quantum_synchrotron_kernel.f90`: quantum synchrotron helper kernel
- `src/Interpolation/SED_interpolation.f90`: observer-frame EATS/Doppler interpolation
- `src/Interpolation/SED_interpolation_structured.f90`: structured jet interpolation
- `src/Interpolation/interpolation_common.f90`: log-SED accumulation

### Hadronic
`src/Hadronic/FS_hadronic_1d.f90` is the forward-shock f2py entry point, dispatching to:
- `hadronic_transport_kernel.f90`: proton injection, adiabatic/synchrotron loss, log-γ energy advance
- `hadronic_radiation_kernel.f90`: proton synchrotron
- `hadronic_interaction_kernel.f90`: Hummer 2010 photopion operator
- `hadronic_decay_kernel.f90`: π⁰→γ, π±/μ decay, neutrino emissivity
- `hadronic_pair_production_kernel.f90`: γ-γ pair production (Aharonian+1983)
- `hadronic_pair_cascade_kernel.f90`: single-step pair-cascade synchrotron kernel used by the iterative branch
- `hadronic_pp_kernel.f90`: pp δ-channel (Kelner+2006)
- `hadronic_pp_models_kernel.f90`: alternate pp source model helpers
- `hadronic_bethe_heitler_kernel.f90`: Bethe-Heitler pair source + proton loss
- `hadronic_hadronic_ic_kernel.f90`: hadronic inverse Compton
- `hadronic_species_transport_kernel.f90`: neutron/π±/μ± explicit transport (upwind + Strang splitting)
- `hadronic_acceleration_kernel.f90`: species-aware acceleration timescale, injection operator, γ_max estimate
- `hadronic_secondary_radiation_kernel.f90`: pion/muon synchrotron + IC
- `hadronic_common.f90`: shared hadronic constants, grid builders, validation

Reverse-shock hadronic light entry: `src/Hadronic/FS_hadronic_reverse_1d.f90`, proton injection/transport + proton synchrotron only. Full-chain RS hadronic dispatch reuses the formal 1D hadronic kernels from `FS_hadronic_1d` through the Python runtime wrapper, with RS magnetic field, RS seed photons, RS shell energy and RS baryon target density. This is still a formal 1D dispatch path, not 2D / χ-resolved RS hadronic transport.

2D / χ-resolved hadronic transport is intentionally not implemented. The current `chi_grid` belongs to the 2D electron transport contract, while `PhotonFieldState` and `HadronicSolution` are shell-level contracts. The implementation boundary and required physics contract are recorded in `doc/hadronic_chi_transport_decision.md`.

Shell-level entries exposed by `FS_hadronic_1d.f90`: `fs_hadronic_1d`, `fs_hadronic_proton_syn_shell`, `fs_hadronic_syn_polarization_shell`, `fs_hadronic_pgamma_operator_shell`, `fs_hadronic_pair_production_shell`, `fs_hadronic_pp_delta_shell`, `fs_hadronic_bethe_heitler_shell`, `fs_hadronic_hadronic_ic_shell`, `fs_hadronic_species_transport_shell`, `fs_hadronic_acceleration_shell`, `fs_hadronic_secondary_radiation_shell`, `fs_hadronic_decay_operator_shell`, `fs_hadronic_pair_cascade_step`, `fs_hadronic_pp_spectral_source`, `fs_hadronic_quantum_syn_cooling_factor`

## 5. Hadronic Current State

- **Configuration**: `Radiation.epsilon_p`, `.proton_synch`, `.pg`, `.bethe_heitler`, `.hadronic_inverse_compton`, `.pp`, `.neutrino`, `.eta_acc`, `.pgamma_scheme`; `Setups.hadronic_enabled`, `.hadronic_solver`, `.num_gam_p`, `.num_nu_nu`
- **Solver names**: `legacy_1d` (proton transport + proton synchrotron only), `am3_1d` (current formal hadronic main path, forward shock + 1D electron only)
- **`pgamma_scheme`**: `hummer_2010_response` (with transport feedback), `ka2008_reference` (emission-only, no transport feedback), `disabled`
- **Active couplings**: proton injection/cooling, proton synchrotron, photopion (α_p, Q_p^reinj, α_γ^{pγ} as local shell survival), neutrino, Bethe-Heitler (proton cooling + e± feedback), pp (γ/ν/pair/proton-loss), hadronic IC (proton + pion/muon), explicit secondary transport (n/π±/μ±), secondary radiation (pion/muon synch + IC), pair production (observer-side attenuation + pair synchrotron branch)
- **Reverse-shock hadronic**: `FS_hadronic_reverse_1d` remains the light proton-synch path; full-chain RS pγ/BH/pp/secondary/cascade coupling is active when `reverse_epsilon_p > 0` and the corresponding hadronic process flags are enabled.
- **Pair cascade**: `pair_cascade_iterations > 1` now selects a shell-sequence time-dependent γγ pair/synch cascade path; the old single-shell iterative kernel remains available in `hadronic_cascade.py` for diagnostics. IC-mediated electromagnetic cascade is intentionally not part of this branch until a photon/e± source-sink and energy-budget contract is defined in `doc/pair_cascade_extension_boundary.md`.

## 6. Build and Test

Default: WSL Ubuntu + uv, commands prefixed with `rtk`.

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module FS_electron_fullhide_2d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module FS_electron_charint_2d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module FS_hadronic_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module FS_hadronic_reverse_1d --force'
rtk /usr/bin/gfortran --version
```

Smoke tests: `rtk uv run python tests/readme_smoke_bench.py`, `tests/fullhide_2d_smoke_bench.py`, `tests/fullhide_2d_medium_diag.py`

Hadronic regressions: `tests/hadronic_1d_smoke.py`, `tests/hadronic_species_transport_smoke.py`, `tests/hadronic_secondary_radiation_smoke.py`, `tests/hadronic_acceleration_smoke.py`, `tests/hadronic_bethe_heitler_smoke.py`, `tests/hadronic_hadronic_ic_smoke.py`, `tests/hadronic_pair_production_smoke.py`, `tests/hadronic_pp_smoke.py`, `tests/hadronic_pg_neutrino_1d_diag.py`, `tests/hadronic_proton_synch_1d_diag.py`, `tests/hadronic_pgamma_benchmark_report.py`

Benchmark: `tests/vegas_afterglow_comparison.py`, `tests/sed_electron_compare.py`. RS benchmark refresh:
`rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --extra compare python tests/vegas_afterglow_comparison.py --scenario baseline --only reverse_shock_lc reverse_shock_thermal'`

Benchmark refresh protocol: `doc/benchmark_refresh_protocol.md`

## 7. Known Boundaries

Not yet complete:
- 2D / χ-resolved hadronic transport; current decision is to keep the formal 1D path until χ-local photon, hadron, secondary-feedback and observer-projection contracts are defined
- inverse-Compton-mediated electromagnetic pair cascade beyond the current γγ pair/synch contract
- reverse-shock `nu_m` is a diagnostic break derived from local `gamma34` injection and post-crossing `U3/M3`; it does not replace the transported electron spectrum
- jet spreading in the current backend
- user-defined `Medium` kernel dispatch and wind profiles with `k != 2`
- thermal electron branch outside `fullhide_1d`

Architectural boundaries:
- ASGARD = shell-evolving blast-wave + observer projection
- AM3 = microphysics/kernel reference only; must not replace ASGARD's dynamics/electron/observer chain
- Final AM3-derived microphysics → `src/Hadronic/*.f90`; Python → orchestration only
- Non-smooth physical time-evolution → treat as bug first
- Public/backend unsupported boundaries are fixed in `doc/public_backend_limits.md`
