# ASGARD Source Tree

精简源码树索引。代码总览见 `doc/code_overview.md`。

## Python Runtime

- `ASGARD/` — 外部用户 API: `api_model.py`, `api_observe.py`, `api_fit.py`, `api_adaptive.py`
- `asgard_core/` — 内部编排与物理耦合: `asgard_runtime.py`, `asgard_state.py`, `asgard_setup.py`, `asgard_ssc.py`, `asgard_types.py`
- `asgard_core/hadronic_*.py` — hadronic Python wrappers/benchmarks (orchestration only; final microphysics → `src/Hadronic/`)
- `asgard_core/hadronic_reverse.py` — RS hadronic proton transport + synchrotron wrapper
- `asgard_core/hadronic_cascade.py` — iterative pair-production synch branch wrapper, not full time-dependent cascade PDE

## Fortran Runtime

- `src/Constants.f90`
- `src/Dynamics/`: `dynamics_common.f90`, `Dynamics_forward.f90`, `Dynamics_reverse.f90`
- `src/Electron/`:
  - 1D/2D entries: `FS_electron_fullhide_1d.f90`, `FS_electron_fullhide_2d.f90`, `FS_electron_charint_1d.f90`, `FS_electron_charint_2d.f90`, `FS_electron_slc1_1d.f90`, `FS_electron_t2g1_1d.f90`, `FS_electron_weno5_1d.f90`
  - Kernels: `electron_common.f90`, `electron_radiation_kernel.f90`, `electron_cooling_kernel.f90`, `electron_seed_history_kernel.f90`, `electron_transport_2d_kernel.f90`, `electron_injection_profiles.f90`, `electron_transport_common.f90`, `electron_reverse_kernel.f90`, `adaptive_resampling_mod.f90`
- `src/Radiation/`: `radiation_common.f90`, `SSC_spec.f90`, `Annihilation.f90`, `Seed_reverse.f90`, `synchrotron_polarization_kernel.f90`, `quantum_synchrotron_kernel.f90`
- `src/Hadronic/`:
  - Entries: `FS_hadronic_1d.f90`, `FS_hadronic_reverse_1d.f90`
  - Kernels: `hadronic_common.f90`, `hadronic_transport_kernel.f90`, `hadronic_radiation_kernel.f90`, `hadronic_interaction_kernel.f90`, `hadronic_decay_kernel.f90`, `hadronic_pair_production_kernel.f90`, `hadronic_pair_cascade_kernel.f90`, `hadronic_pp_kernel.f90`, `hadronic_pp_models_kernel.f90`, `hadronic_bethe_heitler_kernel.f90`, `hadronic_hadronic_ic_kernel.f90`, `hadronic_species_transport_kernel.f90`, `hadronic_acceleration_kernel.f90`, `hadronic_secondary_radiation_kernel.f90`
- `src/Interpolation/`: `SED_interpolation.f90`, `SED_interpolation_structured.f90`, `interpolation_common.f90`

## Build Entrypoint

- `build_extensions.py` — f2py 编译入口。登记的 active hadronic extensions: `FS_hadronic_1d`, `FS_hadronic_reverse_1d`

## Tests / Benchmarks

- `tests/readme_smoke_bench.py`, `tests/reverse_shock_smoke.py`
- `tests/fullhide_2d_smoke_bench.py`, `tests/fullhide_2d_medium_diag.py`
- `tests/vegas_afterglow_comparison.py`, `tests/sed_electron_compare.py`
- Hadronic: `tests/hadronic_1d_smoke.py`, `tests/hadronic_species_transport_smoke.py`, `tests/hadronic_secondary_radiation_smoke.py`, `tests/hadronic_acceleration_smoke.py`, `tests/hadronic_bethe_heitler_smoke.py`, `tests/hadronic_hadronic_ic_smoke.py`, `tests/hadronic_pair_production_smoke.py`, `tests/hadronic_pp_smoke.py`, `tests/hadronic_pg_neutrino_1d_diag.py`, `tests/hadronic_proton_synch_1d_diag.py`, `tests/hadronic_pgamma_benchmark_report.py`
- Pair/polarization/RS: `tests/hadronic_pair_cascade_smoke.py`, `tests/hadronic_pair_branch_smoke.py`, `tests/hadronic_reverse_shock_smoke.py`, `tests/polarization_smoke.py`, `tests/polarization_baseline_bench.py`
- AM3/reference comparisons: `tests/hadronic_am3_solver_smoke.py`, `tests/hadronic_am3_reference_compare_smoke.py`, `tests/hadronic_am3_acceleration_compare_smoke.py`, `tests/hadronic_am3_bethe_heitler_compare_smoke.py`, `tests/hadronic_pgamma_benchmark_compare.py`
