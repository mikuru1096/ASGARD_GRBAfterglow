# Numerical Methods

本文档按实现模块整理当前数值方法和验证重点。

## Dynamics

Forward dynamics:

- Source: `src/Dynamics/Dynamics_forward.f90`
- Shared helpers: `src/Dynamics/dynamics_common.f90`
- Supports ISM and wind `k=2` backend, density jump and energy injection paths.

Reverse dynamics:

- Source: `src/Dynamics/Dynamics_reverse.f90`
- Returns RS shell state, `gamma34`, `U3/V3`, turbulent and ordered magnetic field diagnostics.
- Magnetized jump follows the current RS baseline contract; `sigma -> 0` is a required regression check.

## Electron Energy Grid

Electron solvers use log-gamma grids. Common helpers live in:

- `src/Electron/electron_common.f90`
- `src/Electron/electron_injection_profiles.f90`
- `src/Electron/adaptive_resampling_mod.f90`

Key numerical concerns:

- injection normalization
- high-energy cutoff
- cooling face speeds
- conservative or implicit transport update
- stable synchrotron/SSA seed recomputation

## 1D Electron Solvers

### `fullhide_1d`

Default public baseline:

- Source: `src/Electron/FS_electron_fullhide_1d.f90`
- Transport helper: `src/Electron/electron_transport_common.f90`
- Cooling: `src/Electron/electron_cooling_kernel.f90`
- Radiation: `src/Electron/electron_radiation_kernel.f90`

Strength:

- stable for stiff cooling and dense environments
- fast enough for fitting
- current baseline for public comparisons

Validation:

- compile `FS_electron_fullhide_1d`
- `-Wline-truncation` on source closure
- `tests/polarization_smoke.py` or relevant electron smoke
- benchmark path when performance changes

### `weno5_1d`

High-order spectrum-resolving solver:

- Source: `src/Electron/FS_electron_weno5_1d.f90`
- Useful for resolving spectral shape, not necessarily default fitting path.

### `charint_1d`, `slc1_1d`, `t2g1_1d`

Alternative transport families retained for comparison, diagnostics and algorithmic experiments. They must remain buildable if touched.

## 2D Electron Solvers

2D path includes energy and chi-resolved electron transport:

- `src/Electron/FS_electron_fullhide_2d.f90`
- `src/Electron/FS_electron_charint_2d.f90`
- `src/Electron/electron_transport_2d_kernel.f90`
- `src/Electron/electron_seed_history_kernel.f90`

2D electron transport does not imply 2D hadronic transport. Hadronic state remains shell-level until a separate contract exists.

Validation:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module FS_electron_fullhide_2d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/fullhide_2d_smoke_bench.py'
```

## Cooling Assembly

Cooling kernel:

- `src/Electron/electron_cooling_kernel.f90`

Responsibilities:

- synchrotron cooling
- IC / Compton auxiliary terms
- SSA cooling
- Nakar/Fan-style Y paths
- batch cooling for 2D chi cells

The cooling kernel separates photon-field preparation from final cooling assembly. This avoids recomputing IC auxiliary terms when SSA seed and IC seed are shared.

## Radiation Integrals

Synchrotron:

- `src/Electron/electron_radiation_kernel.f90`
- `src/Radiation/radiation_common.f90`

SSC:

- `src/Radiation/SSC_spec.f90`

Gamma-gamma:

- `src/Radiation/Annihilation.f90`

Reverse seed:

- `src/Radiation/Seed_reverse.f90`

Numerical notes:

- Fixed-grid synchrotron path is the public fast path.
- Adaptive synchrotron integration is diagnostic.
- SSA transfer must remain continuous across frequency cells.
- Do not add floor/guard terms unless they represent physical boundaries already present in the model.

## Hadronic Kernels

Fortran sources:

- `src/Hadronic/FS_hadronic_1d.f90`
- `src/Hadronic/FS_hadronic_reverse_1d.f90`
- `src/Hadronic/hadronic_transport_kernel.f90`
- `src/Hadronic/hadronic_radiation_kernel.f90`
- `src/Hadronic/hadronic_interaction_kernel.f90`
- `src/Hadronic/hadronic_pair_production_kernel.f90`
- `src/Hadronic/hadronic_pp_kernel.f90`
- `src/Hadronic/hadronic_bethe_heitler_kernel.f90`
- `src/Hadronic/hadronic_hadronic_ic_kernel.f90`
- `src/Hadronic/hadronic_species_transport_kernel.f90`
- `src/Hadronic/hadronic_acceleration_kernel.f90`
- `src/Hadronic/hadronic_secondary_radiation_kernel.f90`
- `src/Hadronic/hadronic_decay_kernel.f90`
- `src/Hadronic/hadronic_pair_cascade_kernel.f90`

Validation should be process-specific:

- proton synch
- p-gamma
- Bethe-Heitler
- pp
- hadronic IC
- species transport
- secondary radiation
- pair production/cascade

## Observer Interpolation

Interpolation sources:

- `src/Interpolation/SED_interpolation.f90`
- `src/Interpolation/SED_interpolation_structured.f90`
- `src/Interpolation/interpolation_common.f90`

Projection integrates local shell radiation over observer time/frequency grids. Projection-layer fixes must not be used to hide dynamics or transport bugs.

## Line-Truncation Check

For touched Fortran source closures, run from `/tmp` with a temporary module directory. Example for the `FS_electron_fullhide_1d` closure:

```bash
rtk bash -lc 'source ~/.wsl_env && rm -rf /tmp/asgard_linecheck && mkdir -p /tmp/asgard_linecheck && cd /tmp && gfortran -cpp -fopenmp -Wall -Wline-truncation -fsyntax-only -J /tmp/asgard_linecheck -I /tmp/asgard_linecheck "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Constants.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Dynamics/dynamics_common.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Radiation/radiation_common.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Radiation/synchrotron_polarization_kernel.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/electron_transport_common.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/adaptive_resampling_mod.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/electron_radiation_kernel.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/electron_injection_profiles.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/electron_common.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/electron_cooling_kernel.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/electron_seed_history_kernel.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/FS_electron_fullhide_1d.f90"'
```

Do not run this from the repository root, where stale `.mod` files can be picked up.

## Performance Notes

- `fullhide_1d` is the default fitting path.
- Cold solve time and cached query time must be reported separately.
- Small grids can expose serial hot spots that are hidden in larger benchmark runs.
- Benchmark scratch files belong under `/tmp`; only reproducible documentation artifacts should be committed.
