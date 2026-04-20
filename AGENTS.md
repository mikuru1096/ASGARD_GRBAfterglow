# ASGARD Agent Notes

## Working Rules

- Reply to the user in Chinese.
- Prefix shell commands with `rtk`.
- Default direct Fortran check environment is WSL Ubuntu:
  - `wsl.exe -d Ubuntu -- bash -lc '...'`
- Use native Windows Python only for runtime extension builds that must produce Windows `.pyd` artifacts.
- Do not add numerical guards to non-simulation code.
- Do not hide discontinuities or non-smooth physical results with plotting or post-processing changes.
- If a physical time-evolution quantity is not continuous and smooth, treat it as a likely numerical or physics bug until checked.
- Avoid adding diagnostic scripts unless they will remain useful as regression tests.
- After important Fortran changes, run:
  - the affected `build_extensions.py` target
  - a `-Wline-truncation` check
  - the smallest relevant smoke or diagnostic test

## Build Commands

Windows runtime extension build:

```powershell
rtk powershell -Command '& "C:\Users\jia\AppData\Local\Programs\Python\Python312\python.exe" build_extensions.py --module FS_electron_fullhide_2d --force'
```

WSL toolchain check:

```powershell
rtk powershell -Command "wsl.exe -d Ubuntu -- bash -lc 'cd /mnt/c/Users/jia/Documents/New\ project/ASGARD_GRBAfterglow && /usr/bin/gfortran --version'"
```

## Public Runtime Status

- Public API:
  - `observe(model, config=...)`
  - `run_fit(config)`
- Main electron solvers:
  - `fullhide_1d`
  - `slc1_1d`
  - `charint_1d`
  - `t2g1_1d`
  - `weno5_1d`
  - `fullhide_2d`
- Runtime aliases remain:
  - `fullhide -> fullhide_1d`
  - `slc1 -> slc1_1d`
  - `charint -> charint_1d`
  - `t2g1 -> t2g1_1d`
  - `weno5 -> weno5_1d`
- `num_chi` semantics:
  - `*_1d -> 1`
  - `*_2d -> explicit value or default 64`

## Electron Kernel Layout

- `src/Electron/electron_common.f90`
- `src/Electron/electron_radiation_kernel.f90`
- `src/Electron/electron_cooling_kernel.f90`
- `src/Electron/electron_seed_history_kernel.f90`
- `src/Electron/electron_transport_2d_kernel.f90`
- `src/Electron/calling_modules.f90` remains the aggregate export layer.

## Current 2D State

- `fullhide_2d` uses:
  - chi-dependent downstream velocity
  - implicit eta/log-chi transport
  - shared density and injection helpers
  - chi-resolved historical photon fields
  - shell-reduced public `P_syn / Seed_syn`
  - shell-level diagnostic `nu_m / nu_c / nu_a`
- High-energy xi/log-gamma boundary now uses the verified 1d implicit coefficient path.
- Internal cooling currently uses a reduced 4-band cooling grid.
- Shell-level cooling state reuse is enabled:
  - shell cooling assembly once per shell
  - substep cooling assembly disabled
  - substeps still update dynamics interpolation, density, injection, eta transport, and xi advance

## Current Test Entrypoints

- `tests/readme_smoke_bench.py`
- `tests/fullhide_2d_smoke_bench.py`
- `tests/fullhide_2d_medium_diag.py`
- `tests/vegas_afterglow_comparison.py`

## Cleanup Policy

- Do not commit `.buildcache/`.
- Do not commit `output/benchmark_exp_tail/`.
- Keep generated documentation plots in `output/asgard_doc/vegas_afterglow_compare/` only when they are current benchmark artifacts.
- Remove one-off debug scripts after their findings have been folded into code, tests, or this note.
