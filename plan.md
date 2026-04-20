# ASGARD Current Plan

## Current Status

- Electron solver names are normalized to:
  - `fullhide_1d`
  - `slc1_1d`
  - `charint_1d`
  - `t2g1_1d`
  - `weno5_1d`
  - `fullhide_2d`
- Legacy runtime aliases remain available:
  - `fullhide -> fullhide_1d`
  - `slc1 -> slc1_1d`
  - `charint -> charint_1d`
  - `t2g1 -> t2g1_1d`
  - `weno5 -> weno5_1d`
- `fullhide_2d` is connected to the Python runtime and exposes:
  - `chi_grid`
  - `d_n_gam_e_chi`
  - shell-reduced `P_syn / Seed_syn`
  - shell-level `nu_m / nu_c / nu_a`
- `details()` returns electron spectra through `gamma_e` and shell-total `dN_dgamma_e`.

## Completed 2D Fixes

- `fullhide_2d` uses chi-dependent downstream transport speed.
- `fullhide_2d` uses implicit `eta/log chi` transport.
- `fullhide_2d` reuses shared density and injection helpers.
- The early shell-total electron deficit was traced to real numerical bugs and fixed.
- The high-energy boundary in `advance_energy_loggamma_chi(...)` now reuses the verified 1d implicit coefficient build and backward sweep.
- The internal cooling field currently uses a reduced 4-band grid; public full-frequency radiation outputs are still produced on the public frequency grid.
- Shell-level cooling state reuse is enabled:
  - expensive cooling assembly is done once per shell
  - shell-internal substeps reuse `dEL_mean_chi`, `kappa2_chi`, and reduced-band cooling auxiliaries
  - substeps still update dynamics interpolation, density, injection, eta transport, and xi advance

## Current Performance Snapshot

- `tests/fullhide_2d_smoke_bench.py` currently passes.
- Recent non-profile smoke timings:
  - `basic_smoke`: about `0.35s`
  - `electron_grid`: about `0.08s`
- `ASGARD_PROFILE_2D=1` shows:
  - `shell_cooling_calls = 7`
  - `substep_cooling_calls = 0`
- This meets the requested 10% runtime target relative to the original smoke baseline.

## Current Diagnostic And Benchmark Entrypoints

- Minimal 2d smoke:
  - `tests/fullhide_2d_smoke_bench.py`
- Medium 2d diagnostic:
  - `tests/fullhide_2d_medium_diag.py`
- Public benchmark/plot comparison:
  - `tests/vegas_afterglow_comparison.py`
- README/public smoke:
  - `tests/readme_smoke_bench.py`

## Build And Check Policy

- Windows runtime extension build remains:

```powershell
rtk powershell -Command '& "C:\Users\jia\AppData\Local\Programs\Python\Python312\python.exe" build_extensions.py --module FS_electron_fullhide_2d --force'
```

- WSL Ubuntu is the default direct Fortran check environment:

```powershell
rtk powershell -Command "wsl.exe -d Ubuntu -- bash -lc 'cd /mnt/c/Users/jia/Documents/New\ project/ASGARD_GRBAfterglow && /usr/bin/gfortran --version'"
```

- Fortran changes require:
  - affected extension build check
  - `-Wline-truncation` check
  - relevant smoke or diagnostic run

## Current Open Items

- Magnetic-decay 2d plots are exploratory benchmark artifacts, not regression criteria.
- If high-resolution physics checks regress, compare shell-level cooling refresh against substep-level cooling refresh on the same grid before changing SSA formulas.
- Avoid adding new diagnostic scripts unless they are expected to remain as regression tests.
