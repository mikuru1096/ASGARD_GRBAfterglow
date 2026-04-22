# ASGARD Agent Notes

## Working Rules

- Reply to the user in Chinese.
- Prefix shell commands with `rtk`.
- All builds, compile checks, and runtime commands must run inside WSL Ubuntu.
- Default environment is `WSL + uv`.
- `rtk` is available inside WSL and must be used from WSL only.
- Do not use `wsl.exe`, PowerShell, Windows CMD, native Windows Python, or any other Windows-side command for this repository unless the user explicitly reverts this policy.
- Do not add numerical guards to non-simulation code.
- Do not hide discontinuities or non-smooth physical results with plotting or post-processing changes.
- If a physical time-evolution quantity is not continuous and smooth, treat it as a likely numerical or physics bug until checked.
- Avoid adding diagnostic scripts unless they will remain useful as regression tests.
- After important Fortran changes, run:
  - the affected `build_extensions.py` target
  - a `-Wline-truncation` check
  - the smallest relevant smoke or diagnostic test

## Build Commands

WSL runtime extension build:

```bash
rtk uv run python build_extensions.py --module FS_electron_fullhide_2d --force
```

WSL runtime extension build (`charint_2d`):

```bash
rtk uv run python build_extensions.py --module FS_electron_charint_2d --force
```

WSL toolchain check:

```bash
rtk /usr/bin/gfortran --version
```

WSL import smoke:

```bash
rtk uv run python - <<'PY'
import ASGARD
print("import-ok")
PY
```

## Public Runtime Status

- Public API:
  - `observe(model, config=...)`
  - `run_fit(config)`
- Main electron solvers:
  - `fullhide_1d`
  - `slc1_1d`
  - `charint_1d`
  - `charint_2d`
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
- `charint_2d` uses:
  - shared `fullhide_2d` outer physics, history, and shell diagnostics
  - characteristic eta/log-chi advection
  - implicit eta/log-chi diffusion
  - characteristic xi/log-gamma advance

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
