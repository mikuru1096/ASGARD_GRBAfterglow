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

Current optimization pass keeps the WSL + uv build surface below fixed; do not switch toolchains mid-pass.
Kernel2 optimization continues under the same fixed WSL + uv build surface and must remain behavior-preserving.
Electron-kernel optimization pass covers `src/Electron/*kernel*.f90` plus shared electron core files, with no public API or physics changes.
Current electron pass also collapses duplicated seed-history mapping and reuses the shared 2D transport core for `charint_2d` while preserving the existing charint timestep cap and eta-window logic.

Cleanup baseline for this pass:

```bash
rtk uv run python build_extensions.py --module FS_electron_fullhide_2d --force
rtk uv run python build_extensions.py --module FS_electron_charint_2d --force
rtk /usr/bin/gfortran --version
rtk uv run python tests/readme_smoke_bench.py
rtk uv run python tests/reverse_shock_smoke.py
rtk uv run python tests/hadronic_1d_smoke.py
rtk uv run python tests/hadronic_species_transport_smoke.py
```

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
- Hadronic `pgamma_scheme` canonical names:
  - `hummer_2010_response`
  - `ka2008_reference`
  - `disabled`
- Legacy `pgamma_scheme` aliases remain accepted (compat only):
  - `hummer_2010` / `hummer2010` / `am3_reference` / `am3_numeric` / `am3_numerical` / `am3` -> `hummer_2010_response`
  - `ka2008` / `aharonian_2008` / `kelner_aharonian_2008` -> `ka2008_reference`
- Current hadronic coupling status:
  - `hummer_2010_response` now feeds `alpha_p` and `Q_p^reinj` back into the 1D proton transport update.
  - `ka2008_reference` remains a reference emission backend only; it does not provide transport feedback terms.
  - `alpha_gamma^{pgamma}` is now evaluated in the hadronic stage as a local shell photon-loss closure
    - `tau_pg(nu, r) = alpha_gamma^{pgamma}(nu, r) * R / (12 Gamma c)`
    - the local photon survival factor is applied before observer projection
    - `tau_pg` is kept as a diagnostic and is no longer passed as observer-only `tau_extra`.
  - The photon-loss transfer factor remains
    - `(1 - exp(-tau_total)) / tau_total`
    - not `exp(-tau_total)`.
- `BetheHeitler` now has two active couplings in the 1D hadronic path:
  - proton continuous cooling from `ABHHadP` is folded into the proton energy-loss operator
  - BH secondary `e±` are advanced as an independent lepton distribution and then merged into the forward electron distribution before recomputing `seed_syn`
- `PairProduction` is now wired through a formal observer-side branch:
  - it computes `alpha_{γγ}` and `Q_{e^\pm,γγ}`
  - converts `alpha_{γγ}` to `tau_pair(nu,r) = alpha_{γγ}(nu,r) * R / (12 Gamma c)`
  - and injects the resulting pair synchrotron back into the hadronic emission chain
- Hadronic photons now use the same observer transfer口径 as leptonic photons:
  - local synchrotron self-absorption uses the Fortran electron-kernel escape ratio `P_syn/P_emit`
  - hadronic photon seeds are added to the `γγ` target field before pair-production attenuation
  - do not recompute SSA transfer in Python from an independent `tau` formula
- Proton synchrotron shell emission now evaluates
  - `∫ N_p(γ) F(x) dγ`
  - on the logarithmic proton grid as `∫ N_p(γ) F(x) γ d lnγ`
  - using geometric-mean cell centers, replacing the earlier linear-`dγ` midpoint rule
- The empty `FS_hadronic_am3_1d` Fortran scaffold and its helper modules have been removed.
- Current reality:
  - active final hadronic microphysics kernels are in `src/Hadronic/`
  - Python hadronic modules are wrappers, orchestration, state assembly, and benchmark glue
  - AM3-derived benchmark code remains reference-only and must not become runtime physics

## Open Hadronic Gaps vs Done-When

- The current formal hadronic path is complete only for:
  - forward shock
  - 1D electron solvers
  - proton synchrotron
  - photopion (`hummer_2010_response`)
  - neutrino output
  - Bethe-Heitler proton cooling and BH-pair feedback into the forward electron synchrotron/seed chain
  - hadronic inverse Compton (currently proton channel active in formal runtime path)
  - proton-proton (`pp`) gamma / neutrino / proton-loss channels
  - pair production observer-side attenuation and pair synchrotron branch
  - explicit 1D secondary-species transport on a shared Lorentz-factor grid for:
    - neutron
    - pion plus / minus
    - muon minus left / right
    - muon plus left / right
  - explicit pion / muon secondary radiation in the formal runtime path:
    - pion synchrotron
    - muon synchrotron
    - pion inverse Compton
    - muon inverse Compton
  - hadronic acceleration / injection operators now active in the formal runtime path for proton injection and `gamma_p,max`
- Still not complete relative to the user's strict done-when:
  - reverse-shock hadronic processes are not implemented
  - 2D / chi-resolved hadronic transport is not implemented
  - pair production is formally coupled as an observer-side attenuation + secondary-pair branch, but not yet as a full time-dependent photon kinetic cascade PDE
  - full-resolution Vegas benchmark figures with the complete hadronic chain must still be regenerated and checked for smoothness

## Unresolved Thread Issues

- The current hadronic runtime must keep Python out of final core microphysics.
  - Final AM3-derived microphysics lives in `src/Hadronic/*.f90`
  - Python may only remain as orchestration, wrapping, benchmarking, and API glue
- The full `tests/vegas_afterglow_comparison.py` hadronic benchmark must be rerun after the explicit-species transport CFL fix.
  - Previous benchmark output included a failed placeholder image and cannot be trusted
  - Placeholder plots are now forbidden
- The hadronic-only VHE SED previously showed a non-smooth structure.
  - Until the regenerated Vegas figures are checked, this must be treated as an unresolved physics/numerics issue
- Proton-synchrotron peak placement still needs a stricter AM3-side physical cross-check.
  - The current diagnosis is that the peak is being limited mainly by the dynamical acceleration ceiling `gamma_dyn`
  - This has not yet been closed as a final validated result
- Pair production is not yet implemented as a full time-dependent cascade PDE.
  - The current formal closure is limited to observer-side attenuation plus secondary-pair feedback

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
- `tests/reverse_shock_smoke.py`
- `tests/fullhide_2d_smoke_bench.py`
- `tests/fullhide_2d_medium_diag.py`
- `tests/vegas_afterglow_comparison.py`
- Active hadronic runtime regressions:
  - `tests/hadronic_1d_smoke.py`
  - `tests/hadronic_proton_synch_1d_diag.py`
  - `tests/hadronic_pg_neutrino_1d_diag.py`
  - `tests/hadronic_species_transport_smoke.py`
  - `tests/hadronic_secondary_radiation_smoke.py`
  - `tests/hadronic_acceleration_smoke.py`
  - `tests/hadronic_bethe_heitler_smoke.py`
  - `tests/hadronic_hadronic_ic_smoke.py`
  - `tests/hadronic_pair_production_smoke.py`
  - `tests/hadronic_pp_smoke.py`
  - `tests/hadronic_pgamma_benchmark_report.py`

## Current Documentation Entrypoints

- `doc/code_overview.md`
  - 当前代码内容总览，覆盖公开 API、Python 编排层、Fortran 内核层、hadronic 边界、构建入口和测试入口。
- `doc/source_tree.md`
  - 精简源码树索引。
- `doc/call_chain.md`
  - 当前 Python / Fortran 调用链图。
- `doc/electron_solver_algorithms.md`
  - 前向激波电子求解器算法细节。
- `doc/hadronic_pgamma_notes.md`
  - hadronic / `pγ` 当前状态与边界。

## Cleanup Policy

- Do not commit `.buildcache/`.
- Do not commit `output/benchmark_exp_tail/`.
- Keep generated documentation plots in `output/asgard_doc/vegas_afterglow_compare/` only when they are current benchmark artifacts.
- Remove one-off debug scripts after their findings have been folded into code, tests, or this note.
- Never generate placeholder plots on failure.
  - If a figure build fails, raise the real exception and stop.
  - Do not save error text, notes, or fallback images in place of the requested plot.

## External Reference Code

- A local AM3 reference checkout is available at:
  - `/mnt/c/Users/jia/Documents/New project/_external/am3_reference`
- Clone provenance:
  - `https://gitlab.desy.de/am3/am3.git`
  - current local reference HEAD:
    - `7aba970b230e24d6cfe327826522084b70ed406e`
- This checkout is reference-only. Do not mix its build system, Python package layout, or solver grid assumptions into ASGARD unless the user explicitly asks for a direct port.

## AM3 / ASGARD Coexistence Rules

- AM3 should be used as a microphysics and process-structure reference, not as a template for replacing ASGARD's existing electron solver chain.
- Preserve ASGARD's current forward-shock flow:
  - `dynamics -> electron solver -> photon target field -> hadronic add-on -> observer projection`
- Do not replace ASGARD's 1D/2D electron transport, shell history handling, or observer-side projection with AM3's homogeneous one-zone PDE structure.
- When borrowing from AM3, only port the physically targeted modules needed for hadronic coupling:
  - photopion interaction kernel
  - pion decay
  - muon decay
  - optional photon-loss term due to photopion interactions
- Keep the coexistence boundary explicit:
  - ASGARD remains the driver for blast-wave dynamics, shell evolution, and observer mapping
  - AM3-derived code only supplies source/loss operators on top of ASGARD's local hadronic state
- Any AM3-derived backend must be benchmarked against the local reference checkout above before it is described as `AM3-consistent`.
- New hard constraint from the user:
  - every AM3-derived microphysics implementation must be written in Fortran under `src/Hadronic/`
  - Python may orchestrate, wrap, benchmark, and expose state, but must not remain the home of the final AM3-transcribed physics kernels
  - existing Python AM3/Hümmer/KA2008 ports are transitional references only and must be replaced by Fortran backends before the work can be considered done
