# AM3 Migration Plan (Historical Reference)

本文档记录 ASGARD 中 AM3 微物理迁移的历史规划。大部分 P0/P1 项已完成。

## Completed Migration

- AM3-derived hadronic core kernels → `src/Hadronic/*.f90` (photopion, decay, Bethe-Heitler, pair production, pp, hadronic IC, secondary species transport, secondary radiation, acceleration/injection)
- Local photopion photon-loss closure: `α_γ^{pγ}` → `τ_{pγ} = α_γ · R/(12 Γ c)` → shell survival factor `(1-e^{-τ})/τ`
- Explicit secondary species transport: n, π±, μ± (left/right)
- Secondary radiation: pion/muon synchrotron + IC
- Hadronic acceleration/injection operators

## Remaining Gaps

- Full time-dependent pair cascade PDE
- Reverse-shock pγ/BH/pp/secondary/cascade coupling
- 2D / χ-resolved hadronic transport
- Regenerated Vegas benchmark figures with complete hadronic chain

RS hadronic proton injection/transport + proton synchrotron now exists in `FS_hadronic_reverse_1d`. An iterative pair-production synch branch exists, but it is not the full time-dependent pair cascade PDE.

## Architecture Boundary

- ASGARD = shell-evolving blast-wave + observer projection
- AM3 = microphysics/kernel reference only
- Do NOT replace ASGARD's dynamics/electron/observer chain with AM3's one-zone driver
- Final AM3-derived microphysics must be Fortran under `src/Hadronic/`

## AM3 Reference

- `/mnt/c/Users/jia/Documents/New project/_external/am3_reference` (HEAD: `7aba970b`)
- `~/projects/_external/am3_reference` is a WSL-home mirror at the same HEAD.
