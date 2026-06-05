# AM3 迁移计划（历史参考）

本文档记录 ASGARD 中 AM3 微物理迁移的历史规划。大部分 P0/P1 项已完成；当前权威 TODO 以 `PLAN.md` 为准。

## 已完成迁移

- AM3-derived hadronic core kernels 已迁移到 `src/Hadronic/*.f90`，覆盖 photopion、decay、Bethe-Heitler、pair production、pp、hadronic IC、secondary species transport、secondary radiation 和 acceleration/injection。
- Local photopion photon-loss closure 已实现：`α_γ^{pγ}` → `τ_{pγ} = α_γ · R/(12 Γ c)` → shell survival factor `(1-e^{-τ})/τ`。
- Explicit secondary species transport 已覆盖 n、π±、μ± 的左右向输运。
- Secondary radiation 已覆盖 pion/muon synchrotron + IC。
- Hadronic acceleration/injection operators 已接入。
- Shell-sequence time-dependent γγ pair/synch cascade path 已接入。
- Reverse-shock full-chain hadronic dispatch 已通过 formal 1D kernels 接入。

## 保留缺口

- 2D / chi-resolved hadronic transport。
- 含完整 hadronic chain 的 Vegas benchmark figures 刷新。
- 超出当前 γγ pair/synch contract 的 inverse-Compton-mediated electromagnetic pair cascade。

RS hadronic proton injection/transport + proton synchrotron 已由 `hadronic_reverse_1d` 覆盖。Full-chain RS pγ/BH/pp/secondary/cascade dispatch 通过 runtime wrapper 复用 `hadronic_forward_1d` formal kernels。Legacy iterative pair-production synch branch 只保留作诊断；主 cascade path 是 shell-sequence time-dependent γγ pair/synch cascade。

反激波动力学基线使用局部 `gamma34` 注入和显式 region-3 `U3/V3` thermal-state evolution。VegasAfterglow 是 comparison backend，不是物理目标。

## 架构边界

- ASGARD = shell-evolving blast-wave + observer projection。
- AM3 = microphysics/kernel reference。
- 不用 AM3 one-zone driver 替代 ASGARD 的 dynamics/electron/observer chain。
- 最终 AM3-derived microphysics 必须写入 `src/Hadronic/` 下的 Fortran。

## AM3 参考代码

- `/mnt/c/Users/jia/Documents/New project/_external/am3_reference`，HEAD 为 `7aba970b`。
- `~/projects/_external/am3_reference` 是 WSL home mirror，HEAD 相同。
