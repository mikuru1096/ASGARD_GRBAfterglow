# ASGARD Work Plan

本文档是当前工作树的后续工作基线。实现状态细节见 `doc/code_overview.md`，电子算法细节见 `doc/electron_solver_algorithms.md`，AM3 迁移历史见 `doc/am3_migration_plan.md`。

## Current Baseline

- Public runtime: forward-shock afterglow API、电子同步/SSC/SSA、`gamma-gamma` 吸收、observer projection 和 benchmark workflow 可用。
- Electron solvers: 1D `fullhide/slc1/charint/t2g1/weno5` 与 2D `fullhide_2d/charint_2d` 均在当前工作树登记；`fullhide_1d` 是默认稳定基线。
- Reverse shock: electron synchrotron、RS SSC、FS/RS cross-zone IC 已接入；RS 注入能标使用 shock-front `gamma34`，区域 3 turbulent field 和 post-crossing 热演化使用显式 `U3/V3` thermal state；可选 upstream magnetization `sigma` 已接入 Python/Fortran 动力学契约，使用 Vegas MHD jump 的 `sigma` 依赖和 upstream ordered-field 公式，且以 ASGARD 非磁化 jump condition 固定 `sigma -> 0` 极限；RS hadronic light path 覆盖 1D proton injection/transport + proton synchrotron，full-chain path 复用正式 1D hadronic kernels 覆盖 RS pγ/BH/pp/secondary/cascade coupling。
- Forward-shock hadronic: `legacy_1d` 覆盖 proton transport + proton synchrotron；`am3_1d` 是当前正式研究路径，覆盖 p-gamma、BH、pp、hadronic IC、secondary species transport、secondary radiation、pair production branch、neutrino 等过程。`Radiation.pair_production=True` 且 `Hadronic.pair_cascade_iterations>1` 时，pair branch 使用 shell-sequence time-dependent γγ pair/synch cascade。
- Polarization: `Model.polarization(...)` 已实现同步辐射 Stokes 路径，覆盖 FS/RS electron synch 与 FS/RS hadronic synch；非同步分支不混入偏振 Stokes。
- AM3: 只作为 hadronic 微物理核和 benchmark 参考；ASGARD 的 dynamics/electron/observer 主链不由 AM3 替代。

## Completed RS Work

- Magnetized RS baseline 已完成：`ReverseShock.upstream_sigma` 控制 upstream magnetization；`Dynamics_reverse` 接收 `sigma_r`，Python 暴露 total `B3` 与 `ReverseShockDynamics.ordered_magnetic_cross_g`；`sigma -> 0` 固定回到当前非磁化 jump baseline。
- RS thermal-state baseline 已完成：动力学显式输出 `U3_th/V3_comoving/Gamma34_inst`；`B3` 由 turbulent `sqrt(8 pi epsilon_B,r U3/V3)` 加 ordered upstream component 给出；`details().rvs.nu_m` 是 diagnostic break，不替代输运后的电子谱峰。
- RS hadronic baseline 已完成：light path 覆盖 RS proton injection/transport + proton synchrotron；full-chain path 复用正式 1D hadronic kernels，使用 RS seed photons、RS `B3`、shell energy 和 baryon target density。
- 旧 Python benchmark/comparison 脚本已移除。新的正式刷新入口必须先进入 `tests/` 并给出明确物理验收口径。VegasAfterglow 是 comparison backend，不是光变目标或物理基准。

## Completed Non-RS Decisions

- 2D / χ-resolved hadronic transport: 决策记录见 `doc/hadronic_chi_transport_decision.md`。当前不实现；正式 hadronic path 保持 1D shell 契约，直到 χ-local photon、hadron、secondary feedback 和 observer projection 契约完成。
- Pair cascade extension: 决策记录见 `doc/pair_cascade_extension_boundary.md`。当前主链停在 shell-sequence time-dependent γγ pair/synch cascade；不进入 IC-mediated electromagnetic cascade，除非先给出 photon/e± source-sink 方程、IC kernel 契约和能量守恒 benchmark。
- Polarization timing diagnostic: Lan 2023 overlay 的峰时偏早当前归入 dynamics/jet-evolution 模型边界；旧 overlay 诊断脚本和文档已删除。
- Benchmark baseline refresh protocol: 重新生成 Vegas、literature 或 hadronic benchmark 时必须记录命令、HEAD、tracked diff 状态、扩展 build 状态、输出路径和物理验收口径；旧 Python benchmark 脚本不再作为正式入口。
- Public API/backend limits: 边界记录见 `doc/public_backend_limits.md`。Jet spreading、自定义 `Medium` kernel dispatch、wind `k != 2`、thermal electrons outside `fullhide_1d` 都是明确未支持边界，只有在目标观测或物理问题需要时才进入实现。

## Todo

当前唯一 TODO / 未完成项入口是 `TODO.md`。本文档只记录基线和已完成决策，不再维护分散待办列表。

## Validation Baseline

文档改动只需跑 `git diff --check`。Fortran 或物理路径变更按受影响范围执行：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_charint_2d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_forward_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/hadronic_1d_smoke.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/hadronic_reverse_shock_smoke.py'
```

当前 Vegas baseline benchmark artifacts:
- `output/asgard_doc/vegas_afterglow_compare/compare_reverse_shock_lc.png`
- `output/asgard_doc/vegas_afterglow_compare/compare_reverse_shock_thermal_benchmark.png`
- `output/asgard_doc/vegas_afterglow_compare/compare_speed_profile.png`
- `output/asgard_doc/vegas_afterglow_compare/compare_basic_lc_spec.png`
