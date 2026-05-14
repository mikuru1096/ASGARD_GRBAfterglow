# ASGARD Work Plan

本文档是当前工作树的后续工作基线。实现状态细节见 `doc/code_overview.md`，电子算法细节见 `doc/electron_solver_algorithms.md`，AM3 迁移历史见 `doc/am3_migration_plan.md`。

## Current Baseline

- Public runtime: forward-shock afterglow API、电子同步/SSC/SSA、`gamma-gamma` 吸收、observer projection 和 benchmark workflow 可用。
- Electron solvers: 1D `fullhide/slc1/charint/t2g1/weno5` 与 2D `fullhide_2d/charint_2d` 均在当前工作树登记；`fullhide_1d` 是默认稳定基线。
- Reverse shock: electron synchrotron、RS SSC、FS/RS cross-zone IC 已接入；RS 注入能标使用 shock-front `gamma34`，区域 3 turbulent field 和 post-crossing 热演化使用显式 `U3/V3` thermal state；可选 upstream magnetization `sigma` 已接入 Python/Fortran 动力学契约，使用 Vegas MHD jump 的 `sigma` 依赖和 upstream ordered-field 公式，且以 ASGARD 非磁化 jump condition 固定 `sigma -> 0` 极限；RS hadronic light path 覆盖 1D proton injection/transport + proton synchrotron，full-chain path 复用正式 1D hadronic kernels 覆盖 RS pγ/BH/pp/secondary/cascade coupling。
- Forward-shock hadronic: `legacy_1d` 覆盖 proton transport + proton synchrotron；`am3_1d` 是当前正式研究路径，覆盖 p-gamma、BH、pp、hadronic IC、secondary species transport、secondary radiation、pair production branch、neutrino 等过程。`Radiation.pair_production=True` 且 `Setups.pair_cascade_iterations>1` 时，pair branch 使用 shell-sequence time-dependent γγ pair/synch cascade。
- Polarization: `Model.polarization(...)` 已实现同步辐射 Stokes 路径，覆盖 FS/RS electron synch 与 FS/RS hadronic synch；非同步分支不混入偏振 Stokes。
- AM3: 只作为 hadronic 微物理核和 benchmark 参考；ASGARD 的 dynamics/electron/observer 主链不由 AM3 替代。

## Todo

1. Magnetized-RS full-MHD upgrade decision
   - 当前已完成基线：`ReverseShockConfig.sigma` / `Setups.reverse_sigma` 控制 upstream magnetization；`Dynamics_reverse` 在 `f_e_r` 后接收 `sigma_r`，并在 `gamma_m_cross` 后输出 `B3_ordered_cross`，Python 暴露为 `ReverseShockDynamics.ordered_magnetic_cross_g`。
   - 当前闭合：Vegas MHD jump condition 只提供压缩比的 `sigma` 依赖；ASGARD 用 `(4 gamma34 + 3)` 固定零磁化极限；`B3 = sqrt(8 pi epsilon_B,r U3/V3) + B3_ordered`，crossing 后有序场随 `V3_cross/V3` 膨胀。
   - 后续若要修改 `dU3_shock=(gamma34-1) dm3 c^2`，必须先给出完整 MHD 能量方程和文献 benchmark；不要用光变匹配或 Vegas 平均热态作为替代理由。
   - 增加磁化 RS benchmark 图时只报告 `sigma` 扫描下 `B3/gamma34/U3/V3/nu_m/nu_c/nu_a` 平滑性和 `sigma -> 0` 极限，不设 “match Vegas light curve” 为验收条件。

2. RS thermal-state diagnostics hardening
   - 当前 RS 动力学显式输出 `U3_th/V3_comoving/Gamma34_inst`；`B3=sqrt(8 pi epsilon_B,r U3/V3)`。
   - `details().rvs.nu_m` 是 diagnostic break：crossing 前用局部 `gamma34` 注入能标，crossing 后按 `U3/M3` 的绝热能量比演化，不替代输运后的电子谱峰。
   - 后续若要暴露更多 public detail，需要把 `U3/V3/gamma34` 作为诊断状态传出，而不是重新从 Vegas-style 平均热态反推。

3. RS hadronic full-chain validation
   - 为 RS pγ/BH/pp/secondary/cascade 分别补最小能量预算诊断，确认 RS seed photons、`B3`、shell energy、baryon target density 的单位和 comoving/lab frame 契约。
   - 先验证过程级功率与粒子数演化的平滑性，再考虑更高维 transport。

4. 2D / chi-resolved hadronic transport decision
   - 先决定是否需要和 `fullhide_2d/charint_2d` 共用 `chi_grid` 与历史 photon field。
   - 若不能提供可检验的时空平滑诊断和能量守恒 benchmark，不进入实现。

5. Pair cascade extension boundary
   - 当前主链已有 shell-sequence time-dependent γγ pair/synch cascade；legacy iterative pair-production synch kernel 只保留作诊断。
   - 若要继续扩展到 IC-mediated electromagnetic cascade，必须先给出闭合的 photon/e± source-sink 方程和能量守恒 benchmark。

6. Polarization timing diagnostic
   - Lan, Wu & Dai 2023 对照中峰值幅度基本一致，峰时偏早仍需定位。
   - 优先检查 dynamics/EATS 时间映射和真实 jet 面元权重，不做后处理平滑。

7. Benchmark baseline refresh protocol
   - 重新生成 Vegas / literature benchmark 时，必须同时记录命令、git 基线、模块 build 状态和图像生成脚本。
   - RS benchmark 图必须标明 ASGARD 使用 “local gamma34 injection + explicit U3/V3 thermal state”，Vegas 使用平均 thermal-state closure。
   - VegasAfterglow 只作为 comparison backend；若 ASGARD 与 Vegas 差异来自平均 `Gamma_th` 与局部 shock-front injection 的模型差异，不作为失败。
   - 不提交失败占位图、临时 debug 脚本或 `.buildcache/`。

8. Public API/backend limits
   - Jet spreading、自定义 `Medium` kernel dispatch、wind `k != 2`、thermal electrons outside `fullhide_1d` 都是明确未支持边界。
   - 这些边界只有在目标观测或物理问题需要时才进入实现。

## Validation Baseline

文档改动只需跑 `git diff --check`。Fortran 或物理路径变更按受影响范围执行：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module FS_electron_fullhide_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module FS_hadronic_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/hadronic_1d_smoke.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/hadronic_reverse_shock_smoke.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/polarization_smoke.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --extra compare python tests/vegas_afterglow_comparison.py --scenario baseline --only reverse_shock_lc reverse_shock_thermal'
```

当前 RS benchmark artifacts:
- `output/asgard_doc/vegas_afterglow_compare/compare_reverse_shock_lc.png`
- `output/asgard_doc/vegas_afterglow_compare/compare_reverse_shock_thermal_benchmark.png`
