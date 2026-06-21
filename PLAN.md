# ASGARD Work Plan

本文档是当前工作树的后续工作基线。实现状态细节见 `doc/code_overview.md`，电子算法细节见 `doc/electron_solver_algorithms.md`，AM3-derived hadronic 说明见 `doc/hadronic_pgamma_notes.md`。

## Current Baseline

- Public runtime: forward-shock afterglow API、电子同步/SSC/SSA、`gamma-gamma` 吸收、observer projection 和 benchmark workflow 可用。
- Electron solvers: 1D `fullhide/slc1/charint/t2g1/weno5/dg_1d` 与 2D `fullhide_2d/charint_2d` 均在当前工作树登记；`fullhide_1d` 是默认稳定基线，`dg_1d` 是 opt-in 高阶 1D 谱元路径。`dg_1d` 默认启用 troubled-cell positive-kernel 滤波，只在高阶模态能量或负值暴露的局部谱元及相邻谱元上衰减高阶 Legendre 模态，保留 cell average；`ASGARD_DG1D_POSITIVE_KERNEL=0` 仅作为诊断关闭开关。`charint_2d` 的 downstream chi 能量列使用公共无源特征线冷却 primitive，shock-front 注入和 eta 输运仍保持 2D transport 方程。
- Reverse shock: electron synchrotron、RS SSC、FS/RS cross-zone IC 已接入；RS 注入能标使用 shock-front `gamma34`，区域 3 turbulent field 和 post-crossing 热演化使用显式 `U3/V3` thermal state；可选 upstream magnetization `sigma` 已接入 Python/Fortran 动力学契约，使用 Vegas MHD jump 的 `sigma` 依赖和 upstream ordered-field 公式，且以 ASGARD 非磁化 jump condition 固定 `sigma -> 0` 极限；RS hadronic light path 覆盖 1D proton injection/transport + proton synchrotron，full-chain path 复用正式 1D hadronic kernels 覆盖 RS pγ/BH/pp/secondary/cascade coupling。
- Forward-shock hadronic: `legacy_1d` 覆盖 proton transport + proton synchrotron；`am3_1d` 是当前正式研究路径，覆盖 p-gamma、BH、pp、hadronic IC、secondary species transport、secondary radiation、pair production branch、neutrino 等过程。`Radiation.pair_production=True` 且 `Hadronic.pair_cascade_iterations>1` 时，pair branch 使用 shell-sequence time-dependent γγ pair/synch cascade。
- Polarization: `Model.polarization(...)` 已实现同步辐射 Stokes 路径，覆盖 FS/RS electron synch 与 FS/RS hadronic synch；非同步分支不混入偏振 Stokes。
- AM3: 只作为 hadronic 微物理核和 benchmark 参考；ASGARD 的 dynamics/electron/observer 主链不由 AM3 替代。

## Completed RS Work

- Magnetized RS baseline 已完成：`ReverseShock.upstream_sigma` 控制 upstream magnetization；`Dynamics_reverse` 接收 `sigma_r`，Python 暴露 total `B3` 与 `ReverseShockDynamics.ordered_magnetic_cross_g`；`sigma -> 0` 固定回到当前非磁化 jump baseline。
- RS thermal-state baseline 已完成：动力学显式输出 `U3_th/V3_comoving/Gamma34_inst`；`B3` 由 turbulent `sqrt(8 pi epsilon_B,r U3/V3)` 加 ordered upstream component 给出；断频只通过 `nu_callback` 临时收集，不写入常规 `details()`。
- FS/RS shared 1D electron transport baseline 已完成：`electron_solver="dg_1d"` 同时驱动 FS 与 primary RS 电子输运；FS/RS 共用 solver id 和 shell transport primitive，RS kinetic 注入保持 \((\gamma-1)^{-p}\)，FS 注入保持 \(\gamma^{-p}\)。RS DG 状态跨 shell 持久保存，并在输出到固定 \(\log\gamma\) 网格时使用守恒 cell 投影。默认 troubled positive-kernel 只抑制局部 Gibbs 振荡，不把物理尖 break 或窄 cutoff 当作失败。
- RS DG 验收基线已固定：`num_gam_e=121` 下 DG 高能尾比 fullhide 窄；`num_gam_e=241` fullhide 高能尾同步收缩，说明 121 格 fullhide 宽尾主要是低阶隐式上风数值扩散。当前不把 121 格 fullhide 高能尾作为 DG 强贴合目标；正式判断以高分辨率 fullhide、post-crossing direct characteristic map 和守恒谱形诊断为准。post-crossing direct map 使用闭合低能边界，避免纯冷却电子数从网格低能端流失。最新保留图见 `output/asgard_doc/dg_radiation_stability_scan/postcross_direct_map_effective_support/`，早期分辨率对照图保留在 `output/asgard_doc/reverse_dg_fullhide_electron_spectrum/`。
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
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_dg_1d --module electron_reverse_kernel --module structured_jet_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_charint_2d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_forward_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git diff --check'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --with "mkdocs<2" --with "mkdocs-material>=9.5" mkdocs build --strict --site-dir /tmp/asgard_mkdocs_site'
```

当前已知验证阻塞不计入普通文档基线：RS full-chain hadronic 和 joint feedback 会在 formal hadronic electron-energy 网格契约处失败，报 `electron_energy_gev must be logarithmically uniform`；RS DG 谱形诊断当前暴露 sawtooth-turn 问题。这些是待修真实问题，不能用删除断言、调宽阈值或 fallback 掩盖。

当前 Vegas baseline benchmark artifacts:
- `output/asgard_doc/vegas_afterglow_compare/compare_reverse_shock_lc.png`
- `output/asgard_doc/vegas_afterglow_compare/compare_reverse_shock_thermal_benchmark.png`
- `output/asgard_doc/vegas_afterglow_compare/compare_speed_profile.png`
- `output/asgard_doc/vegas_afterglow_compare/compare_basic_lc_spec.png`

当前 DG/RS electron baseline artifacts:
- `output/asgard_doc/dg_1d_benchmark_troubled/dg_troubled_fullhide_lightcurves.png`
- `output/asgard_doc/dg_1d_benchmark_troubled/dg_troubled_fullhide_lightcurve_ratios.png`
- `output/asgard_doc/dg_1d_benchmark_troubled/dg_troubled_fullhide_electron_spectra.png`
- `output/asgard_doc/dg_1d_benchmark_troubled/dg_troubled_fullhide_radiation_seds.png`
- `output/asgard_doc/dg_1d_sigma_benchmark_troubled/dg_troubled_sigma_fullhide_lightcurves.png`
- `output/asgard_doc/dg_1d_sigma_benchmark_troubled/dg_troubled_sigma_fullhide_lightcurve_ratios.png`
- `output/asgard_doc/dg_1d_sigma_benchmark_troubled/dg_troubled_sigma_fullhide_electron_spectra.png`
- `output/asgard_doc/dg_1d_sigma_benchmark_troubled/dg_troubled_sigma_fullhide_radiation_seds.png`
- `output/asgard_doc/dg_radiation_stability_scan/postcross_direct_map_effective_support/dynamic_upper_problem_case_electron_spectra.png`
- `output/asgard_doc/dg_radiation_stability_scan/postcross_direct_map_effective_support/fixed_previous_shell_dynamic_upper_electron_spectra.png`
- `output/asgard_doc/dg_radiation_stability_scan/postcross_direct_map_effective_support/dynamic_upper_problem_case_lightcurve_slices.png`
- `output/asgard_doc/dg_radiation_stability_scan/postcross_direct_map_effective_support/dynamic_upper_problem_case_summary.csv`
