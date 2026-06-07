# ASGARD 文档总入口

本文档是当前工作树的统一入口。它服务三类读者：

- 使用者：安装项目、运行光变/谱/偏振/天图、读取输出。
- 研究者：理解物理模块、开关、边界和验收口径。
- 开发者：修改 Fortran/Python 主链、重建扩展、刷新 benchmark、提交可追溯结果。

## 推荐阅读路径

首次使用：

1. `README.md`
2. `doc/installation.md`
3. `doc/user_guide.md`
4. `doc/public_api.md`

理解物理和数值主链：

1. `doc/physics_model.md`
2. `doc/numerical_methods.md`
3. `doc/code_overview.md`
4. `doc/call_chain.md`

开发或刷新基准图：

1. `AGENTS.md`
2. `PLAN.md`
3. `TODO.md`
4. `doc/developer_guide.md`
5. `doc/validation_and_benchmarks.md`
6. `doc/benchmark_refresh_protocol.md`

专题决策记录：

- `doc/public_backend_limits.md`：public API 与 backend 的不支持/部分支持边界。
- `TODO.md`：唯一 TODO / 未完成项入口。
- `doc/hadronic_chi_transport_decision.md`：当前不实现 2D / chi-resolved hadronic transport 的理由和前置物理契约。
- `doc/pair_cascade_extension_boundary.md`：当前 gamma-gamma pair/synch cascade 与 IC-mediated electromagnetic cascade 的边界。
- `doc/polarization_timing_diagnostic.md`：Lan 2023 偏振峰时诊断。
- `doc/hadronic_pgamma_notes.md`：p-gamma 微物理和基准说明。
- `doc/am3_migration_plan.md`：AM3 共存、迁移和引用边界。
- `doc/electron_solver_algorithms.md`：电子输运算法说明。

## 当前能力摘要

ASGARD 当前主线是 GRB afterglow 的 shell-evolving blast-wave + observer projection 模型。Python 层负责 public API、配置、状态机、后处理、拟合和 benchmark；Fortran 层负责高代价数值核。

已登记并可用的主功能：

- Forward-shock dynamics、电子输运、synchrotron、SSC、SSA、gamma-gamma absorption、observer projection。
- 1D 电子求解器：`fullhide_1d`, `slc1_1d`, `charint_1d`, `t2g1_1d`, `weno5_1d`。
- 2D 电子求解器：`fullhide_2d`, `charint_2d`。
- Reverse shock electron synchrotron、RS SSC、FS/RS cross-zone IC。
- Reverse shock thermal/magnetized baseline：使用 shock-front `gamma34` 注入能标、显式 `U3/V3` thermal state、可选 upstream `sigma` 和 ordered magnetic component。
- Forward-shock hadronic `legacy_1d` 与 formal `am3_1d` research path。
- Reverse-shock hadronic light proton-synch path，以及复用 formal 1D hadronic kernels 的 full-chain dispatch。
- Gamma-gamma pair production 与 shell-sequence time-dependent pair/synch cascade。
- 同步辐射偏振 Stokes 投影，覆盖 FS/RS electron synch 与 FS/RS hadronic synch。
- `Model` public API、`Fitter` 拟合 API、benchmark/report 脚本。

未完成项、不支持 backend 边界和实现准入条件集中维护在根目录 `TODO.md` 与 `doc/public_backend_limits.md`，避免文档间分散待办列表。

## 核心入口

Public Python API：

- `ASGARD/api_model.py`：`Model`, `ISM`, `Wind`, jet classes, `Observer`, `Radiation`, `Setups`。
- `ASGARD/api_observe.py`：`observe`, `run_fit` 以及 patch/observer 投影内部实现。
- `ASGARD/api_fit.py`：`Fitter`, `Param`, `FitResult`。

Fortran 构建入口：

- `build_extensions.py`

默认强制编译命令：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
```

最小文档/格式检查：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git diff --check'
```

## 文档维护原则

- 文档必须描述当前工作树，而不是理想设计。
- 物理边界必须写明原因、开关、验收口径和不可用范围。
- Benchmark 图像和 CSV 只有在能由脚本复现且命令已记录时才进入版本库。
- Fortran 或物理路径改动后的验证必须记录编译命令、line-truncation 检查和最小 smoke test。
- 不使用经验 smoothing、fallback、后处理补丁掩盖非连续或非光滑结果；物理量随时间/空间不光滑时优先查 bug。
