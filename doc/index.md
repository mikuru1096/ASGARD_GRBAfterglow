# ASGARD 文档总入口

<p class="asgard-home-logo">
  <img src="assets/logo.png" alt="ASGARD logo">
</p>

<p class="asgard-home-tagline-cn">
  ASGARD 得名于北欧神话中的阿斯加德彩虹桥。它的目标是沿着爆波动力学、粒子输运和观测者投影这条主链，给出可追溯的伽马射线暴余辉计算。
</p>

本文档是当前工作树的网页入口。它服务三类读者：

- 使用者：安装项目、运行光变/谱/偏振/天图、读取输出。
- 研究者：理解物理模块、开关、边界和验收口径。
- 开发者：修改 Fortran/Python 主链、重建扩展、刷新基准测试、提交可追溯结果。

公开入口：

- 项目 README：<https://github.com/mikuru1096/ASGARD_private#readme>
- 网页文档：<https://hetools.cn/asgard-doc/>

## 推荐阅读路径

网页导航按用途分为五组。

入门：

1. `doc/terminology.md`：术语与中英混排规则。
2. `doc/installation.md`：安装和构建。
3. `doc/quickstart.md`：最小光变、谱和频段积分。
4. `doc/command_line.md`：构建、作图、文档和基准测试命令。
5. `doc/troubleshooting.md`：常见失败模式。

建模与接口：

1. `doc/examples.md`：多频光变、宽频谱、辐射分量、内部量和观测预测教程。
2. `doc/parameter_reference.md`：参数含义和选择建议。
3. `doc/public_api.md`：公开 API 选择手册。
4. `doc/fitting_workflow.md`：从数据到似然函数和拟合结果的完整路径。
5. `doc/mcmc_fitting.md`：emcee 与 PyMultiNest 专题。
6. `doc/external_inference.md`：外部采样器边界和包装方式。
7. `doc/magnetized_rs_dg1d_tutorial.md`：磁化反向激波、`dg_1d` 高阶电子输运和当前收敛阶。
8. `doc/prompt_internal_shock_tutorial.md`：prompt internal-shock snapshot 的物理推导、代码入口、formal plotting 和边界。

物理模型：

1. `doc/project_physics_design.md`：全项目物理设计总纲。
2. `doc/physical_processes.md`：动力学、电子、辐射、强子、级联和投影的过程说明。
3. `doc/physics_model.md`：公开模型能力和边界。
4. `doc/joint_secondary_feedback_physics.md`：二级粒子联合反馈物理。
5. `doc/pair_cascade_extension_boundary.md`：对级联边界。
6. `doc/hadronic_chi_transport_decision.md`：强子 \(\chi\) 分辨输运暂不实现的决策。
7. `doc/hadronic_pgamma_notes.md`：\(p\gamma\) 微物理和基准说明。

算法与数值：

1. `doc/project_algorithm_design.md`：全项目算法设计总纲。
2. `doc/algorithm_workflow.md`：数组维度、离散方程、缓存和验证矩阵。
3. `doc/physics_algorithm_crosswalk.md`：物理问题、离散变量、Fortran 程序单元和验收指纹的交叉指南。
4. `doc/shock_shell_adaptive_algorithms.md`：有限厚壳层、反向激波和自适应网格。
5. `doc/numerical_methods.md`：数值方法和构建检查。
6. `doc/electron_solver_algorithms.md`：电子输运算法。
7. `doc/fullhide2d_pwn_cr_transport.md`：`fullhide2d_transport_model="pwn_cr_v1"` 的物理契约。
8. `doc/joint_secondary_feedback_algorithm.md`：`electron_photon_coupling="joint"` 的状态机和数组契约。

开发、验证与发布：

1. `AGENTS.md`、`PLAN.md`、`TODO.md`：工作规则、当前计划和唯一未完成项入口。
2. `doc/developer_guide.md`：开发流程。
3. `doc/validation_and_benchmarks.md`：构建门禁、冒烟测试和基准刷新。
4. `doc/public_backend_limits.md`：公开 API 与后端能力边界。
5. `doc/code_overview.md`、`doc/source_tree.md`、`doc/call_chain.md`：代码结构和调用链。
6. `doc/fortran_kernel_index.md`：按源文件逐个列出 Fortran `module`、`subroutine` 和 `function` 的程序单元索引。
7. `doc/web_docs.md`：网页文档发布和 HEtools 托管流程。

除 `README.md`、`AGENTS.md`、`PLAN.md` 和 `TODO.md` 这类根目录入口/开发状态文件外，用户说明、物理说明、算法设计、API 参考、教程、拟合说明和验证说明都应维护在本 `doc/` 网页文档树中。专题计划若已经沉淀为当前能力或明确边界，应合并到对应网页章节，不再保留根目录副本文档。

## 当前能力摘要

ASGARD 当前主线是 GRB 余辉的壳层演化爆波和观测者投影模型。Python 层负责公开 API、配置、状态机、后处理、拟合和基准测试；Fortran 层负责高代价数值核。

已登记并可用的主功能：

- 正向激波动力学、电子输运、同步辐射、同步自康普顿、同步自吸收、\(\gamma\gamma\) 吸收和观测者投影。
- 1D 电子求解器：`fullhide_1d`, `slc1_1d`, `charint_1d`, `dg_1d`, `t2g1_1d`, `weno5_1d`；`dg_1d` 是正向/反向激波共享、需要显式启用的 P12 LGL-DG 路径，默认使用问题单元正性核。
- 2D 电子求解器：`fullhide_2d`, `charint_2d`。
- `chi_eats_2d` 观测者投影：FS 同步辐射 + SSA 使用 finite \(q\)-shell 的局域半径、bulk Lorentz factor 和权重；`projection_kind="lightcurve"` 走专用光变投影，`projection_kind="sed"` 走通用 SED 插值器。
- 反向激波电子同步辐射、RS SSC、FS/RS 跨区逆康普顿。
- 反向激波热/磁化基线：使用激波前沿 `gamma34` 注入能标、显式 `U3/V3` 热状态、可选 `ReverseShock.upstream_sigma`、有限强度 MHD 跳跃条件和有序磁场分量。
- 正向激波强子 `legacy_1d` 与正式研究路径 `am3_1d`。
- `electron_photon_coupling="joint"` 是需要显式启用的壳层级联合闭合：正向激波电子、光子场、正式强子输运、BH/pp/\(\gamma\gamma\) 二级 \(e^\pm\) 和光子存活率在同一 \(R\) 网格上迭代；当前在 formal hadronic electron-energy grid contract 处失败，修复入口见 `TODO.md`。
- 反向激波强子轻量质子同步辐射路径可执行；复用正式 1D 强子核的完整链路调度已接入，但当前 RS full-chain 在 formal hadronic electron-energy grid contract 处失败。
- \(\gamma\gamma\) 对产生与壳层序列含时对产生/同步辐射级联。
- 同步辐射偏振 Stokes 投影，覆盖 FS/RS 电子同步辐射与 FS/RS 强子同步辐射。
- `prompt/` 内部激波快照诊断，覆盖两壳碰撞、磁化 FS/RS jump、同步/SSC/\(\gamma\gamma\) 和 EATS 投影；它不是 `asgard_core` 顶层 public API 或拟合入口。
- `Model` 公开 API、`Fitter` 拟合 API、基准测试和报告脚本。

未完成项、不支持的后端边界和实现准入条件集中维护在根目录 `TODO.md` 与 `doc/public_backend_limits.md`，避免文档间分散待办列表。

## 核心入口

Public Python API：

- `asgard_core/api_model.py`：`Model`, `UniformMedium`, `WindMedium`, `TabulatedMedium`, `top_hat_jet`, `gaussian_jet`, `power_law_jet`, `Observer`, `Radiation`, `Numerics`, `ObserverGrid`, `SolverOptions`, `ReverseShock`, `Hadronic`，以及 `Model` 查询调度和 direct/patch solve 入口。
- `asgard_core/api_observe.py`：内部/旧配置观测工具，以及 sky image / polarization / observation dataset helpers；`observe` 和 `run_fit` 不从 `asgard_core` 顶层导出，不作为新教程入口。
- `asgard_core/api_fit.py`：`Fitter`, `Param`, `FitResult`。

Fortran 构建入口：

- `build_extensions.py`
- Fortran 程序单元级指南：`doc/fortran_kernel_index.md`；物理到算法的交叉索引：`doc/physics_algorithm_crosswalk.md`。

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
- 基准测试图像和 CSV 只有在能由脚本复现且命令已记录时才进入版本库。
- Fortran 或物理路径改动后的验证必须记录编译命令、行截断检查和最小冒烟测试。
- 不使用经验平滑、兜底路径或后处理补丁掩盖非连续或非光滑结果；物理量随时间/空间不光滑时优先查 bug。
