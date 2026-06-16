# ASGARD 源码树

本文档是精简源码树索引。文档阅读入口见 `doc/index.md`；代码总览见 `doc/code_overview.md`。

## Python 运行层

- `asgard_core/`：公开 Python API 和运行时 glue，包括 `api_model.py`, `api_observe.py`, `api_fit.py`, `api_adaptive.py`。`api_model.py` 承担 `Model` 查询调度、`Model -> RuntimeConfig` 适配、direct/patch solve 入口和 details 打包；`api_observe.py` 保留 `observe` / `run_fit` 兼容入口、sky image、polarization 和观测数据集 helpers；`api_adaptive.py` 保留 shared observe packing、数组签名缓存和曝光时间平均工具。
- `asgard_core/`：内部编排与物理耦合，包括 `asgard_config.py`, `asgard_runtime.py`, `asgard_state.py`, `asgard_setup.py`, `asgard_types.py`。
- `asgard_core/structured_jet_kernel.py`：结构化喷流 Fortran backend 薄中间层。
- `asgard_core/hadronic_*.py`：hadronic Python wrappers 和 glue，只做 orchestration；最终微物理写入 `src/Hadronic/`。
- RS hadronic light proton transport + synchrotron wrapper 位于 `asgard_core/asgard_runtime.py`。
- `asgard_core/hadronic_cascade.py`：shell-sequence time-dependent γγ pair/synch cascade。

## 文档

- `doc/index.md`：文档总入口。
- `TODO.md`：唯一 TODO / 未完成项入口。
- `doc/installation.md`：环境、安装和 native extension 构建。
- `doc/public_api.md`：public API 当前契约。
- `doc/project_physics_design.md`：全项目物理设计总纲。
- `doc/project_algorithm_design.md`：全项目算法设计总纲。
- `doc/physics_model.md`：已实现物理模块和边界。
- `doc/numerical_methods.md`：数值核、求解器和 line-truncation 检查。
- `doc/joint_secondary_feedback_physics.md`：含时 BH / 二级 e± / 光子场 joint feedback 的物理方程、源汇预算和支持边界。
- `doc/joint_secondary_feedback_algorithm.md`：joint feedback 的 public switch、状态机、数组契约、函数入口、测试和 benchmark。
- `doc/validation_and_benchmarks.md`：build gate、smoke tests、benchmark refresh 和 artifact policy。
- `doc/developer_guide.md`：开发工作流、提交前检查和 review checklist。
- `doc/call_chain.md`：public call chain 到 runtime/Fortran kernel 的路径。
- `doc/code_overview.md`：代码结构、运行主链和关键边界。
- `doc/electron_solver_algorithms.md`：电子输运算法说明。
- `doc/hadronic_pgamma_notes.md`：p-gamma 微物理和基准说明。
- `doc/am3_migration_plan.md`：AM3 共存、迁移和引用边界。
- `doc/hadronic_chi_transport_decision.md`：2D / \(\chi\) 分辨 hadronic transport 的当前决策边界。
- `doc/pair_cascade_extension_boundary.md`：IC-mediated electromagnetic cascade 的扩展边界。
- `doc/fullhide2d_pwn_cr_transport.md`：`fullhide2d_transport_model="pwn_cr_v1"` 的物理契约。
- `doc/public_backend_limits.md`：public API/backend 未支持和部分支持边界。
- `doc/web_docs.md` 与 `mkdocs.yml`：网页文档发布配置和导航入口；新增或改名文档必须同步 `nav` 并跑 strict build。

## Fortran 运行层

- `src/Constants.f90`
- `src/Dynamics/`：`dynamics_common.f90`, `Dynamics_forward.f90`, `Dynamics_reverse.f90`
- `src/Electron/`：
  - 1D/2D entries：`electron_forward_fullhide_1d.f90`, `electron_forward_transport_2d.f90`, `electron_forward_charint_1d.f90`, `electron_forward_slc1_1d.f90`, `electron_forward_t2g1_1d.f90`, `electron_forward_weno5_1d.f90`；`electron_forward_charint_2d` extension 由 `electron_forward_transport_2d.f90` 的 `fs_electron_transport_2d_core` 构建
  - Kernels：`electron_common.f90`, `electron_radiation_kernel.f90`, `electron_cooling_kernel.f90`, `electron_seed_history_kernel.f90`, `electron_transport_2d_kernel.f90`, `electron_injection_profiles.f90`, `electron_transport_common.f90`, `electron_reverse_kernel.f90`, `adaptive_resampling_mod.f90`
- `src/Radiation/`：`radiation_common.f90`, `radiation_ssc_spectrum.f90`, `radiation_gamma_gamma_absorption.f90`, `synchrotron_polarization_kernel.f90`, `quantum_synchrotron_kernel.f90`
- `src/Hadronic/`：
  - Entries：`hadronic_forward_1d.f90`, `hadronic_reverse_1d.f90`
  - Kernels：`hadronic_common.f90`, `hadronic_transport_kernel.f90`, `hadronic_transport_remap_kernel.f90`, `hadronic_radiation_kernel.f90`, `hadronic_interaction_kernel.f90`, `hadronic_pgamma_hummer_1d.f90`, `hadronic_decay_kernel.f90`, `hadronic_pair_production_kernel.f90`, `hadronic_pair_cascade_kernel.f90`, `hadronic_pp_kernel.f90`, `hadronic_pp_models_kernel.f90`, `hadronic_bethe_heitler_kernel.f90`, `hadronic_hadronic_ic_kernel.f90`, `hadronic_species_transport_kernel.f90`, `hadronic_acceleration_kernel.f90`, `hadronic_secondary_radiation_kernel.f90`
- `src/Structured/`：`structured_jet_1d.f90` 聚合结构化喷流 theta/theta-phi 网格调度，复用现有 Fortran 动力学、电子、辐射、强子和 SED 插值核。
- `src/Interpolation/`：`SED_interpolation.f90`, `SED_interpolation_structured.f90`, `interpolation_common.f90`

## 构建入口

- `build_extensions.py`：f2py 编译入口。已登记 active hadronic extensions：`hadronic_forward_1d`, `hadronic_reverse_1d`。

## 测试与基准

- 基础 smoke：`tests/readme_smoke_bench.py`, `tests/reverse_shock_smoke.py`
- 2D electron：`tests/fullhide_2d_smoke_bench.py`
- Hadronic：`tests/hadronic_1d_smoke.py`, `tests/hadronic_species_transport_smoke.py`, `tests/hadronic_secondary_radiation_smoke.py`, `tests/hadronic_acceleration_smoke.py`, `tests/hadronic_bethe_heitler_smoke.py`, `tests/hadronic_hadronic_ic_smoke.py`, `tests/hadronic_pair_production_smoke.py`, `tests/hadronic_pp_smoke.py`
- Pair / RS：`tests/hadronic_pair_cascade_smoke.py`, `tests/hadronic_pair_branch_smoke.py`, `tests/hadronic_reverse_shock_smoke.py`
- AM3-derived process smoke：`tests/hadronic_am3_solver_smoke.py`
- Electron-photon joint feedback：`tests/electron_photon_coupling_config_smoke.py`, `tests/electron_photon_joint_smoke.py`, `tests/electron_photon_joint_secondary_feedback_smoke.py`, `tests/electron_photon_separated_regression_smoke.py`, `tests/electron_photon_ic_consistency_smoke.py`, `tests/hadronic_r_coordinate_smoke.py`

## 生成产物

- RS-specific artifacts：`compare_reverse_shock_lc.png`, `compare_reverse_shock_thermal_benchmark.png`。
