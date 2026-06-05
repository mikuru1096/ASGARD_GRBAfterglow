# ASGARD 源码树

本文档是精简源码树索引。代码总览见 `doc/code_overview.md`。

## Python 运行层

- `ASGARD/`：外部用户 API，包括 `api_model.py`, `api_observe.py`, `api_fit.py`, `api_adaptive.py`。
- `asgard_core/`：内部编排与物理耦合，包括 `asgard_config.py`, `asgard_runtime.py`, `asgard_state.py`, `asgard_setup.py`, `asgard_ssc.py`, `asgard_types.py`。
- `asgard_core/hadronic_*.py`：hadronic Python wrappers/benchmarks，只做 orchestration；最终微物理写入 `src/Hadronic/`。
- `asgard_core/hadronic_reverse.py`：RS hadronic light proton transport + synchrotron wrapper。
- `asgard_core/hadronic_cascade.py`：pair-production synch diagnostics 和 shell-sequence time-dependent γγ pair/synch cascade。

## 文档

- `doc/index.md`：文档总入口。
- `doc/installation.md`：环境、安装和 native extension 构建。
- `doc/user_guide.md`：常用 public API 工作流。
- `doc/public_api.md`：public API 当前契约。
- `doc/physics_model.md`：已实现物理模块和边界。
- `doc/numerical_methods.md`：数值核、求解器和 line-truncation 检查。
- `doc/validation_and_benchmarks.md`：build gate、smoke tests、benchmark refresh 和 artifact policy。
- `doc/developer_guide.md`：开发工作流、提交前检查和 review checklist。
- `doc/hadronic_chi_transport_decision.md`：2D / chi-resolved hadronic transport 的当前决策边界。
- `doc/pair_cascade_extension_boundary.md`：IC-mediated electromagnetic cascade 的扩展边界。
- `doc/polarization_timing_diagnostic.md`：Lan 2023 polarization 峰时偏移诊断。
- `doc/benchmark_refresh_protocol.md`：benchmark 重新生成的命令、build、artifact 和物理验收协议。
- `doc/public_backend_limits.md`：public API/backend 未支持和部分支持边界。

## Fortran 运行层

- `src/Constants.f90`
- `src/Dynamics/`：`dynamics_common.f90`, `Dynamics_forward.f90`, `Dynamics_reverse.f90`
- `src/Electron/`：
  - 1D/2D entries：`electron_forward_fullhide_1d.f90`, `electron_forward_transport_2d.f90`, `electron_forward_charint_1d.f90`, `electron_forward_slc1_1d.f90`, `electron_forward_t2g1_1d.f90`, `electron_forward_weno5_1d.f90`；`electron_forward_charint_2d` extension 由 `electron_forward_transport_2d.f90` 的 `fs_electron_transport_2d_core` 构建
  - Kernels：`electron_common.f90`, `electron_radiation_kernel.f90`, `electron_cooling_kernel.f90`, `electron_seed_history_kernel.f90`, `electron_transport_2d_kernel.f90`, `electron_injection_profiles.f90`, `electron_transport_common.f90`, `electron_reverse_kernel.f90`, `adaptive_resampling_mod.f90`
- `src/Radiation/`：`radiation_common.f90`, `radiation_ssc_spectrum.f90`, `radiation_gamma_gamma_absorption.f90`, `radiation_reverse_seed.f90`, `synchrotron_polarization_kernel.f90`, `quantum_synchrotron_kernel.f90`
- `src/Hadronic/`：
  - Entries：`hadronic_forward_1d.f90`, `hadronic_reverse_1d.f90`
  - Kernels：`hadronic_common.f90`, `hadronic_transport_kernel.f90`, `hadronic_radiation_kernel.f90`, `hadronic_interaction_kernel.f90`, `hadronic_decay_kernel.f90`, `hadronic_pair_production_kernel.f90`, `hadronic_pair_cascade_kernel.f90`, `hadronic_pp_kernel.f90`, `hadronic_pp_models_kernel.f90`, `hadronic_bethe_heitler_kernel.f90`, `hadronic_hadronic_ic_kernel.f90`, `hadronic_species_transport_kernel.f90`, `hadronic_acceleration_kernel.f90`, `hadronic_secondary_radiation_kernel.f90`
- `src/Interpolation/`：`SED_interpolation.f90`, `SED_interpolation_structured.f90`, `interpolation_common.f90`

## 构建入口

- `build_extensions.py`：f2py 编译入口。已登记 active hadronic extensions：`hadronic_forward_1d`, `hadronic_reverse_1d`。

## 测试与基准

- 基础 smoke：`tests/readme_smoke_bench.py`, `tests/reverse_shock_smoke.py`
- 2D electron：`tests/fullhide_2d_smoke_bench.py`, `tests/fullhide_2d_medium_diag.py`
- 比较与谱：`tests/vegas_afterglow_comparison.py`, `tests/sed_electron_compare.py`
- Hadronic：`tests/hadronic_1d_smoke.py`, `tests/hadronic_species_transport_smoke.py`, `tests/hadronic_secondary_radiation_smoke.py`, `tests/hadronic_acceleration_smoke.py`, `tests/hadronic_bethe_heitler_smoke.py`, `tests/hadronic_hadronic_ic_smoke.py`, `tests/hadronic_pair_production_smoke.py`, `tests/hadronic_pp_smoke.py`, `tests/hadronic_pgamma_benchmark_report.py`
- Pair / polarization / RS：`tests/hadronic_pair_cascade_smoke.py`, `tests/hadronic_pair_branch_smoke.py`, `tests/hadronic_reverse_shock_smoke.py`, `tests/polarization_smoke.py`, `tests/polarization_baseline_bench.py`
- AM3/reference comparisons：`tests/hadronic_am3_solver_smoke.py`, `tests/hadronic_am3_reference_compare_smoke.py`, `tests/hadronic_am3_acceleration_compare_smoke.py`, `tests/hadronic_am3_bethe_heitler_compare_smoke.py`, `tests/hadronic_pgamma_benchmark_compare.py`

## 生成产物

- `output/asgard_doc/vegas_afterglow_compare/`：由 `tests/vegas_afterglow_comparison.py` 生成的 Vegas/ASGARD comparison figures。
- RS refresh artifacts：`compare_reverse_shock_lc.png`, `compare_reverse_shock_thermal_benchmark.png`。
