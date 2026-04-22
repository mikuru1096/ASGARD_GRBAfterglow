# ASGARD Source Tree

本文档按当前工作树实际文件生成，反映 `2026-04-17` 这次检查时的源码组织，而不是旧版 `AGENTS.md` 中已经过期的门面结构。

## Top-Level Tree

```text
ASGARD_GRBAfterglow/
├── ASGARD/
│   ├── __init__.py
│   ├── api.py
│   ├── api_adaptive.py
│   ├── api_fit.py
│   ├── api_model.py
│   ├── api_observe.py
│   └── units.py
├── asgard_core/
│   ├── __init__.py
│   ├── asgard_config.py
│   ├── asgard_coupling.py
│   ├── asgard_fit.py
│   ├── asgard_models.py
│   ├── asgard_observables.py
│   ├── asgard_paths.py
│   ├── asgard_physics_utils.py
│   ├── asgard_plot.py
│   ├── asgard_postprocess.py
│   ├── asgard_presets.py
│   ├── asgard_runtime.py
│   ├── asgard_setup.py
│   ├── asgard_ssc.py
│   ├── asgard_state.py
│   ├── asgard_types.py
│   └── support/
│       ├── chi2_light_curve.py
│       ├── chi2_spectrum.py
│       ├── chi2_utils.py
│       ├── data_light_curve/
│       ├── data_spectrum/
│       ├── extinction_cur.py
│       └── extinc.py
├── doc/
│   ├── call_chain.md
│   ├── plan.md
│   └── source_tree.md
├── scripts/
│   └── fitting/
│       ├── mcmc_fit.py
│       ├── mpi_run.sh
│       ├── multinest_fit.py
│       └── multinest_marginals_corner.py
├── src/
│   ├── __init__.py
│   ├── Constants.f90
│   ├── Dynamics/
│   │   ├── __init__.py
│   │   ├── dynamics_common.f90
│   │   ├── Dynamics_forward.f90
│   │   └── Dynamics_reverse.f90
│   ├── Electron/
│   │   ├── __init__.py
│   │   ├── adaptive_resampling_mod.f90
│   │   ├── calling_modules.f90
│   │   ├── electron_common.f90
│   │   ├── FS_electron_fullhide_1d.f90
│   │   ├── FS_electron_weno5_1d.f90
│   │   ├── FS_electron_t2g1_1d.f90
│   │   ├── FS_electron_slc1_1d.f90
│   │   └── FS_electron_charint_1d.f90
│   ├── Interpolation/
│   │   ├── __init__.py
│   │   ├── interpolation_common.f90
│   │   ├── SED_interpolation.f90
│   │   └── SED_interpolation_structured.f90
│   └── Radiation/
│       ├── __init__.py
│       ├── Cal_ebl.py
│       ├── radiation_common.f90
│       ├── SSC_spec.f90
│       ├── Annihilation.f90
│       └── Seed_reverse.f90
├── tests/
│   ├── doc_figures.py
│   ├── fullhide_2d_medium_diag.py
│   ├── fullhide_2d_smoke_bench.py
│   ├── ic_doc_plots.py
│   ├── readme_smoke_bench.py
│   ├── sed_electron_compare.py
│   ├── slc1_vs_fullhide_bench.py
│   ├── vegas_afterglow_comparison.py
│   └── wind_vs_fullhide_bench.py
├── build_extensions.py
├── install.py
├── install.ps1
├── install.sh
├── lc_spec_demo.py
├── pyproject.toml
├── setup.py
└── README.md
```

## Directory Roles

- `ASGARD/`
  - 面向外部用户的高层 API。
  - `api.py` 现在是薄门面，负责重导出公开 API。
  - `api_model.py` 放对象模型与 `Model`。
  - `api_fit.py` 放 `Param/Fitter/FitResult`。
  - `api_adaptive.py` 放自适应采样、缓存和插值 helper。
  - `api_observe.py` 负责把 `Jet/Medium/Observer` 等对象模型映射到 `FitConfig`，再调 `asgard_core/asgard_state.py`。

- `asgard_core/`
  - 放内部 Python 编排层、运行时绑定、物理辅助和后处理模块。
  - `support/` 子目录放 `chi^2`、extinction 和观测表辅助函数。

- `src/`
  - 真正的数值核与 f2py 编译产物所在目录。
  - `__init__.py` 和各子目录的 `__init__.py` 负责把 `.pyd/.so` 暴露成 Python 可导入接口。

- `tests/`
  - 当前主要承担 benchmark、图像生成、对比检查和 README/doc 配图脚本。
  - 不是传统单元测试仓库布局，更接近“诊断脚本 + 产图脚本”。

- `scripts/fitting/`
  - 放采样、MultiNest 后处理和批量运行脚本。
  - 不再占用根目录。

## Build Graph

`build_extensions.py` 当前登记的 f2py 编译目标如下：

```text
Constants
Dynamics_reverse
Dynamics_forward
FS_electron_weno5
FS_electron_slc1
FS_electron_charint
FS_electron_fullhide
FS_electron_t2g1
SED_interpolation
SED_interpolation_structured
Annihilation
Seed_reverse
SSC_spec
```

对应源码依赖主干：

- `Constants` 由 `src/Constants.f90` 单独构建。
- `Dynamics_*` 依赖 `Constants.f90`，其中 forward 额外依赖 `dynamics_common.f90`。
- 所有 `FS_electron_*` 依赖：
  - `Constants.f90`
  - `adaptive_resampling_mod.f90`
  - `electron_common.f90`
  - `calling_modules.f90`
- `SED_interpolation*` 依赖：
  - `Constants.f90`
  - `interpolation_common.f90`
- `SSC_spec / Annihilation / Seed_reverse` 依赖：
  - `Constants.f90`
  - `radiation_common.f90`

## Current Architectural Center

按当前代码实际使用频率，最值得优先看的核心文件是：

1. `ASGARD/api_model.py`
2. `asgard_core/asgard_fit.py`
3. `asgard_core/asgard_state.py`
4. `asgard_core/asgard_runtime.py`
5. `build_extensions.py`
6. `src/Electron/electron_common.f90`
7. `src/Radiation/SSC_spec.f90`
8. `src/Interpolation/SED_interpolation.f90`
