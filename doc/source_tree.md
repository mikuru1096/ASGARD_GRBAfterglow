# ASGARD Source Tree

本文档按当前工作树实际文件生成，反映 `2026-04-17` 这次检查时的源码组织，而不是旧版 `AGENTS.md` 中已经过期的门面结构。

## Top-Level Tree

```text
ASGARD_GRBAfterglow/
├── ASGARD/
│   ├── __init__.py
│   ├── api.py
│   └── units.py
├── doc/
│   ├── source_tree.md
│   └── call_chain.md
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
│   │   ├── FS_electron_fullhide.f90
│   │   ├── FS_electron_weno5.f90
│   │   ├── FS_electron_t2g1.f90
│   │   ├── FS_electron_slc1.f90
│   │   └── FS_electron_charint.f90
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
│   ├── ic_doc_plots.py
│   ├── readme_smoke_bench.py
│   ├── sed_electron_compare.py
│   ├── slc1_vs_fullhide_bench.py
│   ├── vegas_afterglow_comparison.py
│   └── wind_vs_fullhide_bench.py
├── asgard_coupling.py
├── asgard_fit.py
├── asgard_models.py
├── asgard_observables.py
├── asgard_paths.py
├── asgard_plot.py
├── asgard_postprocess.py
├── asgard_presets.py
├── asgard_runtime.py
├── asgard_setup.py
├── asgard_ssc.py
├── asgard_state.py
├── build_extensions.py
├── cal_chi2_light_curve.py
├── cal_chi2_spectrum.py
├── lc_spec_demo.py
└── README.md
```

## Directory Roles

- `ASGARD/`
  - 面向外部用户的高层 API。
  - `api.py` 负责把 `Jet/Medium/Observer` 等对象模型映射到 `FitConfig`，再调 `asgard_state.py`。

- `src/`
  - 真正的数值核与 f2py 编译产物所在目录。
  - `__init__.py` 和各子目录的 `__init__.py` 负责把 `.pyd/.so` 暴露成 Python 可导入接口。

- `tests/`
  - 当前主要承担 benchmark、图像生成、对比检查和 README/doc 配图脚本。
  - 不是传统单元测试仓库布局，更接近“诊断脚本 + 产图脚本”。

- 顶层 `asgard_*.py`
  - 当前实际的 Python 编排层。
  - `asgard_fit.py` 负责拟合问题编译与 loglike。
  - `asgard_state.py` 负责主物理链装配。
  - `asgard_runtime.py` 负责具体挑选并调用哪个 Fortran 绑定。
  - `asgard_postprocess.py` 负责观测者系插值和卡方。
  - `asgard_ssc.py` 与 `asgard_coupling.py` 负责 SSC 辅助网格和跨区 seed 场。

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

1. `ASGARD/api.py`
2. `asgard_fit.py`
3. `asgard_state.py`
4. `asgard_runtime.py`
5. `build_extensions.py`
6. `src/Electron/electron_common.f90`
7. `src/Radiation/SSC_spec.f90`
8. `src/Interpolation/SED_interpolation.f90`
