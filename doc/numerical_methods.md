# 数值方法

本文档按实现模块整理当前数值方法和验证重点。

## 动力学

正激波动力学：

- Source：`src/Dynamics/Dynamics_forward.f90`
- Shared helpers：`src/Dynamics/dynamics_common.f90`
- 支持 ISM、wind `k=2` backend、density jump 和 energy injection paths。

反激波动力学：

- Source：`src/Dynamics/Dynamics_reverse.f90`
- 返回 RS shell state、`gamma34`、`U3/V3`、turbulent 和 ordered magnetic field diagnostics。
- Magnetized jump 遵循当前 RS baseline 契约；`sigma -> 0` 是必须检查的回归极限。

## 电子能量网格

电子求解器使用 log-gamma 网格。公共 helper 位于：

- `src/Electron/electron_common.f90`
- `src/Electron/electron_injection_profiles.f90`
- `src/Electron/adaptive_resampling_mod.f90`

关键数值问题：

- injection normalization
- high-energy cutoff
- cooling face speeds
- conservative or implicit transport update
- stable synchrotron/SSA seed recomputation

## 1D 电子求解器

### `fullhide_1d`

默认 public baseline：

- Source：`src/Electron/FS_electron_fullhide_1d.f90`
- Transport helper：`src/Electron/electron_transport_common.f90`
- Cooling：`src/Electron/electron_cooling_kernel.f90`
- Radiation：`src/Electron/electron_radiation_kernel.f90`

优势：

- stiff cooling 和高密度环境下稳定。
- 足够快，适合拟合。
- 当前 public comparison 的默认基线。

验证：

- 编译 `FS_electron_fullhide_1d`。
- 对 source closure 执行 `-Wline-truncation`。
- 跑 `tests/polarization_smoke.py` 或相关 electron smoke。
- 性能改动需要跑 benchmark path。

### `weno5_1d`

高阶电子谱解析求解器：

- Source：`src/Electron/FS_electron_weno5_1d.f90`
- 适合解析谱形，不一定是默认拟合路径。

### `charint_1d`, `slc1_1d`, `t2g1_1d`

这些是保留的比较、诊断和算法实验路径。若触及相关代码，必须保持可编译。

## 2D 电子求解器

2D path 包含 energy 和 chi-resolved electron transport：

- `src/Electron/FS_electron_fullhide_2d.f90`
- `src/Electron/FS_electron_charint_2d.f90`
- `src/Electron/electron_transport_2d_kernel.f90`
- `src/Electron/electron_seed_history_kernel.f90`

2D electron transport 不等于 2D hadronic transport。Hadronic state 在定义新契约前仍保持 shell-level。

验证：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module FS_electron_fullhide_2d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/fullhide_2d_smoke_bench.py'
```

## 冷却组装

Cooling kernel：

- `src/Electron/electron_cooling_kernel.f90`

职责：

- synchrotron cooling
- IC / Compton auxiliary terms
- SSA cooling
- Nakar/Fan-style Y paths
- 2D chi cells 的 batch cooling

Cooling kernel 把 photon-field preparation 与最终 cooling assembly 分开，避免共享 seed 时重复计算 IC auxiliary terms。

## 辐射积分

Synchrotron：

- `src/Electron/electron_radiation_kernel.f90`
- `src/Radiation/radiation_common.f90`

SSC：

- `src/Radiation/SSC_spec.f90`

Gamma-gamma：

- `src/Radiation/Annihilation.f90`

Reverse seed：

- `src/Radiation/Seed_reverse.f90`

数值注意事项：

- Fixed-grid synchrotron path 是 public fast path。
- Adaptive synchrotron integration 是 diagnostic path。
- SSA transfer 必须在频率 cell 间保持连续。
- 不添加没有物理边界支撑的 floor 或 guard terms。

## 强子核

Fortran 源文件：

- `src/Hadronic/FS_hadronic_1d.f90`
- `src/Hadronic/FS_hadronic_reverse_1d.f90`
- `src/Hadronic/hadronic_transport_kernel.f90`
- `src/Hadronic/hadronic_radiation_kernel.f90`
- `src/Hadronic/hadronic_interaction_kernel.f90`
- `src/Hadronic/hadronic_pair_production_kernel.f90`
- `src/Hadronic/hadronic_pp_kernel.f90`
- `src/Hadronic/hadronic_bethe_heitler_kernel.f90`
- `src/Hadronic/hadronic_hadronic_ic_kernel.f90`
- `src/Hadronic/hadronic_species_transport_kernel.f90`
- `src/Hadronic/hadronic_acceleration_kernel.f90`
- `src/Hadronic/hadronic_secondary_radiation_kernel.f90`
- `src/Hadronic/hadronic_decay_kernel.f90`
- `src/Hadronic/hadronic_pair_cascade_kernel.f90`

验证应按过程拆分：

- proton synch
- p-gamma
- Bethe-Heitler
- pp
- hadronic IC
- species transport
- secondary radiation
- pair production/cascade

## 观测者插值

插值源文件：

- `src/Interpolation/SED_interpolation.f90`
- `src/Interpolation/SED_interpolation_structured.f90`
- `src/Interpolation/interpolation_common.f90`

Projection 对本地 shell radiation 做 observer time/frequency grid 上的投影。Projection-layer fix 不能用来掩盖 dynamics 或 transport bug。

## 行截断检查

Fortran source closure 必须从 `/tmp` 执行，并指定临时 module 目录。`FS_electron_fullhide_1d` 示例：

```bash
rtk bash -lc 'source ~/.wsl_env && rm -rf /tmp/asgard_linecheck && mkdir -p /tmp/asgard_linecheck && cd /tmp && gfortran -cpp -fopenmp -Wall -Wline-truncation -fsyntax-only -J /tmp/asgard_linecheck -I /tmp/asgard_linecheck "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Constants.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Dynamics/dynamics_common.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Radiation/radiation_common.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Radiation/synchrotron_polarization_kernel.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/electron_transport_common.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/adaptive_resampling_mod.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/electron_radiation_kernel.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/electron_injection_profiles.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/electron_common.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/electron_cooling_kernel.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/electron_seed_history_kernel.f90" "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow/src/Electron/FS_electron_fullhide_1d.f90"'
```

不要从仓库根目录跑这个检查；仓库根目录里的旧 `.mod` 文件可能被误读。

## 性能说明

- `fullhide_1d` 是默认拟合路径。
- Cold solve time 和 cached query time 必须分开报告。
- 小网格可能暴露大网格 benchmark 隐藏的串行热点。
- Benchmark scratch files 放在 `/tmp`；只有可复现的文档 artifact 才应提交。
