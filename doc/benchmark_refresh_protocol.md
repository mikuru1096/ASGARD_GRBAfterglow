# 基准刷新协议

本文档固定 ASGARD benchmark 重新生成的最低基线。目标是让每张图、每个 CSV 和每个 comparison report 都能回溯到代码状态、命令、build 状态和物理验收口径。

## 必须记录的元数据

每次刷新 benchmark 前后都记录：

- Git HEAD：`git rev-parse --short HEAD`
- 工作树状态：`git status --short --branch`
- Tracked diff：`git diff --stat`
- 命令原文，包括 `rtk bash -lc 'source ~/.wsl_env && cd "...ASGARD_GRBAfterglow" && ...'`
- 是否重编译 Fortran extension；若 benchmark 依赖 Fortran 改动，记录具体 `build_extensions.py --module ... --force` 命令。
- 输出路径；图像和 CSV 必须由脚本复现，不手工编辑。
- 验收口径：物理平滑性、过程开关、comparison backend 的角色，以及不作为失败条件的模型差异。

## 标准命令

Vegas baseline full comparison：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --extra compare python tests/vegas_afterglow_comparison.py --scenario baseline'
```

该命令运行 `tests/vegas_afterglow_comparison.py` 的全部 builder。当前 baseline 会刷新 `output/asgard_doc/vegas_afterglow_compare/` 下由脚本生成的 comparison figures；`magnetic_decay_2d` 会额外输出 broadband spectrum 图。

RS / Vegas comparison：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --extra compare python tests/vegas_afterglow_comparison.py --scenario baseline --only reverse_shock_lc reverse_shock_thermal'
```

Magnetized RS sigma scan：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/magnetized_rs_sigma_benchmark.py --medium ism --mode quick'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/magnetized_rs_sigma_benchmark.py --medium wind --mode quick'
```

Lan 2023 polarization overlay：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/polarization_literature_overlay.py'
```

Hadronic pγ report：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/hadronic_pgamma_benchmark_report.py'
```

## 构建门槛

文档或 plotting-only 脚本改动只需：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git diff --check'
```

Fortran 或物理路径改动必须按受影响范围执行：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module Dynamics_forward --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module Dynamics_reverse --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_charint_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_charint_2d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_slc1_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_forward_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_reverse_kernel --force'
```

Fortran line-truncation 检查使用 WSL gfortran，并在 touched source group 上启用 `-Wline-truncation`。如果相关 extension 已改动但未重建，不应把 benchmark 视为已刷新。

## 物理验收

- RS 图必须说明 ASGARD 使用局部 `gamma34` 注入和显式 `U3/V3` thermal state。VegasAfterglow 是 comparison backend，不是目标答案。
- 如果 ASGARD 与 Vegas 的差异来自 Vegas 使用 averaged `Gamma_th` 而 ASGARD 使用局部 shock-front injection，应记录为模型差异。
- Magnetized RS sigma scan 必须检查 `B3`, `gamma34`, `U3/V3`, `nu_m`, `nu_c`, `nu_a` 的平滑性和 `sigma -> 0` 极限。
- Polarization overlay 必须分开报告 peak time 和 peak amplitude。不要用经验 time shift 或 smoothing 修正 timing mismatch。
- Hadronic report 必须写明启用过程，以及路径是 formal 1D shell transport 还是已记录边界。

## 产物策略

- 提交用于复现 benchmark 的脚本和文档。
- 只有当 benchmark figures 是预期文档资产并可由记录命令复现时才提交。
- 不提交 `.buildcache/`、临时 debug 脚本、partial CSV、失败占位图或一次性 notebook export。
- Untracked benchmark output directories 可以保留在本地，直到有意提升为文档资产。
