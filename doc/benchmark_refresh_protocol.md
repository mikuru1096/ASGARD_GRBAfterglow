# Benchmark Refresh Protocol

本文固定 ASGARD benchmark 重新生成的最低基线。目标是让每张图、每个 CSV 和每个 comparison report 都能回溯到代码状态、命令、build 状态和物理验收口径。

## Required Metadata

每次刷新 benchmark 前后都记录：

- Git HEAD: `git rev-parse --short HEAD`
- 工作树状态: `git status --short --branch`
- Tracked diff: `git diff --stat`
- 命令原文，包括 `rtk bash -lc 'source ~/.wsl_env && cd "...ASGARD_GRBAfterglow" && ...'`
- 是否重编译 Fortran extension；若 benchmark 依赖 Fortran 改动，记录具体 `build_extensions.py --module ... --force` 命令。
- 输出路径；图像和 CSV 必须由脚本复现，不手工编辑。
- 验收口径：物理平滑性、过程开关、comparison backend 的角色，以及不作为失败条件的模型差异。

## Standard Commands

RS / Vegas comparison:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --extra compare python tests/vegas_afterglow_comparison.py --scenario baseline --only reverse_shock_lc reverse_shock_thermal'
```

Magnetized RS sigma scan:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/magnetized_rs_sigma_benchmark.py --medium ism --mode quick'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/magnetized_rs_sigma_benchmark.py --medium wind --mode quick'
```

Lan 2023 polarization overlay:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/polarization_literature_overlay.py'
```

Hadronic pγ report:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/hadronic_pgamma_benchmark_report.py'
```

## Build Gate

文档或 plotting-only 脚本改动只需：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git diff --check'
```

Fortran 或物理路径改动必须按受影响范围执行：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module FS_electron_fullhide_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module FS_hadronic_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module Dynamics_reverse --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_reverse_kernel --force'
```

Fortran line-truncation checks use WSL gfortran with `-Wline-truncation` on the touched source group. Do not treat a benchmark as refreshed if the relevant extension was changed but not rebuilt.

## Physical Acceptance

- RS plots must state that ASGARD uses local `gamma34` injection plus explicit `U3/V3` thermal state. VegasAfterglow is a comparison backend, not the target answer.
- If ASGARD and Vegas differ because Vegas uses averaged `Gamma_th` while ASGARD uses local shock-front injection, record that as a model difference.
- Magnetized RS sigma scans must check smoothness of `B3`, `gamma34`, `U3/V3`, `nu_m`, `nu_c`, `nu_a`, and the `sigma -> 0` limit.
- Polarization overlays must report peak time and peak amplitude separately. Do not correct a timing mismatch with empirical time shifts or smoothing.
- Hadronic reports must name enabled processes and whether the path is formal 1D shell transport or a documented boundary.

## Artifact Policy

- Commit scripts and documentation needed to reproduce the benchmark.
- Commit benchmark figures only when they are intended documentation assets and are reproducible from the recorded command.
- Do not commit `.buildcache/`, temporary debug scripts, partial CSVs, failed placeholder plots, or one-off notebook exports.
- Untracked benchmark output directories can remain local until a refresh is intentionally promoted.
