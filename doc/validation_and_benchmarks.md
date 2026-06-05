# 验证与基准

本文档固定 ASGARD 的验证层级。核心原则：验证应回答明确假设或支持决策，不做低信息增益的穷举。

## 验证层级

### 文档或 Python-only 轻量改动

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git diff --check'
```

如果改动涉及 Markdown 链接或示例命令，还应人工确认路径和命令仍对应当前工作树。

### Fortran 语法和行截断

Fortran 重要改动必须跑：

- 受影响 `build_extensions.py --module ... --force`
- `gfortran -Wall -Wline-truncation -fsyntax-only`
- 最小相关 smoke test

Line check 必须从 `/tmp` 执行，并指定临时 `-J` / `-I` 目录。

### 冒烟测试

Forward/electron：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/readme_smoke_bench.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/fullhide_2d_smoke_bench.py'
```

Reverse shock：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/reverse_shock_smoke.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/hadronic_reverse_shock_smoke.py'
```

Hadronic：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/hadronic_1d_smoke.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/hadronic_pair_cascade_smoke.py'
```

Polarization：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/polarization_smoke.py'
```

## 按区域划分的构建门槛

Electron 1D：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
```

Electron 2D：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_transport_2d --force'
```

反激波电子：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_reverse_kernel --force'
```

Hadronic：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_forward_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_reverse_1d --force'
```

Dynamics：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module Dynamics_forward --module Dynamics_reverse --force'
```

## 基准刷新

Benchmark refresh 协议固定在 `doc/benchmark_refresh_protocol.md`。

最少元数据：

- 刷新前后 Git HEAD。
- `git status --short --branch`。
- `git diff --stat`。
- 完整命令文本。
- 是否重建 Fortran extension。
- 输出路径。
- 物理验收口径。

RS / Vegas comparison：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --extra compare python tests/vegas_afterglow_comparison.py --scenario baseline --only reverse_shock_lc reverse_shock_thermal'
```

Speed profile：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --extra compare python tests/vegas_afterglow_comparison.py --scenario baseline --only speed_compare'
```

Magnetized RS sigma scan：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/magnetized_rs_sigma_benchmark.py --medium ism --mode quick'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/magnetized_rs_sigma_benchmark.py --medium wind --mode quick'
```

Runtime breakdown：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python scripts/benchmarks/runtime_breakdown_benchmark.py'
```

## 产物策略

可以提交：

- 作为文档资产的可复现 benchmark figures。
- 由已提交脚本生成且物理验收通过的 CSV summaries。
- 复现图像所需的脚本和文档。

不要提交：

- `.buildcache/`
- `__pycache__/`
- 临时 debug 脚本
- 失败占位图
- 不完整 CSV
- 一次性 notebook export
- 本地编译 `.so`, `.o`, `.mod`

## 物理验收

Forward-shock：

- Light curves 应平滑，除非物理 density jump 或 injection event 产生已记录特征。
- Characteristic frequencies 应连续演化。
- SSA breaks 不应出现 grid-cell discontinuity。

Reverse-shock：

- `sigma -> 0` 必须恢复 unmagnetized baseline。
- `B3`, `gamma34`, `U3/V3`, `nu_m`, `nu_c`, `nu_a` 应平滑。
- VegasAfterglow 是 comparison backend，不是目标真值。

Hadronic：

- 启用过程必须明确列出。
- Formal path 在定义 chi-local contract 前保持 shell-level 1D。
- Pair cascade 必须保持 gamma-gamma pair/synch branch 的解释，不扩展成 IC-mediated cascade。

Polarization：

- Peak amplitude 和 peak time 是分开的诊断。
- 不用经验 time shift 或 smoothing 修正 timing mismatch。

## 失败处理

验证失败时：

1. 记录准确命令和 artifact。
2. 判断失败属于 compile、API、numerical 还是 physical。
3. 检查最小 source closure。
4. 除非失败发生在真实系统边界，不添加 fallback 或 guard code。
5. 只重跑能覆盖修复的最窄检查。
