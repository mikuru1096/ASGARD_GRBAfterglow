# 命令行与脚本

ASGARD 当前没有独立的物理模型命令行前端；正式入口是 Python public API、Fortran build 命令和少量测试入口。本文按新手工作流列出应使用的命令，避免把构建、作图、benchmark 和网页发布混在一起。

## 1. 构建默认数值核

首次运行正向激波 `fullhide_1d` 教程前，先构建 Fortran 扩展：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
```

该命令只构建指定模块。修改其他 Fortran 数值核后，应构建受影响模块，而不是用无关模块代替验证。

## 2. 构建网页文档

本地构建私有网页文档：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --with "mkdocs<2" --with "mkdocs-material>=9.5" mkdocs build --strict'
```

输出目录是 `site/`，该目录不进入版本库。合作者可见的 GitHub Pages 发布流程见 `doc/web_docs.md`。

## 3. 文档文本检查

跨 Windows/WSL 修改文档后，先检查编码和行尾：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
```

再检查空白错误：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git diff --check'
```

## 4. Benchmark 与正式产物

当前不保留一次性 compare-plot CLI。Benchmark 图像和 CSV 只有在先明确假设、决策价值和物理验收口径后，才新建正式 `tests/` 入口或任务明确的 `/tmp` 临时脚本生成。刷新前后记录：

- git HEAD 和 `git status --short --branch`。
- 完整命令。
- 受影响 Fortran build 状态。
- 输出路径。
- 物理验收口径。

当前保留的 q-shell 诊断 benchmark 入口：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/benchmark_theta_j_multiples_magnetic_decay.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/benchmark_skymap_centroid_motion.py'
```

当前 prompt snapshot formal plotting 入口：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python prompt/run_formal_results.py'
```

这些入口都是诊断或图件生成，不是通用命令行物理模型前端。

## 5. 何时不用命令行

如果任务是交互式建模、拟合 likelihood、检查内部物理量或组合多个观测数据集，优先写一个短 Python 脚本调用 `Model` 和 `Fitter`。命令行适合构建、benchmark 刷新和构建文档，不适合承载复杂研究逻辑。
