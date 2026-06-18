- 用中文回复。
- Shell 命令通过 `wsl -e bash -c` 在 WSL Ubuntu 内执行，非交互式 shell 需 `source ~/.wsl_env` 初始化 PATH 和代理。
- 默认环境 `WSL + uv`，使用`rtk`命令降低token使用量。
- 当前工作基线见 `PLAN.md`；代码结构见 `doc/code_overview.md` 与 `doc/source_tree.md`。
- 禁止加数值保护。
- 物理量随时间演化若不连续不平滑，优先怀疑数值或物理 bug。
- Fortran 代码保持 idiomatic Fortran 风格。
- Fortran 重要改动后必须跑：受影响的 `build_extensions.py --force`、`-Wline-truncation` 检查、最小相关 smoke test。
- 不要提交 `.buildcache/`、临时 debug 脚本、失败占位图。
- 文本源文件统一 UTF-8 无 BOM + LF；跨 Windows/WSL 操作前运行 `rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'`。禁止用 PowerShell 重定向或未显式 UTF-8 的 `Set-Content` 写源码/文档；Python 文本 IO 必须显式 `encoding="utf-8"`；Fortran/gfortran 构建通过 `LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONUTF8=1 PYTHONIOENCODING=utf-8` 固定运行环境，不使用对 Fortran 前端无效的 C/C++ charset flags。

## Build Commands

### ApJ/AASTeX 论文工作流
- 论文主目录为 `paper/`；正式正文从 `paper/main.tex` 编译，引用库为 `paper/references.bib`，本地固定 AASTeX 文件为 `paper/aastex702.cls`、`paper/aastex7.cls` 和 `paper/aasjournalv7.1.bst`。
- 正式主图只由 `tests/generate_paper_figures.py` 生成；输出进入 `paper/figures/` 和 `paper/source_data/`。未被该入口复制到 `paper/source_data/` 的 `output/` 图件不作为论文主证据。
- 图件生成命令：`rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/generate_paper_figures.py'`。
- Windows TeX Live 编译命令：`rtk latexmk -cd -pdf -interaction=nonstopmode -halt-on-error -outdir=build paper/main.tex`。`paper/build/` 是编译产物，不提交。

### 完整 f2py 编译
```bash
rtk bash -lc "source ~/.wsl_env && cd \"/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow\" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force"
```
**声明块压缩规则**: 同类型声明合并到一行，语义相近的量分组。B 类子程序声明块 ≤15 行。禁止每行只声明一个变量。

- **Reverse-shock hadronic**: light backend `hadronic_reverse_1d` covers RS proton injection/transport + proton synchrotron; when RS hadronic full-chain flags are enabled, the runtime reuses the formal 1D hadronic kernels for RS pγ/BH/pp/secondary/cascade coupling with RS seed photons and RS shell targets.
- **Pair cascade**: `pair_cascade_iterations > 1` now uses a shell-sequence time-dependent γγ pair/synch cascade path; the legacy single-shell iterative kernel remains for low-level diagnostics. IC-mediated electromagnetic cascade boundary is tracked in `TODO.md`.
- **Reverse-shock thermal/magnetized baseline**: RS 注入能标使用 shock-front `gamma34`；区域 3 turbulent field 和 post-crossing 热演化使用显式 `U3/V3` thermal state；`E_iso` 是总 ejecta 能量，有限 upstream `sigma` 下 baryonic ejecta mass 为 `E_iso/[(1+sigma) Gamma0 c^2]`，MHD jump 同时给出压缩比、有序磁场和热比内能，有序磁场焓进入主 RS 动力学惯性，返回的 `B3` 是 turbulent + ordered total field；局部压缩比的 `sigma -> 0` 极限必须使用有限强度 hydrodynamic jump，不能使用 ultra-relativistic `4*gamma43+3` 近似。VegasAfterglow 作为 jump-condition comparison backend，不作为光变目标或全局动力学基准。
- **Vegas benchmark refresh**: 旧 Vegas Python comparison 入口已删除；需要重新建立时先说明物理假设和验收口径，再把可复用入口放入 `tests/`。
- **Benchmark refresh protocol**: 旧 `scripts/benchmarks/` Python 入口已删除。重新生成 benchmark 前必须先说明假设和决策价值；刷新前后记录 git HEAD、`git status --short --branch`、`git diff --stat`、完整命令、受影响 Fortran build 状态、输出路径和物理验收口径。正式可复用入口应进入 `tests/`，临时研究脚本放 `/tmp`。
- **2D/chi-resolved hadronic decision**: 当前不实现 2D/χ hadronic transport；正式 hadronic path 保持 1D shell 契约，直到 χ-local photon field、hadron density、secondary feedback 和 observer projection 的物理契约完成。
- **Polarization timing diagnostic**: Lan 2023 overlay 的峰值幅度已匹配，峰时偏早主要指向 dynamics/jet-evolution benchmark；禁止在 polarization projection 层用经验时间因子或 smoothing 修正。
- **Public/backend limits**: Jet spreading、自定义 `Medium` kernel dispatch、wind `k != 2`、thermal electrons outside `fullhide_1d` 是明确未支持边界；详见 `doc/public_backend_limits.md`。
- **TODO index**: 当前未完成项集中维护在根目录 `TODO.md`；不要在其他文档新增分散待办列表。

## AM3 / ASGARD Coexistence

- AM3 本地参考: `/mnt/c/Users/jia/Documents/New project/_external/am3_reference` (HEAD: `7aba970b`); WSL home mirror `~/projects/_external/am3_reference` 同 HEAD。
- Python 只做 orchestration、wrapping、benchmark、API glue
