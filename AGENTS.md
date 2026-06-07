- 用中文回复。
- Shell 命令通过 `wsl -e bash -c` 在 WSL Ubuntu 内执行，非交互式 shell 需 `source ~/.wsl_env` 初始化 PATH 和代理。
- 默认环境 `WSL + uv`，使用`rtk`命令降低token使用量。
- 当前工作基线见 `PLAN.md`；代码结构见 `doc/code_overview.md` 与 `doc/source_tree.md`。
- 禁止加数值保护。
- 物理量随时间演化若不连续不平滑，优先怀疑数值或物理 bug。
- Fortran 代码保持 idiomatic Fortran 风格。
- Fortran 重要改动后必须跑：受影响的 `build_extensions.py --force`、`-Wline-truncation` 检查、最小相关 smoke test。
- 不要提交 `.buildcache/`、临时 debug 脚本、失败占位图。

## Build Commands

### 完整 f2py 编译
```bash
rtk bash -lc "source ~/.wsl_env && cd \"/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow\" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force"
```
**声明块压缩规则**: 同类型声明合并到一行，语义相近的量分组。B 类子程序声明块 ≤15 行。禁止每行只声明一个变量。

- **Reverse-shock hadronic**: light backend `hadronic_reverse_1d` covers RS proton injection/transport + proton synchrotron; when RS hadronic full-chain flags are enabled, the runtime reuses the formal 1D hadronic kernels for RS pγ/BH/pp/secondary/cascade coupling with RS seed photons and RS shell targets.
- **Pair cascade**: `pair_cascade_iterations > 1` now uses a shell-sequence time-dependent γγ pair/synch cascade path; the legacy single-shell iterative kernel remains for low-level diagnostics. IC-mediated electromagnetic cascade boundary is tracked in `TODO.md`.
- **Reverse-shock thermal/magnetized baseline**: RS 注入能标使用 shock-front `gamma34`；区域 3 turbulent field 和 post-crossing 热演化使用显式 `U3/V3` thermal state；可选 upstream `sigma` 引入 MHD jump 压缩比和有序磁场，返回的 `B3` 是 turbulent + ordered total field，`sigma -> 0` 必须保持当前 unmagnetized baseline。VegasAfterglow 作为 jump-condition 来源和 comparison backend，不作为光变目标或物理基准。
- **Vegas benchmark refresh**: `tests/vegas_afterglow_comparison.py --scenario baseline` 生成 baseline 全量 comparison figures；RS 窄刷新可用 `--only reverse_shock_lc reverse_shock_thermal`。图像是 artifact，必须能由脚本复现。
- **Benchmark refresh protocol**: 重新生成 benchmark 前后记录 git HEAD、`git status --short --branch`、`git diff --stat`、完整命令、受影响 Fortran build 状态、输出路径和物理验收口径；详见 `doc/benchmark_refresh_protocol.md`。
- **2D/chi-resolved hadronic decision**: 当前不实现 2D/χ hadronic transport；正式 hadronic path 保持 1D shell 契约，直到 χ-local photon field、hadron density、secondary feedback 和 observer projection 的物理契约完成。
- **Polarization timing diagnostic**: Lan 2023 overlay 的峰值幅度已匹配，峰时偏早主要指向 dynamics/jet-evolution benchmark；禁止在 polarization projection 层用经验时间因子或 smoothing 修正。
- **Public/backend limits**: Jet spreading、自定义 `Medium` kernel dispatch、wind `k != 2`、thermal electrons outside `fullhide_1d` 是明确未支持边界；详见 `doc/public_backend_limits.md`。
- **TODO index**: 当前未完成项集中维护在根目录 `TODO.md`；不要在其他文档新增分散待办列表。

## AM3 / ASGARD Coexistence

- AM3 本地参考: `/mnt/c/Users/jia/Documents/New project/_external/am3_reference` (HEAD: `7aba970b`); WSL home mirror `~/projects/_external/am3_reference` 同 HEAD。
- Python 只做 orchestration、wrapping、benchmark、API glue
