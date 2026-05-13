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
rtk bash -lc "source ~/.wsl_env && cd \"/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow\" && TMPDIR=/tmp uv run python build_extensions.py --module FS_electron_fullhide_1d --force"
```
**声明块压缩规则**: 同类型声明合并到一行，语义相近的量分组。B 类子程序声明块 ≤15 行。禁止每行只声明一个变量。

- **Reverse-shock hadronic**: light backend `FS_hadronic_reverse_1d` covers RS proton injection/transport + proton synchrotron; when RS hadronic full-chain flags are enabled, the runtime reuses the formal 1D hadronic kernels for RS pγ/BH/pp/secondary/cascade coupling with RS seed photons and RS shell targets.
- **Pair cascade**: `pair_cascade_iterations > 1` now uses a shell-sequence time-dependent γγ pair/synch cascade path; the legacy single-shell iterative kernel remains for low-level diagnostics.
- **Not yet implemented**: 2D/chi-resolved hadronic transport; inverse-Compton-mediated electromagnetic pair cascade beyond the current γγ pair/synch contract.

## AM3 / ASGARD Coexistence

- AM3 本地参考: `/mnt/c/Users/jia/Documents/New project/_external/am3_reference` (HEAD: `7aba970b`); WSL home mirror `~/projects/_external/am3_reference` 同 HEAD。
- Python 只做 orchestration、wrapping、benchmark、API glue
