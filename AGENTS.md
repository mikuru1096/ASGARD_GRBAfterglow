- 用中文回复。
- Shell 命令通过 `wsl -e bash -c` 在 WSL Ubuntu 内执行，非交互式 shell 需 `source ~/.wsl_env` 初始化 PATH 和代理。
- 默认环境 `WSL + uv`，使用`rtk`命令降低token使用量。
- 禁止加数值保护。
- 物理量随时间演化若不连续不平滑，优先怀疑数值或物理 bug。
- Fortran 代码保持 idiomatic Fortran 风格。
- Fortran 重要改动后必须跑：受影响的 `build_extensions.py --force`、`-Wline-truncation` 检查、最小相关 smoke test。
- 不要提交 `.buildcache/`、临时 debug 脚本、失败占位图。

## Build Commands

### 完整 f2py 编译
```bash
wsl -e bash -c "source ~/.wsl_env && cd ~/projects/ASGARD_GRBAfterglow && uv run python build_extensions.py --module FS_electron_fullhide_1d --force"
```
**声明块压缩规则**: 同类型声明合并到一行，语义相近的量分组。B 类子程序声明块 ≤15 行。禁止每行只声明一个变量。

- **Reverse-shock hadronic**: formal 1D reverse-shock proton injection/transport + proton synchrotron branch is implemented; reverse-shock pγ/BH/pp/secondary species/cascade are not implemented.
- **Not yet implemented**: 2D/chi-resolved hadronic transport, full time-dependent pair cascade PDE

## AM3 / ASGARD Coexistence

- AM3 本地参考: `~/projects/_external/am3_reference` (HEAD: `7aba970b`)
- Python 只做 orchestration、wrapping、benchmark、API glue
