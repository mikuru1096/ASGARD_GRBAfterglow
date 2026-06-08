# ASGARD 代码规模统计

统计时间戳：2026-06-08 19:38:31 CST

## 统计口径

- 统计对象：当前工作树中 Git 跟踪的源码文件。
- 总代码行数扩展名：`.py`, `.f90`, `.F90`, `.f`, `.for`, `.f95`, `.c`, `.h`, `.hpp`, `.cpp`, `.sh`。
- Fortran 核代码行数：`src/` 下 Git 跟踪的 Fortran 源文件，扩展名为 `.f90`, `.F90`, `.f`, `.for`, `.f95`。
- 行数口径：物理行数，包含空行和注释；不统计文档、数据、输出图、缓存和未跟踪附件。

## 当前统计

| 指标 | 文件数 | 行数 |
| --- | ---: | ---: |
| 项目总代码 | 173 | 39,588 |
| `src/` Fortran 核代码 | 86 | 19,770 |

## Fortran 核分目录

| 目录 | 文件数 | 行数 |
| --- | ---: | ---: |
| `src` | 1 | 69 |
| `src/Dynamics` | 3 | 1,007 |
| `src/Electron` | 17 | 6,968 |
| `src/Electron/slatec` | 37 | 5,647 |
| `src/Hadronic` | 18 | 4,325 |
| `src/Interpolation` | 3 | 396 |
| `src/Radiation` | 6 | 1,027 |
| `src/Structured` | 1 | 331 |

## 本次统计命令

```bash
rtk wsl -e bash -c 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && python3 - <<'"'"'PY'"'"'
from pathlib import Path
import subprocess

code_exts = {".py", ".f90", ".F90", ".f", ".for", ".f95", ".c", ".h", ".hpp", ".cpp", ".sh"}
files = subprocess.check_output(["git", "ls-files"], text=True).splitlines()
code_files = [Path(p) for p in files if Path(p).exists() and Path(p).suffix in code_exts]
fortran_kernel_files = [
    p for p in code_files
    if p.parts and p.parts[0] == "src" and p.suffix.lower() in {".f90", ".f", ".for", ".f95"}
]

def count(paths):
    return sum(len(p.read_text(encoding="utf-8", errors="replace").splitlines()) for p in paths)

print("total_code_files", len(code_files))
print("total_code_lines", count(code_files))
print("fortran_kernel_files", len(fortran_kernel_files))
print("fortran_kernel_lines", count(fortran_kernel_files))
PY'
```
