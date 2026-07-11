# 安装与构建

推荐 WSL Ubuntu 或 Linux。Windows 原生环境只用于明确要求的 `.pyd` 构建，不是默认开发路径。

## 依赖

- Python 3.12 与 `uv`
- GCC、G++、GNU Fortran
- Meson/f2py 相关 Python 依赖由 `uv.lock` 固定
- Git

Ubuntu/Debian 可先安装系统工具：

```bash
sudo apt update
sudo apt install -y build-essential gfortran python3-dev git
```

## 克隆与 Python 环境

```bash
git clone https://github.com/mikuru1096/ASGARD_GRBAfterglow.git
cd ASGARD_GRBAfterglow
uv sync
uv run python -c "import asgard_core; print('Python API ready')"
```

不要从 `asgard-private` URL 克隆公开仓库。依赖变化通过 `uv.lock` 管理，不在系统 Python 中逐个补包。

## 构建 Fortran 扩展

默认电子路径：

```bash
TMPDIR=/tmp uv run python build_extensions.py \
  --module electron_forward_fullhide_1d --force
```

常用模块：

```bash
TMPDIR=/tmp uv run python build_extensions.py \
  --module Dynamics_forward --module Dynamics_reverse --force
TMPDIR=/tmp uv run python build_extensions.py \
  --module hadronic_forward_1d --module hadronic_reverse_1d --force
TMPDIR=/tmp uv run python build_extensions.py \
  --module structured_jet_1d --force
```

可用模块名由 `build_extensions.py` 定义。`/mnt/c` 上构建必须令 `TMPDIR=/tmp`，避免 Meson/f2py 在 Windows temp 路径混用文件语义。

## 验证

```bash
uv run python - <<'PY'
from src.Electron.electron_forward_fullhide_1d import fs_electron_fullhide_1d
print("Fortran extension ready")
PY
uv run python tools/check_text_encoding.py
git diff --check
```

修改 Fortran 后，用相同 source closure 加 `-Wline-truncation` 编译检查；只导入 Python 包不能替代真实扩展构建。

## 失败定位

- `ModuleNotFoundError`：先确认目标模块已经由 `build_extensions.py` 构建。
- 编译器缺失：检查 `gfortran --version`，不要切换到另一个 ABI 的预编译产物。
- `/mnt/c` 构建异常：确认 `TMPDIR=/tmp` 且命令在 WSL 内执行。
- 依赖导入失败：运行 `uv sync --frozen`，不要手动改 `.venv`。

正式运行前转到 [快速开始](quickstart.md)；构建矩阵和验证层级见 [开发指南](developer_guide.md)。
