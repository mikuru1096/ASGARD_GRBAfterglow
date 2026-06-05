# 安装与构建

本文档固定 ASGARD 当前推荐安装、构建和本地运行方式。默认开发环境是 Windows 主机上的 WSL Ubuntu，Python 包管理使用 `uv`，shell 命令需要 `rtk` 前缀。

## 环境要求

最低要求：

- Python `>=3.9`
- GNU 工具链：`gcc`, `g++`, `gfortran`
- OpenMP runtime：`libgomp`
- Python 依赖：`numpy`, `scipy`, `matplotlib`, `astropy`, `extinction`, `h5py`, `tqdm`
- 构建依赖：`setuptools`, `wheel`, `meson`, `ninja`

推荐开发环境：

- WSL Ubuntu
- `uv`
- 非交互式 shell 中先执行 `source ~/.wsl_env`
- 从 Windows 路径访问项目：`/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow`

## 克隆和安装

```bash
git clone https://github.com/mikuru1096/ASGARD_GRBAfterglow
cd ASGARD_GRBAfterglow
pip install -r Requirements.txt
python install.py
```

当前机器推荐命令：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv sync'
```

本地源码安装：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv pip install -e .'
```

## 构建 Fortran 扩展

ASGARD 的高代价数值核由 Fortran + f2py 构建。构建入口是 `build_extensions.py`。

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
```

`TMPDIR=/tmp` 是从 `/mnt/c` 构建时的推荐设置，可避免 Windows temp 目录时间戳和 Meson/f2py 缓存问题。

常用模块：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_transport_2d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_reverse_kernel --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_forward_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_reverse_1d --force'
```

多个模块可以一次构建：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --module electron_forward_transport_2d --force'
```

## 编译产物

构建后会生成平台相关扩展：

- `src/*.so`
- `src/Electron/*.so`
- `src/Radiation/*.so`
- `src/Hadronic/*.so`
- `src/Dynamics/*.so`
- `src/Interpolation/*.so`

这些扩展是本地构建产物，不应把 `.so`, `.o`, `.mod`, `.smod`, `.buildcache/` 作为源码提交。

## 运行 demo

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python lc_spec_demo.py'
```

输出：

- `Radiation_Lightcurves.pdf`
- `Radiation_Spectra.pdf`

## 常见问题

### 找不到 Fortran 扩展

先确认对应模块已构建：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && find src -name "*.so" | sort'
```

若缺失，按上面的 `build_extensions.py --module ... --force` 重建。

### f2py 或 Meson 从 Windows temp 目录失败

使用 `TMPDIR=/tmp`，并从 WSL shell 调用构建命令。

### `gh` 或 GitHub CLI 不存在

项目构建不依赖 `gh`。推送或创建远端仓库时可使用 Windows Git Credential Manager，或使用 Codex GitHub 连接器；构建和测试仍默认走 WSL。

### Windows 原生编译

默认不使用 Windows/MSYS2 编译。只有当目标是 Windows `f2py` / `.pyd` 运行时才进入 Windows 原生工具链。
