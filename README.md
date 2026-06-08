# ASGARD：伽马射线暴余辉分析工具

ASGARD 是面向伽马射线暴余辉的数值模拟、辐射计算和拟合工具。代码以 Fortran 数值核为核心，用 Python 提供公开 API、运行编排、观测投影、拟合接口和基准脚本。

ASGARD 的目标是从爆波动力学、粒子谱演化和辐射转移出发，直接计算多波段余辉光变、谱、天图、偏振以及强子辐射相关诊断。当前主线覆盖正激波电子同步辐射、SSC、SSA、gamma-gamma 吸收、观测者等到达时间投影；反激波电子辐射、反激波 SSC、FS/RS cross-zone IC；以及正式 1D shell 契约下的强子过程。

## 主要特性

1. 正激波电子连续性方程的数值 PDE 求解器，包括多个 1D 求解器和已登记的 2D 电子输运路径。
2. 自洽同步辐射、SSC、SSA 和 gamma-gamma 衰减。
3. 反激波电子同步辐射、RS SSC、FS/RS cross-zone IC，以及反激波强子路径。
4. 正激波强子 `legacy_1d` 与 formal `am3_1d` research path。
5. 同步辐射偏振 Stokes 投影，覆盖 FS/RS electron synch 与 FS/RS hadronic synch。
6. `Model` public API、`Fitter` 拟合 API、benchmark/report 脚本和可复现 artifact 协议。

## 文档入口

完整文档从 `doc/index.md` 开始。

推荐阅读：

- `doc/installation.md`：环境、安装、本地 Fortran 扩展构建。
- `doc/user_guide.md`：常用 `Model` 工作流，包括光变、谱、RS、hadronic、偏振、天图和拟合。
- `doc/public_api.md`：当前 public API 契约。
- `doc/physics_model.md`：已实现物理模块和明确边界。
- `doc/numerical_methods.md`：Fortran 数值核、求解器族和数值验证目标。
- `doc/validation_and_benchmarks.md`：build gate、smoke test、benchmark refresh 和 artifact policy。
- `doc/developer_guide.md`：开发工作流和 review checklist。
- `doc/code_overview.md`、`doc/source_tree.md`、`doc/call_chain.md`：实现结构图。

当前开发基线：

- `AGENTS.md`
- `PLAN.md`
- `TODO.md`
- `doc/code_overview.md`
- `doc/public_backend_limits.md`

Benchmark 和 comparison 脚本位于 `tests/` 与 `scripts/benchmarks/`。`output/` 下的图像和 CSV 是 artifact；只有在能由已提交脚本复现且符合 `doc/benchmark_refresh_protocol.md` 时才应提交。

## 快速开始

确保系统中有 GNU 编译器、Python、NumPy、SciPy、Astropy 和 Matplotlib。Ubuntu/Debian 可先安装编译器：

```shell
sudo apt install gcc g++ gfortran
```

克隆仓库：

```shell
git clone https://github.com/mikuru1096/ASGARD_GRBAfterglow
cd ASGARD_GRBAfterglow
```

安装依赖：

```shell
pip install -r Requirements.txt
```

运行自动安装脚本：

```shell
python install.py
```

Linux 也可以使用：

```shell
bash install.sh
```

Windows PowerShell 可使用：

```powershell
.\install.ps1
```

编译完成后运行 demo：

```shell
python lc_spec_demo.py
```

当前开发机器推荐使用 WSL Ubuntu + uv：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python lc_spec_demo.py'
```

## 当前状态

Public runtime 可用于正激波余辉计算和 benchmark 工作流。反激波使用局部 shock-front `gamma34` 作为新注入电子能标，区域 3 磁场和 crossing 后热演化由显式 `U3/V3` thermal state 闭合。VegasAfterglow 是 comparison backend，不是 ASGARD 的 RS 物理目标。

正激波强子分支当前是正式 1D research path。反激波强子包含 light proton-synch path，并在开启 full-chain 过程时复用正式 1D hadronic kernels 处理 p-gamma、BH、pp、secondary 和 cascade coupling。当前 pair cascade 是 shell-sequence time-dependent gamma-gamma pair/synch cascade；未完成项和 public/backend 边界见 `TODO.md` 与 `doc/public_backend_limits.md`。

## 许可

Copyright (c) 2025 Jia Ren

本项目使用 BSD 3-Clause License。

## 引用要求

如果在其他软件项目、论文或公开材料中使用、改写或参考本项目核心算法，请在文档、说明页或论文中明确注明 ASGARD 项目来源。

推荐引用：

```bibtex
@ARTICLE{2024ApJ...962..115R,
       author = {{Ren}, Jia and {Wang}, Yun and {Dai}, Zi-Gao},
       title = "{Jet Structure and Burst Environment of GRB 221009A}",
       journal = {\apj},
       keywords = {Gamma-ray bursts, 629, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2024,
        month = feb,
       volume = {962},
       number = {2},
          eid = {115},
        pages = {115},
          doi = {10.3847/1538-4357/ad1bcd},
       archivePrefix = {arXiv},
       eprint = {2310.15886},
       primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024ApJ...962..115R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

项目名称：ASGARD

项目地址：<https://github.com/mikuru1096/ASGARD_GRBAfterglow>

## Web 界面

在线界面位于 <https://hetools.cn>，可用于比较 ASGARD 和 jetsimpy 的结果。
