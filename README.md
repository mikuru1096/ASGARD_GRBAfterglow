<p align="center">
  <img src="doc/assets/logo.png" alt="ASGARD logo" width="75%">
</p>

# ASGARD

ASGARD 是面向伽马射线暴余辉的数值建模、辐射计算和参数拟合工具。项目把高代价的动力学、粒子演化和辐射核放在 Fortran/f2py 数值层中，用 Python 提供公开 API、观测者投影、拟合接口、文档和验证脚本。

当前公开工作流以 `Model` 构建物理模型，以 `flux_density_grid`、`flux_density`、`spectrum`、`sky_image`、`polarization` 查询观测量，以 `Fitter` 连接 `emcee` 或 PyMultiNest 做参数推断。完整导出类型、构造器字段和选项含义见 `doc/public_api.md`。

## 当前能力

- 正向激波电子同步辐射、SSC、SSA、gamma-gamma 吸收和观测者等到达时间投影。
- 多个 1D 电子求解器，以及登记的 `fullhide_2d` / `charint_2d` 输运路径。
- 反向激波电子同步辐射、反向激波 SSC、FS/RS cross-zone IC。
- 正向激波和反向激波的 1D hadronic research path，包括 proton synchrotron、p-gamma、BH、pp、secondary 和 pair-cascade 相关诊断。
- 结构化喷流 patch 投影、偏振 Stokes 投影、天图渲染和频段积分。
- `Fitter.fit(..., sampler="emcee")` 与 `sampler="pymultinest"` 的公开拟合入口。

## 明确边界

- 重要数值物理在 Fortran 核中实现；Python 主要做 API glue、编排、验证和文档示例。
- `solver_options.geometry_projection="chi_eats_2d"` 当前只替换 FS synchrotron+SSA 的观测者投影；SSC、hadronic 和 pair cascade 仍按 shell-level 契约处理。
- Hadronic 正式路径保持 1D shell 契约，直到 chi-local photon field、hadron density、secondary feedback 和 observer projection 的物理契约完成。
- Jet spreading、自定义 `Medium` kernel dispatch、wind `k != 2`、`fullhide_1d` 外的 thermal-electron public runtime 是明确未支持边界。见 `doc/public_backend_limits.md`。
- 光变、谱断频和反映真实物理演化的参数应连续平滑；如果出现孤立跳变，应回到动力学、输运、源项或投影查 bug，不用 smoothing 或经验补丁掩盖。

## 快速开始

推荐环境是 Linux 或 WSL Ubuntu，Python 依赖通过 `uv` 管理，Fortran 扩展通过 `build_extensions.py` 构建。

```bash
git clone https://github.com/mikuru1096/ASGARD_private
cd ASGARD_private
uv sync
TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force
```

Ubuntu/Debian 若缺少系统编译工具，先安装 `gcc`、`g++`、`gfortran` 和 Python 开发头文件。

构建后先确认默认扩展能加载：

```bash
uv run python - <<'PY'
from src.Electron.electron_forward_fullhide_1d import fs_electron_fullhide_1d
print("ASGARD Fortran extension loaded")
PY
```

第一次建模请直接读 `doc/quickstart.md`。该文档给出完整可运行的 `Model(...)` 构造器，包括 `Jet`、`Medium`、`Observer`、`Radiation`、`Numerics`、`SolverOptions`、`ReverseShock` 和 `Hadronic` 的显式字段。

模型构建后，常用查询形状如下：

```python
import numpy as np

times_s = np.logspace(2.0, 7.0, 80)
frequencies_hz = np.array([1.0e9, 1.0e14, 1.0e18])

result = model.flux_density_grid(times_s, frequencies_hz)
print(result.total.shape)  # (num_frequency, num_time)
```

逐点观测数据使用 `model.flux_density(times_s, frequencies_hz)`；固定时刻宽频谱使用 `model.spectrum(time_s, nu_hz)`；频段积分使用 `model.flux(...)`。返回类型如 `FluxResult`、`SkyImage`、`PolarizationResult` 的字段说明见 `doc/public_api.md`。

## 验证

README 面向的公开查询路径由 smoke test 覆盖：

```bash
uv run python tests/readme_smoke_bench.py
```

文档和文本编码检查：

```bash
uv run python tools/check_text_encoding.py
git diff --check
```

修改 Fortran 后还必须重建受影响模块，并运行 `-Wline-truncation` 语法检查和最窄相关 smoke test。验证分层见 `doc/validation_and_benchmarks.md`。

## 文档地图

- `doc/quickstart.md`：第一次运行 ASGARD 的完整路径。
- `doc/public_api.md`：公开 API 选择手册，解释每个关键字的可选值、物理意义、效果和注意事项。
- `doc/parameter_reference.md`：参数路径、单位和拟合参数建议。
- `doc/fitting_workflow.md`：从观测数据到 `emcee` 拟合的完整教程。
- `doc/mcmc_fitting.md`：`emcee` 和 PyMultiNest 专题。
- `doc/physics_model.md`、`doc/physical_processes.md`：已实现物理模块和过程说明。
- `doc/algorithm_workflow.md`、`doc/numerical_methods.md`：数值链路、离散方程和求解器族。
- `doc/public_backend_limits.md`：公开 backend 的已知边界。
- `doc/developer_guide.md`、`doc/source_tree.md`、`doc/call_chain.md`：开发者入口。
- `doc/web_docs.md`：通过 `asgard-private` 发布合作者文档站。

当前未完成项集中维护在 `TODO.md`；当前计划见 `PLAN.md`。

## 开发约束

- 不添加物理或数值 fallback 来绕过失败；边界外输入应暴露问题。
- 非数值模拟 Python 代码不添加数值保护。
- benchmark 和 comparison 入口保留在 `tests/`；临时研究脚本放 `/tmp`。
- `output/` 下图像和 CSV 只有在可由已提交脚本复现并完成物理验收时才提交。

## 许可

Copyright (c) 2025 Jia Ren

本项目使用 BSD 3-Clause License。

## 引用

如果在论文、报告、网页或其他软件中使用、改写或参考 ASGARD 的核心算法，请注明项目来源。

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

项目地址：<https://github.com/mikuru1096/ASGARD_private>

网页文档：<https://hetools.cn/asgard-doc/>

在线界面：<https://hetools.cn>
