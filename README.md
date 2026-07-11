# ASGARD

ASGARD 是伽马射线暴余辉的动力学、粒子输运、辐射和参数推断工具。高代价数值核由 Fortran/f2py 实现，Python 公开层负责物理配置、观测者投影和拟合。

## 能力与边界

- 正向激波同步辐射、SSC、SSA、γγ 吸收及等到达时间投影。
- 多种 1D 电子求解器和有限 q-shell 的 2D 输运。
- 反向激波同步、SSC、跨区 IC 与磁化激波诊断。
- 1D 强子研究路径：质子同步、pγ、BH、pp、次级粒子与级联诊断。
- top-hat、Gaussian、power-law 和表格角结构喷流；天图与偏振。
- `Fitter` 的 `emcee` 与 PyMultiNest 接口。

`chi_eats_2d` 当前只替换 FS synchrotron+SSA 投影；强子、SSC 和 pair cascade 仍使用 shell-level 契约。Jet spreading、自定义介质 kernel、wind `k != 2` 及若干组合尚未进入公开稳定范围，详见 [后端边界](doc/public_backend_limits.md)。

## 安装

推荐 WSL Ubuntu 或 Linux，需 Python 3.12、`uv`、GCC/G++ 和 gfortran。

```bash
git clone https://github.com/mikuru1096/ASGARD_GRBAfterglow.git
cd ASGARD_GRBAfterglow
uv sync
TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force
```

验证公开包和扩展：

```bash
uv run python - <<'PY'
import asgard_core
from src.Electron.electron_forward_fullhide_1d import fs_fullhide_1d
print("ASGARD ready")
PY
```

完整环境与按模块构建方法见 [安装指南](doc/installation.md)。

## 最小查询

完整 `Model` 构造器见 [快速开始](doc/quickstart.md)。创建 `model` 后：

```python
import numpy as np

times = np.logspace(2, 7, 80)
frequencies = np.array([1e9, 1e14, 1e18])
result = model.flux_density_grid(times, frequencies)
print(result.total.shape)  # (3, 80)

points = model.flux_density(times[:3], frequencies)
sed = model.spectrum(1e4, np.logspace(8, 25, 200))
band = model.flux(time_s=1e4, nu_min_hz=1e14, nu_max_hz=1e18)
```

公开量使用 cgs：时间 s、频率 Hz、距离 cm、角度 rad、能量 erg。输出形状和分量见 [公开 API](doc/public_api.md)。

## 文档

- [文档入口](doc/index.md)：按用户、物理和开发任务导航。
- [快速开始](doc/quickstart.md)：首个可运行模型。
- [示例](doc/examples.md)：查询、反向激波、2D 和 prompt。
- [公开 API](doc/public_api.md)：构造器、selector 和返回类型。
- [拟合工作流](doc/fitting_workflow.md)：数据、参数、likelihood 和采样。
- [物理模型](doc/physics_model.md)与[数值方法](doc/numerical_methods.md)。
- [开发指南](doc/developer_guide.md)：代码修改、构建与验证。

活动工作只记在 [TODO](TODO.md)，已知未修缺陷只记在 [BUG](BUG.md)。光变、断频及连续物理参数若出现无物理来源的跳变，应回到动力学、输运、源项或投影查错，不以平滑或经验后处理掩盖。

## 开发验证

```bash
uv run python tools/check_text_encoding.py
uv run --with "mkdocs<2" --with "mkdocs-material>=9.5" mkdocs build --strict --site-dir /tmp/asgard-site
git diff --check
```

修改 Fortran 后必须重建受影响扩展，并用 `-Wline-truncation` 检查干净 source closure；具体命令见[开发指南](doc/developer_guide.md)。临时 benchmark 放 `/tmp`，不提交生成物。

## 许可与引用

Copyright (c) 2025 Jia Ren。代码采用 BSD 3-Clause License。

使用 ASGARD 核心算法时请引用项目，并引用：Ren, Wang & Dai (2024), *ApJ*, 962, 115, DOI [10.3847/1538-4357/ad1bcd](https://doi.org/10.3847/1538-4357/ad1bcd)。

项目：<https://github.com/mikuru1096/ASGARD_GRBAfterglow>
