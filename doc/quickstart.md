# 新手快速开始

本文给出从零运行 ASGARD 的最短路径。完整安装细节见 `doc/installation.md`，公开接口细节见 `doc/public_api.md`。

## 1. 构建默认数值核

当前开发环境固定为 WSL Ubuntu + `uv`。首次运行正向激波电子模型前，先构建默认 Fortran 扩展：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
```

只做文档或纯 Python 示例时不需要重复构建。若修改了 Fortran，按 `doc/validation_and_benchmarks.md` 跑受影响模块。

## 2. 运行最小光变

```python
import numpy as np

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet

model = Model(
    jet=TophatJet(E_iso=1.0e52, Gamma0=300.0, theta_j=0.1),
    medium=ISM(n_ism=1.0),
    observer=Observer(z=0.1, theta_obs=0.0),
    fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.3, xi_N=0.1, ssc=True),
    setups=Setups(
        electron_solver="fullhide_1d",
        num_r=120,
        num_tobs=120,
        num_gam_e=121,
        num_nu=121,
        observer_time_min_s=1.0e2,
        observer_time_max_s=1.0e7,
    ),
)

times = np.logspace(2, 7, 80)
freqs = np.array([1.0e9, 1.0e14, 1.0e18])
result = model.flux_density_grid(times, freqs)
print(result.total.shape)
```

`result.total` 的形状是 `(num_frequency, num_time)`。常用分量：

- `result.fwd.sync`：正向激波同步辐射。
- `result.fwd.ssc`：正向激波同步自康普顿。
- `result.rev.sync`：反向激波同步辐射，未启用反向激波时为空或为零贡献。
- `result.cross_ic`：跨区逆康普顿，未启用时为 `None`。

对应输出图：

![ASGARD 多频光变](assets/tutorials/quick_light_curves.png)

## 3. 单点、光谱和频段积分

观测数据常是逐点的 \((t_i,\nu_i)\)。这种情况下使用 `flux_density`：

```python
times = np.array([1.0e3, 3.0e4, 1.0e5])
freqs = np.array([1.0e14, 1.0e9, 1.0e18])
points = model.flux_density(times, freqs)
```

固定时刻光谱使用：

```python
nu = np.logspace(8, 25, 200)
sed = model.spectrum(1.0e4, nu)
```

对应谱图：

![ASGARD 宽频谱](assets/tutorials/quick_spectra.png)

频段积分使用：

```python
band_flux = model.flux(
    time_s=1.0e4,
    nu_min_hz=1.0e14,
    nu_max_hz=1.0e18,
    num_points=96,
)
```

## 4. 单位

公开 API 使用 cgs 和观测者频率/时间：

| 量 | 单位 |
| --- | --- |
| 时间 | s |
| 频率 | Hz |
| 距离 | cm |
| 角度 | rad |
| 能量 | erg |

红移 \(z\) 通过 `Observer` 输入。如果没有显式给出 luminosity distance，运行时会由红移和宇宙学工具确定。

## 5. 结果是否可信

最基本的物理检查：

- 光变随时间应连续，除非存在明确的密度跳变、注入事件或 shock crossing。
- \(\nu_m\)、\(\nu_c\)、\(\nu_a\) 应平滑演化。
- 平滑改变参数时，光变峰时和峰值不应出现孤立跳变。
- 反向激波的 `reverse_sigma -> 0` 必须回到非磁化基线。

若这些检查失败，应回到动力学、电子输运、辐射源项和观测投影查 bug，不做 smoothing 或经验修补。
