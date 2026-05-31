# User Guide

本文档给出 ASGARD public API 的常用运行方式。所有示例都以当前 `ASGARD` 包入口为准。

## 最小正激波光变

```python
import numpy as np

from ASGARD import ISM, Model, Observer, Radiation, Setups, TophatJet

model = Model(
    jet=TophatJet(E_iso=1.0e52, Gamma0=300.0, theta_j=0.1),
    medium=ISM(n_ism=1.0),
    observer=Observer(z=0.1, theta_obs=0.0),
    fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.3, xi_N=0.1, ssc=True),
    setups=Setups(
        num_r=120,
        num_tobs=120,
        num_gam_e=121,
        num_nu=121,
        observer_time_min_s=1.0e2,
        observer_time_max_s=1.0e7,
        electron_solver="fullhide_1d",
        ssc_cooling=True,
    ),
)

times = np.logspace(2, 7, 80)
freqs = np.array([1.0e9, 1.0e14, 1.0e18])
result = model.flux_density_grid(times, freqs)

print(result.total.shape)  # (num_freq, num_time)
```

`flux_density_grid(times_s, nu_hz)` 返回 `ModelFluxResult`：

- `total`: 总 flux density。
- `fwd.sync`, `fwd.ssc`: 正激波同步辐射和 SSC。
- `rev.sync`, `rev.ssc`: 反激波同步辐射和 SSC。
- `cross_ic`: FS/RS cross-zone IC，未启用时为 `None`。

## 点对点 flux

如果观测数据是 `(time_i, frequency_i)` 点对，可以使用：

```python
times = np.array([1.0e3, 3.0e4, 1.0e5])
freqs = np.array([1.0e14, 1.0e9, 1.0e18])
points = model.flux_density(times, freqs)
```

当 `times` 和 `freqs` 形状相同时，返回的是对应点的 flux，而不是完整二维网格。

## 谱

```python
nu = np.logspace(8, 25, 200)
sed = model.spectrum(1.0e4, nu)
```

`spectrum(time_s, nu_hz)` 是 `flux_density_grid([time_s], nu_hz).total[:, 0]` 的便捷接口。

## 频段积分 flux

```python
band_flux = model.flux(time_s=1.0e4, nu_min_hz=1.0e14, nu_max_hz=1.0e18, num_points=96)
```

`flux` 对频率网格积分，适合近似宽能段 flux。需要严格仪器响应时，应在外部使用响应矩阵或带权积分。

## 有曝光时间的观测点

```python
times = np.array([1.0e4, 2.0e4])
freqs = np.array([1.0e14, 1.0e14])
exposures = np.array([100.0, 500.0])
avg = model.flux_density_exposures(times, freqs, exposures, num_subsamples=4)
```

该接口对子曝光采样做平均，适合曝光时间不可忽略的光变点。

## Reverse shock

开启反激波：

```python
model = Model(
    jet=TophatJet(E_iso=1.0e52, Gamma0=300.0, theta_j=0.1, duration=20.0),
    medium=ISM(n_ism=1.0),
    observer=Observer(z=0.1, theta_obs=0.0),
    fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-4, p=2.3, xi_N=0.1, ssc=True),
    rvs_rad=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.3, xi_N=0.1, ssc=True),
    setups=Setups(
        rvs_shock=True,
        rvs_ssc=True,
        include_cross_zone_ic=True,
        reverse_delta_t_s=20.0,
        reverse_sigma=0.0,
    ),
)
```

当前 RS baseline：

- 新注入电子能标使用 shock-front `gamma34`。
- 区域 3 磁场和 crossing 后热演化使用显式 `U3/V3` thermal state。
- `reverse_sigma` 控制 upstream magnetization；`sigma -> 0` 必须回到非磁化 baseline。
- `details().rvs.nu_m` 是 diagnostic break，不替代输运后的电子谱峰。

## Hadronic

Forward-shock hadronic formal path：

```python
model = Model(
    jet=TophatJet(E_iso=1.0e52, Gamma0=300.0, theta_j=0.1),
    medium=ISM(n_ism=1.0),
    observer=Observer(z=0.1, theta_obs=0.0),
    fwd_rad=Radiation(
        eps_e=0.1,
        eps_B=1.0e-3,
        p=2.3,
        xi_N=0.1,
        epsilon_p=0.2,
        proton_synch=True,
        pg=True,
        bethe_heitler=True,
        pp=True,
        neutrino=True,
        pair_production=True,
        pgamma_scheme="hummer_2010_response",
    ),
    setups=Setups(
        hadronic_enabled=True,
        hadronic_solver="am3_1d",
        num_gam_p=161,
        pair_cascade_iterations=2,
    ),
)
```

Hadronic 约束：

- 当前 formal path 是 1D shell contract。
- 2D / chi-resolved hadronic transport 未实现。
- IC-mediated electromagnetic cascade 未实现。
- Python hadronic 模块只做 orchestration、wrapping、benchmark；最终微物理核在 `src/Hadronic/*.f90`。

## Polarization

```python
times = np.logspace(3, 6, 40)
freqs = np.array([1.0e14])
pol = model.polarization(times, freqs, magnetic_geometry="shock_random")

print(pol.pi)      # polarization fraction
print(pol.angle)   # polarization angle
print(pol.q, pol.u)
```

当前偏振路径只混合同步辐射 Stokes 贡献。非同步分支不进入偏振 Stokes。

## Sky image

```python
from ASGARD import units

image = model.sky_image(t_obs=1.0e5, nu_obs=1.0e9, fov=500.0 * units.uas, npixel=128)
print(image.image.shape)
print(image.x_uas, image.y_uas)
```

`sky_image` 使用 observer projection 结果渲染天空图。对高分辨图像，应先用较小 `npixel` 做 smoke run。

## Details

```python
details = model.details(t_min=1.0e2, t_max=1.0e7)
print(details.fwd.nu_m)
print(details.fwd.nu_c)
print(details.fwd.nu_a)
```

`details()` 返回内部诊断轨道，例如 dynamics、break frequencies、RS thermal state、hadronic branches。它用于理解物理状态，不应替代主观测量 projection。

## Fitting

基础拟合对象：

```python
from ASGARD import Fitter, Param, Scale

fitter = Fitter(model=model)
fitter.add_flux_density(times_s, nu_hz, flux, flux_err)

fitter.params = [
    Param("logE", "jet.E_iso", 50.0, 54.0, Scale.LOG10),
    Param("logn", "medium.n_ism", -4.0, 2.0, Scale.LOG10),
    Param("p", "fwd_rad.p", 2.0, 2.8),
]

ll = fitter.loglike({"logE": 52.0, "logn": 0.0, "p": 2.3})
```

采样脚本：

- `scripts/fitting/mcmc_fit.py`
- `scripts/fitting/multinest_fit.py`
- `scripts/fitting/mpi_run.sh`

正式拟合的 CPU/GPU 网格、采样和 posterior 产品应按项目当前 fitting skill 或对应任务说明执行。

## 输出单位

ASGARD public API 使用 cgs 和 Hz/s：

- time: second
- frequency: Hz
- distance: cm
- angle: radian，`ASGARD.units` 提供常用角单位转换。
- flux density: 当前 API 输出沿用项目内部 cgs/observer projection 约定，绘图脚本负责转换显示。

## 使用边界

不要把下面能力当作已支持：

- 自定义 `Medium` 直接进入 Fortran kernel dispatch。
- wind `k != 2`。
- jet spreading 动力学 backend。
- thermal electrons outside `fullhide_1d`。
- 2D/chi-resolved hadronic transport。
- IC-mediated electromagnetic cascade。
