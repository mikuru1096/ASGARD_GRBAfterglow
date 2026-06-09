# 使用指南

本文档给出 ASGARD public API 的常用运行方式。所有示例都以当前 `ASGARD` 包入口为准；类名、函数名、物理缩写和命令参数保持英文标识。

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
print(result.total.shape)
```

`flux_density_grid(times_s, nu_hz)` 返回 `FluxResult`。其中：

- `total`：总 flux density。
- `fwd.sync`, `fwd.ssc`：正激波同步辐射和 SSC。
- `rev.sync`, `rev.ssc`：反激波同步辐射和 SSC。
- `cross_ic`：FS/RS cross-zone IC；未启用时为 `None`。

默认 `projection_kind="lightcurve"`，适合光变、拟合和多频段时间序列。若启用 `geometry_kernel="chi_eats_2d"`，该路径对 FS synchrotron+SSA 使用 χ 分辨专用 EATS 投影；SSC、hadronic 和 pair cascade 仍保持 shell-level contract。

## 点对点流量

观测数据若是 `(time_i, frequency_i)` 点对，可使用：

```python
times = np.array([1.0e3, 3.0e4, 1.0e5])
freqs = np.array([1.0e14, 1.0e9, 1.0e18])
points = model.flux_density(times, freqs)
```

当 `times` 和 `freqs` 形状相同时，返回对应点的 flux，而不是完整二维网格。

## 谱

```python
nu = np.logspace(8, 25, 200)
sed = model.spectrum(1.0e4, nu)
```

`spectrum(time_s, nu_hz)` 是 `projection_kind="sed"` 的便捷接口，等价于 `flux_density_grid([time_s], nu_hz, projection_kind="sed").total[:, 0]`。该路径使用通用 shell SED 插值器，适合固定时刻扫频率；如果需要强制比较光变专用投影，可显式传入 `projection_kind="lightcurve"`。

## 频段积分流量

```python
band_flux = model.flux(time_s=1.0e4, nu_min_hz=1.0e14, nu_max_hz=1.0e18, num_points=96)
```

`flux` 对频率网格积分，默认也使用 `projection_kind="sed"`，适合近似宽能段 flux。若需要严格仪器响应，应在外部使用响应矩阵或带权积分。

## 2D χ-EATS 投影

```python
model = Model(
    jet=TophatJet(E_iso=1.0e52, Gamma0=300.0, theta_j=0.1),
    medium=ISM(n_ism=1.0),
    observer=Observer(z=0.1, theta_obs=0.05),
    fwd_rad=Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.3, xi_N=0.1, ssc=True),
    setups=Setups(
        electron_solver="fullhide_2d",
        geometry_kernel="chi_eats_2d",
        num_chi=24,
        num_phi=30,
    ),
)

lc = model.flux_density_grid(times, freqs, projection_kind="lightcurve")
sed = model.spectrum(1.0e5, np.logspace(8, 25, 200))
```

`chi_eats_2d` 只改变 FS synchrotron+SSA 的 observer projection。对轴观测可使用轴对称 φ 折叠；偏轴观测必须设置 `num_phi >= 2`。光变默认走 χ 分辨专用投影，单时刻 SED 默认走通用 SED 插值器。

## 曝光时间平均

```python
times = np.array([1.0e4, 2.0e4])
freqs = np.array([1.0e14, 1.0e14])
exposures = np.array([100.0, 500.0])
avg = model.flux_density_exposures(times, freqs, exposures, num_subsamples=4)
```

该接口对子曝光采样取平均，适合曝光时间不可忽略的光变点。

## 反激波

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

## 强子过程

正激波强子 formal path 示例：

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

强子路径约束：

- 当前 formal path 是 1D shell contract。
- 未完成的 hadronic transport / cascade 扩展集中记录在根目录 `TODO.md`。
- Python hadronic 模块只做 orchestration、wrapping、benchmark；最终微物理核在 `src/Hadronic/*.f90`。

## 偏振

```python
times = np.logspace(3, 6, 40)
freqs = np.array([1.0e14])
pol = model.polarization(times, freqs, magnetic_geometry="shock_random")
print(pol.pi, pol.angle, pol.q, pol.u)
```

当前偏振路径只混合同步辐射 Stokes 贡献；非同步分支不进入偏振 Stokes。

## 天图

```python
from ASGARD import units

image = model.sky_image(t_obs=1.0e5, nu_obs=1.0e9, fov=500.0 * units.uas, npixel=128)
print(image.image.shape)
```

`sky_image` 使用 observer projection 结果渲染天空图。高分辨图像应先用较小 `npixel` 做 smoke run。

## 诊断详情

```python
details = model.details(t_min=1.0e2, t_max=1.0e7)
print(details.fwd.nu_m)
print(details.fwd.nu_c)
print(details.fwd.nu_a)
```

`details()` 返回内部诊断轨道，例如 dynamics、break frequencies、RS thermal state 和 hadronic branches。它用于理解物理状态，不应替代主观测量 projection。

## 拟合

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

正式拟合使用 `ASGARD.Fitter` API 组织 likelihood；CPU/GPU 网格、采样和 posterior 产品应按当前任务的 fitting 规则执行。

## 输出单位

ASGARD public API 使用 cgs 和 Hz/s：

- time：second
- frequency：Hz
- distance：cm
- angle：radian；`ASGARD.units` 提供常用角单位转换。
- flux density：沿用项目内部 cgs / observer projection 约定，绘图脚本负责转换显示。

## 使用边界

不要把 public API 暴露但 backend 尚未支持的能力当作已支持。完整列表集中维护在根目录 `TODO.md` 和 `doc/public_backend_limits.md`。
