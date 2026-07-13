# 公开 API

本页描述 `asgard_core` 顶层稳定入口。签名的最终事实来源是 `asgard_core/__init__.py` 与 `api_model.py`。

## 导入

```python
from asgard_core import (
    AdaptiveFluxResult, CharTrack, FitResult, Fitter, FluxPair, FluxResult,
    Hadronic, Magnetar, Medium, Model, Numerics, Observer, ObserverGrid,
    Param, PolarizationResult, Radiation, ReverseShock, Scale, SkyImage,
    SolverOptions, TabulatedMedium, TrackBundle, UniformMedium, WindMedium,
    gaussian_jet, power_law_jet, tabulated_angular_jet, top_hat_jet,
)
```

公开数值默认使用 cgs；字段名后缀明确单位。

## Model

```text
Model(
    *, jet, medium, observer, fwd_rad, numerics, observer_grid,
    solver_options, reverse_shock, hadronic, rvs_rad=None,
)
```

`Model` 聚合物理配置并提供观测查询。`rvs_rad` 只在反向激波启用时需要独立设置；省略时不应假定任意 RS 组合有效。

主要方法：

| 方法 | 用途 | 主要返回 |
| --- | --- | --- |
| `flux_density_grid(times, nu)` | 时间×频率笛卡尔积 | `FluxResult` |
| `flux_density(times, nu)` | 配对观测点 | `FluxResult` |
| `flux_density_grid_adaptive(times, nu)` | 自适应投影 | `AdaptiveFluxResult` |
| `flux_density_exposures(...)` | 曝光时间平均 | `FluxResult` |
| `spectrum(time, nu)` | 固定时刻频谱 | ndarray |
| `flux(...)` | 频段积分 | ndarray/float |
| `details(t_min, t_max)` | 内部状态 | `TrackBundle` |
| `sky_image(...)` | 天图 | `SkyImage` |
| `polarization(...)` | Stokes 投影 | `PolarizationResult` |

`projection_kind` 的默认值由各方法定义；不要把内部 EATS selector 当成任意字符串传入。

## Medium

### UniformMedium

```python
UniformMedium(number_density_cm3)
```

均匀数密度介质，最稳妥的公开基线。

### WindMedium

```text
WindMedium(a_star, density_floor_cm3, density_cap_cm3)
```

公开 wind 固定为 `k=2`。`density_floor_cm3` 给出远处 ISM 下限；可选 `density_cap_cm3` 限制内侧最大数密度。

### TabulatedMedium

```python
TabulatedMedium(radius_cm, density_cm3, label)
```

半径必须递增，密度为正。表格 profile 代表真实外部输入，不应通过内部 fallback 补点。

## Jet

### top_hat_jet

```python
top_hat_jet(
    energy_iso_erg, initial_lorentz_factor, opening_angle_rad,
    shell_duration_s, magnetar, spreading,
)
```

### gaussian_jet

```python
gaussian_jet(
    energy_iso_erg, initial_lorentz_factor, core_angle_rad,
    outer_angle_rad, shell_duration_s, magnetar, spreading,
)
```

### power_law_jet

```python
power_law_jet(
    energy_iso_erg, initial_lorentz_factor, core_angle_rad,
    energy_index, lorentz_index, outer_angle_rad,
    shell_duration_s, magnetar, spreading,
)
```

`lorentz_index=None` 时复用能量指数。`tabulated_angular_jet` 接收角度、能量与 Lorentz-factor 表格，适合外部结构模型。

结构化 profile 还需要 `Numerics.structured_num_*` 与 `SolverOptions.structured_*`；只更换 profile 不会自动增加 patch 分辨率。

## Observer

```python
Observer(z, viewing_angle_rad, viewing_azimuth_rad, luminosity_distance_cm)
```

`luminosity_distance_cm=None` 时由红移推得距离。两个角均为 rad；非轴对称结构必须设置 azimuth。

## Radiation

```text
Radiation(
    epsilon_e, epsilon_B, p, *, proton_energy_fraction,
    epsilon_b_floor=None, magnetic_decay_alpha_t,
    magnetic_decay_t0_s, accelerated_electron_fraction,
    thermal_electrons, include_ssc, proton_synch, include_pgamma,
    bethe_heitler, hadronic_inverse_compton, pp, neutrino,
    acceleration_efficiency, reverse_proton_energy_fraction,
    pgamma_scheme, pair_production,
)
```

关键字段：

| 字段 | 含义 |
| --- | --- |
| `epsilon_e`, `epsilon_B`, `p` | 电子、磁场能量分数及非热谱指数 |
| `accelerated_electron_fraction` | 加速电子数分数 |
| `include_ssc` | 发射 SSC；冷却模式另由 solver selector 控制 |
| `thermal_electrons` | 热电子路径；并非所有 solver 支持 |
| `proton_energy_fraction` | FS 非热质子能量分数 |
| `reverse_proton_energy_fraction` | RS 非热质子能量分数 |
| `pgamma_scheme` | 辐射配置侧的 pγ 选择，需与 hadronic 配置匹配 |

强子 bool 只启用过程；`Hadronic.enabled` 和 solver 仍需同时满足。未使用的 KN bool 已从公开 API 删除，因为 KN 核不由该旧字段控制。

## Numerics

```python
Numerics(
    num_radius,
    structured_num_theta, structured_num_phi,
    eats_num_theta, eats_num_phi,
    downstream_num_chi,
    num_observer_time, num_electron_gamma, num_photon_frequency,
    num_threads,
    electron_adaptive_substeps, electron_substep_rtol,
    electron_substep_min, electron_substep_max,
    initial_radius_cm,
)
```

三套几何网格不可混用：

- `structured_num_*`：jet patch。
- `eats_num_*`：观测者角积分。
- `downstream_num_chi`：有限下游坐标；1D 路径可为 `None`。

电子与光子频率格点分别控制粒子输运和辐射积分。网格更密不保证更正确，正式结果必须做收敛比较。

## ObserverGrid

```python
ObserverGrid(time_min_s, time_max_s)
```

这是内部计算覆盖范围；查询时间超出范围不应依赖外推。

## SolverOptions

前四个 selector（电子求解器、动力学、投影和电子--光子耦合）必须显式提供；
其余字段采用项目规范默认值。其中 `ssc_cooling_mode` 默认采用
`nakar_y_thomson`，也可显式设为 `none` 或 `numeric_ic_kn`：

| 字段 | 常用值 | 作用 |
| --- | --- | --- |
| `electron_solver` | `fullhide_1d`, `charint_1d`, `fullhide_2d` 等 | 电子输运核 |
| `dynamics_solver` | `forward_legacy` | 动力学核 |
| `geometry_projection` | `sed_legacy`, `sed_adaptive_theta`, `chi_eats_2d` | 观测投影 |
| `electron_photon_coupling` | `separated`, `joint` | 次级反馈 |
| `ssc_cooling_mode` | `nakar_y_thomson`（默认）, `none`, `numeric_ic_kn` | SSC 冷却 |
| `synchrotron_integration` | `fixed_grid`, `cyclotron` | 同步核 |
| `cooling_kernel` | `legacy` | 冷却实现 |
| `structured_backend` | `fortran_1d` | patch 后端 |
| `patch_sampling` | `uniform` 或受支持自适应值 | patch 采样 |
| `structured_parallel_mode` | `outer`, `inner`, `nested` | 并行所有权 |
| `fullhide2d_transport_model` | `legacy`, `pwn_cr_v1` | 2D 输运 |
| `fullhide2d_escape_mode` | `closed`, `free_outer` | 2D 逃逸边界 |

patch pilot、beaming、线程和 2D stochastic acceleration 参数只对相应 selector 生效。额外默认字段包括：

```text
adaptive_grid="manual"
projection_adaptive_rtol=2e-2
projection_adaptive_max_depth=4
structured_adaptive_rtol=0
structured_adaptive_max_depth=4
nu_callback=None
```

字符串到 Fortran ID 只在公开边界映射一次；内部核只接收命名常量对应的整数 ID。

## ReverseShock

```python
ReverseShock(
    enabled, shell_duration_s, upstream_sigma,
    include_cross_zone_ic, include_ssc,
)
```

`upstream_sigma` 是未激波 ejecta 的磁化参数。`include_cross_zone_ic` 与 RS/FS photon ownership 绑定；不能视作独立于 RS 的通用 IC 开关。

## Hadronic

```python
Hadronic(
    enabled, solver, num_proton_gamma, num_neutrino_frequency,
    pgamma_scheme, pair_cascade_iterations, pp_gamma_model="delta",
)
```

`solver` 支持正式登记的 `legacy_1d` 与 `am3_1d` 路径。pγ selector 要与 `Radiation.include_pgamma/pgamma_scheme` 共同解释。

`pp_gamma_model`：

| 值 | 行为 |
| --- | --- |
| `delta` | 默认，保持旧 pp gamma 路径 |
| `geant4` | Kafexhiu/AM3 Geant4 π0 gamma 谱 |
| `sibyll` | SIBYLL π0 gamma 谱 |
| `qgsjet` | QGSJET π0 gamma 谱 |
| `pythia8` | Pythia8 π0 gamma 谱 |

它只改变 π0 gamma；质子损失、neutrino 与 secondary pair 保持现有链。详细模型具有上游参数化固有分段跳变，属于 opt-in 研究功能。

## 返回类型

`FluxResult` 的核心字段：

- `total`：总 flux-density 网格。
- `fwd`、`rev`：各为 `FluxPair(sync, ssc)`。
- `cross_ic`：跨区 IC，可为 `None`。
- 强子与 pair 分量仅在对应过程运行后存在。

`TrackBundle` 聚合 FS/RS `CharTrack`；`CharTrack` 包含 `radius`、`Gamma`、`B_comv`、粒子网格和可选过程数组。可选数组为 `None` 表示该路径没有产生它，不表示零谱。

`SkyImage` 保存图像、坐标和矩；`PolarizationResult` 保存 Stokes 与派生量；`AdaptiveFluxResult` 同时保存 flux 和自适应诊断。

## Fitter

```python
Fitter(model, data=None, params=None, num_workers=1)
Param(name, path, lower, upper, scale=Scale.LINEAR)
```

`Scale` 包括 `LINEAR`、`LOG10` 和 `FIXED`。主要数据入口：

```python
fitter.add_flux_density(times, frequencies, flux, flux_err)
fitter.add_spectrum(time, frequencies, flux, flux_err)
fitter.add_flux(nu_min, nu_max, time, flux, flux_err, num_points=64)
fitter.loglike(values)
result = fitter.fit(params, sampler="emcee")
```

`sampler` 支持 `emcee` 与 `pymultinest`。PyMultiNest 当前要求所有参数 active；固定参数只由 emcee 路径处理。完整流程见 [拟合工作流](fitting_workflow.md)。

## 支持边界

公开字段存在不代表任意笛卡尔组合都已验收。下列组合尤其需要先查 [公开后端边界](public_backend_limits.md)：2D hadronic、thermal electron 与非默认 solver、structured+RS/hadronic、spreading、自定义介质以及 reverse cross-zone transfer。遇到边界外组合应暴露问题，不添加 fallback。
