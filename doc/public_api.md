# 公开 API 参考

本文档记录当前 public API 的稳定入口和语义边界。实现文件主要是 `ASGARD/api_model.py`、`ASGARD/api_observe.py` 和 `ASGARD/api_fit.py`。

## 导入入口

```python
from ASGARD import Model, ISM, Wind, TophatJet, GaussianJet, PowerLawJet
from ASGARD import TwoComponentJet, StepPowerLawJet, Ejecta
from ASGARD import Observer, Radiation, Setups, Fitter, Param, Scale
from ASGARD import observe, run_fit, units
```

## 介质

`ISM` 和 `Wind` 是 public constructor aliases，返回 `Medium` 对象。backend dispatch 读取 `medium.kind` 与 `to_kernel_params()`，而不是依赖 Python subclass。

### `ISM`

```python
ISM(n_ism=1.0)
ISM(n0=1.0)
```

Kernel 参数：

- `d_ne = n_ism`
- `a_star = -1`

### `Wind`

```python
Wind(A_star=1.0, n_ism=0.1, n0=None, k=2.0)
```

当前 backend 只支持 `k=2.0`。如果设置 `n0`，wind density 会由相应 transition radius 截断。

### `Medium`

`Medium(rho=callable)` 可以在 Python 层评估密度，但当前 Fortran kernel dispatch 不支持任意用户自定义介质。边界见 `doc/public_backend_limits.md`。

## 喷流

`TophatJet`, `GaussianJet`, `PowerLawJet`, `TwoComponentJet`, `StepPowerLawJet` 和 `Ejecta` 是 public constructor aliases，返回 `JetProfile` 对象。runtime 通过 `jet.kind` 选择直接 tophat kernel 或 patch projection。

### `TophatJet`

```python
TophatJet(E_iso=1.0e52, Gamma0=300.0, theta_j=0.1)
TophatJet(E_iso=1.0e52, lf=300.0, theta_j=0.1)
```

常用字段：

- `E_iso`
- `lf` / `Gamma0`
- `theta_j` / `theta_c`
- `duration`
- `magnetar`
- `spreading`

`spreading=True` 在对象层可接受，但当前 backend 没有实现 jet spreading dynamics。

### 结构化喷流

```python
GaussianJet(E_iso=1.0e52, Gamma0=300.0, theta_c=0.05, theta_max=0.5)
PowerLawJet(E_iso=1.0e52, Gamma0=300.0, theta_c=0.05, k=2.0, theta_max=0.5)
TwoComponentJet(E_iso_c=1.0e52, Gamma0_c=300.0, E_iso_outer=1.0e50, Gamma0_outer=30.0, theta_c=0.05, theta_o=0.2)
```

Public model 通过 `energy_iso(phi, theta)` 和 `gamma0(phi, theta)` 评估 patch。

## 观测者

```python
Observer(z=0.1, theta_obs=0.0, phi_obs=0.0)
Observer(lumi_dist_cm=1.0e28, z=0.1, theta_obs=0.0)
```

如果未提供 luminosity distance，runtime 会从 redshift 和宇宙学工具确定。

## 辐射参数

```python
Radiation(eps_e=0.1, eps_B=1.0e-3, p=2.3, xi_N=0.1, ssc=True)
```

主要电子开关：

- `eps_e`
- `eps_B`
- `p`
- `xi_N`
- `thermal_electrons`
- `ssc`
- `kn`
- `epsilon_b_floor`
- `magnetic_decay_alpha_t`
- `magnetic_decay_t0_s`

强子开关：

- `epsilon_p`
- `proton_synch`
- `pg`
- `bethe_heitler`
- `hadronic_inverse_compton`
- `pp`
- `neutrino`
- `eta_acc`
- `pgamma_scheme`
- `pair_production`
- `reverse_epsilon_p`

## 数值设置

`Setups` 控制数值网格、求解器选择、RS flags、hadronic flags 和观测者时间范围。

重要字段：

- `num_r`, `num_theta`, `num_phi`, `num_tobs`
- `num_gam_e`, `num_nu`, `num_chi`
- `observer_time_min_s`, `observer_time_max_s`
- `electron_solver`
- `index_y`
- `index_syn_integr`
- `ssc_cooling`
- `rvs_shock`, `rvs_ssc`, `include_cross_zone_ic`
- `reverse_delta_t_s`, `reverse_sigma`
- `hadronic_enabled`, `hadronic_solver`, `num_gam_p`, `num_nu_nu`
- `pair_cascade_iterations`
- `num_threads`

Solver aliases：

- `fullhide` -> `fullhide_1d`
- `slc1` -> `slc1_1d`
- `charint` -> `charint_1d`
- `t2g1` -> `t2g1_1d`
- `weno5` -> `weno5_1d`

已登记 solver names：

- `fullhide_1d`
- `slc1_1d`
- `charint_1d`
- `t2g1_1d`
- `weno5_1d`
- `fullhide_2d`
- `charint_2d`

## `Model`

构造方式：

```python
model = Model(jet=jet, medium=medium, observer=observer, fwd_rad=fwd_rad, rvs_rad=rvs_rad, setups=setups)
```

兼容的位置参数形式：

- `Model(jet, medium, observer, radiation)`
- `Model(medium, jet, observer, radiation)`

主要方法：

- `flux_density_grid(times_s, nu_hz, projection_kind="lightcurve")`：完整二维网格投影，输出 `total` 形状为 `(num_nu, num_time)`。
- `flux_density(times_s, nu_hz, projection_kind="lightcurve")`：点对投影，适合 matched `(time_i, nu_i)` 数组。
- `flux_density_exposures(times_s, nu_hz, exposures_s, num_subsamples=4)`：曝光平均 flux density。
- `spectrum(time_s, nu_hz, projection_kind="sed")`：单时刻谱。
- `flux(time_s, nu_min_hz, nu_max_hz, num_points=64, projection_kind="sed")`：频段积分 flux。
- `sky_image(t_obs, nu_obs, fov, npixel=128)`：观测者平面天图。
- `polarization(times_s, nu_hz, magnetic_geometry="shock_random", local_emissivity="analytic_then_kernel")`：同步辐射 Stokes 和偏振诊断。
- `details(t_min=None, t_max=None)`：内部动力学、电子、强子和观测者诊断状态。

`projection_kind` 只接受：

- `"lightcurve"`：光变、拟合和多频段时间序列默认路径。对 `geometry_kernel="chi_eats_2d"`，FS synchrotron+SSA 使用 χ 分辨专用 EATS 投影；非 χ 分量保持 shell-level projection。
- `"sed"`：固定时刻扫频率和频段积分默认路径。该路径使用通用 shell SED 插值器，避免为 SED-only 图走光变专用投影热路。

`geometry_kernel="chi_eats_2d"` 是 opt-in 几何核，只支持 2D electron solver。它当前只使 FS synchrotron+SSA 的 lightcurve projection 使用 χ 分辨有限厚壳层；SSC、hadronic 和 pair cascade 仍保持 shell-level contract。

## `Fitter`

```python
from ASGARD import Fitter, Param, Scale

fitter = Fitter(model=model)
fitter.add_flux_density(times_s, nu_hz, flux, err)
fitter.params = [Param("logE", "jet.E_iso", 50.0, 54.0, Scale.LOG10)]
```

`Param` 形式：

- `Param(name, lower, upper, scale)`
- `Param(name, path, lower, upper, scale)`

`Scale.LOG` 和 `Scale.LOG10` 通过 `10**value` 转换；`Scale.FIXED` 固定使用 `lower`。

采样辅助：

- `Fitter.run_emcee(initial, nwalkers, nsteps)`
- `Fitter.run_multinest(...)`

## `observe` / `run_fit`

`observe(model, config=..., spectrum_output=...)` 是 demo 和脚本使用的较低层执行入口。普通交互使用优先选择 `Model` 方法。

`run_fit(config)` 保留给 config-style 工作流兼容使用。

## 返回类型

常见结果对象：

- `FluxResult`：`total`, `fwd`, `rev`, `cross_ic`
- `FluxPair`：`sync`, `ssc`
- `SkyImage`：image array 和 coordinate grids
- `PolarizationResult`：Stokes 和偏振诊断
- `TrackBundle`：forward/reverse dynamics, electron, hadronic, observer diagnostics

## 公开边界

Public API 接受的部分选项先于 backend 完成。未支持选项应显式失败，而不是静默 fallback。当前固定边界见 `doc/public_backend_limits.md`。
