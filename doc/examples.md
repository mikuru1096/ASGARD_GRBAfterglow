# 示例

以下片段假设已按 [快速开始](quickstart.md) 创建 `model`，并已导入 NumPy。

## 多频光变

```python
times = np.logspace(2, 7, 120)
nu = np.array([3e9, 5e14, 2.42e17])
result = model.flux_density_grid(times, nu)

for frequency, curve in zip(nu, result.total):
    print(frequency, curve.max())
```

`total`、`fwd.sync`、`fwd.ssc` 的轴顺序都是 `(frequency, time)`。RS 未启用时，反向分量不贡献总光变。

## 逐点数据与宽频谱

```python
time_obs = np.array([1e3, 3e4, 2e5])
nu_obs = np.array([5e14, 3e9, 2.42e17])
prediction = model.flux_density(time_obs, nu_obs)

nu_sed = np.logspace(8, 25, 240)
sed = model.spectrum(1e4, nu_sed)
```

逐点接口要求时间和频率可广播为同一形状；它不是笛卡尔积查询。

## 频段积分

```python
xray = model.flux(
    time_s=1e4,
    nu_min_hz=0.3 * 2.417989e17,
    nu_max_hz=10.0 * 2.417989e17,
    num_points=96,
)
```

积分在频率上执行；输入仍是 Hz，输出单位取决于 flux-density 契约，参见 [公开 API](public_api.md)。

## 内部状态

```python
tracks = model.details(t_min=1e2, t_max=1e6)
forward = tracks.fwd
print(forward.radius.shape, forward.Gamma.shape)
print(forward.gamma_e.shape, forward.dN_dgamma_e.shape)
```

`details()` 用于诊断而不是拟合热路径。可选数组只在对应 solver/process 生成时存在。

## 自适应投影

```python
adaptive = model.flux_density_grid_adaptive(
    np.logspace(2, 6, 60),
    np.array([1e9, 1e14]),
)
print(adaptive.flux.total.shape)
```

自适应路径的容差和最大深度由 `SolverOptions` 控制；正式运行应与固定网格做一次收敛比较。

## 天图与偏振

```python
image = model.sky_image(t_obs=1e5, nu_obs=3e9, fov=2e-9, npixel=128)
pol = model.polarization(
    np.logspace(3, 6, 40),
    np.array([3e9]),
    magnetic_geometry="shock_random",
)
```

天图需要明确的角视场；偏振结果包含 Stokes 量及派生偏振度。可选磁场几何以运行时签名为准。

## Gaussian 结构化喷流

```python
from asgard_core import gaussian_jet

jet = gaussian_jet(
    energy_iso_erg=1e53,
    initial_lorentz_factor=300,
    core_angle_rad=0.05,
    outer_angle_rad=0.35,
    shell_duration_s=None,
    magnetar=None,
    spreading=False,
)
```

把 `jet` 传给一个新 `Model`，并将 `structured_num_theta/phi` 设为实际 patch 数、选择已支持的 `structured_backend`。不要修改已构造模型的字段，也不要只更换 profile 而保留单 patch 网格。

## 磁化反向激波

```python
from asgard_core import ReverseShock

model.reverse_shock = ReverseShock(
    enabled=True,
    shell_duration_s=10.0,
    upstream_sigma=1e-2,
    include_cross_zone_ic=False,
    include_ssc=True,
)
model.rvs_rad = model.fwd_rad
```

改变已建模型字段后，推荐重新构造 `Model`，确保 runtime config 与公开配置同步。RS 微物理可以通过独立 `rvs_rad` 设置。

## Hadronic pp gamma selector

```python
from asgard_core import Hadronic

hadronic = Hadronic(
    enabled=True,
    solver="am3_1d",
    num_proton_gamma=161,
    num_neutrino_frequency=121,
    pgamma_scheme="disabled",
    pair_cascade_iterations=1,
    pp_gamma_model="geant4",
)
```

默认 `delta` 保持旧路径。`geant4`、`sibyll`、`qgsjet`、`pythia8` 只改变 pp 的 π0 gamma 谱；它们是存在固有分段跳变的研究选项，见根目录 `BUG.md`。

## 最小拟合数据

```python
from asgard_core import Fitter, Param, Scale

fitter = Fitter(model)
fitter.add_flux_density(time_obs, nu_obs, flux_obs, flux_err)
params = [
    Param("energy", "jet.energy_iso_erg", 51.0, 54.0, Scale.LOG10),
    Param("density", "medium.number_density_cm3", -4.0, 2.0, Scale.LOG10),
]
print(fitter.loglike({"energy": 52.0, "density": 0.0}))
```

参数路径会经过公开配置解析；正式采样、支持路径与结果解释见 [拟合工作流](fitting_workflow.md)。

## 运行前检查

```python
assert np.all(np.isfinite(result.total))
assert np.all(result.total >= 0.0)
```

还应检查网格收敛、能量预算和随参数连续变化。示例成功运行只说明入口有效，不等于该物理组合已完成验收。
