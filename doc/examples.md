# 教程示例

本文按“先得到光变，再拆分物理分量，再检查内部量，再进入拟合”的顺序组织。图像是文档 artifact；示例代码是当前可执行的 public API 入口。

## 1. 建立一个最小 ASGARD 模型

本页示例沿用 `doc/quickstart.md` 第 2 节的完整 `model = Model(...)` 构造。先运行 quickstart 中的模型构造代码，再执行本页后续代码块。

这样做是为了避免在多个教程中复制同一套 90 行构造器。所有 API 字段的含义、可选值和注意事项统一维护在 `doc/public_api.md`；本页只展示如何查询光变、谱、分量和内部状态。

quickstart 的模型包含正向激波同步辐射与 SSC。没有启用反向激波时，`result.rev.sync` 和 `result.rev.ssc` 不应被解释为物理反向激波信号。

## 2. 多频光变

```python
import numpy as np

times = np.logspace(2, 7, 80)
freqs = np.array([1.0e9, 1.0e14, 1.0e18])
result = model.flux_density_grid(times, freqs)
```

输出矩阵 `result.total` 的形状为 `(num_frequency, num_time)`。

![ASGARD 多频光变](assets/tutorials/quick_light_curves.png)

光变应随 \(t_{\rm obs}\) 平滑演化。若某个频段出现孤立尖峰，首先检查动力学、电子输运、频率网格和投影路径，不使用 smoothing 掩盖。

## 3. 宽频谱

```python
times = np.array([1.0e3, 1.0e5, 1.0e7])
freqs = np.logspace(8, 22, 110)
result = model.flux_density_grid(times, freqs, projection_kind="sed")
```

![ASGARD 宽频谱](assets/tutorials/quick_spectra.png)

固定时刻的谱用于检查 \(\nu_m\)、\(\nu_c\)、\(\nu_a\) 的相对位置，以及 SSC 是否进入观测频段。`projection_kind="sed"` 明确表示固定时刻扫频率。

## 4. 辐射分量拆分

```python
times = np.logspace(2, 7, 80)
result = model.flux_density_grid(times, np.array([1.0e14]))

total = result.total[0]
fwd_sync = result.fwd.sync[0]
fwd_ssc = result.fwd.ssc[0]
```

![ASGARD 辐射分量](assets/tutorials/component_breakdown.png)

拆分分量时要先确认对应开关已经启用。例如 `Radiation(..., include_ssc=True)` 控制正向激波 SSC 光子输出；`SolverOptions(ssc_cooling_mode="nakar_y_thomson")` 表示电子冷却方程用 Nakar \(Y\) 参数近似加入 SSC/IC 冷却。

## 5. 内部量演化

```python
details = model.details(1.0e2, 1.0e7)
track = details.fwd

radius = track.radius
gamma_bulk = track.Gamma
B_comv = track.B_comv
```

![ASGARD 内部量演化](assets/tutorials/internal_quantities.png)

内部量检查是光变验收的必要步骤。若 \(\Gamma(R)\)、\(B'(R)\)、电子谱或最终光变不连续，应回到产生该量的数值核，而不是只看最终光变。

## 6. 逐点观测预测

观测数据通常不是完整网格，而是逐点的 \((t_i,\nu_i)\)。此时使用 `flux_density`：

```python
times = np.array([1.0e3, 3.0e4, 1.0e5])
freqs = np.array([1.0e14, 1.0e9, 1.0e18])
points = model.flux_density(times, freqs)
```

`points.total` 是一维数组，与输入点一一对应。

## 7. 频段积分

```python
band_flux = model.flux(
    time_s=1.0e4,
    nu_min_hz=1.0e14,
    nu_max_hz=1.0e18,
    num_points=96,
)
```

该接口对 \(F_\nu\) 做频率积分：

\[
F_{\rm band}(t)
=
\int_{\nu_{\min}}^{\nu_{\max}} F_\nu(t,\nu)\,{\rm d}\nu .
\]

## 8. 模型配置自检

```python
print(model.jet.kind)
print(model.medium.kind)
print(model.observer.z)
print(model.fwd_rad.epsilon_e, model.fwd_rad.epsilon_B, model.fwd_rad.p)
```

拟合或 benchmark 前应记录模型配置。若修改了喷流、介质或辐射开关，必须同时更新图像脚本和文档说明。

## 9. 天图入口

ASGARD 提供 `model.sky_image(t_obs, nu_obs, fov, npixel)`。天图依赖观测投影设置，适合放在专题 benchmark 中验证。新手教程先使用光变、谱和内部量完成最小闭环。

## 10. 磁化反向激波

磁化反向激波不是光变后处理凸起。它会改变 ejecta baryonic mass、反向激波触发、MHD jump compression、区域 3 热状态、ordered field 和磁压焓惯性。完整公式见 `magnetized_rs_dg1d_tutorial.md`。

在 quickstart 的完整构造器中替换反向激波对象：

```python
reverse_shock = ReverseShock(
    enabled=True,
    shell_duration_s=30.0,
    upstream_sigma=0.1,
    include_cross_zone_ic=False,
    include_ssc=True,
)
```

运行后检查：

```python
times = np.logspace(-1, 6, 120)
freqs = np.array([1.0e9, 1.0e14, 1.0e18])
result = model.flux_density_grid(times, freqs)

rs_sync = result.rev.sync
rs_ssc = result.rev.ssc
```

验收时先看 `upstream_sigma -> 0` 是否回到非磁化基线，再看穿越附近的光变和内部量是否连续。

## 11. DG1D 与 fullhide 对照

`dg_1d` 是 FS/RS 共享、需要显式启用的高阶 1D 电子路径。它默认使用 P12 LGL-DG、对数四速度坐标和问题单元正性核；默认 `fullhide_1d` 不变。

在 quickstart 的 `SolverOptions` 中只改：

```python
solver_options = SolverOptions(
    electron_solver="dg_1d",
    dynamics_solver="forward_legacy",
    geometry_projection="sed_legacy",
    electron_photon_coupling="separated",
    ssc_cooling_mode="nakar_y_thomson",
    synchrotron_integration="fixed_grid",
    cooling_kernel="legacy",
    radiation_kernel="legacy",
    structured_backend="fortran_1d",
    patch_sampling="uniform",
    patch_projection="auto",
    patch_sampling_pilot_theta=0,
    patch_sampling_num_times=12,
    patch_sampling_beaming_factor=3.0,
    patch_sampling_beaming_resolution=8.0,
    structured_parallel_mode="outer",
    structured_outer_threads=None,
    structured_inner_threads=None,
    fullhide2d_transport_model="legacy",
    fullhide2d_stochastic_accel_norm=0.0,
    fullhide2d_escape_mode="closed",
)
```

对照测试只改变 `electron_solver`，其它物理参数保持不变。DG 谱峰低能侧若有小幅局部振荡，但活动支撑连续、无零洞、辐射光变平滑，可以接受；若光变或断频出现孤立尖跳，应回到电子输运、冷却项或动力学检查。

## 12. 有限 q-shell 2D 投影

`fullhide_2d` / `charint_2d` 的厚壳层输运主坐标是有限 \(q\)-mass shell。`chi_grid` 是每个 \(q\) cell 的 BM 等效诊断坐标；实际投影几何由 `chi_radius_cm`、`chi_gamma_bulk` 和 `chi_dvolume_weight` 决定。

在 quickstart 的 `Numerics` 和 `SolverOptions` 中改成：

```python
numerics = Numerics(
    num_radius=88,
    num_theta=36,
    num_phi=24,
    num_observer_time=88,
    num_electron_gamma=61,
    num_photon_frequency=72,
    num_chi=12,
    num_threads=4,
    electron_adaptive_substeps=True,
    electron_substep_rtol=0.02,
    electron_substep_min=16,
    electron_substep_max=256,
    initial_radius_cm=1.0e14,
)

solver_options = SolverOptions(
    electron_solver="fullhide_2d",
    dynamics_solver="forward_legacy",
    geometry_projection="chi_eats_2d",
    electron_photon_coupling="separated",
    ssc_cooling_mode="none",
    synchrotron_integration="fixed_grid",
    cooling_kernel="legacy",
    radiation_kernel="legacy",
    structured_backend="fortran_1d",
    patch_sampling="uniform",
    patch_projection="auto",
    patch_sampling_pilot_theta=0,
    patch_sampling_num_times=12,
    patch_sampling_beaming_factor=3.0,
    patch_sampling_beaming_resolution=8.0,
    structured_parallel_mode="outer",
    structured_outer_threads=None,
    structured_inner_threads=None,
    fullhide2d_transport_model="legacy",
    fullhide2d_stochastic_accel_norm=0.0,
    fullhide2d_escape_mode="closed",
)
```

随后使用普通 `Model` 查询：

```python
times = np.geomspace(1.0e2, 1.0e7, 88)
freqs = np.array([1.0e10, 1.0e14, 1.0e18])
result = model.flux_density_grid(times, freqs, projection_kind="lightcurve")
```

`projection_kind="lightcurve"` 才会让 `chi_eats_2d` 替换 FS synch+SSA 的投影；`projection_kind="sed"` 仍走 shell-level SED 插值。若 off-axis 观测角非零，`num_phi` 必须至少为 2。

多观测角 benchmark 可以复用同一个 2D solve state，只重跑 `project_flux`，条件是底层物理参数、时间网格、频率 seed 网格、电子输运状态不变。保留入口是：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/benchmark_theta_j_multiples_magnetic_decay.py'
```

输出的 1D thin shell 与 2D q-shell 差异是有限厚度、磁场衰减和离轴 EATS 几何共同作用的诊断，不应解释为强子、SSC 或 pair cascade 已经 \(q\)-local。

## 13. Prompt 内部激波快照

Prompt internal-shock snapshot 使用 `prompt/` 包，不从 `asgard_core` 顶层导出。最小代码：

```python
import numpy as np

from prompt.eats import EATSNumerics
from prompt.internal_shock import InternalShockNumerics, InternalShockShell, simulate_internal_shock
from prompt.radiation import InternalShockMicrophysics, RadiationNumerics, compute_prompt_observed_flux

slow = InternalShockShell(gamma=200.0, total_energy_iso_erg=4.0e52, duration_s=0.2, sigma=0.01)
fast = InternalShockShell(gamma=600.0, total_energy_iso_erg=4.0e52, duration_s=0.2, sigma=0.03)
microphysics = InternalShockMicrophysics(epsilon_e=0.1, epsilon_b=0.01, electron_index_p=2.3)

solution = simulate_internal_shock(
    slow,
    fast,
    engine_gap_s=0.2,
    redshift=0.5,
    luminosity_distance_cm=1.0e28,
    opening_angle_rad=0.1,
    epsilon_e=microphysics.epsilon_e,
    epsilon_b=microphysics.epsilon_b,
    numerics=InternalShockNumerics(num_branch_steps=64),
)

flux = compute_prompt_observed_flux(
    solution,
    observer_frequency_hz=np.logspace(16.0, 24.0, 64),
    observer_time_s=np.linspace(1.0e-4, 2.0, 128),
    microphysics=microphysics,
    radiation_numerics=RadiationNumerics(num_electron_gamma=121, num_photon_frequency=161, num_threads=4),
    eats_numerics=EATSNumerics(num_theta=32, num_phi=1, num_threads=4),
)
```

完整物理推导、formal plotting 命令和边界见 `prompt_internal_shock_tutorial.md`。当前 prompt snapshot 是诊断工作流，不是余辉 `Model` 的拟合分支。
