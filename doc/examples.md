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
nu_m = track.nu_m
nu_c = track.nu_c
nu_a = track.nu_a
```

![ASGARD 内部量演化](assets/tutorials/internal_quantities.png)

内部量检查是光变验收的必要步骤。若 \(\Gamma(R)\)、\(B'(R)\)、\(\nu_m(R)\)、\(\nu_c(R)\) 或 \(\nu_a(R)\) 不连续，应回到产生该量的数值核，而不是只看最终光变。

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
