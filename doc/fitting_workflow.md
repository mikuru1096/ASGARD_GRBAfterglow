# 拟合工作流

公开拟合层把观测数据编译成固定的查询计划，并在每组参数下重建模型、计算预测和 Gaussian log-likelihood。先确认单次物理模型可靠，再启动采样。

## 1. 准备数据

flux-density 数据需要配对的一维数组：

```python
import numpy as np

time_s = np.array([1e3, 3e3, 1e4, 3e4])
nu_hz = np.array([3e9, 3e9, 5e14, 2.42e17])
flux = np.array([0.12, 0.20, 0.08, 0.015])
flux_err = np.array([0.01, 0.02, 0.01, 0.003])
```

时间、频率、flux 和误差必须逐点对应。单位必须与公开模型输出一致；不要在 likelihood 内混用 mJy、Jy 和 cgs。

固定时刻频谱的数据结构是一个时间加频率数组。频段数据则由时间、上下限频率、积分 flux 与误差组成。

## 2. 建立并验证基线模型

先按 [快速开始](quickstart.md) 创建 `model`，随后直接运行观测网格：

```python
baseline = model.flux_density(time_s, nu_hz)
assert baseline.shape == flux.shape
assert np.all(np.isfinite(baseline))
assert np.all(baseline >= 0.0)
```

采样前至少确认：

- 时间与频率覆盖观测范围；
- 网格对峰时、峰值和谱形已收敛；
- 光变随时间和参数连续；
- 可选 RS/hadronic 分量确实来自所选 backend；
- 单次运行成本适合计划的样本数。

采样器不能修复错误的物理配置或未收敛数值核。

## 3. 添加观测

```python
from asgard_core import Fitter

fitter = Fitter(model, num_workers=1)
fitter.add_flux_density(time_s, nu_hz, flux, flux_err)
```

还可添加：

```python
fitter.add_spectrum(
    time_s=1e4,
    frequencies_hz=np.logspace(9, 18, 12),
    flux=sed_flux,
    flux_err=sed_err,
)

fitter.add_flux(
    nu_min_hz=0.3 * 2.417989e17,
    nu_max_hz=10.0 * 2.417989e17,
    time_s=1e4,
    flux=xray_flux,
    flux_err=xray_err,
    num_points=96,
)
```

多个数据块会合并到同一 observation plan；相同查询点可被复用。

## 4. 定义参数

```python
from asgard_core import Param, Scale

params = [
    Param("energy", "jet.energy_iso_erg", 51.0, 54.0, Scale.LOG10),
    Param("density", "medium.number_density_cm3", -4.0, 2.0, Scale.LOG10),
    Param("epsilon_e", "fwd_rad.epsilon_e", -3.0, -0.3, Scale.LOG10),
    Param("epsilon_B", "fwd_rad.epsilon_B", -6.0, -0.3, Scale.LOG10),
    Param("p", "fwd_rad.p", 2.01, 2.8, Scale.LINEAR),
]
```

`LOG10` 表示采样变量在给定区间线性变化、写入模型前取 `10**value`。`FIXED` 使用 lower 值。显式路径比名称推断更清楚。

常用可推断名称包括：

| 名称 | 路径 |
| --- | --- |
| `energy_iso_erg` | `jet.energy_iso_erg` |
| `initial_lorentz_factor` | `jet.initial_lorentz_factor` |
| `viewing_angle_rad` | `observer.viewing_angle_rad` |
| `epsilon_e`, `epsilon_B`, `p` | `fwd_rad.*` |
| `number_density_cm3` | `medium.number_density_cm3` |
| `reverse_shock_upstream_sigma` | `reverse_shock.upstream_sigma` |

无法推断的名称必须传 `path`。selector、网格规模及改变后端拓扑的字段通常不适合作为连续采样参数。

## 5. 先检查 likelihood

```python
fitter.params = params

trial = {
    "energy": 52.0,
    "density": 0.0,
    "epsilon_e": -1.0,
    "epsilon_B": -3.0,
    "p": 2.3,
}
print(fitter.loglike(trial))
```

在参数上下界、中点和一个物理上可解释的偏移点分别运行。log-likelihood 应有限，并对小参数变化连续。若单点失败，先修模型或参数范围，不在 likelihood 外包异常捕获或罚函数。

可以显式检查重建模型：

```python
trial_model = fitter.build_model(trial)
trial_flux = trial_model.flux_density(time_s, nu_hz)
```

## 6. emcee

```python
result = fitter.fit(
    params,
    sampler="emcee",
    total_steps=512,
    burn_frac=0.5,
    thin=2,
    nwalkers=32,
)

print(result.best_params)
print(result.best_loglike)
print(result.samples.shape)
```

默认 walker 初值在各 prior 内确定性生成。正式采样的步数和 walker 数需由自相关、链混合与重复运行决定；短链只能 smoke test。

固定参数可写为：

```python
Param("p", "fwd_rad.p", 2.3, 2.3, Scale.FIXED)
```

## 7. PyMultiNest

```python
result = fitter.fit(
    params,
    sampler="pymultinest",
    outputfiles_basename="chains/asgard-",
    verbose=True,
)
print(result.logz, result.logzerr)
```

PyMultiNest 需要系统库和 Python binding，且当前公开实现要求所有参数 active；不要混入 `FIXED` 参数。evidence 只有在先验、likelihood 与模型族都明确时才可比较。

## 8. 读取结果

`FitResult` 提供：

| 字段 | 内容 |
| --- | --- |
| `best_params` | 最佳参数字典 |
| `best_loglike` | 最大 log-likelihood |
| `samples` | 后验样本 |
| `log_prob` | 样本 log-likelihood/probability |
| `logz`, `logzerr` | nested-sampling evidence |
| `labels` | 样本列标签 |
| `top_k_params` | 候选最佳解 |

用最佳参数重建模型，而不是直接把采样空间中的 log10 值写入物理字段：

```python
best_model = fitter.build_model(result.best_params)
best_curve = best_model.flux_density_grid(
    np.logspace(2, 7, 200),
    np.unique(nu_hz),
)
```

## 9. 后验物理检查

- 绘制每个波段的数据、最佳解和可信区间。
- 检查 posterior 是否贴住先验边界或存在未探索多峰。
- 对代表性 posterior 样本复查网格收敛，而非只检查最佳解。
- 检查能量分数、Lorentz factor、密度、RS 磁化和强子预算是否物理可行。
- 检查模型随时间、频率和参数连续，避免 sampler 跨越数值跳变。

只有当简模型在残差中留下有物理结构的系统偏差，才增加 RS、结构化喷流、能量注入或强子过程。复杂模型必须回答明确假设，不能仅靠额外自由度降低残差。

## 10. 性能与复现

固定 Git commit、`uv.lock`、线程、网格、随机种子和采样器配置。记录单次 likelihood 的至少三次 median，再估算总成本。并行所有权必须唯一：若外层 walker 并行，通常不要让每个模型继续占满所有 OpenMP 线程。

正式产物应保存参数定义、数据单位、模型配置、结果数组和重建命令。临时链与 profile 放 `/tmp`；可发表链只有在复现与物理验收后进入 artifact 目录。
