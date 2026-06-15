# MCMC 拟合专题

本文专门说明 ASGARD 中两个内置采样入口：`emcee` 和 `pymultinest`。完整的数据准备和模型构建见 `doc/fitting_workflow.md`。

## 1. 先理解采样变量

`Param` 定义的是采样器变量，不一定是直接写入模型的物理量：

```python
Param("logE", "jet.energy_iso_erg", 50.0, 54.0, Scale.LOG10)
Param("p", "fwd_rad.p", 2.01, 2.8, Scale.LINEAR)
```

对 `Scale.LOG10`，采样器看到的是

\[
x=\log_{10}E_{\rm iso},
\]

写入模型的是

\[
E_{\rm iso}=10^x.
\]

对 `Scale.LINEAR`，采样值直接写入模型。`result.samples` 存储的也是采样变量，所以 `logE` 样本仍是 \(\log_{10}E\)，不是 \(E\)。

## 2. emcee：默认起步选择

### 2.1 什么时候用

用 `emcee` 当第一轮采样器：

- 想快速看参数范围是否合理。
- 想看 posterior 是否贴边、是否有明显退化。
- 模型维度不高，likelihood 单次计算还能接受。
- 不需要 Bayesian evidence。

### 2.2 最小运行

```python
result = fitter.fit(
    sampler="emcee",
    total_steps=512,
    burn_frac=0.5,
    thin=2,
    nwalkers=24,
)
```

`sampler` 只接受 `emcee` 和 `pymultinest`；其它字符串会直接报错。

### 2.3 参数怎么选

| 参数 | 意思 | 效果 | 注意事项 |
| --- | --- | --- | --- |
| `nwalkers` | walker 数量 | 越多越能探索退化方向，但每步计算更多模型。 | 至少 \(2N_{\rm dim}+2\)。初学者可用 \(4N_{\rm dim}\) 到 \(8N_{\rm dim}\)。 |
| `total_steps` | 每个 walker 步数 | 决定链长。 | quick check 用 128 到 512；正式结果不能只靠这个数字判断收敛。 |
| `burn_frac` | 丢弃前多少比例 flattened samples | 去掉初始过渡段。 | trace 未稳定时，增大 burn 不能让结果自动可信。 |
| `thin` | 抽稀因子 | 减少输出样本数。 | 只影响保存样本，不解决相关性。 |

### 2.4 初始化和边界

当前实现只对 active 参数采样。固定参数不会进入 emcee 维度。若参数范围为 \([a,b]\)，walker 初值在

\[
a+0.2(b-a)
\le x \le
a+0.8(b-a)
\]

中均匀生成。采样时若

\[
x<a
\quad{\rm or}\quad
x>b,
\]

log probability 返回 \(-\infty\)。

这意味着边界就是先验边界。边界应该来自物理假设和数据约束能力，而不是为了让链“好看”。

### 2.5 判断 emcee 是否可用

至少检查：

- 每个 walker 的 trace 是否离开初始区间并混合。
- `best_params` 是否贴着边界。
- posterior 是否出现明显多峰。
- 不同 quick grid 与 formal grid 的结论是否一致。
- best-fit 物理轨道是否平滑。

如果 posterior 多峰、退化很强，或者需要 evidence，再用 PyMultiNest。

## 3. PyMultiNest：evidence 和多峰 posterior

### 3.1 什么时候用

用 `pymultinest` 的场景：

- 需要比较模型证据 \(Z\)。
- 怀疑 posterior 多峰。
- 参数退化导致 emcee 混合慢。
- 想要 nested sampling 的全局探索，而不只是局部链。

不需要 evidence 时，不要因为“更正式”就默认用 MultiNest。它更重，运行环境也更挑剔。

### 3.2 最小运行

先确认 Python 包和底层动态库能加载：

```python
try:
    from pymultinest.solve import solve
except ImportError as exc:
    raise RuntimeError(
        "当前环境不能导入 pymultinest；先安装 pymultinest 并确认 MultiNest 动态库在加载路径中。"
    ) from exc
else:
    print("pymultinest OK")
```

这个检查只验证导入；正式运行前还要确认 `outputfiles_basename` 指向的目录可写。

```python
result = fitter.fit(
    sampler="pymultinest",
    outputfiles_basename="chains/asgard-",
    verbose=True,
)
```

`sampler` 只接受 `emcee` 和 `pymultinest`。PyMultiNest 路径必须写 `sampler="pymultinest"`。

当前实现内部调用：

```python
from pymultinest.solve import solve
```

所以运行环境必须同时满足：

- Python 能 import `pymultinest`。
- 底层 MultiNest 库能被动态链接器找到。
- `outputfiles_basename` 指向的目录可写。

### 3.3 先验如何定义

`sampler="pymultinest"` 使用单位 cube 到参数边界的线性映射：

\[
x_i
=
x_{i,\min}
+
u_i(x_{i,\max}-x_{i,\min}),
\qquad
u_i\in[0,1].
\]

如果参数是 `Scale.LOG10`，这里的 \(x_i\) 是 log 参数。例如

\[
x=\log_{10}E_{\rm iso}\in[50,54]
\]

对应物理先验

\[
E_{\rm iso}\in[10^{50},10^{54}]\,{\rm erg},
\]

且在 \(\log_{10}E\) 上均匀。

### 3.4 当前限制

PyMultiNest 路径当前要求所有 `Param` 都是 active 参数。下面这种固定参数会触发错误：

```python
Param("p", "fwd_rad.p", 2.3, 2.3, Scale.FIXED)
```

若要固定 \(p\)，不要把它放进 `fitter.params`。错误写法和正确写法如下：

```python
# 错误：PyMultiNest 当前不接受 fixed Param
fitter.params = [
    Param("logE", "jet.energy_iso_erg", 50.0, 54.0, Scale.LOG10),
    Param("p", "fwd_rad.p", 2.3, 2.3, Scale.FIXED),
]
```

```python
# 正确：先把固定值写进 model，只把 active 参数交给采样器
model.fwd_rad.p = 2.3
fitter = Fitter(model=model)
fitter.params = [
    Param("logE", "jet.energy_iso_erg", 50.0, 54.0, Scale.LOG10),
]
```

`total_steps`、`burn_frac`、`thin` 和 `nwalkers` 不控制 PyMultiNest。它们只属于 emcee。

### 3.5 结果字段

PyMultiNest 的 `FitResult` 额外包含：

```python
print(result.logz)
print(result.logzerr)
```

\(\log Z\) 只能在同一数据、同一 likelihood、同一先验定义下比较。若两个运行的参数范围不同，evidence 差异同时包含先验体积差异。

比较两个模型 A 和 B 时，先算

\[
\Delta\ln Z = \ln Z_A-\ln Z_B,
\qquad
K_{A/B}=\exp(\Delta\ln Z).
\]

```python
delta_logz = result_A.logz - result_B.logz
bayes_factor = np.exp(delta_logz)
print(delta_logz, bayes_factor)
```

经验读法：\(K\gtrsim3\) 是有一定支持，\(K\gtrsim10\) 是强支持，\(K\gtrsim100\) 才能称为很强支持。这个判断只在两个运行使用同一数据、同一 likelihood 和可比较先验时成立。

## 4. emcee 和 PyMultiNest 怎么选

| 目标 | 推荐 | 原因 |
| --- | --- | --- |
| 第一次跑通拟合 | emcee | 接口简单，容易看 trace 和 corner。 |
| 调参数范围 | emcee | 快速暴露边界和单位问题。 |
| 低维、单峰 posterior | emcee | 足够直接。 |
| 多峰 posterior | PyMultiNest | nested sampling 更适合全局探索。 |
| 比较 top-hat vs structured jet evidence | PyMultiNest | 需要 \(\log Z\)。 |
| 只要 best-fit 曲线 | emcee 或手动优化外部工具 | MultiNest evidence 不必要。 |

## 5. quick grid 和 formal grid

拟合通常分两档：

| 档位 | 目的 | 典型设置 |
| --- | --- | --- |
| quick | 检查单位、边界、退化方向、明显 bug。 | `num_radius=72`, `num_observer_time=72`, `num_electron_gamma=81`, `num_photon_frequency=81`。 |
| formal | 生成最终 posterior、best-fit 图和物理诊断。 | 增加半径、时间、电子和频率网格，记录完整命令。 |

quick 结果只能用于决策，不应作为论文级 posterior。formal 结果必须记录：

- git HEAD。
- `git status --short --branch`。
- 文档和代码 diff 摘要。
- Fortran 构建命令。
- 采样器、参数范围、数据版本。
- 输出目录。
- 物理验收图和结论。

## 6. 拟合失败时先查什么

| 症状 | 优先检查 |
| --- | --- |
| `loglike = -inf` 或 NaN | 数据误差、参数越界、扩展未构建。 |
| 模型比数据差很多数量级 | flux 单位、距离、红移、频率单位。 |
| posterior 贴边 | 边界太窄、模型缺少物理分量、单位错误。 |
| emcee trace 不混合 | 参数退化强、walker 太少、模型太慢、初值区域不合适。 |
| MultiNest evidence 不稳定 | live point 设置需要通过外部 PyMultiNest 配置扩展；当前 `Fitter` 只暴露最小 solve 参数。 |
| 光变有孤立尖峰 | 动力学/输运/投影 bug，不能用 smoothing 修。 |

## 7. 输出不要怎么处理

- 不用 posterior smoothing 让光变变平滑。
- 不用改误差条掩盖模型缺项。
- 不用低信息增益的全组合 ablation 填表。
- 不把 cached warm query 时间当作完整 cold likelihood 时间。
- 不把 `best_params` 当成物理结论；必须结合 posterior、残差和物理轨道。
