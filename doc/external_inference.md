# 外部采样接口

ASGARD 当前没有内置 Redback 适配层。若需要与 Redback、bilby、BlackJAX、DIME 或其他采样器结合，应把 ASGARD 作为确定性 forward model 调用，而不是在外部采样器中重新实现物理。

## 1. 最小包装函数

```python
import numpy as np

from asgard_core import Fitter

fitter = Fitter(model=model)
fitter.add_flux_density(times_s, nu_hz, flux, flux_err)
fitter.params = params

def log_likelihood(theta):
    values = {
        "logE": theta[0],
        "logn": theta[1],
        "logepsB": theta[2],
        "p": theta[3],
    }
    return fitter.loglike(values)
```

外部采样器只负责产生 `theta` 和记录 posterior。参数写入、模型求解和 likelihood 计算仍由 `Fitter` 管理。

## 2. 先验变换

若采样器使用单位立方体 \(u_i\in[0,1]\)，先验变换应显式写出：

\[
x_i
=
x_{i,\min}
+
u_i\left(x_{i,\max}-x_{i,\min}\right).
\]

对 `Scale.LOG10` 参数，采样变量 \(x_i\) 是 \(\log_{10}\) 空间中的值，写入模型时由 `Param.transform` 转成 \(10^{x_i}\)。

## 3. 数据边界

外部观测文件进入 ASGARD 前完成：

- 时间统一为 s。
- 频率统一为 Hz。
- flux density 统一为 cgs。
- 缺测值、非有限值和文件列名错误在读入边界暴露。

ASGARD 内部不为不可能发生的内部状态添加防御式兜底。

## 4. 与 Redback 的当前边界

Redback 页面可以作为“外部采样器如何组织模型、先验和数据”的参考，但 ASGARD 仓库当前没有正式 `redback` model 名称、参数表或注册函数。因此文档不提供伪接口。需要 Redback 支持时，应先实现一个明确的适配层，并加入最小端到端验证。

## 5. 正式采样记录

外部采样正式运行时记录：

- ASGARD git HEAD。
- 数据文件路径和 hash。
- Fortran 构建命令。
- 参数定义和先验边界。
- 采样器名称、版本、随机种子和并行设置。
- best-fit 光变、posterior、trace 或 SMC 诊断图。

这些记录属于拟合产品的一部分，不应在事后凭记忆补写。
