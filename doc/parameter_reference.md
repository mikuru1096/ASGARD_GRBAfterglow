# 参数参考

本文列出新手最常用的公开参数。完整构造函数以 `doc/public_api.md` 为准。

## 1. 喷流

| 字段 | 单位 | 常用范围 | 说明 |
| --- | --- | --- | --- |
| `jet.E_iso` | erg | \(10^{50}\) 到 \(10^{54}\) | 各向同性等效能量。 |
| `jet.lf` / `Gamma0` | 无量纲 | \(30\) 到 \(1000\) | 初始 Lorentz 因子。 |
| `jet.theta_j` | rad | \(0.02\) 到 \(0.3\) | top-hat 半开角。 |
| `jet.duration` | s | 任务相关 | 反向激波 shell 宽度相关参数。 |

结构化喷流使用 `GaussianJet`、`PowerLawJet`、`TwoComponentJet` 等构造器。当前 jet spreading 后端未支持，边界见 `doc/public_backend_limits.md`。

## 2. 外部介质

| 构造器 | 参数 | 说明 |
| --- | --- | --- |
| `ISM(n_ism=...)` | \(n_{\rm ISM}\,[{\rm cm^{-3}}]\) | 均匀介质。 |
| `Wind(A_star=..., n_ism=...)` | \(A_\star\)、低密度 floor | wind 只支持 \(k=2\)。 |

自定义 `Medium(rho=...)` 可以在 Python 层计算密度，但当前 Fortran 后端不支持任意介质 kernel dispatch。

## 3. 辐射参数

| 字段 | 常用范围 | 说明 |
| --- | --- | --- |
| `eps_e` | \(10^{-3}\) 到 \(0.5\) | 电子能量分数。 |
| `eps_B` | \(10^{-6}\) 到 \(10^{-1}\) | 磁场能量分数。 |
| `p` | \(2.0\) 到 \(2.8\) | 非热电子谱指数。 |
| `xi_N` | \(10^{-3}\) 到 \(1\) | 参与加速的电子比例。 |
| `ssc` | `True/False` | 是否计算同步自康普顿。 |

同步辐射特征频率通常满足近似标度：

\[
\nu_m \propto \Gamma B\gamma_m^2,
\qquad
\nu_c \propto \Gamma B\gamma_c^2.
\]

这些标度只用于理解趋势，正式结果由 Fortran 辐射核计算。

## 4. 数值设置

| 字段 | 说明 |
| --- | --- |
| `num_r` | 半径壳层数。 |
| `num_tobs` | 观测者时间网格数。 |
| `num_gam_e` | 电子能量网格数。 |
| `num_nu` | 频率网格数。 |
| `electron_solver` | 默认推荐 `fullhide_1d`；2D 使用 `fullhide_2d` 或 `charint_2d`。 |
| `geometry_kernel` | `chi_eats_2d` 只替换 FS synchrotron+SSA 的 lightcurve projection。 |
| `num_threads` | Fortran/OpenMP 线程数。 |

网格加密应服务明确问题：例如峰时定位、谱断点解析、\(\chi\) 投影收敛或强子能量预算。不要为了填满表格做低信息增益扫描。

## 5. 拟合参数路径

常见 `Param` 绑定：

```python
Param("logE", "jet.E_iso", 50.0, 54.0, Scale.LOG10)
Param("logn", "medium.n_ism", -4.0, 2.0, Scale.LOG10)
Param("logepsB", "fwd_rad.eps_B", -6.0, -1.0, Scale.LOG10)
Param("p", "fwd_rad.p", 2.0, 2.8, Scale.LINEAR)
Param("theta_obs", "observer.theta_obs", 0.0, 0.3, Scale.LINEAR)
```

参数范围是观测问题的一部分，不是通用常数。若 posterior 长时间贴边，应重新检查物理模型、单位转换和数据选择。
