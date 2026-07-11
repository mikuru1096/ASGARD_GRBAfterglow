# 参数参考

本文列出新手最常用的公开参数。完整构造函数以 `doc/public_api.md` 为准。

## 1. 喷流

| 字段 | 单位 | 常用范围 | 说明 |
| --- | --- | --- | --- |
| `jet.energy_iso_erg` | erg | \(10^{50}\) 到 \(10^{54}\) | 各向同性等效能量。 |
| `jet.initial_lorentz_factor` | 无量纲 | \(30\) 到 \(1000\) | 初始 Lorentz 因子。 |
| `jet.opening_angle_rad` | rad | \(0.02\) 到 \(0.3\) | top-hat 半开角。 |
| `jet.shell_duration_s` | s | 任务相关 | 反向激波 shell 宽度相关参数；这是推荐 `Param` 路径。 |

公开结构化喷流使用 `gaussian_jet(...)`、`power_law_jet(...)` 或 `tabulated_angular_jet(...)` 构造器。当前 jet spreading 后端未支持，边界见 `doc/public_backend_limits.md`。内部研究用喷流 helper 不属于稳定 public API，新文档和教程不以这些名称作为用户入口。

## 2. 外部介质

| 构造器 | 参数 | 说明 |
| --- | --- | --- |
| `UniformMedium(number_density_cm3=...)` | \(n_{\rm ISM}\,[{\rm cm^{-3}}]\) | 均匀介质。 |
| `WindMedium(a_star=..., density_floor_cm3=...)` | \(A_\star\)、低密度 floor | wind 只支持 \(k=2\)。 |
| `TabulatedMedium(radius_cm=..., density_cm3=...)` | 显式 \(n(R)\) 表 | 用对数插值描述 density shell / WR bubble。 |

自定义 `Medium(rho=...)` 可以在 Python 层计算密度，但当前 Fortran 后端不支持任意介质 kernel dispatch。

密度剖面必须表达清楚物理来源。平滑 bump 或 shell 可以产生真实的再增亮，但没有密度结构时，光变中的孤立尖峰优先视为数值/投影 bug。

## 3. 辐射参数

| 字段 | 常用范围 | 说明 |
| --- | --- | --- |
| `epsilon_e` | \(10^{-3}\) 到 \(0.5\) | 电子能量分数。 |
| `epsilon_B` | \(10^{-6}\) 到 \(10^{-1}\) | 磁场能量分数。 |
| `p` | \(2.0\) 到 \(2.8\) | 非热电子谱指数。 |
| `accelerated_electron_fraction` | \(10^{-3}\) 到 \(1\) | 参与加速的电子比例。 |
| `include_ssc` | `True/False` | 是否计算同步自康普顿。 |
| `proton_energy_fraction` | \(0\) 到问题相关上限 | 正向激波质子能量分数，非零时才有强子能量预算。 |
| `include_pgamma` | `True/False` | 是否计算 p-gamma 过程。 |
| `bethe_heitler` | `True/False` | 是否计算 BH pair production。 |
| `pp` | `True/False` | 是否计算 pp 二级产物。 |
| `neutrino` | `True/False` | 是否输出中微子谱；中微子不反馈电磁输运。 |
| `pgamma_scheme` | `disabled`, `hummer_2010_response` | formal p-gamma feedback 使用 `hummer_2010_response`。 |
| `pair_production` | `True/False` | 是否启用 gamma-gamma pair branch。 |

同步辐射特征频率通常满足近似标度：

\[
\nu_m \propto \Gamma B\gamma_m^2,
\qquad
\nu_c \propto \Gamma B\gamma_c^2.
\]

这些标度只用于理解趋势，正式结果由 Fortran 辐射核计算。

## 4. 数值设置

| 字段 | 可以选什么 | 效果 | 需要注意什么 |
| --- | --- | --- | --- |
| `num_radius` | 正整数 | 控制半径壳层数。 | 影响动力学、输运和投影的时间解析度；正式图需要做收敛检查。 |
| `num_observer_time` | 正整数 | 控制内部观测者时间网格。 | 过低会让峰时和谱断点插值粗糙。 |
| `num_electron_gamma` | 正整数 | 控制电子能量网格。 | 过低会抹平冷却断点；过高会增加运行时间。 |
| `num_photon_frequency` | 正整数 | 控制光子频率网格。 | SSC、\(\gamma\gamma\) 和强子 target photon 对它敏感。 |
| `electron_solver` | `fullhide_1d`, `fullhide_1d_hz`, `slc1_1d`, `charint_1d`, `dg_1d`, `t2g1_1d`, `weno5_1d`, `fullhide_2d`, `charint_2d` | 选择电子输运核。 | 新手先用 `fullhide_1d`；`dg_1d` 是 FS/RS 共享的 opt-in 高阶 1D 路径；2D solver 只在需要 finite \(q\)-shell 厚度结构时使用，`*_hz` 不是普通拟合默认。`fullhide_2d_pic` 没有跟踪源码和构建登记，运行时映射已删除。 |
| `geometry_projection` | `sed_legacy`, `sed_adaptive_theta`, `chi_eats_2d` | 选择观测者投影核。 | `sed_adaptive_theta` 对壳层级 EATS 做 θ 方向自适应积分；`chi_eats_2d` 要配 2D electron solver，只替换 FS synchrotron+SSA lightcurve projection。 |
| `num_threads` | 正整数 | 控制 Fortran/OpenMP 线程数。 | 拟合时并行层级要和外部采样器并行统一规划。 |
| `electron_photon_coupling` | `separated`, `joint` | `separated` 是默认后处理语义；`joint` 做 shell-level 含时二级反馈闭合。 | `joint` 约束更强，需要 `ssc_cooling_mode="numeric_ic_kn"`、formal hadronic 和相容开关。 |
| `ssc_cooling_mode` | `none`, `numeric_ic_kn`, `nakar_y_thomson` | 控制电子冷却方程中的 IC/SSC 项。 | `numeric_ic_kn` 是含 KN/Jones 核的数值 IC 损失积分；`nakar_y_thomson` 是 Nakar \(Y\) 参数 Thomson 近似；`include_ssc=True` 只控制 SSC 光子输出。 |
| `patch_sampling` | `uniform`, `dominant_region_ioka_v1`, `dominant_region_ioka_time_v1` | 控制结构化喷流角向 patch 的采样策略。 | dominant-region 当前只支持 `structured_backend="python_patch"`。 |
| `structured_adaptive_rtol` / `structured_adaptive_max_depth` | 非负数 / 非负整数 | 控制 axisymmetric structured `chi_eats_2d` theta-ring 自适应窗口。 | 只重排 observer projection ring 采样；不改变动力学、电子输运、强子源项或 seed photon 场。 |
| `fullhide2d_transport_model` | `legacy`, `pwn_cr_v1` | 选择 2D transport 研究路径。 | 普通 afterglow 使用 `legacy`。 |
| `projection_adaptive_rtol` / `projection_adaptive_max_depth` | 正数 / 非负整数 | 控制 `sed_adaptive_theta` 的角向积分误差阈值和最大递归深度。 | 只控制 observer projection quadrature，不改变动力学、电子输运或强子源项。 |
| `pair_cascade_iterations` | 正整数 | `1` 为低层诊断路径；`>1` 为 shell-sequence 含时 pair/synch cascade。 | IC-mediated cascade 边界仍按 TODO 记录，不要过度解释为完整 EM cascade。 |
| `pp_gamma_model` | `delta`, `geant4`, `sibyll`, `qgsjet`, `pythia8` | 选择 pp 的 π⁰ gamma 谱；默认 `delta`。 | 详细模型只改变 gamma 分量，其他 pp 二级产物仍用现有 delta 链。 |

网格加密应服务明确问题：例如峰时定位、谱断点解析、\(\chi\) 投影收敛或强子能量预算。不要为了填满表格做低信息增益扫描。

## 5. 物理开关组合

| 目标 | 推荐最小组合 | 物理图像 | 验收重点 |
| --- | --- | --- | --- |
| 标准 FS afterglow | `top_hat_jet`, `UniformMedium`, `electron_solver="fullhide_1d"` | 爆波扫掠外介质，非热电子同步/SSC 辐射，经 EATS 投影成光变。 | \(\Gamma(R)\)、\(\nu_m,\nu_c,\nu_a\) 和光变平滑。 |
| wind afterglow | `WindMedium`, `forward_legacy` | 外介质 \(n\propto R^{-2}\)，减速标度比 ISM 更浅。 | \(\Gamma(t)\) 是否接近 wind 标度，早期 SSA 是否合理。 |
| 结构化喷流 off-axis | `gaussian_jet`、`power_law_jet` 或 `tabulated_angular_jet`, `structured_backend="fortran_1d"` | 每个角向 patch 有不同能量和 \(\Gamma_0\)，观测峰由 Doppler cone 进入视线决定。 | 角向采样收敛、峰时随 `viewing_angle_rad` 连续变化。 |
| 壳层级自适应 EATS | `geometry_projection="sed_adaptive_theta"` | 普通 shell-level SED 投影在 θ 方向用嵌套中点规则估计积分误差并细分。 | 适合检查 off-axis / 窄 beaming 的角向积分收敛；φ 方向仍由 `num_phi` 控制。 |
| 反向激波 | `ReverseShock(enabled=True, upstream_sigma=...)` | ejecta 被反向激波加热后形成 region 3，同步辐射叠加到 `rev.sync`；有限 `upstream_sigma` 同时改变 baryonic mass、MHD jump、ordered field 和磁压焓惯性。 | `upstream_sigma -> 0` 回到非磁化基线，crossing 前后状态连续。 |
| formal 强子 | `Hadronic(enabled=True, solver="am3_1d", pgamma_scheme="hummer_2010_response")` | 质子输运与 p-gamma/BH/pp 二级产物在 shell-level 计算。 | proton loss、secondary pair、photon survival 能量预算一致。 |
| 联合二级反馈 | `electron_photon_coupling="joint"` 加 BH、`am3_1d`、`ssc_cooling_mode="numeric_ic_kn"` | 主电子、强子二级 \(e^\pm\)、光子 target/sink 在同一 \(R\) 网格闭合。 | 弱反馈回到 separated baseline；强反馈时谱和光变仍连续。 |
| 有限厚度投影 | `electron_solver="fullhide_2d"`, `geometry_projection="chi_eats_2d"` | FS synchrotron+SSA 在 finite \(q\)-shell 上做 EATS 投影。 | `num_chi` 收敛；不要把强子解释成 chi-local。 |

## 6. 拟合参数路径

常见 `Param` 绑定：

```python
Param("logE", "jet.energy_iso_erg", 50.0, 54.0, Scale.LOG10)
Param("logn", "medium.number_density_cm3", -4.0, 2.0, Scale.LOG10)
Param("logepsB", "fwd_rad.epsilon_B", -6.0, -1.0, Scale.LOG10)
Param("p", "fwd_rad.p", 2.0, 2.8, Scale.LINEAR)
Param("viewing_angle_rad", "observer.viewing_angle_rad", 0.0, 0.3, Scale.LINEAR)
```

参数范围是观测问题的一部分，不是通用常数。若 posterior 长时间贴边，应重新检查物理模型、单位转换和数据选择。
