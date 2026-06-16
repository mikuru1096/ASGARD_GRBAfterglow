# 公开 API 选择手册

本文回答四个问题：可以选什么、选项是什么意思、选择后模型会发生什么、需要注意什么。所有示例都按当前公开入口写；不要在新脚本里使用旧 alias 或未从 `asgard_core` 导出的内部 helper。

## 1. 导入入口

```python
from asgard_core import (
    Fitter,
    Hadronic,
    Model,
    Numerics,
    Observer,
    ObserverGrid,
    Param,
    Radiation,
    ReverseShock,
    Scale,
    SolverOptions,
    TabulatedMedium,
    UniformMedium,
    WindMedium,
    gaussian_jet,
    power_law_jet,
    top_hat_jet,
    units,
)
```

`from asgard_core import units` 导入的是轻量单位常量模块，便于写可读的时间、角度和频率：

| 常量 | 数值含义 | 例子 |
| --- | --- | --- |
| `units.day`, `units.hour`, `units.year` | 以秒为基准。 | `3.0 * units.day`。 |
| `units.deg`, `units.arcsec`, `units.mas` | 以弧度为基准。 | `5.0 * units.deg`。 |
| `units.GHz`, `units.keV`, `units.GeV` | 以 Hz 为基准。 | `1.0 * units.keV`。 |
| `units.mJy`, `units.uJy` | 以 Jy 为基准。 | 转成 cgs flux density 仍需乘 \(10^{-23}\)。 |

## 2. 一个完整 Model 必须给什么

`Model` 当前应显式给全套对象：

```python
model = Model(
    jet=...,
    medium=...,
    observer=...,
    fwd_rad=...,
    rvs_rad=None,
    numerics=...,
    observer_grid=...,
    solver_options=...,
    reverse_shock=...,
    hadronic=...,
)
```

不要写 `Model(jet, medium, observer, radiation)` 这类短构造作为教程。当前 API 已经要求用户明确数值网格、求解器、反向激波和强子开关。

`rvs_rad=None` 是反向激波辐射的独立 `Radiation` 对象入口。新手和普通公开路径通常保持 `None`，表示反向激波不单独给另一套辐射参数；只有明确需要 FS/RS 不同微物理参数的研究路径才传入第二个 `Radiation`。

## 3. Medium：外部介质怎么选

### 3.1 可选入口

| 入口 | 必填参数 | 物理意义 | 选择后效果 | 注意事项 |
| --- | --- | --- | --- | --- |
| `UniformMedium` | `number_density_cm3` | 常数 ISM 数密度。 | blast wave 在均匀介质中减速。 | 默认教学和第一轮拟合优先选它。 |
| `WindMedium` | `a_star`, `density_floor_cm3`, `density_cap_cm3` | \(n\propto R^{-2}\) 恒星风。 | 早期减速和光变斜率变成 wind-like。 | 当前公开 kernel 只支持 wind \(k=2\)。`density_cap_cm3=None` 表示不设 cap。 |
| `TabulatedMedium` | `radius_cm`, `density_cm3`, `label` | 显式半径-密度表。 | 动力学按表格密度插值。 | 半径必须严格递增，密度必须为正；不是任意 Python callable kernel。 |

内置 WR/Mesler 密度剖面可用于物理情景测试，但拟合前要说明 progenitor 假设。不要把任意 callable `Medium` 写成已支持 public runtime。

### 3.2 选择建议

- 不知道环境时：先用 `UniformMedium`。
- 观测斜率明确支持恒星风时：用 `WindMedium`。
- 有外部模拟或文献密度剖面时：用 `TabulatedMedium`，并保存原始密度表来源。

## 4. Jet：喷流结构怎么选

### 4.1 `top_hat_jet`

```python
jet = top_hat_jet(
    energy_iso_erg=1.0e52,
    initial_lorentz_factor=300.0,
    opening_angle_rad=0.1,
    shell_duration_s=None,
    magnetar=None,
    spreading=False,
)
```

| 参数 | 意思 | 效果 | 注意事项 |
| --- | --- | --- | --- |
| `energy_iso_erg` | isotropic-equivalent 能量。 | 控制整体能量和通量量级。 | 拟合通常用 `Scale.LOG10`。 |
| `initial_lorentz_factor` | 初始 Lorentz 因子。 | 控制早期 coasting 和减速。 | 数据没有早期点时常退化。 |
| `opening_angle_rad` | top-hat 半开角。 | 控制 jet edge 和角向面积。 | off-axis/jet break 数据才强约束。 |
| `shell_duration_s` | ejecta shell duration。 | 与反向激波 shell 宽度相关。 | 不启用反向激波时通常不是第一轮参数；拟合路径用 `jet.shell_duration_s`。 |
| `magnetar` | 能量注入对象。 | 可引入 engine injection。 | 只有明确 plateau/注入假设时打开。 |
| `spreading` | 横向扩张开关字段。 | 当前不要当成已验证 jet spreading。 | 文档和拟合不应依赖它解释光变。 |

### 4.2 `gaussian_jet`

```python
jet = gaussian_jet(
    energy_iso_erg=1.0e52,
    initial_lorentz_factor=300.0,
    core_angle_rad=0.05,
    outer_angle_rad=0.5,
    shell_duration_s=None,
    magnetar=None,
    spreading=False,
)
```

能量和 \(\Gamma_0-1\) 按角度高斯下降。它适合 off-axis 结构化喷流问题。选择后必须检查角向 patch 收敛，特别是 `viewing_angle_rad` 附近的峰时和峰值是否连续。

### 4.3 `power_law_jet`

```python
jet = power_law_jet(
    energy_iso_erg=1.0e52,
    initial_lorentz_factor=300.0,
    core_angle_rad=0.05,
    energy_index=2.0,
    lorentz_index=None,
    outer_angle_rad=0.5,
    shell_duration_s=None,
    magnetar=None,
    spreading=False,
)
```

核内常数，核外按幂律下降。`lorentz_index=None` 时沿用 `energy_index`。大指数会制造更尖的角向边界，需要更认真地检查 patch sampling。

### 4.4 `Magnetar`

```python
magnetar = Magnetar(L0=1.0e47, t0=1.0e5, q=2.0)
```

| 字段 | 意思 | 选择后效果 | 注意事项 |
| --- | --- | --- | --- |
| `L0` | 初始注入 luminosity。 | 给 jet 增加 engine energy-injection 项。 | 只有 plateau 或明确注入假设时使用。 |
| `t0` | 注入特征时间。 | 控制注入衰减开始的时间尺度。 | 要和观测时间单位 s 保持一致。 |
| `q` | 注入衰减指数。 | 控制晚期注入强弱。 | 不要用它替代缺失的介质或喷流结构物理。 |

## 5. Observer：观测者怎么填

```python
observer = Observer(
    z=0.1,
    viewing_angle_rad=0.0,
    viewing_azimuth_rad=0.0,
    luminosity_distance_cm=None,
)
```

| 字段 | 意思 | 效果 | 注意事项 |
| --- | --- | --- | --- |
| `z` | 红移。 | 影响时间膨胀、频率红移和距离。 | 必须显式给。 |
| `viewing_angle_rad` | 视线与喷流轴夹角。 | 控制 Doppler cone、峰时和 off-axis 亮度。 | on-axis 用 0。 |
| `viewing_azimuth_rad` | 视线方位角。 | 结构化/非轴对称投影时有意义。 | axisymmetric 模型通常给 0。 |
| `luminosity_distance_cm` | 光度距离。 | 控制通量归一化。 | 可显式给 `None`，由红移推导。 |

不要写 `Observer(z=0.1, viewing_angle_rad=0.0)` 作为完整示例；当前构造器要求方位角和距离字段也明确。

## 6. Radiation：辐射和微物理怎么选

`Radiation` 必须显式写全当前字段：

```python
fwd_rad = Radiation(
    epsilon_e=0.1,
    epsilon_B=1.0e-3,
    p=2.3,
    proton_energy_fraction=0.0,
    epsilon_b_floor=None,
    magnetic_decay_alpha_t=0.0,
    magnetic_decay_t0_s=1.0,
    accelerated_electron_fraction=0.1,
    thermal_electrons=False,
    include_ssc=True,
    include_kn_correction=False,
    proton_synch=True,
    include_pgamma=False,
    bethe_heitler=False,
    hadronic_inverse_compton=False,
    pp=False,
    neutrino=False,
    acceleration_efficiency=1.0,
    reverse_proton_energy_fraction=0.0,
    pgamma_scheme="disabled",
    pair_production=False,
)
```

### 6.1 电子和磁场

| 字段 | 意思 | 效果 | 注意事项 |
| --- | --- | --- | --- |
| `epsilon_e` | shock 能量给电子的比例。 | 控制 \(\gamma_m\) 和通量。 | 通常 log/linear 都可能退化，先用物理范围。 |
| `epsilon_B` | shock 能量给磁场的比例。 | 控制 \(B'\)、同步辐射、冷却频率。 | 通常用 `Scale.LOG10`。 |
| `p` | 非热电子谱指数。 | 控制谱斜率和衰减斜率。 | 典型要求 \(p>2\)，次级 RS 分支也要求 \(p>2\)。 |
| `accelerated_electron_fraction` | 非热电子数分数 \(\xi_N\)。 | 改变 \(\gamma_m\) 和归一化。 | 与 \(\epsilon_e\)、\(E\)、\(n\) 强退化。 |
| `thermal_electrons` | thermal electron branch。 | 改变电子分布。 | 当前要求 `electron_solver="fullhide_1d"` 或兼容路径。 |
| `epsilon_b_floor` | 磁场 floor。 | 约束磁场衰减下限。 | 只有明确物理理由时设置。 |
| `magnetic_decay_alpha_t`, `magnetic_decay_t0_s` | 磁场衰减参数。 | 影响冷却和谱断点。 | 不是默认 afterglow 拟合第一轮参数。 |

### 6.2 SSC 和冷却

| 字段 | 意思 | 效果 | 注意事项 |
| --- | --- | --- | --- |
| `include_ssc` | 输出 SSC 分量。 | `result.fwd.ssc` 可非零。 | 是否参与冷却还由 `SolverOptions.ssc_cooling_mode` 控制。 |
| `include_kn_correction` | KN 相关开关字段。 | 影响兼容路径。 | 公开拟合优先通过 `ssc_cooling_mode` 说明冷却模式。 |

### 6.3 强子相关

| 字段 | 意思 | 效果 | 前置条件和注意事项 |
| --- | --- | --- | --- |
| `proton_energy_fraction` | shock 能量给质子的比例。 | 打开强子时控制质子注入能量。 | `Hadronic.enabled=False` 时不应解释为强子输出。 |
| `proton_synch` | 质子同步辐射。 | 控制 hadronic proton synch 分量。 | `Hadronic.solver="legacy_1d"` 覆盖 proton transport + proton synch；formal 强子问题用 `am3_1d`。 |
| `include_pgamma` | \(p\gamma\) 过程。 | 进入 formal strong-interaction path。 | 需要 `Hadronic.enabled=True`、`Hadronic.solver="am3_1d"`；joint feedback 还要求 `Hadronic.pgamma_scheme="hummer_2010_response"`。 |
| `bethe_heitler` | BH pair production。 | 产生 secondary \(e^\pm\) 和 photon sink。 | joint feedback 要求 `electron_photon_coupling="joint"`、`Hadronic.solver="am3_1d"`、`pgamma_scheme="hummer_2010_response"`、`ssc_cooling_mode="numeric_ic_kn"` 且固定子步。 |
| `pp` | pp 相互作用。 | 产生 secondary 和 neutrino 相关输出。 | 需要 baryon target density 解释；不要在缺少靶密度物理时打开。 |
| `neutrino` | neutrino 输出。 | 输出逃逸 neutrino luminosity。 | 需要 formal 强子路径；neutrino 不反馈到电子或光子方程。 |
| `pgamma_scheme` | p-gamma 核。 | 选择 `disabled`、`hummer_2010_response`。 | joint feedback 要求 `hummer_2010_response`。 |
| `pair_production` | pair production branch。 | 可进入 \(\gamma\gamma\) pair/synch cascade。 | `pair_cascade_iterations>1` 使用 shell-sequence time-dependent pair/synch cascade；IC-mediated electromagnetic cascade 未实现。 |

## 7. Numerics：网格怎么选

```python
numerics = Numerics(
    num_radius=300,
    num_theta=300,
    num_phi=1,
    num_observer_time=200,
    num_electron_gamma=201,
    num_photon_frequency=201,
    num_chi=None,
    num_threads=8,
    electron_adaptive_substeps=False,
    electron_substep_rtol=0.02,
    electron_substep_min=100,
    electron_substep_max=1000,
    initial_radius_cm=1.0e14,
)
```

| 字段 | 意思 | 增大后的效果 | 注意事项 |
| --- | --- | --- | --- |
| `num_radius` | 半径壳层数。 | 动力学和输运更细。 | 影响冷启动成本。 |
| `num_theta`, `num_phi` | 角向投影网格。 | off-axis/结构化喷流更准。 | on-axis axisymmetric 可用 `num_phi=1`。 |
| `num_observer_time` | 内部观测时间网格。 | 光变插值更细。 | 不等于用户查询点数量。 |
| `num_electron_gamma` | 电子能量网格。 | 谱断点和冷却更稳。 | 太低会产生非物理谱形。 |
| `num_photon_frequency` | photon/辐射频率网格。 | SSC、吸收、强子 target 更稳。 | 强子和 feedback 对它敏感。 |
| `num_chi` | \(\chi\) 网格数。 | 2D 厚壳层解析。 | 只有 2D electron solver / `chi_eats_2d` 需要。 |
| `num_threads` | OpenMP/并行线程。 | 可能加速。 | 小网格可能被线程开销支配。 |
| `electron_adaptive_substeps` | 自适应电子子步。 | 根据误差调步长。 | joint feedback 当前要求固定子步。 |
| `electron_substep_rtol/min/max` | 自适应子步控制。 | 控制误差和步数范围。 | 不用于后处理 smoothing。 |
| `initial_radius_cm` | 初始半径。 | 控制积分起点。 | 早期数据敏感时需要检查。 |

quick grid 用于调试，formal grid 用于最终结果。不要把 quick grid posterior 当最终结论。

## 8. ObserverGrid：内部时间范围

```python
observer_grid = ObserverGrid(time_min_s=1.0e2, time_max_s=1.0e8)
```

它定义内部默认观测时间范围。若数据早于 `time_min_s` 或晚于 `time_max_s`，先扩大这里，而不是让模型在边界外插值。

## 9. SolverOptions：求解器和投影怎么选

```python
solver_options = SolverOptions(
    electron_solver="fullhide_1d",
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

需要临时检查断频时，可传 `nu_callback(label, nu_m, nu_c, nu_a)`；它只接收底层核已经返回的数组，不写入 `details()`、`FitResult` 或常规 runtime state。

### 9.1 电子求解器

| `electron_solver` | 意思 | 选择后效果 | 注意事项 |
| --- | --- | --- | --- |
| `fullhide_1d` | 默认 1D 隐式电子输运。 | 最稳健 public baseline。 | 拟合第一选择。 |
| `fullhide_1d_hz` | hybrid/hz 变体。 | 特定实验路径。 | 不作为教学默认。 |
| `slc1_1d` | semi-Lagrangian family。 | 用于方法比较。 | 结果应与 baseline 物理诊断一致。 |
| `charint_1d` | characteristic integration。 | 检查冷却断点/重映射。 | 非默认拟合路径。 |
| `t2g1_1d` | legacy implicit path。 | 回归比较。 | 不作为新工作默认。 |
| `weno5_1d` | 高阶 WENO5。 | 诊断谱形和数值耗散。 | 需要额外平滑性检查。 |
| `fullhide_2d` | \((\gamma,\chi)\) 电子输运。 | 给厚壳层电子历史和 `chi_eats_2d`。 | 不代表强子/SSC 也 \(\chi\) 局域。 |
| `charint_2d` | 2D characteristic 混合路径。 | 加速/比较 2D 输运。 | 同样只服务已定义的 2D 电子契约。 |
| `fullhide_2d_pic` | 2D PIC/实验路径。 | 特定研究路径。 | 不写入新手拟合默认。 |

### 9.2 投影和耦合

| 字段 | 可选项 | 意思 | 注意事项 |
| --- | --- | --- | --- |
| `dynamics_solver` | `forward_legacy` | 正向激波动力学主链。 | 反向激波由 `ReverseShock` 接入，不是另一个 public dynamics solver。 |
| `geometry_projection` | `sed_legacy`, `chi_eats_2d` | 观测者投影核。 | `chi_eats_2d` 要配 2D electron solver；只替换 FS synch+SSA lightcurve projection。 |
| `electron_photon_coupling` | `separated`, `joint` | 是否做壳层级 electron/photon/hadronic 联合闭合。 | `joint` 只支持 `fullhide_1d`、`am3_1d`、BH enabled、`pgamma_scheme="hummer_2010_response"`、`ssc_cooling_mode="numeric_ic_kn"` 和固定子步。 |
| `ssc_cooling_mode` | `none`, `numeric_ic_kn`, `nakar_y_thomson` | 电子冷却方程中如何加入 IC/SSC 冷却。 | `none` 不把 SSC/IC 写入电子冷却；`numeric_ic_kn` 对 seed photon field 做含 KN/Jones 核的数值 IC 损失积分；`nakar_y_thomson` 用 Nakar 型 \(Y\) 参数近似放大同步冷却。 |
| `synchrotron_integration` | `fixed_grid` | 同步辐射积分核。 | 当前公开只支持这个值。 |
| `cooling_kernel` | `legacy` | 电子冷却核族。 | 目前 public path 只接受 `legacy`。 |
| `radiation_kernel` | `legacy` | 辐射核族。 | 目前 public path 只接受 `legacy`。 |

`include_ssc=True` 决定是否输出 SSC 光子分量；`ssc_cooling_mode` 决定这些光子场是否以及如何反馈到电子冷却。两者不能混为一个开关。

### 9.3 结构化喷流 patch

| 字段 | 可选项 | 效果 | 注意事项 |
| --- | --- | --- | --- |
| `structured_backend` | `fortran_1d`, `python_patch` | 选择结构化喷流 patch 求解后端。 | `fortran_1d` 要求 `electron_solver="fullhide_1d"`。 |
| `patch_sampling` | `uniform`, `dominant_region_ioka_v1`, `dominant_region_ioka_time_v1` | 角向采样策略。 | dominant-region 当前只支持 `structured_backend="python_patch"`。 |
| `patch_projection` | `auto`, `tophat_cell`, `surface_element` | patch 面元投影方式。 | `auto` 会按采样策略选择。 |
| `patch_sampling_pilot_theta` | int | dominant-region pilot 中的 theta 采样控制。 | `uniform` 下不改变物理求解。 |
| `patch_sampling_num_times` | int | time-aware dominant-region 采样使用的时间节点数。 | 只在 `dominant_region_ioka_time_v1` 有决策意义。 |
| `patch_sampling_beaming_factor` | float | dominant-region 采样中 Doppler cone 宽度因子。 | 过小会漏掉有效角区，需做角向收敛。 |
| `patch_sampling_beaming_resolution` | float | beaming 区域解析度因子。 | 增大通常增加采样点和成本。 |
| `structured_parallel_mode` | `outer`, `inner`, `nested` | 结构化喷流并行模式：`outer` 并行 patch，`inner` 并行单个 patch 内核，`nested` 同时分配外层和内层线程。 | `nested` 必须同时设置 `structured_outer_threads` 和 `structured_inner_threads`；线程乘积不能超过 `num_threads` 和机器 CPU 数。 |
| `structured_outer_threads` | int 或 `None` | 外层 patch 并行线程数。 | `None` 表示由 runtime 默认分配。 |
| `structured_inner_threads` | int 或 `None` | 单个 patch 内核线程数。 | 小网格不一定随线程数加速。 |

不要在文档中承诺结构化 Fortran 后端支持所有 RS SSC、cross-zone IC、BH、pp、pair cascade 组合。复杂组合需要看 `doc/public_backend_limits.md` 和验证结果。

### 9.4 2D 附加选项

| 字段 | 可选项 | 意思 | 注意事项 |
| --- | --- | --- | --- |
| `fullhide2d_transport_model` | `legacy`, `pwn_cr_v1` | 2D electron transport 附加物理。 | `pwn_cr_v1` 是专题研究路径，不是默认 afterglow 拟合。 |
| `fullhide2d_stochastic_accel_norm` | float | stochastic acceleration 强度。 | 只有对应物理假设时使用。 |
| `fullhide2d_escape_mode` | `closed`, `free_outer` | \(\chi\) 方向逃逸边界；`closed` 不加外边界逃逸，`free_outer` 把外侧 ghost density 置零并允许最外层向外逃逸 sink。 | `free_outer` 属于 2D/PWN-CR 研究路径，使用前要检查电子数/能量预算和光变平滑性。 |

## 10. ReverseShock：反向激波怎么选

```python
reverse_shock = ReverseShock(
    enabled=False,
    shell_duration_s=10.0,
    upstream_sigma=0.0,
    include_cross_zone_ic=False,
    include_ssc=False,
)
```

| 字段 | 意思 | 打开后的效果 | 注意事项 |
| --- | --- | --- | --- |
| `enabled` | 是否计算反向激波。 | `result.rev.sync` 可有物理贡献。 | 只在早期 excess 或科学问题要求时打开。 |
| `shell_duration_s` | ejecta shell duration。 | 影响 RS crossing 和 shell 宽度。 | 与 jet `shell_duration_s` 语义相关，需保持物理一致。 |
| `upstream_sigma` | 上游磁化 \(\sigma\)。 | 增加 ordered field 和磁化 jump。 | \(\sigma\to0\) 必须回到非磁化基线。 |
| `include_cross_zone_ic` | FS/RS 跨区 IC。 | 增加 `cross_ic`。 | 只在 seed field 和数据需要时打开。 |
| `include_ssc` | RS SSC。 | `result.rev.ssc` 可非零。 | 不等于 RS 强子 full-chain。 |

密度增强触发的次级反向激波由 medium density jump 数组进入动力学，详见 `doc/shock_shell_adaptive_algorithms.md`。

## 11. Hadronic：强子路径怎么选

```python
hadronic = Hadronic(
    enabled=False,
    solver="legacy_1d",
    num_proton_gamma=161,
    num_neutrino_frequency=121,
    pgamma_scheme="disabled",
    pair_cascade_iterations=1,
)
```

| 字段 | 可选项 | 意思 | 注意事项 |
| --- | --- | --- | --- |
| `enabled` | `False`, `True` | 是否运行强子路径。 | 强子目前是壳层级 1D 契约。 |
| `solver` | `legacy_1d`, `am3_1d` | 强子求解器。 | `legacy_1d` 只适合 proton transport + proton synch；p-gamma/neutrino 用 `am3_1d`。 |
| `num_proton_gamma` | int | 质子能量网格。 | formal 强子结果需做网格收敛。 |
| `num_neutrino_frequency` | int | neutrino 频率网格。 | 只影响 neutrino 输出。 |
| `pgamma_scheme` | `disabled`, `hummer_2010_response` | p-gamma 过程核。 | joint feedback 要求 `hummer_2010_response`。 |
| `pair_cascade_iterations` | 正整数 | \(\gamma\gamma\) pair/synch cascade 迭代。 | `>1` 使用 shell-sequence time-dependent pair/synch cascade；IC-mediated cascade 未实现。 |

强子路径不支持 2D/\(\chi\) 局域 hadronic transport。若 `electron_solver` 是 2D，不要把 hadronic 输出解释成 \(\chi\)-local 反馈。

## 12. 查询接口怎么选

| 方法 | 输入 | 输出 | 什么时候用 |
| --- | --- | --- | --- |
| `flux_density_grid(times_s, nu_hz, projection_kind="lightcurve")` | 时间数组和频率数组 | `(num_frequency, num_time)` | 多频光变、画模型曲线。 |
| `flux_density(times_s, nu_hz, projection_kind="lightcurve")` | 同长度时间和频率数组 | 逐点通量 | 观测数据点和 fitting。 |
| `flux_density_grid_adaptive(times_s, nu_hz, projection_kind="lightcurve")` | 查询时间和频率 | 自适应内部时间网格结果 | direct top-hat 且需要自动加密观测时间网格时。 |
| `flux_density_exposures(times_s, nu_hz, exposures_s, num_subsamples=4)` | 曝光中心、频率、曝光时间 | 曝光平均通量 | 长曝光观测；`num_subsamples` 控制曝光内采样。 |
| `spectrum(time_s, nu_hz, projection_kind="sed")` | 一个时刻，多频率 | SED | 固定时刻扫频率。 |
| `flux(time_s, nu_min_hz, nu_max_hz, num_points=64, projection_kind="sed")` | 频段上下界 | 频段积分 flux | X-ray 或高能频段积分。 |
| `details(t_min=None, t_max=None)` | 时间范围 | 内部状态轨道 | 物理验收和 debug。 |
| `sky_image(t_obs, nu_obs, fov, npixel=128)` | 观测时刻、频率、视场、像素数 | `SkyImage` | 天图和 centroid；需要显式 `luminosity_distance_cm`。 |
| `polarization(times_s, nu_hz, magnetic_geometry="shock_random", local_emissivity="analytic_then_kernel")` | 时间、频率、磁场几何和局域发射近似 | `PolarizationResult` | 偏振诊断；不要用经验时间因子修峰时。 |

`projection_kind="lightcurve"` 面向固定频率、多时刻光变和拟合；`projection_kind="sed"` 面向固定时刻扫频率和频段积分。二者共享同一物理 solve state，但投影插值目标不同。

| `projection_kind` | 意思 | 对输出的影响 | 什么时候选 |
| --- | --- | --- | --- |
| `"lightcurve"` | 光变投影。 | `flux_density_grid` 返回 `(num_frequency, num_time)`，插值目标是同一频率随时间演化。 | 光变、拟合、曝光平均默认用它。 |
| `"sed"` | SED 投影。 | 返回形状仍是 `(num_frequency, num_time)`，但插值目标是固定时刻扫频率。 | `spectrum()`、`flux()` 和固定时刻宽频谱用它。 |

## 13. 返回类型怎么读

### 13.1 `FluxResult`

`flux_density_grid`、`flux_density`、`spectrum` 和 `flux` 返回的核心对象都是 `FluxResult` 或可从中取出的数组。`FluxResult.fwd` 和 `FluxResult.rev` 是 `FluxPair`，分别用 `.sync` 和 `.ssc` 保存同步辐射与 SSC 分量。矩阵约定是频率在前、时间在后。

| 字段 | 类型/形状 | 意思 | 什么时候非零或非空 |
| --- | --- | --- | --- |
| `total` | `np.ndarray` | 总通量密度。 | 总是存在。 |
| `fwd.sync` | `np.ndarray` | 正向激波同步辐射。 | 正向激波主分量。 |
| `fwd.ssc` | `np.ndarray` | 正向激波 SSC。 | `Radiation.include_ssc=True` 时可非零。 |
| `rev.sync` | `np.ndarray` | 反向激波同步辐射。 | `ReverseShock.enabled=True` 时才解释。 |
| `rev.ssc` | `np.ndarray` | 反向激波 SSC。 | `ReverseShock.enabled=True` 且 `ReverseShock.include_ssc=True` 时才解释。 |
| `cross_ic` | `np.ndarray` 或 `None` | FS/RS 跨区 IC。 | `ReverseShock.include_cross_zone_ic=True` 时可非空。 |

### 13.2 `TrackBundle` 和 `CharTrack`

`details(t_min=None, t_max=None)` 返回 `TrackBundle`。`details.fwd` 是正向激波轨道，`details.rev` 是反向激波轨道或 `None`，`details.patches` 用于结构化喷流 patch 诊断。

| `CharTrack` 字段 | 意思 | 何时使用 |
| --- | --- | --- |
| `t_obs` | 观测者时间 \([{\rm s}]\)。 | 对齐光变和内部状态。 |
| `radius` | shock 半径 \([{\rm cm}]\)。 | 检查动力学和密度结构。 |
| `Gamma` | bulk Lorentz factor。 | 检查减速是否平滑。 |
| `B_comv` | 共动磁场 \([{\rm G}]\)。 | 检查 \(\epsilon_B\)、冷却和 RS 磁化效应。 |
| `gamma_e`, `dN_dgamma_e` | 电子 Lorentz 因子网格和电子谱。 | 诊断电子输运、冷却和数值振荡。 |
| `N_p`, `Doppler` | 壳层粒子数和 Doppler 因子。 | 诊断归一化和 off-axis 投影。 |

上表列出新手和拟合诊断最常用字段。其余 `CharTrack` 字段服务 hadronic、2D \(\chi\) 分辨、次级反向激波和详细能量预算，只有对应物理模块启用时才非空；完整字段定义见 `asgard_core/api_model.py` 中的 `CharTrack` 类。

### 13.3 其他返回类型

| 类型 | 来源方法 | 关键字段 | 意思 |
| --- | --- | --- | --- |
| `AdaptiveFluxResult` | `flux_density_grid_adaptive(...)` | `time_s`, `flux`, `total`, `fwd`, `rev`, `cross_ic` | 自适应内部观测时间网格和对应 `FluxResult`。 |
| `SkyImage` | `sky_image(...)` | `image`, `extent`, `pixel_solid_angle`, `pixel_size`, `direct_flux`, `rendered_flux`, `normalization_scale`, `x_centroid`, `y_centroid`, `centroid`, `flux_ratio` | 天图、像素尺度、直接通量/渲染通量和质心。 |
| `PolarizationResult` | `polarization(...)` | `I_sync`, `Q`, `U`, `linear_polarization`, `polarization_angle_rad`, `components` | 同步辐射 Stokes 量、线偏振度、偏振角和分量拆分。 |

## 14. Fitter、Param 和采样器

```python
fitter = Fitter(model=model)
fitter.add_flux_density(times_s, nu_hz, flux, flux_err)
fitter.params = [
    Param("logE", "jet.energy_iso_erg", 50.0, 54.0, Scale.LOG10),
    Param("p", "fwd_rad.p", 2.01, 2.8, Scale.LINEAR),
]
```

| 对象 | 可以选什么 | 意思 | 注意事项 |
| --- | --- | --- | --- |
| `Fitter(model, data=None, params=None, num_workers=1)` | 拟合器构造。 | 保存一个确定性 forward model 和观测数据。 | `model` 必须是已构建的 `Model`；`num_workers` 记录 worker 数，当前不自动创建 emcee pool。 |
| `add_flux_density(times_s, frequencies_hz, flux, flux_err)` | 添加逐点通量密度。 | likelihood 加入 \((t_i,\nu_i,F_{\nu,i})\)。 | 单位固定为 s、Hz、cgs flux density。 |
| `add_spectrum(time_s, frequencies_hz, flux, flux_err)` | 添加一个时刻的 SED。 | likelihood 加入固定时刻多频点。 | `time_s` 是单个时刻。 |
| `add_flux(nu_min_hz, nu_max_hz, time_s, flux, flux_err, num_points=64)` | 添加频段积分通量。 | 模型先积分 \(F_\nu\)，再进入 likelihood。 | X-ray/高能 band flux 用这个。 |
| `Param(name, path, lower, upper, scale)` | 显式参数路径。 | 采样器变量写入 `path`。 | 推荐用这一种，路径必须是 `Model` 字段路径。 |
| `Param(name, lower, upper, scale)` | canonical 名字推断路径。 | 只对下方白名单做精确匹配。 | 新手不要依赖推断，优先显式写 path。 |
| `Scale.LINEAR` | 线性 | 采样值直接写入模型。 | 适合 \(p\)、角度等。 |
| `Scale.LOG10` | log10 | 写入模型前做 \(10^x\)。 | 适合能量、密度、\(\epsilon_B\)。 |
| `Scale.FIXED` | 固定 | 不作为 active 参数。 | PyMultiNest 当前不支持 fixed 参数放进 params。 |
| `fit(param_defs=None, sampler="emcee", total_steps=128, burn_frac=0.5, thin=1, nwalkers=None, outputfiles_basename="chains/asgard-", verbose=False)` | 运行采样。 | 返回 `FitResult`。 | `sampler` 只接受 `emcee` 或 `pymultinest`；网格通过 `Numerics` 设置，不在 `fit()` 里临时覆盖。 |
| `sampler="emcee"` | emcee ensemble sampler | 默认 MCMC。 | `nwalkers/total_steps/burn_frac/thin` 生效。 |
| `sampler="pymultinest"` | PyMultiNest nested sampler | 输出 evidence。 | 依赖外部 `pymultinest` 和 MultiNest 库；`total_steps` 等 emcee 参数不生效。 |

自动路径推断只适合少量常用连续参数；正式教程推荐显式写 `Param(name, path, lower, upper, scale)`。当前 `Param(name, lower, upper, scale)` 的精确白名单如下：

| `name` | 等价显式 `path` | 推荐尺度 | 注意事项 |
| --- | --- | --- | --- |
| `energy_iso_erg`, `log10_energy_iso_erg` | `jet.energy_iso_erg` | `Scale.LOG10` | 采样变量是 \(\log_{10}E_{\rm iso}\) 时用 `Scale.LOG10`。 |
| `initial_lorentz_factor`, `log10_initial_lorentz_factor` | `jet.initial_lorentz_factor` | `Scale.LINEAR` 或 `Scale.LOG10` | 早期数据不足时强退化。 |
| `opening_angle_rad` | `jet.opening_angle_rad` | `Scale.LINEAR` | top-hat 半开角。 |
| `core_angle_rad` | top-hat 中为 `jet.opening_angle_rad`，结构化喷流中为 `jet.core_angle_rad` | `Scale.LINEAR` | 为避免歧义，结构化喷流优先显式写 path。 |
| `outer_angle_rad` | `jet.outer_angle_rad` | `Scale.LINEAR` | 结构化喷流外边界。 |
| `energy_index` | `jet.energy_index` | `Scale.LINEAR` | power-law jet 核外能量指数。 |
| `lorentz_index` | `jet.lorentz_index` | `Scale.LINEAR` | power-law jet 核外 Lorentz 因子指数。 |
| `shell_duration_s` | `jet.shell_duration_s` | `Scale.LOG10` 或 `Scale.LINEAR` | 与反向激波 shell 宽度相关。 |
| `viewing_angle_rad` | `observer.viewing_angle_rad` | `Scale.LINEAR` | off-axis 拟合常用参数。 |
| `viewing_azimuth_rad` | `observer.viewing_azimuth_rad` | `Scale.LINEAR` | axisymmetric 模型通常固定为 0。 |
| `z` | `observer.z` | `Scale.LINEAR` | 有观测红移时不应拟合。 |
| `luminosity_distance_cm` | `observer.luminosity_distance_cm` | `Scale.LOG10` | 若由红移推导，通常不拟合。 |
| `epsilon_e` | `fwd_rad.epsilon_e` | `Scale.LOG10` 或 `Scale.LINEAR` | 与能量、密度退化。 |
| `epsilon_B` | `fwd_rad.epsilon_B` | `Scale.LOG10` | 常见拟合参数。 |
| `p` | `fwd_rad.p` | `Scale.LINEAR` | 典型要求 \(p>2\)。 |
| `proton_energy_fraction` | `fwd_rad.proton_energy_fraction` | `Scale.LOG10` 或 `Scale.LINEAR` | 只有强子路径启用时才解释。 |
| `epsilon_b_floor` | `fwd_rad.epsilon_b_floor` | `Scale.LOG10` | 只有明确磁场 floor 假设时使用。 |
| `magnetic_decay_alpha_t` | `fwd_rad.magnetic_decay_alpha_t` | `Scale.LINEAR` | 磁场衰减专题参数。 |
| `magnetic_decay_t0_s` | `fwd_rad.magnetic_decay_t0_s` | `Scale.LOG10` | 磁场衰减专题参数。 |
| `accelerated_electron_fraction` | `fwd_rad.accelerated_electron_fraction` | `Scale.LOG10` 或 `Scale.LINEAR` | 与 \(\epsilon_e\) 强退化。 |
| `acceleration_efficiency` | `fwd_rad.acceleration_efficiency` | `Scale.LOG10` 或 `Scale.LINEAR` | 影响最大能量估计。 |
| `reverse_proton_energy_fraction` | `fwd_rad.reverse_proton_energy_fraction` | `Scale.LOG10` 或 `Scale.LINEAR` | 只在 RS 强子问题中解释。 |
| `number_density_cm3` | `medium.number_density_cm3` | `Scale.LOG10` | `UniformMedium` 常用。 |
| `a_star` | `medium.A_star` | `Scale.LOG10` | `WindMedium` 常用。 |
| `density_floor_cm3` | `medium.density_floor_cm3` | `Scale.LOG10` | wind floor。 |
| `density_cap_cm3` | `medium.density_cap_cm3` | `Scale.LOG10` | `None` 不能作为采样值，拟合时需给数值上界。 |
| `reverse_shock_enabled` | `reverse_shock.enabled` | `Scale.FIXED` | 布尔开关不建议连续采样；按模型假设固定。 |
| `reverse_shock_shell_duration_s` | `reverse_shock.shell_duration_s` | `Scale.LOG10` | RS crossing 相关。 |
| `reverse_shock_upstream_sigma` | `reverse_shock.upstream_sigma` | `Scale.LOG10` 或 `Scale.LINEAR` | \(\sigma\to0\) 应回到非磁化基线。 |
| `hadronic_num_proton_gamma` | `hadronic.num_proton_gamma` | `Scale.FIXED` | 网格数不应作为物理 posterior 参数。 |
| `hadronic_num_neutrino_frequency` | `hadronic.num_neutrino_frequency` | `Scale.FIXED` | 网格数用于收敛检查。 |
| `hadronic_pair_cascade_iterations` | `hadronic.pair_cascade_iterations` | `Scale.FIXED` | 算法设置，不是连续物理参数。 |

完整拟合教程见 `doc/fitting_workflow.md`；emcee/PyMultiNest 对比见 `doc/mcmc_fitting.md`。

## 15. 不要承诺的组合

以下内容不是当前 public API 的稳定承诺：

- 任意 callable medium kernel。
- wind \(k\ne2\)。
- 已验证 jet spreading。
- 2D/\(\chi\) 分辨 hadronic transport。
- 结构化 Fortran 后端覆盖所有 RS SSC、cross-zone IC、BH、pp、pair cascade 组合。
- 未构建 Fortran 扩展时依赖 Python fallback 继续拟合。

遇到不支持组合，应直接暴露错误并回到物理/算法契约，不写兜底代码或后处理修正。
