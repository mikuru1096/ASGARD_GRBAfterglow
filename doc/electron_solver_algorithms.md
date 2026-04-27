# ASGARD 电子求解器算法说明

描述当前工作树中前向激波电子求解器的方程、离散和数值特性。对应源文件在 `src/Electron/`。

## 1. 统一物理问题

### 1.1 1D 变量

所有 1D 求解器在 `x = log10(gamma_e)` 上求解。主守恒变量：
- `dN_x = dN / dlog10(gamma_e) = gamma_e * ln(10) * dN/dgamma_e`

原因：GRB afterglow 电子谱跨多个数量级；源项和冷却断点在对数空间更自然；特征线/半拉格朗日/隐式迎风都更容易在 x 空间写成守恒更新。

### 1.2 2D 变量

2D 求解器展开下游流体结构：
- `x = log10(gamma_e)`, `eta = log10(chi)`, `chi = 1 + 8 * Gamma_sh^2 * x_downstream / R`
- 守恒变量: `U(log10(gamma_e), log10(chi), R) = dN / (dlog10(gamma_e) dlog10(chi))`

### 1.3 共享微观物理

所有求解器每个壳层共享：
- 外部介质密度 `dNe`（ISM/wind/密度跳变/展宽）
- 磁场 `DB = 0.39 * sqrt(epsilon_b * dNe * (Gamma * (Gamma - 1)))`
- `gamma_max = 3 * m_e c^2 / sqrt(8 B e^3)`
- `gamma_m` 由 `electron_gamma_m_exact` 计算（p>2, 1<p<2, p=2 三种分支）
- `gamma_c` 解析估计用于壳层诊断；2D 另有 `electron_gamma_c_from_loss_mean` 从实际损失系数反推
- 注入源项：`electron_injection_prefactor` + `electron_build_source_term_exp_cutoff`（γ^{-p} × 指数截断）

### 1.4 共享冷却模型

由 `electron_cooling_kernel.f90` 组装 `dEl(gamma)`（每单位半径的能量损失系数）：
- 同步辐射损失 + SSC/IC 损失 + SSA 回热修正
- 绝热项通过传输系数中的 `1/R` 项进入

`index_Y` 分支：
- 0: 仅同步辐射主项 `dEl = (f_r - dot_gamma_ssa * scale) * gamma`
- 1: 数值 IC 冷却辅助量
- 2: Nakar 型 Y 参数
- 3: Fan 型 Y 参数

### 1.5 共享辐射输出

主要输出: `gam_e`, `dN_gam_e`, `P_syn`, `Seed_syn`, `nu_m`, `nu_c`, `nu_a`。`weno5_1d` 除外（Fortran 只返回 `gam_e/dN/P_syn/Seed_syn`，Python 运行时再算 `nu_m/nu_c/nu_a`）。

## 2. 统一数值骨架

### 网格

所有求解器调用 `electron_initialize_spectrum` 构造对数 `gamma_e` 网格。`charint_1d/slc1_1d` 还显式构造单元边界 `electron_profile_log_cell_edges`。

### 源项处理

三种写法：
- 中心点近似 `electron_build_source_term_exp_cutoff`
- 单元边界积分 `electron_build_source_term_exp_cutoff_edges`
- 特征线法：固定谱型预处理为线性重构系数，四点求积加入

差异在：是否按单元平均、高能截止是否平滑、与特征线时间积分是否同阶。

### 冷却系数均值

多数非 WENO 求解器构造 `dEL_mean`：相邻单元中心 `dEl` 算术平均除以 `ln(10)`，变成 `x = log10(gamma)` 空间上的面速度。

### 子步策略

四类：固定子步、冷却 CFL 控制子步、全步/半步 Richardson 自适应子步、2D 中 ξ 与 η 双方向最小步长。

## 3. 各算法

### 3.1 `fullhide_1d`

一阶隐式迎风格式。核心更新 `electron_fullhide_step`:
- 速度系数: `face_coeff = dEL_mean + 1/(R ln 10)`
- 单下对角三对角系统隐式更新

优点: 稳定性强，1D 最稳妥基准。缺点: 数值扩散偏强，高能端截止被抹平。

### 3.2 `slc1_1d`

Strang splitting 半拉格朗日: 半步源项 → 半拉格朗日输运 → 半步源项。调用 `electron_semi_lagrangian_step`。

优点: 比一阶隐式更接近"沿流线搬运"。缺点: 大量固定子步(100-1000)，精度靠步数堆出来。

### 3.3 `charint_1d`

转 `u = 1/gamma` 空间沿特征线回溯单元边界，守恒 remap。`index_Y=0` 用仿射特征线，`index_Y≠0` 用分段仿射。

处理: 回溯每个单元边界 → 非均匀网格保守 remap → 四点求积沿特征线积分源项。

优点: 1D 高能截止和冷却断点通常最锋利。缺点: 实现复杂，子步大小对精度极敏感。

### 3.4 `t2g1_1d`

三层时间推进（BDF2），启动退回单步法。`右端 ≈ 2*dN^n - 0.5*dN^{n-1} + source`。

优点: 时间精度高于一阶。缺点: 三层法对系数快速变化敏感，startup 误差会传播。

### 3.5 `weno5_1d`

显式守恒律: WENO5 空间重构 + SSP RK3 时间推进。通量 `fp = dEl1 * dN_x`，按速度符号切换 `fpx/fmx`。常值外推 ghost cells。正性维护: `where(dN_x < 0) dN_x = 0`。

优点: 平滑区高阶。缺点: 显式 CFL 限制，刚性问题步长极小，截断负值影响守恒。

### 3.6 `fullhide_2d`

χ 方向展开后激波区: χ-dependent 下游速度/磁场/历史光子场/SSA/pair opacity。

数值: η 方向隐式平流-扩散-源项三对角; ξ 方向一阶隐式迎风。子步 `dDR_try = min(dDR_xi, dDR_eta, dDD)`。冷却用 reduced 6-band 网格 (`Num_nu_cool = min(6, Num_nu)`)。

优点: 最完整 2D 物理基线。缺点: 代价大，ξ 一阶隐式抹平高能端。

### 3.7 `charint_2d`

当前公开实现: η 隐式 + ξ 特征线 remap（不是全方向特征线）。共享 `fullhide_2d` 外层物理/历史/诊断。

零源项 ξ 步跳过多点求积；`active_chi_hi` 只推进电子占据明显的 χ 列。子步上限 `L1 ≤ 512`。

优点: 比纯隐式更能保高能截止，比 `fullhide_2d` 快。缺点: η 还未切到特征线方案；ξ 子步限制过紧会伤高能端；不是完整 2D characteristic solver。

## 4. 数值对比

| 求解器 | 能量方向方法 | χ 方向 | 子步策略 | 当前角色 |
| --- | --- | --- | --- | --- |
| `fullhide_1d` | 一阶隐式迎风 | 无 | 固定/自适应 | 1D 稳定基线 |
| `slc1_1d` | 半拉格朗日+Strang | 无 | 固定大量子步 | 方法比较 |
| `charint_1d` | 特征线保守 remap | 无 | 固定/自适应 | 1D 高保真方法 |
| `t2g1_1d` | 三层二阶隐式 | 无 | 固定大量子步 | 时间推进比较 |
| `weno5_1d` | WENO5+SSP RK3 | 无 | CFL 显式 | 高阶显式比较 |
| `fullhide_2d` | 一阶隐式迎风 | 隐式平流-扩散 | 双方向 CFL | 2D 物理基线 |
| `charint_2d` | 特征线 remap | 隐式平流-扩散 | ξ CFL+上限 | 2D 加速混合版 |

## 5. 连续方程与离散

### 5.1 1D 守恒方程

`N(gamma, R) = dN/dgamma`, `gamma' = dgamma/dR`:

```text
∂N/∂R + ∂/∂gamma [ gamma' * N ] = Q
```

变到 `x = log10(gamma)`, `dN_x = dN/dx`:

```text
∂dN_x/∂R + ∂/∂x [ A_x * dN_x ] = S_x
A_x = gamma' / (gamma * ln(10)), S_x = gamma * ln(10) * Q
```

代码中 `face_coeff ≈ dEL_mean + 1/(R ln 10)`，第二项来自球对称膨胀的几何项。

### 5.2 1D 离散

- `fullhide_1d`: 一阶隐式有限体积，迎风通量，三对角 backward sweep
- `slc1_1d`: operator splitting 半拉格朗日，沿能量流线回溯 + 保守重映射 + 源项 split
- `charint_1d`: `u = 1/gamma` 空间，回溯单元边界 → 非均匀网格保守 remap → 四点求积源项积分
- `weno5_1d`: 显式守恒律 `∂U/∂R + ∂F(U)/∂x = S`，WENO5 重构，SSP RK3 时间推进

### 5.3 2D 连续方程

`U(x, eta, R) = dN/(dx deta)`:

```text
∂U/∂R + ∂(A_x U)/∂x + ∂(A_eta U)/∂eta = ∂(D_eta ∂U/∂eta)/∂eta + S
```

x 方向是典型一阶平流型；η 方向是平流+扩散+局域源混合型。

### 5.4 2D 离散

`fullhide_2d`:
- η: 隐式平流-扩散-源项三对角 (`advance_eta_logchi_implicit`)
- ξ: 一阶隐式迎风每 χ 列 (`advance_energy_loggamma_chi`)

`charint_2d` (当前公开):
- η: 与 `fullhide_2d` 同
- ξ: 特征线保守 remap (`advance_energy_loggamma_chi_charint`)

`charint_2d` 的 η 方程有扩散项，扩散不能被写成单纯特征线 remap，因此当前设计是对的。

## 6. PDE 项与代码例程映射

| PDE 项 | 负责代码 |
| --- | --- |
| γ 网格 | `electron_common.f90`: `electron_initialize_spectrum`, `electron_profile_log_cell_edges` |
| χ/η 几何 | `electron_transport_2d_kernel.f90`: `compute_log_chi_geometry`, `compute_downstream_comoving_grid` |
| 注入 S | `electron_injection_profiles.f90`: `electron_build_source_term_exp_cutoff` / `_edges` |
| 冷却 A_x | `electron_cooling_kernel.f90` 组装 `dEl`; `electron_transport_common.f90` 转为面系数 |
| 1D ξ 隐式 | `electron_transport_common.f90`: `electron_fullhide_step` |
| 1D 半拉格朗日 | `electron_transport_common.f90`: `electron_semi_lagrangian_step` |
| 1D 特征线 | `electron_transport_common.f90`: `electron_characteristic_update` |
| 2D η 输运 | `electron_transport_2d_kernel.f90`: `advance_eta_logchi_implicit` |
| 2D ξ 隐式 | `electron_transport_2d_kernel.f90`: `advance_energy_loggamma_chi` |
| 2D ξ 特征线 | `electron_transport_2d_kernel.f90`: `advance_energy_loggamma_chi_charint` |
| 辐射谱 | `electron_radiation_kernel.f90`: `get_syn_state`, `get_syn_adaptive`, `get_nu_a` |
| 历史场 | `electron_seed_history_kernel.f90`: `accumulate_comoving_history_fields` |
| 诊断量 | `nu_m`: Υ_m 直接换算; `nu_c`: `electron_gamma_c_from_loss_mean` 反解; `nu_a`: `get_nu_a_2d_cell_path` |

## 7. 关键风险

- **`charint_2d` ξ 子步粗细**: 高能端是否过早截断的第一嫌疑项
- **2D reduced cooling bands**: `Num_nu_cool = min(6, Num_nu)` 是明确近似，窄频带敏感时 `dEl_chi` 可能系统偏差
- **诊断量一致性**: `nu_m/nu_c/nu_a` 必须与主输运同一套局域状态，否则时间不平滑

共性 bug 信号: `nu_m/nu_c/nu_a` 随时间跳变、相邻时刻总谱锯齿、电子总数突变、参数平滑变化时高频端台阶式缩短。
