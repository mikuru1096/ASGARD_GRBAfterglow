# ASGARD 电子求解器算法说明

本文档描述当前工作树中前向激波电子求解器的实际实现状态，重点回答四件事：

1. 各个求解器到底在解什么方程。
2. 共享了哪些物理模块，哪些地方彼此不同。
3. 各方法在数值离散上分别做了什么。
4. 每个方法当前的优缺点和适用场景是什么。

本文以当前源码实现为准，而不是以历史 README、旧版计划或理想化算法目标为准。对应核心文件主要包括：

- `src/Electron/FS_electron_fullhide_1d.f90`
- `src/Electron/FS_electron_slc1_1d.f90`
- `src/Electron/FS_electron_charint_1d.f90`
- `src/Electron/FS_electron_t2g1_1d.f90`
- `src/Electron/FS_electron_weno5_1d.f90`
- `src/Electron/FS_electron_fullhide_2d.f90`
- `src/Electron/FS_electron_charint_2d.f90`
- `src/Electron/electron_common.f90`
- `src/Electron/electron_transport_2d_kernel.f90`
- `src/Electron/electron_cooling_kernel.f90`
- `src/Electron/electron_radiation_kernel.f90`
- `src/Electron/electron_seed_history_kernel.f90`

## 1. 统一物理问题

### 1.1 1D 变量

1D 求解器统一在电子洛伦兹因子方向上求解电子分布。代码中的主守恒变量通常不是直接的 `dN/dgamma`，而是

- `x = log10(gamma_e)`
- `dN_x = dN / dlog10(gamma_e) = gamma_e * ln(10) * dN/dgamma_e`

这样做的原因是：

- GRB afterglow 电子谱跨很多数量级，线性 `gamma` 网格不合适。
- 源项和冷却断点本来就在对数空间里更自然。
- 特征线法、半拉格朗日法和隐式迎风法都更容易在 `x` 空间中写成守恒更新。

### 1.2 2D 变量

2D 求解器把下游流体内部结构也显式展开。当前代码采用：

- `x = log10(gamma_e)`
- `eta = log10(chi)`
- `chi = 1 + 8 * Gamma_sh^2 * x_downstream / R`

内部守恒变量定义为：

- `U(log10(gamma_e), log10(chi), R) = dN / (dlog10(gamma_e) dlog10(chi))`

这对应 `src/Electron/FS_electron_fullhide_2d.f90` 和 `src/Electron/FS_electron_charint_2d.f90` 文件开头的注释。

### 1.3 共享微观物理

所有求解器在每个壳层上都会共享以下物理量构造：

- 外部介质密度 `dNe`
  - 由 `electron_initial_density` 和 `electron_external_density` 计算。
  - 统一支持 ISM、wind、密度跳变与展宽参数。
- 磁场
  - `DB = 0.39 * sqrt(epsilon_b * dNe * (Gamma * (Gamma - 1)))`
- 最大电子能量 `gamma_max`
  - `gamma_max = 3 * m_e c^2 / sqrt(8 B e^3)` 这种代码常数化形式
- 最低注入能量 `gamma_m`
  - 当前主用 `electron_gamma_m_exact`
  - 在 `p > 2`、`1 < p < 2`、`p = 2` 三种情况下分支不同
- 冷却断点 `gamma_c`
  - 解析估计用于壳层诊断
  - 2D 中另有 `electron_gamma_c_from_loss_mean` 从实际损失系数反推诊断值
- 注入源项
  - 注入前因子由 `electron_injection_prefactor` 给出
  - 代码形式为一个壳层体积增量乘以密度与加速电子比例
- 源项谱型
  - 当前主用 `electron_build_source_term_exp_cutoff`
  - 即 `gamma^{-p}` 乘高能指数截止

### 1.4 共享冷却模型

冷却统一由 `electron_cooling_kernel.f90` 组装。核心对象是 `dEl(gamma)`，它表示每单位半径上的能量损失系数。其来源包括：

- 同步辐射损失
- SSC/IC 损失
- SSA 回热修正
- 绝热项不是在 `dEl` 中直接加成到同一个表达式里，而是大多通过传输系数中的 `1/R` 项进入

当前 `index_Y` 分支含义可以直接按代码理解为：

- `index_Y = 0`
  - 只做同步辐射主项，`dEl = (f_r - dot_gamma_ssa * scale) * gamma`
- `index_Y = 1`
  - 使用数值 IC 冷却辅助量 `cooling_aux`
- `index_Y = 2`
  - Nakar 型 `Y`
- `index_Y = 3`
  - Fan 型 `Y`

### 1.5 共享辐射输出

所有主要求解器最终都会输出：

- `gam_e`
- `dN_gam_e`
- `P_syn`
- `Seed_syn`
- `nu_m`
- `nu_c`
- `nu_a`

例外是 `weno5_1d` 的 Fortran 接口本身只返回 `gam_e/dN/P_syn/Seed_syn`，Python 运行时再额外计算 `nu_m/nu_c/nu_a`。

## 2. 统一数值骨架

### 2.1 网格

所有主流求解器都调用 `electron_build_gamma_grid` 构造对数 `gamma_e` 网格。`charint_1d/slc1_1d` 这类方法还会显式构造单元边界：

- `electron_log_cell_edges`

边界网格很重要，因为：

- 特征线法需要在单元面回溯
- 保守 remap 需要旧单元和新单元的边界交叠
- 指数截止源项如果只在中心点采样会更容易截断失真

### 2.2 源项处理

目前代码中常见的三种源项写法：

- 中心点近似 `electron_build_source_term_exp_cutoff`
- 基于单元边界积分的 `electron_build_source_term_exp_cutoff_edges`
- 特征线法中先把固定谱型预处理成线性重构系数，再按四点求积加入

这三种写法的差异主要在：

- 是否按单元平均而不是点值加源
- 是否对高能截止附近做更平滑积分
- 是否与特征线时间积分保持同阶

### 2.3 冷却系数均值

多数非 WENO 求解器都先构造单元面上的平均损失系数：

- `electron_loss_mean`

形式是相邻单元中心 `dEl` 的算术平均，再除以 `ln(10)`。这样做的含义是把 `dEl` 变成 `x = log10(gamma)` 空间上的面速度。

### 2.4 子步策略

当前代码库中有四类子步策略：

- 单步或固定子步
- 由冷却 CFL 控制的固定子步
- 基于全步/半步误差比较的自适应子步
- 2D 中按 `xi` 与 `eta` 两个方向的最大系数取最小步长

这件事对谱形影响极大，特别是：

- 高能截止的平滑程度
- `nu_m/nu_c` 的连续性
- 2D `charint_2d` 的速度/精度权衡

## 3. 各算法逐个说明

## 3.1 `fullhide_1d`

### 数学形式

`fullhide_1d` 本质上是 1D 对数能量空间上的隐式迎风格式。实际更新函数在 `electron_fullhide_step` 中，可写成：

- 速度系数：`face_coeff = dEL_mean + 1 / (R * ln 10)`
- 之后用单下对角三对角系统做隐式更新

代码路径：

- `FS_electron_fullhide_1d.f90`
- `electron_fullhide_step`
- `electron_prepare_implicit_coeffs`
- `electron_backward_sweep`

### 数值处理

- 守恒变量：`dN_x`
- `xi` 方向：一阶隐式迎风
- 源项：显式加到右端
- 子步：
  - 可固定子步
  - 也可用全步/半步 Richardson 误差估计做自适应
- 非均匀外介质：
  - 子步中重新采样密度并按 `dNe_mid/dNe_shell` 重缩放损失系数

### 物理特点

- 壳层辐射和冷却每个主壳层都重新计算
- `gamma_m/gamma_c/gamma_max` 使用当前壳层条件
- 绝热项通过 `+ 1 / (R ln 10)` 进入能量通量

### 优点

- 稳定性强，是当前最稳妥的 1D 基准方法之一。
- 对强冷却区、密度突变和长时间积分比较稳健。
- 自适应子步逻辑明确，适合做参考解或回归基线。

### 缺点

- 数值扩散偏强，特别是高能端截止和尖锐冷却断点会被抹平。
- 一阶隐式迎风在光滑解上不够精细。
- 想要更尖锐的谱形，通常需要更多网格点或更多子步。

## 3.2 `slc1_1d`

### 数学形式

`slc1_1d` 调用 `electron_semi_lagrangian_step`，对应：

1. 半步源项
2. 半拉格朗日输运
3. 半步源项

这是一个典型 Strang splitting 风格的半拉格朗日离散。

### 数值处理

- 守恒变量：`dN_x`
- 传输：`electron_semi_lagrangian_transport`
- 源项：显式对称分裂
- 子步数：
  - 当前代码非常保守，固定在 `100` 到 `1000` 范围内
- 无自适应误差控制

### 物理特点

- 仍沿用与 `fullhide_1d` 相同的微观物理
- 仍是 1D、单区、壳层平均模型

### 优点

- 比纯一阶隐式迎风更接近“沿流线搬运”的直觉。
- 对平移类问题的形状保持通常好于全隐式迎风。

### 缺点

- 当前实现依赖大量固定子步，代价不小。
- 精度主要靠“步数堆出来”，不够自适应。
- 当前在代码体系中更像一个比较方法，而不是主参考算法。

## 3.3 `charint_1d`

### 数学形式

`charint_1d` 是当前 1D 中最有“方法论独立性”的算法。它不是在 `x = log10(gamma)` 上直接做迎风，而是转到：

- `u = 1 / gamma`

然后沿特征线回溯单元边界，再做保守 remap。

### 数值处理

对于 `index_Y = 0`：

- 冷却方程在 `u` 空间是仿射的
- 调用 `electron_characteristic_step_affine_u_prepared_source`

对于 `index_Y != 0`：

- `dEl(gamma)` 不再是简单仿射
- 先在每个单元内构造 `u dgamma/dR` 的分段仿射逼近
- 调用 `electron_characteristic_step_piecewise_u_prepared_source`

源项处理：

- 先把固定谱型源项做成保守 remap 预处理
- 再用 4 点求积在子步内积分

当前代码常数：

- `charint_quad_order = 4`
- `charint_cfl_relax = 32`

### 物理特点

- 对同步-only 分支和 IC/KN 分支分别用不同特征线推进器
- 保守性比简单中心点源项要好
- 对高能指数截止更敏感，也更容易保真

### 优点

- 1D 高能截止、冷却断点和谱峰附近通常最锋利。
- 物理上更接近“沿特征线搬运分布函数”的结构。
- 对源项和输运的守恒处理更细。

### 缺点

- 实现复杂度明显更高。
- 子步大小对精度非常敏感，过粗会直接伤高能端。
- 在复杂冷却律下的局部分段仿射近似需要小心验证。

## 3.4 `t2g1_1d`

### 数学形式

`t2g1_1d` 是三层时间推进格式，目标是得到二阶时间精度。代码里：

- 启动步或前两步使用单步隐式格式
- 之后使用三层公式：
  - 新步右端约为 `2 * dN^n - 0.5 * dN^{n-1} + source`

### 数值处理

- 守恒变量：`dN_x`
- `xi` 方向：仍是隐式三对角体系
- 时间离散：三层二阶
- 子步：当前实现固定在 `100` 到 `1000`

### 物理特点

- 与 `fullhide_1d` 共用同一套冷却和辐射模块
- 主要差别在时间推进格式，而不是微观物理

### 优点

- 理论上时间精度高于一阶隐式。
- 适合检验“误差来自时间离散还是来自物理建模”。

### 缺点

- 三层法在系数快速变化时更容易受启动步和历史步误差影响。
- 仍然需要很多子步，代码里的默认选择并不省。
- 当前更像方法试验件，不是主参考求解器。

## 3.5 `weno5_1d`

### 数学形式

`weno5_1d` 在 `x = log10(gamma)` 上直接把方程写成显式守恒律，然后采用：

- 五阶 WENO 空间重构
- 三阶 SSP Runge-Kutta 时间推进

### 数值处理

- 守恒变量：`dN_x`
- 通量：`fp = dEl1 * dN_x`
- 按局部速度符号切换 `fpx` 或 `fmx`
- 边界：
  - 当前实现用常值外推 ghost cells
- 正性：
  - 每步后直接 `where(dN_x < 0) dN_x = 0`

### 物理特点

- 冷却项和源项仍沿用同一套物理模块
- 但数值核心是显式高阶有限体积，而不是隐式/半拉格朗日/特征线

### 优点

- 对平滑解有高阶空间精度潜力。
- 在不受 CFL 限制时能给出很干净的波形传播。

### 缺点

- 显式格式对 CFL 很敏感，冷却强时步长很小。
- 通过截断负值维持正性，这会影响严格守恒与高阶性。
- 当前在整个项目里不像 `fullhide_1d/charint_1d` 那样是主工作流核心。

## 3.6 `fullhide_2d`

### 物理模型

`fullhide_2d` 把后激波区沿 `chi` 方向展开。它不再把整个壳层当成一个均匀区，而是显式计算：

- `chi` 几何
- 每个 `chi` 单元的下游速度
- 每个 `chi` 单元的本地磁场
- 每个 `chi` 单元的历史光子场
- 每个 `chi` 单元的 SSA 与 pair opacity

当前外层框架中，以下量都是按 `chi` 分辨率保存的：

- `P_hist`
- `Seed_hist`
- `Tau_hist`
- `P_hist_cool`
- `Seed_hist_cool`
- `Tau_prop_hist_cool`

### 当前数值结构

当前公开入口 `fs_electron_fullhide_2d` 传给核心的 `use_charint_transport = .false.`，因此公开默认行为是：

- `eta = log10(chi)` 方向
  - 使用单步隐式三对角求解
  - 把平流、扩散、注入同时装配到一个系统里
- `xi = log10(gamma)` 方向
  - 使用与 `fullhide_1d` 同构的一阶隐式更新

子步策略：

- `dDR_try = min(dDR_xi, dDR_eta, dDD)`
- 其中 `dDR_xi` 由最大能量方向系数控制
- `dDR_eta` 由最大 `chi` 方向系数控制

冷却代价控制：

- 当前 `Num_nu_cool = min(6, Num_nu)`
- 即 2D 历史冷却不是用全频段，而是投影到 6 个对数冷却频带
- 每个壳层只装配一次 cooling auxiliary
- 子步中不重复做完整 cooling assembly

### 物理优点

- 是当前最完整的 2D 物理描述。
- 历史种子光子场、SSA、pair opacity、`chi` 依赖磁场都显式保留。
- 作为 2D 基线非常重要。

### 数值优点

- 结构清楚，稳定性强。
- `eta` 和 `xi` 都走隐式路线时很稳。
- 适合作为 `charint_2d` 的对照参考。

### 缺点

- 代价大，主要瓶颈通常不在单个 PDE 更新，而在
  - `syn_state`
  - `prepare_aux`
  - `chi` 分辨历史辐射装配
- `xi` 一阶隐式会额外抹平高能端。
- `eta` 隐式把平流与扩散混在一起，稳定但更耗散。

## 3.7 `charint_2d`

### 设计目标

`charint_2d` 的目标不是重新写一套 2D 外层物理，而是在保留 `fullhide_2d` 外层框架的前提下，把更适合做特征线的部分替换为特征线推进。

### 当前真实状态

当前代码中的开关是：

- `use_charint_eta = .false.`
- `use_charint_xi = .true.`

这意味着当前公开 `charint_2d` 不是“`eta` 和 `xi` 全部特征线化”的版本，而是一个混合版本：

- `eta` 方向仍然使用 `fullhide_2d` 的隐式更新
- `xi` 方向改用特征线保守 remap

这是当前工作树的真实状态，必须和理想计划区分开。

### 当前 `xi` 数值处理

当前 `xi` 更新调用：

- `advance_energy_loggamma_chi_charint`

其中：

- `index_Y = 0`
  - 用 `u = 1/gamma` 的仿射特征线推进
- `index_Y != 0`
  - 用分段仿射 `u` 特征线推进

并且当前已经做了两个关键优化：

1. 零源项 `xi` 步不再做四点求积，而是只做一次特征回溯加保守 remap。
2. `active_chi_hi` 只推进电子占据明显的 `chi` 列，尾部空列不再无意义地更新。

### 当前子步控制

当前实现中：

- `dDR_xi = 0.4 * d_x_E / max_xi_coeff`
- 因为公开状态下 `eta` 仍是隐式，所以 `dDR_try = min(dDD, dDR_xi)`
- 之后又加入了一个上限：
  - `L1 = min(L1, 512)`

这个上限是当前速度优化的重要来源，但也是精度折中的主要来源。

### 当前优点

- 共享 `fullhide_2d` 的全部外层物理，不需要再复制第二套历史光子场逻辑。
- 比纯隐式 `xi` 更有机会保留高能截止形状。
- 在当前 quick benchmark 里，已经可以比 `fullhide_2d` 更快。

### 当前缺点

- 早期高能端仍可能低于 `fullhide_2d`，尤其在最早时刻更明显。
- 因为 `eta` 还没有切到特征线方案，所以它还不是完全意义上的“2D characteristic solver”。
- 子步上限 `L1 <= 512` 本质上是速度/精度折中，过紧会伤高能端，过松会让速度掉回去。

### 当前应如何理解

当前 `charint_2d` 最准确的描述是：

- 它是一个“共享 `fullhide_2d` 外层物理、在 `xi` 方向使用 characteristic remap 的 2D 混合求解器”
- 它还不是一个数学上把 2D 全部输运项都特征线化的终版

### 当前本地 quick 对比印象

以当前工作树、`num_gam_e = 81`、`num_chi = 8`、ISM、`ssc + kn` 的 quick case 为例，本地最近一次对比大致表现为：

- `charint_2d` 速度已经快于 `fullhide_2d`
- 但最早时刻的高频端仍偏低
- 越到后期，两者高频端越接近

这说明当前版本已经不是“又慢又错”的状态，但它仍然不是一个可以完全替代 `fullhide_2d` 作为 2D 物理基准的终版。

## 4. 物理处理异同

### 4.1 相同点

以下内容所有主求解器共享：

- 同一套动力学输入：`R_Tobs, R_Gamma, R`
- 同一套微观参数：`epsilon_e, epsilon_b, p, f_e`
- 同一套 `gamma_m/gamma_max/gamma_c` 构造
- 同一套同步、SSC、SSA、pair opacity 基本模块
- 同一套前向激波注入逻辑

### 4.2 不同点的本质

真正决定求解器差异的不是“物理公式是否不同”，而是以下三类数值选择：

- 状态变量用什么
  - `dN/dgamma`
  - `dN/dloggamma`
  - `dN/(dloggamma dlogchi)`
- 传输怎么做
  - 隐式迎风
  - 半拉格朗日
  - 特征线 remap
  - WENO5 显式守恒律
- 源项怎么并入
  - 右端显式
  - 对称分裂
  - 与特征线同步高阶积分

## 5. 数值处理对比

| 求解器 | 主变量 | 能量方向方法 | `chi` 方向方法 | 子步策略 | 当前角色 |
| --- | --- | --- | --- | --- | --- |
| `fullhide_1d` | `dN_x` | 一阶隐式迎风 | 无 | 固定/自适应 | 1D 稳定基线 |
| `slc1_1d` | `dN_x` | 半拉格朗日 + Strang split | 无 | 固定大量子步 | 方法比较 |
| `charint_1d` | `dN_x` | 特征线保守 remap | 无 | 固定/自适应 | 1D 高保真方法 |
| `t2g1_1d` | `dN_x` | 三层二阶隐式 | 无 | 固定大量子步 | 时间推进比较 |
| `weno5_1d` | `dN_x` | WENO5 + SSP RK3 | 无 | CFL 显式 | 高阶显式比较 |
| `fullhide_2d` | `U(x,eta)` | 一阶隐式迎风 | 隐式平流-扩散-源项 | 双方向 CFL | 2D 物理基线 |
| `charint_2d` | `U(x,eta)` | 特征线 remap | 当前仍为隐式 | `xi` CFL + 上限 | 2D 加速混合版 |

## 6. 每个方法的适用场景

### 优先追求稳健参考解

推荐：

- `fullhide_1d`
- `fullhide_2d`

原因：

- 稳
- 结构清楚
- 最适合做其他算法的对照基准

### 优先追求 1D 高频端与截止保真

推荐：

- `charint_1d`

原因：

- 特征线守恒 remap 对尖锐前沿更友好

### 优先追求 2D 计算速度

当前可考虑：

- `charint_2d`

但必须记住：

- 它当前的速度来自混合设计和子步上限
- 不是完全无代价的“白赚”

### 作为算法比较样本

可以用：

- `slc1_1d`
- `t2g1_1d`
- `weno5_1d`

它们更适合拿来判断：

- 误差来自空间离散还是时间推进
- 高阶显式/半拉格朗日/三层法各自会产生什么谱形差异

## 7. 当前工作树下的现实判断

从当前代码组织和测试入口来看，可以把求解器分成两层：

### 第一层：当前主干工作流

- `fullhide_1d`
- `charint_1d`
- `fullhide_2d`
- `charint_2d`

这四个方法和当前 runtime、测试脚本、2D 历史光子场框架的耦合最紧。

### 第二层：方法库与比较基线

- `slc1_1d`
- `t2g1_1d`
- `weno5_1d`

它们仍然有价值，但更像“保留下来的不同数值路线”，而不是当前 2D 开发的主轴。

## 8. 当前已知问题与风险

### 8.1 `charint_2d` 还没有彻底完成 `eta` 特征线化

当前公开版本仍是：

- `eta` 隐式
- `xi` 特征线

所以如果有人把它理解成“完整 2D characteristic solver”，那是不准确的。

### 8.2 `charint_2d` 的高能端仍需继续校准

当前工作树下，`charint_2d` 相比 `fullhide_2d` 的早期高频端已经比之前明显改善，但仍可能偏低。根因不是单一 bug，而是：

- `xi` 子步大小
- `L1` 上限
- 高能端特征线步长与谱形保真之间的直接权衡

### 8.3 2D 冷却历史不是全频完整积分

当前 2D 内部使用：

- `Num_nu_cool = min(6, Num_nu)`

这意味着 2D 历史冷却求积是一个 reduced-band 近似，而不是对全 `Num_nu` 频带逐点保留。

这不是 bug，但它是 2D 模型和“理想全频历史求积”之间的明确近似。

## 9. 一句话总结

如果把当前代码库只看成一句话：

- `fullhide_*` 是稳健隐式基线
- `charint_*` 是沿特征线做保守 remap 的高保真路线
- `charint_2d` 当前是一个已经实用、但仍带有明显速度/精度折中的混合版本
- `slc1/t2g1/weno5` 主要用于方法学比较，而不是当前 2D 主线的最终形态

## 10. 连续方程与离散推导

这一节只做两件事：

- 先把代码背后的连续方程写清楚
- 再说明每个 solver 实际离散的是哪一部分

### 10.1 1D 电子输运方程

1D 里真正要解的是电子分布随半径的守恒方程。若用

- `N(gamma, R) = dN / dgamma`
- `Q(gamma, R) = dN_inj / (dR dgamma)`

则连续方程可写成

```text
∂N/∂R + ∂/∂gamma [ gamma' * N ] = Q
```

其中

```text
gamma' = dgamma/dR
```

包含：

- 绝热冷却
- 同步辐射冷却
- IC/KN 冷却

代码中普遍不直接在 `gamma` 上做，而是改到

- `x = log10(gamma)`
- `dN_x = dN / dx = gamma * ln(10) * N`

因为

```text
dgamma = gamma * ln(10) * dx
```

所以守恒方程变成

```text
∂dN_x/∂R + ∂/∂x [ A_x(x,R) * dN_x ] = S_x(x,R)
```

这里

```text
A_x = gamma' / (gamma * ln(10))
S_x = gamma * ln(10) * Q
```

这就是 `fullhide_1d`、`slc1_1d`、`t2g1_1d`、`weno5_1d`、`charint_1d` 共享的基本守恒形式。

### 10.2 为什么代码里会出现 `1 / (R ln 10)`

在当前实现里，损失项不是只写成纯微观冷却 `dEL/dR`，还把球对称膨胀带来的守恒 Jacobian 一并并入面对流系数。所以代码会构造

```text
face_coeff ≈ <dEL/dR> + 1 / (R ln 10)
```

这里第二项不是额外拍脑袋加的“保护项”，而是变量变换后从

- 径向膨胀
- 对数网格体积元
- `dN_x = gamma ln(10) N`

一起整理出来的几何项。

因此，对 1D solver 来说，真正被离散的是“冷却导致的能量漂移 + 几何膨胀导致的有效漂移”的合成通量。

### 10.3 `fullhide_1d` 的离散式

`fullhide_1d` 在 `dN_x` 上做一阶隐式有限体积更新。把网格记作 `x_{i-1/2}, x_{i+1/2}`，单元平均记作 `U_i^n`，则它的离散思想是

```text
(U_i^{n+1} - U_i^n) / ΔR
+ [ F_{i+1/2}^{n+1} - F_{i-1/2}^{n+1} ] / Δx_i
= S_i^{n+1}
```

其中面对流通量采用单向迎风

```text
F_{i+1/2}^{n+1} = A_{i+1/2}^{n+1} * U_upwind^{n+1}
```

结果是一个三对角线性系统，用 backward sweep 解出。

这个格式的意义很直接：

- 无条件稳定性强
- 对很硬的冷却问题可靠
- 但数值耗散也最大

### 10.4 `slc1_1d` 的离散式

`slc1_1d` 把 1D 问题看成

```text
输运 + 注入
```

然后做 split：

1. 沿能量流线回溯到旧时刻
2. 用保守重映射得到平流后的分布
3. 再加注入源

它本质上是半拉格朗日法。和 `charint_1d` 的差别不是“是否回溯特征线”，而是：

- `slc1_1d` 的时间推进更接近 operator splitting
- `charint_1d` 把源项和保守 remap 结合得更紧

所以 `slc1_1d` 直觉上更像“搬运再加源项”，而 `charint_1d` 更像“沿特征线把整个守恒更新一次做完”。

### 10.5 `charint_1d` 的离散式

`charint_1d` 先改变量

```text
u = 1 / gamma
```

原因是很多冷却律下，`u(R)` 比 `gamma(R)` 更接近线性或分段线性。于是特征方程变成

```text
du/dR = G(u, R)
```

当前代码分两种情况：

- `index_Y == 0` 时，按 affine `u` 特征线处理
- `index_Y != 0` 时，按 piecewise-affine `u` 特征线处理

离散上不是“把单元中心沿特征线搬一下”这么简单，而是：

1. 回溯每个单元边界在上一步对应的位置
2. 用非均匀网格保守 remap 计算旧解有多少质量落入新区间
3. 对注入源用四点求积沿特征线积分

写成形式化表达就是

```text
U_i^{n+1}
= Remap[ X_{i-1/2}^*, X_{i+1/2}^* ; U^n ]
+ Integral_of_source_along_characteristics
```

这里 `X^*` 是回溯后的边界位置。

这个方法的核心优点不是“高阶”这三个字，而是：

- 它直接利用了 PDE 的特征结构
- 高能端 cutoff 和尖前沿通常比一阶隐式更容易保住

### 10.6 `t2g1_1d` 的离散式

`t2g1_1d` 是三层时间推进。它用

- `n-1`
- `n`
- `n+1`

三层构造二阶时间离散，启动时再退回单步法。

抽象写法是

```text
(a U^{n+1} + b U^n + c U^{n-1}) / ΔR
+ L_x(U^{n+1}, U^n, U^{n-1})
= S
```

这里 `L_x` 仍是能量方向输运算子。它想得到的是：

- 比一阶隐式更高的时间精度
- 但仍保持隐式框架对刚性冷却的稳定性

代价是：

- 启动复杂
- 系数耦合更重
- 调参与维护都比 `fullhide_1d` 难

### 10.7 `weno5_1d` 的离散式

`weno5_1d` 把 1D 方程当作显式守恒律来做：

```text
∂U/∂R + ∂F(U)/∂x = S
```

其中空间通量用 WENO5 重构，时间推进用 SSP RK3。它的强项是：

- 光滑区高阶
- 平移型问题色散/耗散特性通常更好

但这里有一个非常关键的物理-数值矛盾：

- 电子冷却问题往往是刚性的
- 显式 RK 的稳定步长必须非常小

所以它更适合做数值方法比较，而不是当前仓库的主力生产 solver。

### 10.8 2D 连续方程的自然写法

2D 里公共变量是

- `x = log10(gamma)`
- `eta = log10(chi)`
- `U(x, eta, R) = dN / (dx deta)`

那么连续方程自然可写成

```text
∂U/∂R
+ ∂/∂x   [ A_x   U ]
+ ∂/∂eta [ A_eta U ]
= ∂/∂eta [ D_eta ∂U/∂eta ]
+ S(x, eta, R)
```

这里

- `A_x` 是能量方向漂移，来自冷却与几何项
- `A_eta` 是壳后空间坐标 `chi` 的平流
- `D_eta` 是 `chi` 方向扩散或展宽项
- `S` 是注入项，通常集中在靠近 shock 的 `chi ~ 1`

关键结论只有一个：

- `x` 方向是典型一阶平流型
- `eta` 方向是“平流 + 扩散 + 局域源项”混合型

因此数学上不能把整个 2D 问题都说成“纯特征线问题”。

### 10.9 `fullhide_2d` 的分裂与离散

`fullhide_2d` 对上面的 2D 方程采取的是稳定优先路线。按当前实现理解，它每个子步主要做两件事：

1. 在 `eta` 方向求解

```text
∂U/∂R + ∂(A_eta U)/∂eta = ∂/∂eta ( D_eta ∂U/∂eta ) + S_eta
```

对应代码中隐式三对角更新。

2. 在 `x` 方向求解

```text
∂U/∂R + ∂(A_x U)/∂x = S_x
```

对应每个 `chi` 列上与 `fullhide_1d` 同构的一阶隐式迎风。

所以 `fullhide_2d` 的真实数值结构是：

- `eta` 方向：隐式平流-扩散-源项
- `x` 方向：隐式迎风

这也是它稳、但高能端更容易被抹平的根本原因。

### 10.10 `charint_2d` 当前版本的真实离散

当前工作树下，`charint_2d` 公开实现不是

```text
eta-charint + xi-charint
```

而是

```text
eta-implicit + xi-charint
```

也就是说它真正做的是：

1. `eta` 方向仍沿用 `fullhide_2d` 的隐式更新

```text
∂U/∂R + ∂(A_eta U)/∂eta = ∂/∂eta ( D_eta ∂U/∂eta ) + S_eta
```

2. `x` 方向对每个 `chi` 列改用 `charint_1d` 同构的特征线保守 remap

```text
∂U/∂R + ∂(A_x U)/∂x = S_x
```

的 characteristic 版本

因此当前 `charint_2d` 的准确说法只能是：

- 它在 2D 外层物理框架里，把最耗散的 `xi/log-gamma` 方向替换成特征线推进
- 它没有把 `eta` 方向的扩散问题伪装成特征线问题

这也是它在数学上仍然自洽的原因。

### 10.11 为什么 `charint_2d` 不能简单做成“全方向特征线”

如果只看

```text
∂U/∂R + ∂(A_eta U)/∂eta = 0
```

那当然可以沿 `eta` 特征线回溯。

但当前 `eta` 方程不是这个，而是

```text
∂U/∂R + ∂(A_eta U)/∂eta = ∂/∂eta ( D_eta ∂U/∂eta ) + S_eta
```

右边有扩散算子。扩散项的本质是二阶抹平，不对应一族单值特征线。于是：

- 平流项可以特征线化
- 扩散项不能被严格写成单纯 remap

所以如果把整个 `eta` 子问题都改叫“charint”，那就不是严格的 PDE 求解陈述，而是近似算法包装。当前代码没有这样做，这是对的。

### 10.12 当前代码中的“理想 PDE”与“真实实现”对应关系

如果只从数学结构出发，最自然的 2D Strang-split 版本其实应当是

```text
eta-advection  ->  eta-diffusion  ->  xi-charint
```

其中：

- `eta-advection` 用保守特征线 remap
- `eta-diffusion` 用隐式三对角
- `xi-charint` 用 `charint_1d` 同构推进

这就是“概念上最干净”的 `charint_2d`。

但当前工作树下真正公开的实现还停在：

```text
eta-implicit  ->  xi-charint
```

所以阅读和比较结果时必须分清：

- “这个 solver 想做成什么”
- “这个 solver 现在实际上做了什么”

这两个问题在当前仓库里不是同一个答案。

## 11. PDE 项与代码例程映射

这一节不再讲抽象算法，而是直接回答：

- 连续方程里的每一项在代码里由谁负责
- 不同 solver 实际调用了哪条数值路径

### 11.1 网格、单元边界和几何量

1D `gamma` 网格与对数单元边界主要在：

- `src/Electron/electron_common.f90`
  - `electron_build_gamma_grid`
  - `electron_log_cell_edges`
  - `electron_cell_geometry`

它们负责：

- 构造 `gam_e`
- 构造 `x = log10(gamma)` 的单元边界
- 给出单元中心与单元宽度

2D `chi/eta` 几何主要在：

- `src/Electron/electron_transport_2d_kernel.f90`
  - `compute_log_chi_geometry`
  - `compute_downstream_comoving_grid`

它们负责：

- 构造 `eta_face / eta_grid`
- 构造 `chi_face / chi_grid`
- 构造下游共动路径长度与单元宽度

这部分对应 PDE 里的：

- `x` 网格
- `eta` 网格
- 2D 下游传播路径上的积分几何权重

### 11.2 注入项 `S`

初始电子谱构造主要在：

- `src/Electron/electron_common.f90`
  - `electron_initial_powerlaw`
  - `electron_initial_powerlaw_exp_cutoff`
  - `electron_initial_powerlaw_exp_cutoff_edges`

步进时的注入源构造主要在：

- `src/Electron/electron_common.f90`
  - `electron_injection_prefactor`
  - `electron_build_source_term`
  - `electron_build_source_term_exp_cutoff`
  - `electron_build_source_term_exp_cutoff_edges`
  - `electron_build_source_term_profile`

当前各 solver 主要调用路径是：

- `fullhide_1d`：`electron_injection_prefactor` + `electron_build_source_term_exp_cutoff`
- `charint_1d`：`electron_injection_prefactor` + `electron_build_source_term_exp_cutoff_edges`
- `t2g1_1d`：`electron_injection_prefactor` + `electron_build_source_term_profile`
- `fullhide_2d / charint_2d`：每个 `chi` 列上仍以 `electron_build_source_term_exp_cutoff` 为主

这里很关键的一点是：

- `charint_1d` 用 edge-based 源项版本，不是装饰，而是因为它要和特征线边界回溯一起做守恒求积

### 11.3 冷却项与损失系数 `A_x`

连续方程中能量方向漂移的核心输入，是每个电子 Lorentz 因子的总损失率。当前实现里它主要由两层代码负责。

第一层是冷却物理项组装：

- `src/Electron/electron_cooling_kernel.f90`
  - `prepare_forward_cooling_aux`
  - `prepare_forward_cooling_aux_batch`
  - `assemble_forward_cooling`
  - `assemble_forward_cooling_split`
  - `assemble_forward_cooling_split_batch`
  - `assemble_forward_cooling_with_ssa`
  - `assemble_forward_cooling_from_terms`

这层把：

- synchrotron
- SSA
- IC/KN
- Compton `Y`

装配成 `dEl`。

第二层是把 `dEl` 变成数值输运需要的平均漂移系数：

- `src/Electron/electron_common.f90`
  - `electron_loss_mean`
  - `electron_gamma_c_from_loss_mean`
  - `electron_prepare_implicit_coeffs`

调用路径上：

- 所有 1D solver 都会先构造 `dEl`
- 再把它压成 `dEL_mean`
- 再喂给各自的输运器

2D 里则是：

- 先按 `chi` 计算 `dEl_chi(:, I_chi)`
- 再得到 `dEL_mean_chi(:, I_chi)`

所以 PDE 中的 `A_x(x,eta,R)` 在代码里不是单独一个函数，而是：

- 冷却核给出局域损失
- 公共模块把损失转换成面系数
- 具体 solver 再按自己的离散方式消费它

### 11.4 1D `x` 方向隐式迎风

这条路径的核心更新器在：

- `src/Electron/electron_common.f90`
  - `electron_fullhide_step`
  - `electron_prepare_implicit_coeffs`
  - `electron_backward_sweep`

对应主入口：

- `src/Electron/FS_electron_fullhide_1d.f90::fs_electron_fullhide_1d`

数学上对应：

```text
∂U/∂R + ∂(A_x U)/∂x = S
```

的一阶隐式有限体积更新。

2D `xi` 方向的隐式版本也复用同一思路，只是外层变成“对每个 `chi` 列逐列求解”。

### 11.5 1D 半拉格朗日路径

这条路径的核心在：

- `src/Electron/electron_common.f90`
  - `electron_semi_lagrangian_step`
  - `electron_semi_lagrangian_step_nonuniform`

对应主入口：

- `src/Electron/FS_electron_slc1_1d.f90::fs_electron_slc1_1d`

这里它承担的是：

- 回溯旧位置
- 做保守重映射
- 与源项 split 后组合

### 11.6 1D 特征线保守 remap 路径

这条路径的核心在：

- `src/Electron/electron_common.f90`
  - `electron_characteristic_step_prepared_core`
  - `electron_conservative_remap_nonuniform_prepared`
  - `electron_trace_affine_u_edges_batch`
  - `electron_trace_piecewise_affine_u_edges_batch`
  - `electron_characteristic_transport_affine_u`
  - `electron_characteristic_transport_piecewise_u`

对应主入口：

- `src/Electron/FS_electron_charint_1d.f90::fs_electron_charint_1d`

这里的调用关系是：

- 有源项时，主入口先准备 edge source
- 然后走 `electron_characteristic_step_*_prepared_source`
- 最终在 `electron_characteristic_step_prepared_core` 中完成边界回溯、保守 remap 和源项积分

所以 1D `charint` 的数值核心并不是“一个黑箱 solver”，而是一套：

- 边界追踪
- 非均匀守恒 remap
- 沿特征线源项求积

的组合。

### 11.7 2D `eta` 方向平流与扩散

2D `eta/log-chi` 方向相关例程集中在：

- `src/Electron/electron_transport_2d_kernel.f90`
  - `eta_face_transport_coeffs`
  - `advance_eta_logchi_implicit`
  - `advance_eta_logchi_advection_charint`
  - `advance_eta_logchi_diffusion_implicit`

其中：

- `eta_face_transport_coeffs` 负责构造 `A_eta_face`
- `advance_eta_logchi_implicit` 是当前生产路径上的主更新器
- `advance_eta_logchi_advection_charint` 是未来可用的纯平流特征线部分
- `advance_eta_logchi_diffusion_implicit` 是 `eta` 扩散的独立隐式步

当前公开实现里：

- `fullhide_2d` 用 `advance_eta_logchi_implicit`
- `charint_2d` 也仍用 `advance_eta_logchi_implicit`

原因不是“还没来得及切”，而是因为当前 `eta` 子问题本来就不是单纯特征线型。

### 11.8 2D `xi/log-gamma` 方向两条路径

2D `xi` 方向隐式路径在：

- `src/Electron/electron_transport_2d_kernel.f90`
  - `advance_energy_loggamma_chi`

它本质上就是：

- 对每个 `chi` 列计算 `dEL_mean_chi`
- 逐列做与 `fullhide_1d` 同构的隐式迎风

2D `xi` 方向 characteristic 路径在：

- `src/Electron/electron_transport_2d_kernel.f90`
  - `advance_energy_loggamma_chi_charint`

它内部会：

- 调 `electron_log_cell_edges`
- 再按 `index_Y` 分支到
  - `electron_characteristic_transport_affine_u`
  - `electron_characteristic_transport_piecewise_u`

当前公开实现里：

- `fullhide_2d` 走 `advance_energy_loggamma_chi`
- `charint_2d` 走 `advance_energy_loggamma_chi_charint`

所以 `charint_2d` 真正替换掉的，是 2D 中最直接影响高能端谱形的那一条 `xi` 漂移路径。

### 11.9 辐射输出、历史光子场和诊断量

辐射谱核心在：

- `src/Electron/electron_radiation_kernel.f90`
  - `get_syn`
  - `get_syn_state`
  - `get_syn_selected`
  - `reduce_syn_shell_from_chi`
  - `get_nu_a`
  - `get_nu_a_nonuniform`
  - `get_nu_a_2d_path`
  - `get_nu_a_2d_cell_path`

历史光子场在：

- `src/Electron/electron_seed_history_kernel.f90`
  - `propagate_downstream_history_fields`
  - `accumulate_comoving_history_fields`
  - `history_transfer_weight`

2D cooling 预处理与 reduced-band 冷却在：

- `src/Electron/electron_cooling_kernel.f90`
  - `prepare_forward_cooling_aux_batch`
  - `assemble_forward_cooling_split_batch`

当前 `nu_m / nu_c / nu_a` 的逻辑可以概括成：

- `nu_m`
  - 1D 中主要由 `Gam_e_m` 直接换算
  - 2D 中先得到每个 `chi` 单元的局域峰，再做壳层汇总
- `nu_c`
  - 由 `dEL_mean` 反解 `Gam_e_c_diag`
  - 再换算成对应冷却断点频率
- `nu_a`
  - 1D 用 `get_nu_a` 或 `get_nu_a_nonuniform`
  - 2D 用 `get_nu_a_2d_cell_path`

也就是说，当前公共输出量并不是“输运器自己顺手带出来的副产物”，而是：

- 先做电子输运
- 再做同步辐射/光子场投影
- 最后从投影结果或平均冷却率中抽取诊断量

### 11.10 主入口到内核的最短调用图

如果只抓主干路径，可以压缩成下面这几条。

`fullhide_1d`：

```text
fs_electron_fullhide_1d
  -> get_forward_cooling
  -> electron_loss_mean
  -> electron_fullhide_step
  -> get_syn_selected / get_nu_a
```

`charint_1d`：

```text
fs_electron_charint_1d
  -> get_forward_cooling
  -> electron_build_source_term_exp_cutoff_edges
  -> electron_characteristic_step_*_prepared_source
  -> get_syn_selected / get_nu_a
```

`fullhide_2d`：

```text
fs_electron_fullhide_2d
  -> reduce_syn_shell_from_chi + accumulate_comoving_history_fields
  -> prepare_forward_cooling_aux_batch + assemble_forward_cooling_split_batch
  -> advance_eta_logchi_implicit
  -> advance_energy_loggamma_chi
  -> get_syn_state + get_nu_a_2d_cell_path
```

`charint_2d`：

```text
fs_electron_charint_2d
  -> reduce_syn_shell_from_chi + accumulate_comoving_history_fields
  -> prepare_forward_cooling_aux_batch + assemble_forward_cooling_split_batch
  -> advance_eta_logchi_implicit
  -> advance_energy_loggamma_chi_charint
  -> get_syn_state + get_nu_a_2d_cell_path
```

## 12. 误差来源与物理风险

这一节只讨论两类东西：

- 方法本身不可避免的误差结构
- 当前实现中特别容易把误差放大成物理假象的地方

如果一个量是“真实物理趋势”，它应该在参数或时间上连续、平滑；如果不是，就优先怀疑这里列出的数值源。

### 12.1 `fullhide_1d`

主要误差来源：

- 一阶隐式迎风的数值耗散
- 在强冷却区对高能端前沿的额外抹平
- 子步过大时断点位置向低能偏移
- `active_hi` 支持度截断会让极弱高能尾更早退出主更新区

最容易受影响的物理量：

- 高频端 cutoff
- `nu_c`
- 电子谱尖锐折点

典型症状：

- 高能端比理论预期更圆
- 早期快冷阶段尾部更“短”

判断上要分清：

- 这是方法本性，不一定是 bug
- 但如果 `nu_c` 随时间出现跳变或锯齿，那就不是单纯的一阶耗散能解释的

### 12.2 `slc1_1d`

主要误差来源：

- split error
- 半拉格朗日 remap 的插值误差
- 子步数固定时，某些阶段过粗、某些阶段过密

最容易受影响的物理量：

- 注入附近谱形
- 断点附近曲率

典型症状：

- 整体平移看起来合理
- 但注入与冷却耦合强时，局域斜率会有系统偏差

它通常不会像显式 WENO 那样炸掉，但也不会像真正沿特征线积分的 `charint_1d` 那样稳地保住截止。

### 12.3 `charint_1d`

主要误差来源：

- 特征线追踪误差
- 源项沿特征线积分的时间离散误差
- 子步控制过粗时的高能端质量流失
- 固定子步分支中的 `L1 <= 64` 硬上限

最容易受影响的物理量：

- cutoff 位置
- 高频同步辐射尾部
- `nu_c` 早期演化

典型症状：

- 子步太粗时，高能端会被“提前吃掉”
- 如果 `index_Y != 0` 分支的 piecewise-affine 特征线没有充分分辨，IC/KN 情形下偏差会更明显
- `dEl` 若在单元内弯曲明显，分段仿射回溯会直接把 cutoff 位置带偏

它的关键风险不是低阶扩散，而是：

- 表面上很“sharp”
- 但如果特征线时间步不够细，错误会直接落在截止位置上

### 12.4 `t2g1_1d`

主要误差来源：

- 三层法启动步误差
- 多时间层耦合对初始两步特别敏感
- 系数不一致时容易出现非物理振荡或相位偏差
- 启动段与后续三层段不是同一离散阶

最容易受影响的物理量：

- 早期时间的 `nu_m / nu_c`
- 刚从初值演化出的谱形曲率

典型症状：

- 中后期还行
- 最前几步容易与其他方法分开

所以比较 `t2g1_1d` 时，不能只看长时间段，也要单独盯前几步。

### 12.5 `weno5_1d`

主要误差来源：

- 显式 CFL 限制
- 刚性冷却下时间步太小，否则空间高阶没有意义
- positivity clipping 会主动改动解
- 常值 ghost cell 外推会把边界误差持续送回内部

最容易受影响的物理量：

- 高频尾部
- 注入低端附近的小尺度结构
- 在极强冷却下的总电子数守恒

典型症状：

- 如果步长足够小，光滑区很漂亮
- 如果步长不够小，会先出现时间推进误差，而不是空间重构误差
- 一旦触发负值截断，谱尾可能被人为削平
- 当前实现的 `where(dN_x < 0) dN_x = 0` 不是中性操作，它确实会改动守恒与谱尾形状

所以 `weno5_1d` 的风险不是“高阶会振荡”这么简单，而是：

- 对刚性问题，时间推进往往先成为瓶颈

### 12.6 `fullhide_2d`

主要误差来源：

- `eta` 方向把平流、扩散、源项混在一个隐式步里
- `xi` 方向仍是一阶隐式
- 2D 冷却历史使用 reduced cooling bands
- 历史回灌与投影走 reduced log-band 映射，而不是全频逐点保留

最容易受影响的物理量：

- 高频端截断
- `chi` 方向展宽后的局域谱形
- `nu_a` 路径积分结果

典型症状：

- 壳层总谱通常平滑、稳
- 但高能端偏保守
- `chi` 方向结构会比真实解更钝

其中 reduced cooling bands 的具体风险是：

- 当前实现用 `Num_nu_cool = min(6, Num_nu)`
- 这会把 2D 历史冷却耦合压缩到少量频带
- 如果某个物理设置对窄频段结构很敏感，`dEl_chi` 就可能出现系统偏差

这不是实现 bug，但它确实是明确近似。

### 12.7 `charint_2d`

主要误差来源：

- `eta` 仍走隐式，`xi` 走 characteristic，误差结构不对称
- `xi` 子步控制直接决定高能端保真
- `active_chi_hi` 和 `L1` 上限会影响真实参与推进的列数与时间分辨率
- 同样继承了 2D reduced cooling bands 近似
- `active_hi` 也会影响每个 `chi` 列实际推进到的高能支持范围

最容易受影响的物理量：

- 高频同步/SSC 端
- 最早时刻的 cutoff
- `chi` 汇总后的壳层总谱

典型症状：

- 总体上比 `fullhide_2d` 更能保高能端
- 但如果 `xi` 子步限制过强，反而会出现更明显的早期高能截断
- 如果 `eta` 方向本身仍偏耗散，最终壳层总谱仍可能和 `fullhide_2d` 拉不开足够差距

当前最需要警惕的实现细节有五个：

- `use_charint_eta = .false.`
  - 所以它现在不是“全 2D charint”
- `dDR_try` 的选择
  - 直接决定 `xi` 回溯是否过粗
- `L1` 上限
  - 限得太狠会把高能端提前削掉
- `active_hi / active_chi_hi` 阈值
  - 会让弱占据高能列更早退出更新
- `Num_nu_cool = min(6, Num_nu)`
  - 会把冷却历史的谱结构进一步简化

换句话说，`charint_2d` 当前最可能出问题的地方，不是某一个神秘的公式常数，而是：

- 子步策略
- 混合分裂结构
- reduced-band 冷却近似

三者叠加后的系统偏差。

### 12.8 跨 solver 的共性风险

不管用哪条路线，下面这些现象一旦出现，都更像 bug 或实现不一致，而不是方法差异。

- `nu_m / nu_c / nu_a` 随时间明显跳变
- 壳层总谱在相邻时间点出现非物理锯齿
- 同一物理设置下，总电子数出现无理由突变
- 高频端在参数平滑变化时出现台阶式缩短
- 2D 各 `chi` 列积分到总谱后反而比单列更不连续

这些症状通常说明：

- 子步控制不连续
- 某个 diagnostic 使用了与主推进不一致的系数
- 历史光子场映射或 `tau` 路径积分有离散跳变
- 某条 solver 路径在边界单元处理上和其他路径不一致

### 12.9 如何区分“方法误差”和“物理错误”

一个实用判断标准是：

如果差异主要表现为

- 更平
- 更圆
- 更保守

那常常是数值耗散。

如果差异主要表现为

- 断点位置错位很多
- 高频端突然断掉
- 时间曲线不连续

那就不能只说“这是方法差异”，必须继续查：

- 子步
- 边界回溯
- 源项积分
- 历史场投影
- 诊断量提取

### 12.10 当前仓库下最值得盯住的三件事

如果只保留三个最重要的核查点，就是：

1. `charint_2d` 的 `xi` 子步是否足够细  
   这是当前高能端是否过早截断的第一嫌疑项。

2. 2D reduced cooling bands 是否在目标物理设置下足够  
   如果不够，`dEl_chi` 和随后 `nu_c`、高频谱尾都会偏。

3. `nu_m / nu_c / nu_a` 的提取是否与主输运保持同一套局域状态  
   如果诊断从不同近似里拼出来，就会表现为时间不平滑。
