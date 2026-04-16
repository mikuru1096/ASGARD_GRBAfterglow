# 解析特征积分电子求解器方案

## Summary

- 你描述的方案我理解了。它在数学上叫 **特征线法**（method of characteristics, MOC）；如果写成单元平均、保守更新的形式，在数值输运文献里属于 **conservative / cell-integrated semi-Lagrangian** 一类；在高能电子动力学文献里也常被表述为 **沿冷却轨迹的半解析解** 或 **Green’s-function / Duhamel 解法**。
- 对当前项目，这条路线应实现为一个**新 solver**，不替换现有 `slc1`。公开名称固定为 `electron_solver="charint"`。
- 物理范围按已确认口径执行：
  - `index_y=0`：走**严格解析**特征映射；
  - `index_y=1/2/3`：走**子步冻结 cooling 系数后的近似解析**特征映射；
  - 这版不再继续基于当前临时 `SL+IMEX` 原型迭代，先回到最近接受的 `slc1` 基线，再新增 `charint`。

## Key Changes

### 1. 新 solver 的数学形式固定

- 保持现有守恒量 `dN_x` 为内部主变量，不改辐射层接口。
- 在每个半径子步，把电子方程写成
  - 累积电子项：旧时层 `dN_x` 沿特征线搬运到新时层；
  - 注入电子项：本子步内源项 `dF1` 通过 Duhamel 积分沿同一特征线映射到新时层。
- 特征参数化固定在 `u = 1 / gamma` 空间完成，因为：
  - 线性冷却 + 绝热项时，`u` 满足仿射 ODE；
  - 对 `index_y=0` 可得到全子步闭式映射；
  - 对 `index_y>0` 可把冻结后的 cooling 写成分段仿射 `du/dR = a_i + b_i u`，逐单元精确过界。

### 2. `index_y=0` 的严格解析模式

- 在子步中点冻结流体量 `R_mid, Gamma_mid, B_mid, dNe_mid`。
- 用物理公式直接构造全局系数：
  - radiative coefficient `a_rad`
  - adiabatic coefficient `b_ad = 1 / R_mid`
- 用闭式特征映射追踪每个目标单元的左右边界，从而得到旧时层 departure interval。
- 累积电子项：
  - 对旧 `dN_x` 做单元内一次重构；
  - 对 departure interval 上的积分解析计算，得到 cell-average 更新。
- 注入电子项：
  - 源项以当前已存在的 cell-average `dF1` 为基础；
  - 视每个源单元在本子步内为常值源；
  - 对“源单元 × 注入时间”落到目标单元的时空重叠体积做解析积分。
- 这一版的“高精度”定义是：**对冻结子步方程精确**，误差只来自子步内物理量冻结，而不来自输运离散。

### 3. `index_y=1/2/3` 的近似解析模式

- 不要求全局闭式冷却律。
- 在子步中点先调用现有 cooling 例程，得到网格上的 `dEl(gamma)`。
- 将其转换为 `u` 空间的 characteristic speed，并在每个 `u` 单元上构造分段仿射系数 `a_i, b_i`。
- 目标单元边界 backward tracing 时：
  - 精确计算每个 crossing 到相邻单元边界的时间；
  - 跨单元时切换到下一段仿射系数继续积分；
  - 直到回到子步起点。
- 累积电子和注入电子都沿这条分段解析特征处理。
- 这一路径的定位是“近似全覆盖”，不是严格解析；它的可接受标准是平滑、保正、守恒和不出现前沿锯齿。

### 4. 代码接入方式固定

- 新增 Fortran 求解器文件：`src/Electron/FS_electron_charint.f90`
- 公共特征映射、departure interval、注入时空重叠积分放入：`src/Electron/electron_common.f90`
- 运行时绑定只新增一条分支，不改现有 solver 语义：`asgard_runtime.py`
- 公开接口唯一变化：
  - `FitConfig.electron_solver` 新增 `"charint"`
- `slc1`、`fullhide`、`weno5` 保持原语义，`charint` 第一版不设为默认。

## Test Plan

### 1. 解析正确性单测

- 新增 exact-characteristics 单测，使用人工构造的冻结系数：
  - 纯搬运、无源项；
  - 纯源项、无旧分布；
  - 搬运 + 源项同时存在。
- 这些测试直接比较 `charint` 数值结果和同一冻结方程的解析单元平均解。
- 通过条件：
  - 相对误差达到机器精度量级；
  - 无负值；
  - 总粒子数变化与解析源积分一致。

### 2. 现有主线验收

- 在 `tests/order_convergence_check.py` 中把 `charint` 纳入与 `slc1/fullhide` 同级的电子链和辐射链检查。
- `index_y=0` 的 acceptance：
  - `charint-electron-peak-gamma > 2`
  - `charint-electron-support-low > 2`
  - `charint-electron-shape-aligned > 2`
- 同时要求辐射侧 `observer-total / bands / nu_a` 不低于当前 `slc1` 的已验收水平。

### 3. 近似全覆盖模式验收

- 新增 `charint` 平滑性检查，覆盖 `index_y=1/2/3`。
- 基线配置上检查 shell-by-shell 的
  - `peak_gamma`
  - `g_lo`
  - `nu_a`
- acceptance 固定为：
  - 所有量有限且为正；
  - 相邻 shell 的 `log10` 跳变 `< 0.25 dex`；
  - 不出现锯齿型一正一负交替振荡的 sawtooth 前沿。
- 继续保留谱图输出：
  - `tests/electron_exp_tail_spectrum_compare.py`
  - 每轮 solver 改动后都要更新图并人工核查前沿和平滑性。

## Assumptions And Defaults

- 默认从**最近接受的 `slc1` 基线**开始，不在当前临时 `SL+IMEX` 原型上继续叠改。
- `charint` 第一版不引入新的公开配置项；是否走严格解析或近似解析，由现有 `index_y` 自动决定。
- 方案命名与文献口径：
  - 数学名：**method of characteristics**
  - 数值输运名：**conservative / cell-integrated semi-Lagrangian**
  - 高能电子动力学名：**semi-analytic cooling-trajectory / Green’s-function solution**
- 主要文献依据：
  - Kardashev 1962：电子连续性方程解析框架
  - Vestrand 1983：沿能量损失轨迹对注入核积分
  - Zacharias & Schlickeiser 2017, A&A：将电子方程变到 `x = 1/gamma` 后得到纯输运/Green’s-function 形式
  - Kaas 2008, Tellus：cell-integrated / locally mass-conserving semi-Lagrangian 口径
