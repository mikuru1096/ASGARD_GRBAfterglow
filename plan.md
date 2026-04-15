# ASGARD 当前主线计划

## 0. 使用规则

- 这份文件是当前主线任务的唯一规划落点。
- 后续每当：
  - 新确认一个主问题
  - 新增一个阶段目标
  - 放弃一条路线
  - 决定更新或拒绝更新基准
  都必须同步更新本文件。
- 后续推进前优先查阅本文件，不再依赖对话记忆。

## 1. 当前稳定基线

- 电子默认求解器：`fullhide`
- 当前主链门槛：
  - `python tests/regression_check.py`
  - `python tests/vegasafterglow_comprehensive_validation.py`
  - `python tests/bug_audit_check.py`
  - `python tests/physical_closure_check.py`
- 当前稳定回归结果：
  - `xrt = 0.010657`
  - `optr = 0.000241`
  - `9GHz = 0.000000`
- 当前结论：
  - 主链可用
  - 还**不能**更新为新的数值基准

## 2. 已确认的主问题

### 2.1 电子谱有效收敛阶偏低

当前结果见：
- `output/vegasafterglow_doc/order_convergence.json`

当前实测：
- `fullhide-electron-spectrum = 0.5442`
- `t2g1-electron-spectrum = 0.6154`
- `weno5-electron-spectrum = 0.4973`

这说明电子谱本身的有效收敛阶偏低，主误差在电子谱求解链内部，不是只在后处理层。

### 2.2 电子谱低阶的已确认原因

#### A. `fullhide` 公共主因

文件：
- `src/Electron/FS_electron_fullhide.f90`
- `src/Electron/electron_common.f90`

已确认问题：
- 单边隐式递推是主骨架
- `electron_backward_sweep(...)` 使用 `max(zero, ...)`
- 这会在尾部、截止点和局部负值附近引入非线性截断
- 对整条电子谱形状的有效阶数有直接压制
- 当前 `fullhide/t2g1` 的默认输运核**不是 Chang-Cooper**
- 它更接近单边隐式上风递推：
  - 没有 Chang-Cooper 界面权重 `delta(w)`
  - 没有以界面通量守恒形式组装三对角系统
  - 也没有按漂移/扩散比构造指数拟合权重

#### B. `fullhide` 和 `t2g1` 共享主因

文件：
- `src/Electron/FS_electron_fullhide.f90`
- `src/Electron/FS_electron_t2g1.f90`
- `src/Electron/calling_modules.f90`

已确认问题：
- `P_syn / Seed_syn -> dEl -> dEL_mean` 先按旧谱计算
- 然后在整个外层壳步里冻结使用
- 这是 lagged cooling / operator splitting
- 会系统性降低全链路有效阶数
- 但根据 `electron_cooling_lag.json` 的最新离线量化，在当前默认 `index_Y=2` 主路径下，它**不是峰位锁网格的主限**

#### C. `t2g1` 额外主因

文件：
- `src/Electron/FS_electron_t2g1.f90`

已确认问题：
- 虽然目标是二阶时间精度
- 但启动阶段仍回退到单步低阶格式
- 当前已修正 `L1` 的明显错误
- 但收敛阶只从约 `0.54` 提高到约 `0.62`
- 说明真正主因不是 `L1` 一个变量

#### D. 源项和谱边界的网格硬切

文件：
- `src/Electron/electron_common.f90`

已确认问题：
- `electron_build_source_term(...)` 仍是硬截断
- 初始谱 `electron_initial_powerlaw(...)` 也是格点级分段
- 当 `gamma_m / gamma_c / gamma_max` 移动时，会形成明显边界层误差

#### E. `weno5` 仍未达到应有精度

文件：
- `src/Electron/FS_electron_weno5.f90`

已确认问题：
- 已从“高分辨率会炸”修到“正收敛”
- 但仍然又慢又不够准
- 当前怀疑重点：
  - 运输项离散与方程并不严格等价
  - 边界处理仍过粗
  - 高能端数值耗散过强

#### F. `slc1` 半拉格朗日实验分支已接入，但还不能替代默认链

文件：
- `src/Electron/electron_common.f90`
- `src/Electron/FS_electron_slc1.f90`

当前状态：
- 已接入守恒型半拉格朗日实验分支 `slc1`
- 本轮已把一阶源项后加改成 `Strang` 式半步-输运-半步
- 重构已从更耗散的 `minmod` 提到二阶 TVD `MC` limiter

当前结果：
- `electron_solver_comparison.json`：
  - `slc1 direct rel_d_n_gam_e: 4.6073 -> 1.8037`
  - `slc1 observed rel_bands_flux: 0.9085 -> 0.2593`
  - `slc1 observed seconds ≈ 0.695`
- `order_convergence.json`：
  - `slc1-bands-flux ≈ 0.622`
  - `slc1-nu_a ≈ 0.655`
  - `slc1-electron-spectrum ≈ 0.330`
  - `slc1-electron-peak-gamma ≈ 1.066`
  - `slc1-electron-shape-aligned ≈ 1.475`

结论：
- `slc1` 的峰位项明显优于 `fullhide/t2g1`
- 但整条电子谱误差仍大，当前还不能作为默认链
- 下一步主攻不是继续加线程，而是把守恒 remap 再升一阶，优先考虑 `PPM/parabolic conservative remap`

## 3. 当前观察到的关键现象

### 3.1 电子谱支撑区和峰位随分辨率明显漂移

当前直接诊断表明：
- `g_lo`
- `g_hi`
- `peak_gamma`

都会随 `num_gam_e` 变化明显漂移。

这说明低阶不是“只差一个整体归一化”，而是谱形几何结构本身不稳定。

### 3.2 电子谱分解测阶的最新结果

最新结果见：
- `output/vegasafterglow_doc/order_convergence.json`

当前新增拆分指标后，结论是：

#### `fullhide`
- `electron-peak-gamma ≈ 0`
- `electron-support-low ≈ 0.41`
- `electron-support-high ≈ 0.81`
- `electron-shape-aligned ≈ 0.51`

#### `t2g1`
- `electron-peak-gamma ≈ 0`
- `electron-support-low ≈ 0.41`
- `electron-support-high ≈ 0.81`
- `electron-shape-aligned ≈ 0.47`

#### `weno5`
- `electron-peak-gamma ≈ 1.07`
- `electron-support-low ≈ 0.53`
- `electron-support-high ≈ 0.26`
- `electron-shape-aligned ≈ 0.77`

这组结果说明：
- `fullhide/t2g1` 的主问题首先是 **峰位锁网格**
- `fullhide/t2g1` 的峰位误差几乎不收敛
- 即使把峰位对齐，谱形本身也仍然只有约 `0.5` 阶
- `weno5` 的峰位比前两者好，但高能支撑区 `g_hi` 仍然很差
- 当前三条算法都还没有把电子谱几何结构真正收敛好

### 3.3 当前测阶方法是“有效收敛阶”，不是 formal order

文件：
- `tests/order_convergence_check.py`

当前做法：
- 三层网格
- Richardson / observed order
- 电子谱误差使用 `log(gamma)` 上的相对 `L1` 形状误差

这套方法对“全链路有效收敛阶”是成立的，
但不能直接当作离散格式 formal order。

### 3.4 固定动力学下的 `peak_gamma` 锁网格诊断

最新结果见：
- `output/vegasafterglow_doc/electron_peak_lock_diagnostic.json`

当前直接诊断结论：

#### `fullhide / t2g1`
- 在固定动力学输入下，`peak_gamma` 不会随着 `gamma_m` 最近格点误差单调下降
- 典型表现是相邻格点之间跳变，而不是连续贴近 `gamma_m`
- 代表性壳层：
  - `shell_frac = 0.5`
  - `gamma_m_theory = 60556.56`
  - `n = 61, 121, 241`
  - `nearest_to_gamma_m = 70026.28 -> 58513.84 -> 58513.84`
  - `fullhide peak_gamma = 70026.28 -> 70026.28 -> 64011.77`
  - `t2g1 peak_gamma = 70026.28 -> 58513.84 -> 64011.77`
- 这说明峰位不是简单地“取最接近 `gamma_m` 的格点”，而是被相邻几个格点的离散竞争锁住
- `g_lo` 同时也随分辨率明显上移：
  - `fullhide shell_frac = 0.3: 458.34 -> 2761.95 -> 7417.03`
  - `t2g1 shell_frac = 0.3: 458.34 -> 2761.95 -> 7417.03`
- 因此 `fullhide/t2g1` 的主问题不只是源项开关位置，而是：
  1. 初始谱和注入项的格点投影
  2. 冷却后谱形在相邻格点间的竞争
  3. 冻结冷却系数导致的谱峰位置滞后

#### `weno5`
- `peak_gamma` 相比 `fullhide/t2g1` 有更明显的收敛趋势
- 但 `g_lo` 和高能支撑区仍明显漂移：
  - `shell_frac = 0.5, g_lo = 8.81 -> 44.38 -> 99.58`
- 因此 `weno5` 的主误差仍集中在支撑区和高能端耗散，不在峰位本身

#### 直接推论
- `t2g1` 已经使用 `electron_build_source_term_profile(...)`，但峰位锁网格现象与 `fullhide` 几乎一样
- 所以“源项硬切”不是 `fullhide/t2g1` 共享主因的唯一来源
- 更值得优先修的是：
  1. 初始谱 `electron_initial_powerlaw(...)` 的格点级硬切
  2. `P_syn / Seed_syn -> dEl -> dEL_mean` 的冻结冷却链
  3. 必要时再看默认链源项的亚格点处理

### 3.5 初始谱亚格点投影的单独实验结果

最新结果见：
- `output/vegasafterglow_doc/electron_initial_projection.json`

当前直接结论：
- 单独把初始谱从“格点中心硬切”改成“log-gamma 单元内积分”，**不能单独解决** `peak_gamma` 锁网格
- 它能改善的是低能支撑起点 `g_lo`
- 但 `peak_gamma` 仍然会停在 `gamma_m` 上方一到两个格点

代表性结果：
- `n = 121`
  - `gamma_m_theory = 60604.89`
  - 最近格点 `= 58513.84`
  - `cell_center peak_gamma = 70026.28`
  - `cell_integrated peak_gamma = 70026.28`
  - 但 `g_lo` 从 `70026.28` 改善到 `58513.84`
- `n = 241`
  - 最近格点 `= 58513.84`
  - `cell_center peak_gamma = 64011.77`
  - `cell_integrated peak_gamma = 64011.77`
  - `g_lo` 从 `64011.77` 改善到 `58513.84`

这说明：
- 初始谱投影的硬切确实是误差来源之一
- 但它主要影响的是支撑区起点，而不是谱峰主位置
- `peak_gamma` 锁网格的更强来源仍然是后续演化中的：
  1. 冻结冷却系数
  2. 相邻格点的离散竞争
  3. 单边回代和正性截断

### 3.6 冻结冷却系数的离线量化结果

最新结果见：
- `output/vegasafterglow_doc/electron_cooling_lag.json`

实验定义：
- 不改默认链
- 只在单个壳层步里，把当前 `fullhide` 的冻结冷却更新，和“每个子步都刷新冷却系数”做直接对照
- 仍使用相同动力学输入、相同 `L1`、相同源项和相同隐式步进

当前结果：
- `shell_frac = 0.3`
  - `refresh_vs_frozen ≈ 1e-11 ~ 4e-11`
- `shell_frac = 0.5`
  - `refresh_vs_frozen ≈ 1.8e-7 ~ 2.1e-7`
- `frozen_matches_solver ≈ 3e-5`

这说明：
- 当前 Python 侧离线复现实验与真实 `fullhide` 输出是一致的
- 把冷却系数从“整壳冻结”改成“每个子步刷新”，对电子谱结果几乎没有影响
- 因此在当前默认 `index_Y=2` 主路径下，`lagged cooling / frozen coefficient` **不是电子谱低阶的主限**

直接推论：
- 前面把冻结冷却列为 `fullhide/t2g1` 共享主因，现在要下调其优先级
- 当前更可能的主限应转向：
  1. 电子谱峰位附近的格点竞争
  2. 单边隐式回代 + `max(zero, ...)` 正性截断
  3. 初始谱/源项/支撑区的格点投影与移动边界

### 3.7 正性截断 `max(zero, ...)` 的离线量化结果

最新结果见：
- `output/vegasafterglow_doc/electron_clipping.json`

实验定义：
- 不改默认链
- 只在单个壳层步里，把当前 `electron_backward_sweep(...)` 的正性截断版本，和“去掉 `max(zero, ...)` 的纯回代版本”做直接对照

当前结果：
- `unclipped_vs_clipped ≈ 1e-15`
- `negative_cells_unclipped = 0`
- `clipped_matches_solver ≈ 1e-5 ~ 4e-5`

这说明：
- 在当前默认 `fullhide` 主路径和代表性壳层上，`max(zero, ...)` 基本没有实际触发
- 因此正性截断本身**不是当前峰位锁网格的主因**
- 单边隐式回代的“方向性与格点竞争”仍然可能重要，但不是因为这个硬截断语句单独起作用

### 3.8 源项硬切 / 源项亚格点投影的离线量化结果

最新结果见：
- `output/vegasafterglow_doc/electron_source_projection.json`

实验定义：
- 不改默认链
- 只在单个壳层步里，把当前 `electron_build_source_term(...)` 的硬窗口源项，和 `electron_build_source_term_profile(...)` 的亚格点积分源项做直接对照

当前结果：
- `shell_frac = 0.3`
  - `integrated_vs_hard ≈ 0.16 -> 0.41`
- `shell_frac = 0.5`
  - `integrated_vs_hard ≈ 0.48 -> 0.50`
- 但大多数 case 下：
  - `peak_gamma` 几乎不变
  - `g_lo / g_hi` 也基本不变

这说明：
- 源项硬切对幅度有显著影响，不能随便直接替换进默认链
- 但它对当前观察到的 `peak_gamma` 锁网格并不是主控制因素
- 因此“把源项改成亚格点积分”不是当前默认链的第一修复手段

### 3.9 operator-only 诊断：默认输运核本身会锁峰位

最新结果见：
- `output/vegasafterglow_doc/electron_operator_only.json`

实验定义：
- 完全去掉源项
- 不刷新冷却，不引入额外耦合
- 只保留默认离散输运核
- 初始条件改成 `log(gamma)` 上的平滑单峰 Gaussian 单元平均，不再使用幂律初始谱

当前结果：
- 初始峰位在 `n = 121, 241` 时位于 `58513.84`
- 但经过 operator-only 演化后：
  - `fullhide final peak_gamma = 48894.06`
  - `t2g1 final peak_gamma = 48894.06`
- 且两条算法在这个实验里几乎给出同样的峰位和支撑区结果

这说明：
- 就算拿掉源项和冷却链，默认输运核本身也会把平滑单峰谱锁到相同格点
- 因此当前 `fullhide/t2g1` 的峰位问题，主因已经可以确认在：
  1. 单边隐式输运离散本体
  2. 相邻格点竞争与数值耗散
- 不是源项、不是冷却滞后、也不是 `max(zero, ...)` 触发造成的假象

### 3.10 纯 operator 定阶结果

最新结果见：
- `output/vegasafterglow_doc/electron_operator_scheme.json`

实验定义：
- 完全去掉源项和额外耦合
- 只保留默认输运核
- 初始条件是 `log(gamma)` 上平滑单峰 Gaussian
- 比较数值解与同样平移量下的解析平移解

当前结果：
- `fullhide-operator-only = 0.2681`
- `t2g1-operator-only = 0.3205`

同时峰位轨迹显示：
- 初始峰位：`1584.89`
- 解析峰位：`251.19`
- 数值峰位：
  - `61 -> 316.23`
  - `121 -> 281.84`
  - `241 -> 266.07`

这说明：
- 默认输运核在最干净的平滑平移问题上，收敛阶本身就非常低
- `t2g1` 相比 `fullhide` 略有改善，但仍远没有达到其名义时间阶数
- 因此当前主问题已经可以进一步收敛为：
  1. 单边隐式上风离散本体过于耗散
  2. 该离散的有效误差在空间上支配了整个电子谱链
  3. 继续只抠外围耦合不会把主阶数抬起来

### 3.11 测阶脚本原先过慢，现已改成缓存

文件：
- `tests/order_convergence_check.py`

已完成：
- 对 `dynamics / observed / direct` 三类 case 增加 `npz` 缓存
- 同一 `solver + n + mode` 不再重复整链求解

当前实测：
- 首轮填缓存仍然慢
- 复跑时间已降到约 `2.91 s`

结论：
- 之后讨论电子谱阶数时，可以高频复跑 `order_convergence_check.py`
- 当前慢点不再是测试脚本，而是数值核本身

## 4. 当前阶段目标

### 阶段 1：把电子谱低阶问题拆干净

目标：
- 把当前电子谱误差拆成三个量分别测：
  1. `peak_gamma` 误差
  2. `g_lo / g_hi` 误差
  3. 峰位对齐后的 shape error

验收：
- 新测试能明确回答：
  - 主要误差是峰位漂移
  - 还是支撑区漂移
  - 还是归一化后形状误差

### 阶段 2：量化 lagged cooling / frozen coefficient 的贡献

目标：
- 不直接改默认链
- 先做离线诊断或实验分支
- 定量比较：
  - 当前冻结冷却系数
  - 子步内更新冷却系数
  对电子谱有效阶数的影响

验收：
- 能清楚判断：
  - 时滞冻结是不是主限
  - 贡献量级有多大

### 阶段 3：推进默认链的真实修复

只允许按下面顺序推进：
1. 先修共享低阶源
2. 再修 `t2g1`
3. 最后再继续修 `weno5`

#### 3A. 共享低阶源

候选方向：
- `gamma` 网格的分段/自适应分布
- `gamma_m / gamma_c / gamma_max` 附近的亚格点处理
- 减少谱支撑区在网格上的锁定漂移

#### 3B. `t2g1`

候选方向：
- 减少启动段低阶回退影响
- 评估是否能在不破坏主链的前提下抬高有效阶数

#### 3C. `weno5`

候选方向：
- 严格对照 `WENO-Z`
- 检查通量分裂是否真正与方程一致
- 修边界 stencil
- 再重新测：
  - 电子谱阶数
  - `bands_flux`
  - `nu_a`

## 5. 基准更新规则

当前规则：
- **只要默认链电子谱有效阶数仍明显偏低，就不更新主基准。**

更新基准前必须同时满足：
- `regression_check` 通过
- `vegasafterglow_comprehensive_validation` 通过
- `bug_audit_check` 通过
- `physical_closure_check` 通过
- 默认电子链的电子谱有效阶数达到可接受水平
- 结果比当前稳定基线更好，而不是单纯“也能跑”

当前状态：
- **不更新新的主基准**

## 6. 每轮优化必须执行的检查

若改动 Fortran：
- `python build_extensions.py --module <touched> --force`
- `rg -n ".{133,}" <touched fortran files>`

每轮都必须执行：
- `python tests/regression_check.py`
- `python tests/vegasafterglow_comprehensive_validation.py`
- `python tests/bug_audit_check.py`
- `python tests/physical_closure_check.py`

涉及电子谱或阶数判断时，额外执行：
- `python tests/order_convergence_check.py`
- `python tests/electron_solver_comparison.py`

## 7. 下一步直接行动项

1. 修改 `tests/order_convergence_check.py`
   - 把电子谱误差拆成：
     - `peak_gamma`
     - `g_lo / g_hi`
     - 对齐峰位后的 shape error

2. 为 lagged cooling / frozen coefficient 加离线诊断
   - 不直接并入默认链
   - 先测贡献量级

3. 在证据足够前，不继续改主基准文件

4. 新的主攻顺序改为：
   1. 先解释 `fullhide/t2g1` 的 `peak_gamma ≈ 0`
   2. 再解释对齐峰位后为什么 shape 仍只有 `~0.5`
   3. 最后再回到 `weno5` 的高能支撑区误差

5. 新的默认链修复顺序改为：
   1. 先做 `electron_initial_powerlaw(...)` 的亚格点投影实验诊断
   2. 再做 `lagged cooling / frozen coefficient` 的离线贡献量化
   3. 只有在前两项证据清楚后，才回到默认链的实际修复
   4. 如果冻结冷却是主限，再决定是否进入默认链的 predictor-corrector 或子步中点更新

6. 根据最新证据，新的默认链主攻顺序调整为：
   1. 优先重写或替换默认链的单边隐式输运离散
   2. 重点比较 `fullhide` 和 `t2g1` 的 operator-only 演化差异，确认二者是否仅在时间骨架上略有差别
   3. 再检查初始谱边界在峰位锁定中的贡献
   4. 源项投影只保留为幅度误差来源，不再作为峰位主限
   5. 冻结冷却和 `max(zero, ...)` 都降级为次级因素
## 7. 最新进展（2026-04-15）

### 7.1 `slc1` 的 O(N) 前缀积分 remap 已恢复并保留

文件：
- `src/Electron/electron_common.f90`

本轮改动：
- 恢复了 `slc1` 半拉格朗日输运里的 O(N) 前缀积分式守恒重映射。
- 现在通过
  - `electron_linear_prefix_integral(...)`
  - `electron_prefix_integral_eval(...)`
  先构造分片线性重构的单元积分前缀和，再用前缀差值做回溯区间积分。
- `electron_semi_lagrangian_transport(...)` 已切回这条 O(N) 路径。

结论：
- 这次恢复后，`slc1` 的 benchmark 和 stage-level profile 都明显改善。
- 当前树里保留 O(N) 版本，不再回退到逐区间积分版。

### 7.2 `slc1` 当前 observed order

结果文件：
- `output/vegasafterglow_doc/order_convergence.json`

当前 `slc1`：
- `bands-flux ≈ 0.622`
- `nu_a ≈ 0.655`
- `electron-spectrum ≈ 0.330`
- `electron-peak-gamma ≈ 1.066`
- `electron-shape-aligned ≈ 1.475`

解释：
- 这些是全链路 observed order，不是 formal order。
- `slc1` 在峰位项和对齐后的谱形项上，已经明显优于 `fullhide/t2g1`。
- 但整条电子谱误差仍大，所以还不能替代默认链。

### 7.3 `slc1` 最新 benchmark / profile / 图

结果文件：
- `output/vegasafterglow_doc/slc1_benchmark_compare.json`
- `output/vegasafterglow_doc/slc1_benchmark_compare.png`
- `output/vegasafterglow_doc/slc1_subroutine_profile.json`
- `output/vegasafterglow_doc/slc1_subroutine_profile.png`
- `output/vegasafterglow_doc/slc1_quick-lc.png`
- `output/vegasafterglow_doc/slc1_quick-spec.png`

关键 benchmark：
- `quickstart: 1.497 s -> 1.314 s`
- `lightcurve_grid: 2.587 s -> 1.565 s`
- `spectrum_grid: 2.220 s -> 1.469 s`
- `pair_points: 2.999 s -> 1.576 s`
- `exposures: 2.461 s -> 1.476 s`
- `details: 2.309 s -> 1.404 s`
- `fitter_cfg: 3.584 s -> 3.137 s`
- `sky_image: 7.282 s -> 4.201 s`

stage-level profile（Tophat sync）：
- `fullhide total ≈ 2.583 s`
- `slc1 total ≈ 1.212 s`
- 其中电子核：
  - `Electron.fs_electron_fullhide ≈ 2.501 s`
  - `Electron.fs_electron_slc1 ≈ 1.137 s`

直接结论：
- 当前最快的继续提速路线仍然是 `slc1`
- 当前最需要补的是 `slc1` 的精度，不是继续给 `fullhide` 加线程

### 7.4 `slc1` 当前主精度瓶颈判断

当前确认结论：
- `slc1` 的主精度瓶颈不是 limiter 本身
- 主因是时间耦合和冻结界面速度
- 因此把源项耦合改成 `Strang` 式半步-输运-半步，这个方向是对的

后续含义：
- 下一步不再优先更换 `MC limiter`
- 优先处理：
  1. 界面速度的时间层更新方式
  2. characteristic / departure point 的时间中点或 predictor-corrector
  3. 在此基础上再评估是否有必要把 remap 从线性 TVD 升到 `PPM/parabolic`

### 7.5 `slc1` 固定 `L1=10` 的试验结论

试验做法：
- 在 `src/Electron/FS_electron_slc1.f90` 里把子步数临时固定成 `L1=10`
- 只做 `slc1` 的电子谱对照、benchmark 和 stage-level profile

结果：
- 速度显著提升：
  - `slc1 quickstart ≈ 0.682 s`
  - `slc1 lightcurve_grid ≈ 0.716 s`
  - `slc1 spectrum_grid ≈ 0.711 s`
  - `slc1 sky_image ≈ 2.130 s`
  - `Tophat sync total ≈ 0.711 s`
  - `Electron.fs_electron_slc1 ≈ 0.634 s`
- 但精度直接坏掉：
  - `direct rel_d_n_gam_e ≈ 55.0`
  - `observed rel_bands_flux ≈ 14.46`

结论：
- `L1=10` 不能保留
- 这条路线证明了：当前 `slc1` 真正缺的不是更多并行，而是“用更少子步仍保持精度”的高阶时间耦合
- 已回退此试验，不保留在树里

### 7.6 `slc1` predictor-corrector 时间耦合试验结论

试验做法：
- 在 `slc1` 子步里引入左右时间层的界面速度和源项
- 用 predictor-corrector 的方式校正 characteristic / departure point

结果：
- 没有形成净收益
- 误差几乎不降：
  - `direct rel_d_n_gam_e` 仍约 `1.80`
  - `observed rel_bands_flux` 仍约 `0.26`
- 速度明显回退：
  - `slc1 observed seconds` 从约 `0.65` 回到约 `0.74`
  - `slc1 TopHat sync total` 从约 `1.21 s` 回到约 `2.21 s`

结论：
- 这版 predictor-corrector 不保留，已经回退
- 当前最该继续做的不是简单多做一次 transport，而是：
  1. 更高阶的守恒 remap
  2. 或者更准确但不翻倍 transport 次数的 characteristic 时间中点近似

### 7.7 `slc1` PPM/parabolic remap 当前状态

本轮在 `src/Electron/electron_common.f90` 中把 `slc1` 的守恒 remap 从线性 TVD 升到 `PPM/parabolic`。

已确认问题：
- 第一版编译失败不是数值公式问题，而是 f2py 生成的 wrapper 名过长。
- 根因是 `electron_parabolic_prefix_integral_eval` 这类过程名在 wrapper 中拼接后超长，导致 `FS_electron_slc1-f2pywrappers2.f90` 解析失败。
- 已通过缩短过程名修复：
  - `electron_cell_parabolic_integral -> electron_ppm_cell_int`
  - `electron_parabolic_prefix_integral -> electron_ppm_prefix`
  - `electron_parabolic_prefix_integral_eval -> electron_ppm_prefix_eval`

当前状态：
- `FS_electron_slc1 / FS_electron_fullhide / FS_electron_t2g1` 已重新编译通过。
- line truncate 检查通过。
- 主链门槛全部通过：
  - `regression_check`
  - `vegasafterglow_comprehensive_validation`
  - `bug_audit_check`
  - `physical_closure_check`

当前数值/速度结果：
- `electron_solver_comparison.json`
  - `slc1 direct rel_d_n_gam_e ≈ 1.802`
  - `slc1 observed rel_bands_flux ≈ 0.260`
  - `slc1 observed seconds ≈ 0.672`
- `slc1_benchmark_compare.json`
  - `quickstart ≈ 1.030 s`
  - `lightcurve_grid ≈ 1.229 s`
  - `spectrum_grid ≈ 1.129 s`
  - `pair_points ≈ 1.238 s`
  - `exposures ≈ 1.190 s`
  - `details ≈ 1.090 s`
  - `fitter_cfg ≈ 2.368 s`
  - `sky_image ≈ 3.406 s`
- `slc1_subroutine_profile.json`
  - `Tophat sync total ≈ 1.608 s`
  - `Electron.fs_electron_slc1 ≈ 1.526 s`

当前 observed order：
- `slc1-bands-flux ≈ 0.622`
- `slc1-nu_a ≈ 0.681`
- `slc1-electron-spectrum ≈ 0.331`
- `slc1-electron-peak-gamma ≈ 1.066`
- `slc1-electron-shape-aligned ≈ 1.473`

结论：
- 这版 `PPM/parabolic` remap 已经能稳定编译、运行、过主链。
- 速度明显优于 `fullhide`，并且在 benchmark 上继续占优。
- 但精度没有出现决定性改善，主精度瓶颈仍然不是 limiter 本身，而是时间耦合和冻结界面速度。
- 下一步不再优先继续堆更复杂的 limiter，而是优先审查 `slc1` 的界面速度时间层和 departure-point 时间精度，同时保住这版可编译的 PPM 基线。

补充试验结论：
- 已试过“左右时间层界面速度 + 中点回溯”的廉价修正，不增加第二次 transport。
- 结果是：
  - `electron_solver_comparison` 中 `slc1` 误差几乎不变；
  - `slc1` benchmark 反而变慢；
  - 没有形成净收益。
- 这条路线已回退，不保留在树里。

### 7.8 星风环境 benchmark

新增脚本：
- `tests/vegasafterglow_slc1_wind_benchmark.py`

输出文件：
- `output/vegasafterglow_doc/slc1_wind_benchmark_compare.json`
- `output/vegasafterglow_doc/slc1_wind_benchmark_compare.png`
- `output/vegasafterglow_doc/slc1_wind_benchmark_compare.pdf`

环境：
- `TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0)`
- `Wind(A_star=0.1, n0=1.0)`
- `Observer(1.0e26, 0.1, 0.0)`
- `Radiation(0.1, 1.0e-3, 2.3)`

结果：
- `quickstart`: `fullhide 1.547 s`, `slc1 1.293 s`
- `lightcurve_grid`: `fullhide 2.380 s`, `slc1 1.543 s`
- `spectrum_grid`: `fullhide 2.156 s`, `slc1 1.441 s`
- `pair_points`: `fullhide 2.357 s`, `slc1 1.606 s`
- `exposures`: `fullhide 2.163 s`, `slc1 1.429 s`
- `details`: `fullhide 2.096 s`, `slc1 1.367 s`
- `fitter_cfg`: `fullhide 3.512 s`, `slc1 2.952 s`
- `sky_image`: `fullhide 6.093 s`, `slc1 4.146 s`

结论：
- 在星风环境下，`slc1` 依然全面快于 `fullhide`。
- `sky_image` 的 flux ratio 在两条求解器下都约 `0.826 ~ 0.827`，这更像当前 wind case 的视场覆盖不足，而不是电子求解器差异。

### 7.9 电子谱 support-high 测阶指标已修正

修改：
- `tests/order_convergence_check.py`
- 把 `g_lo / g_hi` 支撑边界从“离散阈值格点”改成“连续阈值交点”估计。

结果：
- 之前 `slc1-electron-support-high` 的 `coarse_error = 0 -> NaN` 已经消失。
- 现在剩下的是实问题，不是指标定义问题。

当前剩余 failed 项：
- `fullhide-electron-peak-gamma`
- `t2g1-electron-peak-gamma`
- `slc1-electron-support-high`
- `weno5-electron-support-high`

新的结论：
- `slc1` 的高能支撑区现在测出来是负阶，说明随着分辨率提高，高能截止附近的误差没有单调下降。
- `weno5` 的高能支撑区也有同样问题。
- 这把当前主线进一步收窄到：
  1. `fullhide/t2g1` 的峰位锁网格；
  2. `slc1/weno5` 的高能截止与高能支撑区处理。

## 2026-04-15 update

- Added a high-energy exponential injection tail only for `fullhide`, `slc1`, and `weno5`.
- Accepted form:
  - original power law for `gamma <= gamma_max`
  - exponential tail for `gamma > gamma_max`
  - tail factor: `exp(1 - gamma / gamma_max)`
- Rejected form:
  - multiplying the whole injected spectrum by `exp(-gamma/gamma_max)`
  - reason: it changed the regression baseline too much
- New common routines in `src/Electron/electron_common.f90`:
  - `electron_initial_powerlaw_exp_cutoff`
  - `electron_build_source_term_exp_cutoff`
- New benchmark output folder:
  - `output/benchmark_exp_tail`
- New benchmark/plot entry:
  - `tests/electron_exp_tail_benchmark_suite.py`
- Current benchmark scope:
  - solvers: `fullhide`, `slc1`
  - media: `ISM`, `Wind`
  - forward-shock only

## 2026-04-15 slc1_mmg2 implementation

- Added experimental solver `slc1_mmg2` in `src/Electron/FS_electron_slc1.f90`.
- Core design:
  - fixed-N moving mesh on `x = log10(gamma)`
  - one mesh update per shell via monitor/equidistribution
  - conservative remap on nonuniform cells
  - Strang source split + semi-Lagrangian transport on working mesh
  - shell-end remap back to public uniform gamma grid for compatibility
- Added nonuniform-grid helpers in `src/Electron/electron_common.f90`:
  - `electron_gamma_from_edges`
  - `electron_mc_slopes_nonuniform`
  - `electron_ppm_interfaces_nonuniform`
  - `electron_ppm_prefix_nonuniform`
  - `electron_ppm_prefix_eval_nonuniform`
  - `electron_conservative_remap_nonuniform`
  - `electron_semi_lagrangian_transport_nonuniform`
  - `electron_semi_lagrangian_step_nonuniform`
  - `electron_cell_geometry`
  - `electron_build_moving_monitor`
  - `electron_equidistribute_edges`
- Runtime wiring:
  - added `electron_solver = "slc1_mmg2"` to `asgard_runtime.py`
  - exported in `src/Electron/__init__.py`
- Safety fix:
  - `get_syn_simpson` now falls back to `get_syn` when `gamma` grid is not uniform in `log(gamma)`.
- Benchmark configuration update:
  - `fullhide` uses `num_gam_e = 121`
  - `slc1` and `slc1_mmg2` use `num_gam_e = 32`
  - patched benchmark/profile/spectrum scripts to compare `fullhide / slc1 / slc1_mmg2`
- Current results:
  - regression passed: `xrt=0.010782`, `optr=0.000243`, `9GHz=0.000001`
  - comprehensive validation passed: `regression_failed = 0/99`
  - bug audit passed: `0/40`
  - physical closure passed: `0/44`
  - stage profile (`Tophat sync`):
    - `fullhide total ~ 2.689 s`
    - `slc1 total ~ 0.919 s`
    - `slc1_mmg2 total ~ 0.779 s`
  - electron kernel:
    - `Electron.fs_electron_fullhide ~ 2.613 s`
    - `Electron.fs_electron_slc1 ~ 0.830 s`
    - `Electron.slc1_mmg2 ~ 0.705 s`
  - benchmark (ISM):
    - `quickstart`: `1.439 -> 0.343 -> 0.330 s`
    - `lightcurve_grid`: `1.930 -> 0.447 -> 0.420 s`
    - `spectrum_grid`: `1.641 -> 0.354 -> 0.350 s`
  - benchmark (Wind):
    - `quickstart`: `1.504 -> 0.452 -> 0.422 s`
    - `lightcurve_grid`: `2.159 -> 0.654 -> 0.592 s`
    - `spectrum_grid`: `1.841 -> 0.510 -> 0.468 s`
  - observed order with `20/32/40` for `slc1_mmg2`:
    - `bands-flux ~ 0.545`
    - `nu_a ~ 1.121`
    - `electron-spectrum ~ 1.624`
    - `peak-gamma ~ 0.843`
    - `support-low ~ 1.856`
    - `support-high ~ -0.559`
    - `shape-aligned ~ 1.824`
- Conclusion:
  - moving mesh improved speed and electron-shape alignment
  - target 10-20% speed of current fullhide is reached for many benchmark queries, but not yet for stage-level `Tophat sync`
  - second-order target is not yet met on `bands-flux`, `peak-gamma`, and `support-high`
  - next priority: high-energy cutoff tracking and high-energy support handling in `slc1_mmg2`

## [2026-04-15 18:14:18] Power-law-grid radiation integrator update
- Replaced recursive cellwise Simpson in calling_modules.f90:get_syn_adaptive with log-gamma cell integration using 2/3-point Gauss-Legendre plus one local split when the cellwise 2-point and 3-point estimates disagree above el_tol=5e-4.
- Switched get_nu_a tau evaluation onto the same x=ln(gamma) cell quadrature, removing the old midpoint-only tau accumulation from the adaptive path.
- Result: main-chain tests still pass (egression_check, ug_audit_check, physical_closure_check, egasafterglow_comprehensive_validation).
- Observed effect:
  - ullhide_adaptive direct el_p_syn ~ 2.82e-2, observed el_bands_flux ~ 2.49e-2.
  - slc1 and slc1_mmg2 speed recovered relative to the previous recursive adaptive integrator, but ands_flux order is still dominated by electron-front/high-energy-cutoff errors rather than the synchrotron quadrature itself.
- Next focus: keep the log-grid cell quadrature, but do not expand expensive global adaptivity further before fixing slc1_mmg2 high-energy support tracking.


- Workflow note (2026-04-15): from now on, external validation focus is reduced to three items only: order-convergence, light-curve plots, and spectrum plots. Fortran compile and line-truncation checks remain mandatory internal gates.

## [2026-04-15 19:00:41] slc1_mmg2 nonuniform cell-average initialization/source update
- Added nonuniform-cell average construction for the exponential-cutoff initial electron spectrum and source term in `src/Electron/electron_common.f90`.
- New routines:
  - `electron_initial_powerlaw_exp_cutoff_edges(...)`
  - `electron_build_source_term_exp_cutoff_edges(...)`
- Implementation detail:
  - the integrals are taken in `x = log10(gamma)` on the actual moving-mesh cell edges;
  - each active piece of the broken power-law is integrated cellwise with 3-point Gauss-Legendre;
  - cells that cross `gamma_max` are split at `x_max = log10(gamma_max)` before quadrature, so the exponential high-energy tail is averaged instead of sampled at the cell center.
- These new routines are wired only into `fs_electron_slc1_mmg2(...)`; legacy solvers remain unchanged.
- Validation policy updated:
  - `order_convergence_check.py` now only measures the new solver `slc1_mmg2`;
  - cache version bumped to `7` to invalidate old convergence caches after the new source/initialization change.
- Current observed-order result for `slc1_mmg2` with `Num_gam_e = 20/32/40`:
  - `bands-flux ≈ 2.295`
  - `nu_a ≈ 3.637`
  - `electron-spectrum ≈ 2.797`
  - `peak-gamma ≈ 2.847`
  - `support-low ≈ 0.965`
  - `support-high ≈ 3.923`
  - `shape-aligned ≈ 2.414`
- Interpretation:
  - the high-energy cutoff/front tracking problem is largely improved;
  - the remaining weak point is the low-energy support location (`support-low`), not the high-energy end anymore.

## [2026-04-15 19:19:00] slc1_mmg2 low-energy front anchoring attempt
- Added low-energy front extraction and moving-mesh control:
  - `electron_find_low_energy_front(...)`
  - `electron_anchor_low_energy_edges(...)`
  - `electron_rescale_lower_boundary(...)`
- Wired these only into `fs_electron_slc1_mmg2(...)`:
  - dynamic lower-boundary rescaling before monitor/equidistribution
  - low-energy front anchoring after equidistribution
  - moving-mesh monitor now includes a Gaussian weight around the detected low-energy front
- Validation outcome after recompilation and rerunning the reduced workflow (`order_convergence`, light-curve plots, spectrum plots):
  - `support-low` stayed at `~0.965`
  - all other `slc1_mmg2` observed-order quantities were effectively unchanged
- Conclusion:
  - the current low-energy error is not primarily caused by missing low-energy edge anchoring or wasted empty cells below the spectrum;
  - the next place to inspect is the low-energy feature definition and the final diagnostic/public-grid remap, not the moving-mesh boundary placement itself.

## [2026-04-15 19:31:00] slc1_mmg2 log-positive conservative remap attempt
- Replaced the `slc1_mmg2`-only nonuniform remap calls with a new positive log-space conservative remap:
  - `electron_conservative_remap_log_nonuniform(...)`
  - based on a monotone linear reconstruction of `log10(1+q)` rather than direct PPM on `q`
- Scope:
  - only `fs_electron_slc1_mmg2(...)`
  - used for work-grid-to-work-grid remap and final work-grid-to-public-grid export
  - legacy solvers unchanged
- Motivation:
  - the remaining weak point was `support-low`, and the working hypothesis was that steep low-energy fronts were being smeared by direct-`q` PPM remap on nonuniform grids
- Outcome after recompilation and rerunning the reduced workflow:
  - `support-low` remained `~0.965`
  - all other `slc1_mmg2` observed-order quantities remained effectively unchanged
- Conclusion:
  - the remaining low-energy support error is not primarily due to the choice between direct-`q` PPM remap and log-positive monotone remap;
  - the next algorithmic target should move away from remap details and toward the low-energy feature definition itself, likely the `gamma_m/gamma_c` transition structure and/or the way low-energy support is diagnosed from the exported spectrum.

## [2026-04-15 19:43:00] characteristic-energy anchoring for slc1_mmg2
- Reverted the `slc1_mmg2`-only log-positive conservative remap experiment after it produced wrong light curves and spectra.
- Restored the standard nonuniform conservative remap and changed strategy to direct characteristic-energy anchoring.
- New routine:
  - `electron_anchor_characteristic_edges(...)`
- Scope:
  - used only inside `fs_electron_slc1_mmg2(...)`
  - after moving-mesh monitor/equidistribution and low/high-front anchoring
  - targets `x_break = log10(min(gamma_m, gamma_c))`
- Rationale:
  - the remaining unresolved issue was low-energy support and low-energy spectral structure;
  - direct anchoring of the characteristic break is more aligned with the physics goal than threshold-front anchoring or alternative remap formulas.
- Validation outcome after recompilation and rerunning the reduced workflow:
  - `bands-flux ≈ 2.295`
  - `nu_a ≈ 2.614`
  - `electron-spectrum ≈ 1.958`
  - `peak-gamma ≈ 1.484`
  - `support-low ≈ 3.386`
  - `support-high ≈ 3.538`
  - `shape-aligned ≈ 2.306`
- Interpretation:
  - the low-energy support problem is now substantially improved;
  - the remaining weaker items are `peak-gamma` and the mixed whole-spectrum metric, not the support edges;
  - future work should continue on characteristic-energy/break tracking, not on remap replacement.

## [2026-04-15 20:28:00] radiation-integral update on logarithmic frequency cells
- Scope:
  - `src/Electron/calling_modules.f90`
  - updated `get_SSA_numerical(...)`, `get_IC_numerical(...)`, and `get_Y_Nakar(...)`
- Change:
  - removed the remaining linear-`dnu` trapezoid bias on the seed-frequency grid;
  - switched the frequency-cell handling to logarithmic-cell quadrature / midpoint form aligned with the power-law frequency mesh;
  - added:
    - `electron_powerlaw_interp(...)`
    - `electron_log_gauss2_interval(...)`
    - `electron_integrate_powerlaw_segment(...)`
    - `electron_ssa_segment(...)`
- Details:
  - `get_Y_Nakar(...)` now builds the cumulative synchrotron integral with logarithmic-cell segment integrals instead of linear-grid trapezoids;
  - `get_SSA_numerical(...)` now integrates each active frequency cell on the log-frequency domain and explicitly splits cells across `nu_low` / `nu_up`;
  - `get_IC_numerical(...)` now uses logarithmic cell midpoints and `Δnu = nu_mid * Δln(nu)` for both outer and inner photon-frequency loops.
- Reduced validation workflow:
  - recompiled `FS_electron_slc1`, `FS_electron_fullhide`, `FS_electron_t2g1`, `FS_electron_weno5`
  - reran only:
    - `tests/order_convergence_check.py`
    - `tests/electron_exp_tail_doc_style_plots.py`
    - `tests/electron_exp_tail_spectrum_compare.py`
- Current observed order after cache refresh:
  - `bands-flux ≈ 3.065`
  - `nu_a ≈ 1.297`
  - `electron-spectrum ≈ 2.089`
  - `peak-gamma ≈ 1.484`
  - `support-low ≈ 3.015`
  - `support-high ≈ 3.789`
  - `shape-aligned ≈ 2.306`
- Interpretation:
  - the radiation-step change did not destabilize the new solver;
  - the mixed spectral metric crossed the two-order target;
  - the remaining weakest quantities are now `nu_a` and `peak-gamma`;
  - the next radiation-side target should be the SSA / `nu_a` chain, not the synchrotron core integral.

## [2026-04-15 21:35:00] nonuniform SSC / gamma-gamma follow-up and wider SED plots
- Scope:
  - `src/Radiation/SSC_spec.f90`
  - `src/Radiation/Annihilation.f90`
  - `src/Radiation/radiation_common.f90`
  - `tests/electron_exp_tail_doc_style_plots.py`
  - `tests/electron_exp_tail_spectrum_compare.py`
  - `tests/electron_exp_tail_benchmark_suite.py`
- Implemented:
  - added `ssc_spec_nonuniform(...)` for electron-work-grid SSC integration on `x = log10(gamma)` cell edges;
  - upgraded the gamma-gamma annihilation frequency-cell handling from linear-width midpoint style to logarithmic-cell midpoint / power-law interpolation;
  - widened the SED plotting range to `10^7 -- 10^30 Hz`;
  - switched the comparison spectra to `SED = nu F_nu`;
  - switched the doc plots to the IC-dominated parameter group requested in this round.
- Important runtime note:
  - the stable high-level runtime still loads the legacy `FS_electron_slc1_v3` binary first;
  - therefore the new work-grid-history radiation path is not yet the default high-level execution path;
  - this is intentional for now because the current source-built `FS_electron_slc1` branch is still experimental and should not replace the stable chain before passing the hard gates.
- Validation:
  - Fortran compile checks passed for:
    - `FS_electron_slc1`
    - `FS_electron_fullhide`
    - `FS_electron_t2g1`
    - `SSC_spec`
    - `Annihilation`
  - line-truncation checks passed on all touched Fortran files.
  - `tests/order_convergence_check.py` was rerun on the stable runtime path and stayed at the previous stable result:
    - `bands-flux ~ 4.751`
    - `nu_a ~ 1.472`
    - `electron-spectrum ~ 2.042`
    - `peak-gamma ~ 1.484`
    - `support-low ~ 3.081`
    - `support-high ~ 3.861`
    - `shape-aligned ~ 2.306`
- Direct SSC-kernel diagnostic:
  - on a synthetic smooth distribution, `ssc_spec_nonuniform(...)` is finite and converges toward the old uniform-grid `ssc_spec(...)` as the grid is refined;
  - measured median relative difference:
    - `81 x 121`: `~2.53e-1`
    - `161 x 241`: `~1.11e-1`
    - `321 x 481`: `~6.11e-2`
  - but the worst core-region discrepancy is still large (`~4.4e-1`), so this kernel is not yet ready to replace the old SSC path.
- Interpretation:
  - the remaining SSC error is no longer dominated by seed-frequency quadrature;
  - the main unresolved term is the piecewise-constant treatment of `dN_x` inside each electron cell;
  - the next correct step is to add intra-cell electron-spectrum reconstruction for the nonuniform SSC kernel, not to keep widening the seed-frequency quadrature alone.

## [2026-04-15 22:05:00] nonuniform SSC kernel accuracy-vs-cost checkpoint
- Follow-up on `ssc_spec_nonuniform(...)`:
  - added limited linear intra-cell reconstruction of `dN_x(x)` inside each nonuniform electron cell;
  - replaced the tail moments with analytic cellwise `gamma^-2` / `gamma^-4` moment integrals of the reconstructed profile.
- Compile status:
  - `SSC_spec` rebuild passed;
  - line-truncation check passed.
- Diagnostic against a high-resolution old-kernel reference on a synthetic smooth distribution:
  - for `N_gamma = 81`, the new nonuniform kernel gives a smaller overall median error than the old coarse uniform kernel;
  - however the core-emission region still does not improve enough, so the present reconstruction is not yet sufficient.
- Micro-benchmark on a synthetic `81 x 121 x 64` problem:
  - old `ssc_spec(...)`: `~3.63e-2 s`
  - current `ssc_spec_nonuniform(...)`: `~1.87e1 s`
  - slowdown: `~5.17e2 x`
- Interpretation:
  - the direct nonuniform SSC path is still far too expensive for the hard speed target;
  - the remaining bottleneck is the low-energy branch, where the current implementation still does a large per-frequency, per-seed-cell, per-electron-cell sweep with expensive transcendental work inside the innermost loops;
  - the next correct optimization target is not more quadrature tuning, but a structural reduction of repeated low-branch work:
    - precompute output-frequency-dependent kernel coefficients;
    - reduce repeated binary searches / logs / exponentials inside the innermost loops;
    - move toward cell-window / prefix-sum style evaluation for the dominant low-branch contribution.

## [2026-04-15 22:42:00] auxiliary-grid SSC scheme test
- Scope:
  - `tests/ssc_aux_grid_scheme_check.py`
- Goal:
  - test the mixed route:
    - keep `slc1_mmg2` on the moving electron work grid;
    - keep synchrotron / `nu_a` on the work grid;
    - project each shell's electron spectrum onto a small auxiliary uniform `log(gamma)` grid for `SSC`, then reuse the old fast `ssc_spec(...)`.
- Tested reference:
  - source-built `src.Electron.FS_electron_slc1.fs_electron_slc1_mmg2(...)`
  - `Num_gam_e(work) = 32`
  - IC-dominated parameter group
  - reference SSC path: direct `ssc_spec_nonuniform(...)`
- Tested auxiliary projections:
  - `conservative`: conservative cell-average projection to the auxiliary grid;
  - `center_log`: interpolate `log(dN_x)` to auxiliary cell centers, then convert to `dN/dgamma`;
  - `center_log_renorm`: same as `center_log`, then renormalize the shellwise `\int dN_x dx`.
- Key timing result:
  - direct `ssc_spec_nonuniform(...)`: `~1.67e1 s`
  - old public uniform SSC on the coarse work-grid output: `~2.87e-2 s`
  - auxiliary-grid SSC:
    - `center_log, N_aux=40`: `~3.57e-2 s`
    - `center_log, N_aux=64`: `~5.58e-2 s`
    - `center_log, N_aux=121`: `~1.02e-1 s`
  - therefore the auxiliary-grid route recovers the old SSC speed scale and stays `~1.6e2 -- 4.7e2 x` faster than the direct nonuniform kernel.
- Accuracy diagnostic against the direct nonuniform SSC reference:
  - conservative projection produces huge tail outliers and is not the right interface to the old `ssc_spec(...)`, which expects center values rather than cell averages;
  - `center_log` sharply reduces the worst outliers:
    - `N_aux=40`: `median_rel ~ 6.84e-1`, `core_max_rel ~ 3.47`
    - `N_aux=64`: `median_rel ~ 6.63e-1`, `core_max_rel ~ 3.50`
    - `N_aux=121`: `median_rel ~ 6.88e-1`, `core_max_rel ~ 4.56`
  - `center_log_renorm` slightly lowers the worst core outlier for some grids, but does not materially improve the median core error.
- Interpretation:
  - the mixed auxiliary-grid idea is viable on speed and should replace the current direct nonuniform SSC branch as the main recovery path;
  - the dominant unresolved error is no longer the auxiliary-grid size;
  - the main mismatch is the interface between moving-mesh cell averages and the old SSC kernel's center-value quadrature assumptions;
  - the next correct step is to implement the auxiliary-grid route with a dedicated shellwise reconstruction that feeds consistent point values into `ssc_spec(...)`, then validate against the high-resolution physical reference instead of the current slow kernel alone.

## [2026-04-15 23:35:00] high-energy observed synchrotron mismatch check
- Scope:
  - `asgard_component_backend.py`
  - high-energy forward-shock SED / light-curve mismatch among `fullhide`, `slc1`, and `slc1_mmg2`
- Key code comparison result:
  - `fullhide` and `slc1` share the same injection backbone:
    - `electron_gamma_m_exact(...)`
    - exponential-cutoff initial spectrum
    - exponential-cutoff source term
  - `slc1_mmg2` keeps the same physics and only changes the source projection to moving nonuniform cells;
  - `t2g1` is the solver with a genuinely different injection definition, so it should not be mixed into the current three-way interpretation.
- Root cause:
  - the order-of-magnitude mismatch at the high-energy end was not caused by different injection definitions inside `fullhide / slc1 / slc1_mmg2`;
  - the first attempted `nu_M`-based hard truncation was a wrong fix and has been reverted;
  - a single-shell `nu_M(t_obs)` cannot be used as a hard observed cutoff on EATS output, because one observed time receives contributions from different emission radii with different `Gamma`, different Doppler factors, and different high-energy cutoffs;
  - in addition, synchrotron turnover above `nu_c(gamma_max)` is smooth, not a hard step, so truncating the source spectrum at `nu_M_source` destroys the Wind high-energy spectral shape.
- Current status after revert:
  - the incorrect source-side hard mask in `asgard_component_backend.py` has been removed;
  - Wind high-energy synchrotron tails are smooth again;
  - the remaining three-solver high-energy mismatch must be traced without any hard `nu_M` clipping.
- Validation after revert:
  - Wind case, `t_obs ~ 1.826e5 s`:
    - `slc1` observed forward synchrotron decreases smoothly from `~3.76e-35` at `1e24 Hz` to `~4.75e-80` at `1e27 Hz`;
    - at `nu_obs ~ 2.418e26 Hz`, the forward synchrotron component is tiny but nonzero:
      - `fullhide ~ 4.31e-50`
      - `slc1 ~ 1.20e-49`
      - `slc1_mmg2 ~ 1.09e-51`
  - therefore the previous clipped shape was artificial, and the revert is physically required.

## [2026-04-16 00:10:00] ISM experimental-solver mismatch note
- New user-visible issue:
  - in `ISM`, the experimental moving-mesh solver `slc1_mmg2` is not only off in `SSC`;
  - its observed forward synchrotron spectrum is also inconsistent with `fullhide` / `slc1`, with a low-energy shape anomaly / bump and a clear peak shift toward lower frequency.
- Current diagnostic status:
  - with `ssc=False`, `slc1_mmg2` synchrotron peaks are still broadly comparable to the old solvers, but the low-energy shape is not yet trusted;
  - with `ssc=True`, the ISM mismatch becomes much larger because `SSC` dominates the total SED.
  - on the same `slc1_mmg2` electron state, `Radiation.ssc_spec_nonuniform(...)` is systematically lower than the old uniform-grid `ssc_spec(...)` path:
    - representative raw-source ratios are `~0.42 -- 0.74` in the core SSC band for the ISM benchmark;
    - therefore the current nonuniform SSC kernel is a confirmed active error source for the ISM total spectrum / light curve.
- Important interpretation:
  - the present ISM experimental mismatch cannot be blamed only on the moving-mesh electron solver itself;
  - at least one major piece of the error is inside the nonuniform SSC radiation kernel path.
- Plotting requirement added:
  - future spectrum-comparison figures must include electron-spectrum comparisons at the same spectrum snapshot times;
  - put the electron-spectrum comparison in the same figure as radiation SEDs, as subplots, so electron-shape errors and radiation-shape errors can be inspected together.
