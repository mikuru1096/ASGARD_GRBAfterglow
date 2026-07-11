# 电子求解器算法

本文比较电子求解器的离散角色。共享物理公式见算法总纲；公开组合限制见
[`public_backend_limits.md`](public_backend_limits.md)。

## 1. 共同方程

所有 1D solver 求解

\[
\partial_R N_e+\partial_\gamma(\dot\gamma_R N_e)=Q_R,
\]

并共享 shock injection、绝热/同步/SSA/IC cooling、同步辐射和 seed photon
定义。能量网格是等距 `log10(gamma)`；离散守恒量不是逐点 `dN/dgamma`，转换时
必须保留 cell Jacobian。

注入归一化由扫掠电子数和电子能量预算共同确定。首半径点建立初态，后续点只使用
相邻动力学状态间的局域步长，不能借用未来壳层。

## 2. `fullhide_1d`

默认生产路径采用隐式有限体积更新。面通量由 cooling velocity 与迎风状态构造，
三对角系统一次推进整个能量网格。它适合强冷却、宽能区和 hadronic/joint 主链。

`fullhide_1d_hz` 使用 Hz 辐射输出契约，输运离散不变。热电子只在已明确支持的
fullhide/DG 路径进入共享注入与辐射核。

## 3. Characteristic 与 SLC

`charint_1d` 沿 cooling characteristic 回溯分布，适合检查平滑损失下的谱形和
峰位。`slc1_1d` 使用 semi-Lagrangian conservative 更新，在大位移时减少显式
CFL 约束。两者是独立离散，不是 fullhide 的 fallback。

验收时比较积分粒子数、能量预算和 cooling break，而不是要求逐数组 bitwise。

## 4. T2G 与 WENO

`t2g1_1d` 使用二阶时间/空间更新；`weno5_1d` 用高阶重构解析平滑谱段。高阶路径
在边界和陡谱处仍必须遵守守恒通量，不允许通过滤波修改物理解。

这些 solver 不支持所有热电子、joint 或 2D 组合；公开请求在配置边界失败。

## 5. DG

`dg_1d` 在对数能量单元内使用局域多项式自由度和数值通量。质量矩阵、导数矩阵和
filter 表只读复用。filter 只抑制离散高阶模态，不承担修复负谱或不连续的职责。

DG 与有限体积输出的采样位置不同。比较前先投影到共同物理量，再检查粒子数、能量
和辐射积分。

## 6. 2D finite-q shell

`charint_2d` 和 `fullhide_2d` 求解 (N_e(\gamma,q,R))。q 是有限主动壳层的
质量坐标；每个 q cell 具有局域半径、bulk Lorentz factor、体积权重和磁场。

`chi_grid` 仅用于 BM 等效诊断。observer 的 χ-resolved 路径实际读取
`chi_radius_cm`、`chi_gamma_bulk` 和 `chi_dvolume_weight`。

能量输运与 q 空间输运分裂推进。默认 legacy 路径保持原 q 算子；`pwn_cr_v1`
可加入局域绝热损失、q 扩散、外边界逃逸和可选随机加速，详见
[`fullhide2d_pwn_cr_transport.md`](fullhide2d_pwn_cr_transport.md)。

## 7. 子步

固定子步用于可重复的基准和 joint 预算。自适应子步根据最大 cooling displacement
决定，但只在相邻输出半径之间细分，不改变输出网格。事件边界由动力学层先切分。

子步压缩只能合并物理输入相同的区间。任何优化都需证明注入总量、cooling integral
和辐射输出没有系统偏差。

## 8. 共享 kernel

- `electron_common.f90`：网格、常量和基础积分；
- `electron_injection_profiles.f90`：nonthermal/thermal injection；
- `electron_cooling_kernel.f90`：cooling 组装；
- `electron_cooling_ssa_kernel.f90`：同步与 SSA cooling；
- `electron_cooling_ic_kernel.f90`：numeric IC/KN；
- `electron_cooling_y_kernel.f90`：Y 近似；
- `electron_radiation_kernel.f90`：synchrotron、SSA 和 seed；
- `electron_transport_common.f90`：1D transport helper；
- `electron_transport_2d_kernel.f90`：q 输运；
- `electron_dg_transport.f90`：DG 算子。

公开 wrapper 位于对应 `electron_forward_*.f90`，构建闭包以
`build_extensions.py` 为准。

## 9. 选择原则

| 需求 | 首选 |
|---|---|
| 默认生产、强冷却、强子链 | `fullhide_1d` |
| characteristic 互证 | `charint_1d` |
| conservative semi-Lagrangian | `slc1_1d` |
| 高阶谱形研究 | `weno5_1d` / `dg_1d` |
| 下游空间分辨 | `fullhide_2d` |
| 2D characteristic 互证 | `charint_2d` |

求解器选择表达算法，不表达“失败时换一个”。禁止自动降级。

## 10. 验收

每个电子改动至少检查：

1. 禁用新项时回到原路径；
2. 注入粒子数和能量预算；
3. 无冷却、纯绝热和强同步极限；
4. 谱和 light curve 随 R 连续；
5. 受影响 extension 强制构建及 line-truncation；
6. 真实 `Model` 入口；
7. 性能改动的三次 median 与峰值临时内存。

末位浮点差异可以由更少临时量或等价求值顺序产生；物理峰位、积分预算和连续性
不能变坏。

## 11. 明确保留的边界

f2py wrapper、二维物理网格生成、frequency callback 和线程执行回调不是可随意
折叠的薄层。精简前先确认它们是否承载 ABI、单位、坐标或并发语义。

这些边界应由调用搜索和真实运行共同证明。
