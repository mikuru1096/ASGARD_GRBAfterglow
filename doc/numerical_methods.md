# 数值方法索引

本文只回答“使用什么离散、实现在哪里、如何验收”。完整公式见
[`project_algorithm_design.md`](project_algorithm_design.md)，电子细节见
[`electron_solver_algorithms.md`](electron_solver_algorithms.md)。

## 1. 坐标与状态量

主推进坐标是激波半径 (R)。动力学在输出网格给出 (R,Gamma,t_{\rm obs})
和扫掠质量；电子、光子和强子状态在同一 (R) 网格推进。局域 proper-time
步长来自

\[
\Delta t'_i=\frac{\Delta R_i}{2c}
\left[(\Gamma_{i-1}^2-1)^{-1/2}+(\Gamma_i^2-1)^{-1/2}\right].
\]

能量方向使用等距对数网格。粒子数组的 Jacobian 由各 kernel contract 决定，
不能把 `dN/dgamma`、每对数格粒子数和 source rate 互换。

## 2. 动力学

FS 和 RS 都使用显式 RK 推进。密度跃变、注能结束和 RS crossing 是事件边界：
步长在事件处切分，而不是跨步后修改结果。RS 有限强度磁化 jump 位于
`reverse_shock_mhd_jump.f90`；状态和 phase 下标位于 `reverse_shock_state.f90`。

验收关注：半径和时间单调、事件两侧状态连续、`sigma -> 0` 回到流体极限，
以及总能量变化只来自明确源项。

## 3. 电子输运

1D 电子方程是能量空间守恒律：

\[
\frac{\partial N_e}{\partial R}
+\frac{\partial(\dot\gamma_R N_e)}{\partial\gamma}=Q_R.
\]

默认 `fullhide_1d` 使用隐式有限体积推进。charint、SLC、T2G、WENO 和 DG
是互证或研究路径，不改变注入、冷却和辐射的共享定义。2D 路径在有限 q-mass
shell 上增加空间输运；`chi_grid` 只是 BM 等效诊断坐标。

自适应子步由冷却跨格尺度控制。它不能跨越动力学事件，也不能改变输出半径点。
热电子和 solver 组合限制见 backend limits。

## 4. 冷却与辐射

冷却由绝热、同步、SSA、IC/KN 和可选 Y 近似组装。同步 emissivity、SSA
transfer 和 seed photon density 共享同一电子状态及磁场。SSC 使用 Jones/KN
响应；electron loss 与 photon emissivity 必须在同一网格和单位契约下闭合。

光子 luminosity 转 local density 需要体积、escape time 和频率 Jacobian。该转换
只能发生一次，observer luminosity 不能直接作为局域 target field。

## 5. 强子

formal 1D kernel 按壳层推进质子，再依序计算 proton synchrotron、pγ、BH、pp、
次级 species、IC 和 pair。pγ 使用 Hummer response。pp 的 `delta` 为默认；
AM3-derived detailed gamma 模型只替换 π0 gamma 谱，不改变 proton loss、neutrino
和 secondary pair。

每个过程必须声明输入是 density、rate 还是 luminosity。能量坐标转换、mbarn 到
cm² 和 photon survival 各自只允许一个所有者。

## 6. Joint feedback

`separated` 在电子阶段完成后计算强子；`joint` 在同一 R 网格上交替更新 electron、
photon field、hadron 和 secondary source。后者目前只支持严格组合，且仍有已登记的
pair state/source 所有权缺陷，见根目录 `BUG.md`。

算法入口和边界见
[`joint_secondary_feedback_algorithm.md`](joint_secondary_feedback_algorithm.md)。

## 7. Pair cascade

当前 cascade 闭合 γγ pair injection 与 pair synchrotron。它不是包含 IC 再处理的
完整电磁级联。单次和 shell-sequence 路径都必须满足无吸收极限、薄/厚极限和
pair 能量预算。

## 8. Observer 投影

shell-level EATS 对本征光度施加 Doppler、红移、距离和到达时延。adaptive-theta
改变角积分精度；`chi_eats_2d` 只替换 FS synchrotron+SSA 投影，SSC、强子和
secondary feedback 保持 shell-level。

结构化喷流复用每个 ring/patch 的局域求解结果。只改变观测角时，应复用 solve
state，仅重算 projection。

## 9. 数值缓存

缓存键必须覆盖所有改变局域状态的公开参数。只缓存纯函数结果或只读表；线程执行
回调不能共享可写 scratch。缓存优化必须比较结果、命中语义和峰值内存。

## 10. 验收顺序

修改后按风险递增检查：

1. Python compile/encoding 和 Fortran syntax；
2. 受影响 f2py extension 强制构建；
3. 干净 `/tmp` source closure 的 `-Wline-truncation`；
4. 最短真实入口；
5. 禁用过程极限、能量预算和连续性；
6. 至少三次 median 的性能比较。

数值结果允许末位舍入差异，但性能收益必须可测，重要物理曲线必须连续。不得用
clamp、smoothing、经验归一化或 fallback 掩盖错误。

## 11. 实现索引

| 阶段 | Python 编排 | Fortran 主入口 |
|---|---|---|
| Dynamics | `asgard_state.py` | `Dynamics_forward.f90`, `reverse_shock.f90` |
| Electron | `asgard_runtime.py` | `electron_forward_*.f90` |
| Radiation | `asgard_state.py` | `electron_radiation_kernel.f90`, `src/Radiation` |
| Hadronic | `hadronic_am3_solver.py` | `hadronic_forward_1d.f90` |
| Observer | `asgard_postprocess.py` | `SED_interpolation.f90` |

源闭包和 ABI 导出以 `build_extensions.py` 为准；不要维护复制源码符号的静态索引。

专题页只在确有独立算法契约时存在；一般实现细节直接从本表进入源码，避免文档与
程序单元重复漂移。
