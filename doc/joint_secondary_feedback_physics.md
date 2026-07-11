# 含时二级反馈的物理契约

本页定义 `coupling="joint"` 的物理所有权。它不是在 `separated` 结果上追加更多分量，而是在同一壳层演化中闭合质子、电子/正电子和光子场。

## 1. 演化变量

联合状态为

\[
N_p(E_p,R),\qquad N_e(\gamma,R),\qquad n_\gamma(\epsilon,R).
\]

代码沿动力学轨迹使用半径 \(R\)；局域共动时间增量由

\[
dt'=\frac{dR}{\Gamma\beta c}
\]

给出。这里没有把所有壳层宽度泛化为 \(R/(12\Gamma)\)。体积、路径长度和逃逸时间来自当前壳层几何状态。

粒子数组是每个能量单元的谱密度；光子数组是局域共动光子谱。所有源项进入方程前必须与该归一化一致。

## 2. 质子方程

\[
\frac{\partial N_p}{\partial t'}+
\frac{\partial}{\partial E_p}(\dot E_pN_p)
=Q_p-\frac{N_p}{t'_{p,\rm esc}}.
\]

\(\dot E_p\) 可包含绝热、proton synchrotron、pγ、BH 和 pp 损失。每个过程的连续损失与对应产物共用一次能量转移率。

pγ 与 pp 可产生 gamma、pair 和 neutrino；BH 产生 pair。neutrino 离开联合状态，但计入能量预算。

## 3. 电子/正电子方程

\[
\frac{\partial N_e}{\partial t'}+
\frac{\partial}{\partial\gamma}(\dot\gamma N_e)
=Q_{e,\rm pri}+Q_{e,\gamma\gamma}+Q_{e,\rm BH}
+Q_{e,p\gamma}+Q_{e,pp}-\frac{N_e}{t'_{e,\rm esc}}.
\]

主电子和二级电子在进入辐射核后服从相同的同步、IC 与绝热损失。pair 注入属于联合 solver；observer 不得再次从同一 gamma-gamma 吸收量生成 pair。

正负电子在不需要电荷符号的辐射核中以总 pair 谱演化。任何按单粒子种类定义的响应表必须在注入前处理其 multiplicity。

## 4. 光子方程

\[
\frac{\partial n_\gamma}{\partial t'}=
Q_{\rm syn}+Q_{\rm IC}+Q_{p,\rm syn}+Q_{p\gamma}
+Q_{pp}+Q_{\rm ann}-
\frac{n_\gamma}{t'_{\rm esc}}-
\dot n_{\gamma\gamma}-\dot n_{\rm abs}.
\]

`Q_ann` 表示当前实现允许的 annihilation 光子源。SSA 或其他吸收项只有在联合状态显式回收其能量时才构成反馈；否则它们是传输层衰减，不得虚构热化源。

光子逃逸时间由壳层路径长度确定。逸出谱进入 observer，留在壳层内的谱继续作为 SSC、pγ、BH 和 gamma-gamma 的目标场。

## 5. 交换预算

每个内部交换必须成对出现：

| 损失端 | 产物端 |
|---|---|
| 电子同步/IC | 同步/IC 光子 |
| proton synchrotron | proton-synchrotron 光子 |
| gamma-gamma 光子 sink | 电子–正电子注入 |
| BH 质子损失 | BH pair 注入 |
| pγ 质子损失 | gamma、pair、neutrino |
| pp 质子损失 | gamma、pair、neutrino |

离散能量预算检查采用各网格上的能量加权积分。网格边界逸出和 neutrino 输出是允许的系统 sink；内部交换不能同时算作 sink 和净损失。

## 6. `joint` 与 `separated`

`separated`：

1. 先求主电子和局域辐射；
2. 再以冻结目标场计算强子与 pair 分量；
3. 二级产物不回写此前的目标场。

`joint`：

1. 在每个动力学步共同更新三类状态；
2. 二级 pair 会改变后续同步、IC 和目标光子场；
3. 强子反应率随更新后的光子场演化。

两种结果不能逐分量相加。联合模式输出已经包含其拥有的二级反馈。

## 7. 当前支持边界

- 联合输运是壳层平均的一维能量问题。
- 它不是 `fullhide_2d` q 单元或 chi-resolved 强子输运。
- `chi_eats_2d` 不投影联合 SSC、强子或 pair 体积分布。
- reverse hadronic 复用壳层目标场与投影契约，不宣称完整的 RS 空间级联。
- pp 默认 gamma 模型为 `delta`；AM3 详细模型只改变 π0 gamma 源。
- 不用经验 clamp、平滑或事后归一化修复能量闭合。

## 8. 验收

联合路径至少检查：

1. 关闭二级过程时退化到对应主输运；
2. 每个源与 sink 有相同单位和网格测度；
3. gamma-gamma 吸收能与 pair 注入能相符；
4. BH、pγ、pp 产物不超过质子损失预算；
5. `separated` 与 `joint` 不重复组装 pair；
6. 谱有限、非负，阈值以下源为零；
7. 随半径和时间的物理量连续平滑；
8. observer 只投影逸出分量，不改变局域联合状态。

算法离散与实际入口见 [联合反馈算法](joint_secondary_feedback_algorithm.md)。
