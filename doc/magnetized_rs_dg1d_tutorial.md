# 磁化反向激波与 DG1D

本页说明主反向激波（RS）的横向场 MHD jump，以及 crossing 后电子谱的 DG1D 输运。二者共享 RS 状态，但不是同一个算法开关。

## 1. 使用范围

磁化 RS 描述冷 ejecta 上游进入主 RS 的局域一维 jump。DG1D 是电子能量空间的高阶输运后端，可用于 RS crossing 前后的注入与冷却。

次级密度跳变 RS 当前仍是流体诊断，不继承主 RS 的有序场压缩、磁焓反馈或 DG1D 状态。

## 2. 上游磁化

横向场磁化定义为

\[
\sigma=\frac{B_u'^2}{4\pi\rho_u'c^2}.
\]

这里 \(\rho_u'\) 是冷 ejecta 共动静质量密度。\(\sigma\) 表示上游有序 Poynting 能与物质静质量能之比，不等于辐射微物理参数 \(\epsilon_B\)。

RS 与 ejecta 的相对洛伦兹因子为

\[
\gamma_{34}=\Gamma_3\Gamma_4(1-\beta_3\beta_4).
\]

区域 3 位于 RS 下游和接触面内侧，区域 4 是未受激 ejecta。

## 3. MHD jump

激波静止系的质量、能量、动量和磁通守恒构成 jump：

\[
\rho_u'u_{u,s}=\rho_d'u_{d,s},
\]

\[
T^{0x}_u=T^{0x}_d,
\qquad T^{xx}_u=T^{xx}_d,
\]

\[
\frac{B_u'}{\rho_u'}=\frac{B_d'}{\rho_d'}
\]

用于一维横向场理想 MHD 情形。求解量是激波系下游四速度 \(u_{d,s}\)，随后得到压缩比、热能与下游场。

代码使用

\[
\hat\gamma=\frac43+\frac{1}{3\gamma_{34}}
\]

跨接相对论与弱相对论状态。对 \(\sigma\le10^{-6}\) 的速度求解，流体支给出

\[
\frac{u_{u,s}}{u_{d,s}}=4\gamma_{34}
\]

在该有效绝热指数下是精确结果，不应改写为超相对论近似 \(4\gamma_{34}+3\)。

只有 \(\sigma\le0\) 的热能分支严格返回 \(\gamma_{34}-1\)。在 \(0<\sigma\le10^{-6}\) 内，速度采用流体支，但热化仍保留磁焓表达式，并连续趋向零磁化极限。

## 4. 接触面闭合

接触面洛伦兹因子由 FS 与 RS 下游总压力平衡确定：

\[
p'_{2,\rm tot}=p'_{3,\rm tot}.
\]

RS 侧总压力包括热压与压缩有序场的磁压。增大 \(\sigma\) 通常降低可转化为随机热能的份额，但具体趋势需由 jump 解而不是经验 suppression 因子决定。

物理解必须满足：

- \(\Gamma_3\) 位于 ejecta 与 shocked external medium 的允许区间；
- 质量通量和压力残差收敛；
- 压缩比大于 1；
- \(\sigma\to0\) 连续回到流体 RS；
- crossing 时状态不跳变。

## 5. Crossing

RS crossing 前，未受激 ejecta 持续流入 shock，电子注入率由质量通量和耗散热能决定。穿越完成后：

\[
Q_e(\gamma,t'>t'_{\times})=0.
\]

已有电子仍经历同步、IC 和绝热冷却。数值状态必须把接受的流体状态与其积分时间配对；不得用输出网格时间覆盖真正的状态时间。

检查量包括半径、RS 洛伦兹因子、swept ejecta mass 和 crossing fraction。它们应有限、连续，扫掠质量应单调且不超过 ejecta 总质量。

## 6. 电子输运方程

DG1D 求解

\[
\frac{\partial N}{\partial t'}+
\frac{\partial}{\partial\gamma}[\dot\gamma N]
=Q-\frac{N}{t'_{\rm esc}}.
\]

能量区间划分为单元 \(I_j\)，每个单元内以局域多项式表示 \(N_h\)。弱形式为

\[
\int_{I_j}\!\phi_k\frac{\partial N_h}{\partial t'}d\gamma
-\int_{I_j}\!\frac{d\phi_k}{d\gamma}
(\dot\gamma N_h)d\gamma
+\left[\phi_k\widehat{\dot\gamma N_h}\right]_{j-1/2}^{j+1/2}
=\int_{I_j}\!\phi_k Q\,d\gamma.
\]

数值通量按 \(\dot\gamma\) 的符号选择上风状态。辐射冷却通常有 \(\dot\gamma<0\)，信息从高能流向低能。

## 7. DG1D 的阶数

当前 DG1D 使用单元内线性表示，即空间二阶精度的 DG 基础形式。整体收敛阶还受时间积分、非均匀网格、边界和 limiter 影响，不能只凭多项式次数宣称全局阶数。

模态系数分别编码单元平均值和斜率。辐射计算需要在正交点或目标能量网格重构物理解，而不是把所有模态都当作点值。

## 8. 正性与 limiter

电子数密度必须非负。高阶重构在陡峭 cutoff 附近可能产生局域负值，因此仅对被识别为问题的单元缩放高阶模态，保留单元平均值：

\[
N_h^{\rm lim}=\bar N+\theta(N_h-\bar N),
\qquad 0\le\theta\le1.
\]

该操作不是对最终谱逐点 clamp。若单元平均值本身为负，说明时间步或守恒更新已失败，不能由 limiter 掩盖。

limiter 应满足：

- 保持单元粒子数；
- 光滑正值区 \(\theta=1\)；
- 只缩放导致负重构的高阶部分；
- 不平滑真实谱折点或 cutoff。

## 9. 边界与 crossing 后演化

高能边界禁止非物理流入；低能端允许冷却通量离开活动网格，或按当前边界条件累计。escape 是独立物理 sink，不能用边界通量替代。

crossing 后设置 \(Q=0\)，其余输运算子不变。若谱在 crossing 时突跳，优先检查注入积分区间、状态时间和边界通量，而不是对输出平滑。

## 10. 最小配置思路

公开配置应分别选择：

1. reverse shock 是否启用及 ejecta 磁化；
2. electron solver 是否使用 DG1D；
3. 同步/SSC 与 observer 分量。

不要把“磁化 RS”编码成 DG solver 的隐式副作用。相同 jump 状态应能供不同电子 solver 使用，DG1D 也应能在明确支持的非磁化状态上运行。

## 11. 验证

真实入口至少运行：

```bash
TMPDIR=/tmp uv run python build_extensions.py --module Dynamics_reverse --force
```

随后用本页配置运行一次公开 `Model.flux_density_grid`，检查 RS crossing 前后状态与光变；仓库不提供单独的测试文件入口。

电子模块改动还需构建对应 DG1D 扩展，并对完整 source closure 执行 `-Wline-truncation`。

物理验收包括：

1. \(\sigma=0\) 与流体 jump 一致；
2. 小正 \(\sigma\) 的热能连续；
3. 接触面两侧总压力相等；
4. swept ejecta mass 单调且有界；
5. crossing 后注入严格为零；
6. 无源冷却时粒子数变化只来自边界与 escape；
7. 光滑谱区保持预期收敛，cutoff 附近无负重构；
8. 光变与谱随时间、\(\sigma\) 和数值分辨率连续。

## 12. 不支持内容

- 三维 oblique MHD shock；
- ejecta 内磁重联或场拓扑演化；
- 次级密度跳变 RS 的磁化 jump；
- 用 limiter 修补负单元平均值；
- 用经验 suppression 或平滑替代守恒 jump 与输运。

总物理定义见[物理模型](physics_model.md#19)，数值后端见[电子求解器算法](electron_solver_algorithms.md)。
