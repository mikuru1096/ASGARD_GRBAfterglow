# 数值方法

本文档说明 ASGARD 当前工作树的数值方法、实现方式和验证重点。它不是代码清单，而是解释为什么这些算法适合 GRB 余辉的多尺度问题，以及每个方法在什么边界内可信。全项目算法设计见 `doc/project_algorithm_design.md`，逐步流程见 `doc/algorithm_workflow.md`。

## 1. 方法总览

ASGARD 的数值主线可以概括为：

```text
R 坐标激波动力学
-> 对数能量网格上的守恒粒子输运
-> 局域辐射和光子场构建
-> 强子、对产生和二级粒子算子
-> 等到达时间面观测者投影
```

核心数值设计包括：

- 用激波半径 \(R\) 作为粒子、光子和强子反馈的共同推进坐标。
- 用对数洛伦兹因子网格上的守恒变量推进电子和质子谱，显式保留雅可比因子。
- 用隐式迎风/特征线/WENO 等多种电子输运器交叉验证谱演化。
- 用壳层级 `PhotonFieldState` 把 SSC、p-gamma、BH、pp、\(\gamma\gamma\) 和二级粒子反馈组织成可追踪的源-汇系统。
- 用等到达时间面与 \(\chi\) 分辨有限厚壳层投影区分本地输运和观测几何。

这套设计避免把观测光变当成唯一状态量。真正需要验收的是

\[
\Gamma(R),\quad
B'(R),\quad
N_e(\gamma,R),\quad
N_p(\gamma,R),\quad
n_\gamma(\epsilon,R),\quad
F_\nu(t_{\rm obs}).
\]

只有这些内部量和最终观测量同时平滑、连续且预算闭合，结果才可信。

## 2. 统一的 \(R\) 坐标

GRB 余辉中不同过程天然使用不同时间：

- 动力学常用实验室时间或观测者时间。
- 粒子冷却常用共动时间 \(t'\)。
- 光变输出使用 \(t_{\rm obs}\)。

ASGARD 在输运层统一使用激波半径 \(R\)。任意共动率 \(\lambda'\) 进入 \(R\) 坐标前都要做

\[
\lambda_R
=
\lambda'
\frac{\mathrm{d}t'}{\mathrm{d}R}
=
\frac{\lambda'}{\beta\Gamma c}.
\]

壳层步长上的共动时间为

\[
\Delta t'_i
=
\frac{\Delta R_i}{\beta_i\Gamma_i c}.
\]

这使电子冷却、强子相互作用、二级粒子注入和光子汇可以在同一壳层上比较预算。禁止用观测者时间步长代替 \(\Delta t'\)，因为等到达时间面投影中的 \(t_{\rm obs}\) 已经包含几何延迟。

## 3. 动力学积分

正向激波动力学入口是 `src/Dynamics/Dynamics_forward.f90`。主变量可抽象为

\[
\mathbf{y}
=
(\Gamma,M_{\rm sw},U,R),
\qquad
\frac{\mathrm{d}\mathbf{y}}{\mathrm{d}T}
=
\mathbf{f}(T,\mathbf{y};\mathbf{p}).
\]

输出网格由观测者时间相关的 `T` 控制。对 \(T>0\) 的加密，代码在 \(\ln T\) 上推进：

\[
\frac{\mathrm{d}\mathbf{y}}{\mathrm{d}\ln T}
=
T
\frac{\mathrm{d}\mathbf{y}}{\mathrm{d}T}.
\]

RK4 更新为

\[
\mathbf{y}_{n+1}
=
\mathbf{y}_n
+
\frac{h}{6}
(\mathbf{k}_1+2\mathbf{k}_2+2\mathbf{k}_3+\mathbf{k}_4).
\]

因此在右端项平滑且没有事件分裂的动力学区间，动力学 ODE 的全局阶数为四阶：

\[
\mathbf{y}(T)-\mathbf{y}_{\Delta T}(T)
=
{\cal O}(\Delta T^4).
\]

反向激波 crossing、waiting-to-shock 和 pressure-balance 触发是物理分支切换；代码用事件端点捕获保持状态连续。事件附近的验收不是对单一 RK4 步做形式阶数外推，而是检查分支变量、热状态和辐射输出在物理事件两侧连续。

ASGARD 的动力学实现不追求单一解析标度，而是在同一右端项框架中处理均匀星际介质、\(k=2\) 恒星风介质、密度跳变和能量注入。解析 BM 标度只作为验收直觉：

\[
\Gamma_{\rm ISM}(R)\propto R^{-3/2},
\qquad
\Gamma_{\rm wind}(R)\propto R^{-1/2}.
\]

如果无注入、无密度跳变的基线偏离这些趋势，应先检查动力学阶段，而不是在电子或投影层修正光变。

## 4. 反向激波动力学的分支处理

反向激波 f2py 入口 `Dynamics_reverse.dynamics_reverse` 位于 `src/Dynamics/reverse_shock.f90`。first RS 和 density-jump multiple RS 使用同一个状态推进循环：每个输出时刻先把 first 分支推进到目标时间，写主 RS 输出，再用相邻两个状态扫描 jump window 并写 multiple 分支诊断。穿越前后仍作为不同物理分支处理，避免一个 RK 步跨过不连续方程组。

穿越前使用

\[
\xi
=
\ln m_{3,{\rm frac}},
\qquad
m_{3,{\rm frac}}
=
\frac{M_3}{M_{\rm ej}},
\]

作为推进变量之一。穿越端点为

\[
m_{3,{\rm frac}}=1.
\]

穿越后固定已激波加热抛射物质量 \(M_3\)，继续推进 \(U_3\)、\(V_3\) 和整体运动状态。磁化反向激波使用

\[
\sigma
=
\frac{B_4^2}{4\pi \Gamma_4^2\rho_4c^2}
\]

控制上游磁化，并输出湍流磁场和有序磁场合成的总磁场：

\[
B_3
=
\sqrt{B_{\rm turb}^2+B_{\rm ord}^2}.
\]

该方法的关键验收是

\[
\sigma\rightarrow0
\quad\Rightarrow\quad
{\rm magnetized\ reverse\ shock}
\rightarrow
{\rm unmagnetized\ reverse\ shock}.
\]

这比只看某条光变是否相近更严格，因为它同时约束 ejecta baryonic mass、pressure-balance 触发、MHD jump、区域 3 热状态、有序磁场和磁压焓惯性。

当前磁化链条的代码映射为：

| 物理量 | 公式角色 | 实现入口 |
| --- | --- | --- |
| \(M_{\rm ej,b}=E_{\rm iso}/[(1+\sigma)\eta_0c^2]\) | 有限磁化下的 baryonic ejecta mass | `Dynamics_reverse`, `structured_jet_1d` 主流程内显式分支 |
| \(u_{3s}, C, \epsilon_{\rm th,3}, \mathcal{S}_{\rm fast}\) | 下游四速度、压缩比、热比内能和快模判据 | `rs_mhd_state` |
| \(B_{4,{\rm ord}}=\sqrt{4\pi c^2\sigma\rho_4}\) | 上游有序磁场 | `reverse_rhs` 主流程内显式分支 |
| \(M_{B,{\rm eff}}\) | 有序场磁压焓对 bulk 方程的惯性贡献 | `compute_ordered_magnetic_inertia` |

### 4.1 穿越捕获而不是跨步修补

反向激波穿越是真实物理分支切换，不是数值异常。ASGARD 不让 RK 步直接跨过穿越点后再用后处理修正，而是在穿越前阶段捕获

\[
T_{\rm cross},\quad
R_{\rm cross},\quad
M_{3,{\rm cross}},\quad
U_{3,{\rm cross}},\quad
V_{3,{\rm cross}}.
\]

算法上，当试探步会让

\[
m_{3,{\rm frac}}
<1
\]

跨到穿越之后时，代码用专门的 \(M_3\) 穿越前 RK 端点处理，把穿越处的状态固定下来，再进入穿越后方程。这样做的数值意义是：

- 穿越前有持续的新已激波加热抛射物注入。
- 穿越后 \(M_3\) 固定，不再有新的反向激波注入质量。
- \(U_3\) 和 \(V_3\) 的演化从激波加热加膨胀变为穿越后的热演化。

若把这两个阶段混在一个右端项里，\(\gamma_{34}\)、\(M_3\)、\(U_3/V_3\) 和 \(B_3\) 会在穿越附近出现非物理尖峰。

### 4.2 `gamma34` 注入能标

反向激波电子注入不是用正向激波的 \(\Gamma-1\) 近似，而是使用激波前沿相对洛伦兹因子：

\[
\gamma_{34}
=
\Gamma_3\Gamma_4
\left(
1-\beta_3\beta_4
\right).
\]

非热电子最小洛伦兹因子的量级关系是

\[
\gamma_{m,3}
\simeq
1+
\frac{\epsilon_{e,3}}{\xi_{N,3}}
\frac{p_3-2}{p_3-1}
\frac{m_p}{m_e}
(\gamma_{34}-1).
\]

这保证反向激波注入能标跟抛射物/区域 3 的局域跳跃条件绑定，而不是跟正向激波下游状态混用。`Gamma34_inst` 是诊断字段，应用来检查穿越前后注入能标是否平滑。

### 4.3 显式 \(U_3/V_3\) 热状态

反向激波区域 3 的热状态用总内能 \(U_3\) 和共动体积 \(V_3\) 表示：

\[
e_3
=
\frac{U_3}{V_3}.
\]

穿越前，\(U_3\) 的变化包含激波加热和绝热项：

\[
\mathrm{d}U_3
=
\mathrm{d}U_{3,{\rm shock}}
+
\mathrm{d}U_{3,{\rm ad}}.
\]

激波加热的量级是

\[
\mathrm{d}U_{3,{\rm shock}}
\sim
(\gamma_{34}-1)c^2\,\mathrm{d}M_3.
\]

体积变化包含新的已激波加热物质加入和膨胀：

\[
\mathrm{d}V_3
=
\mathrm{d}V_{3,{\rm shock}}
+
\mathrm{d}V_{3,{\rm exp}}.
\]

穿越后

\[
\mathrm{d}M_3=0,
\qquad
\mathrm{d}U_{3,{\rm shock}}=0,
\qquad
\mathrm{d}V_{3,{\rm shock}}=0,
\]

只保留膨胀和热演化。这个显式热状态使磁场、电子注入和穿越后冷却都来自同一组 \(U_3,V_3,M_3\)，而不是由光变层经验参数控制。

### 4.4 湍流磁场和有序磁场

反向激波总磁场分成湍流分量和有序分量：

\[
B_{3,{\rm tot}}
=
\sqrt{
B_{3,{\rm turb}}^2
+
B_{3,{\rm ord}}^2
}.
\]

湍流分量来自区域 3 热能密度：

\[
B_{3,{\rm turb}}
=
\sqrt{8\pi\epsilon_{B,3}e_3}.
\]

有序分量来自上游磁化 \(\sigma\)。穿越前，它由上游抛射物磁场和磁流体压缩给出；穿越后不再有新的上游物质注入，有序磁场随壳层膨胀演化。代码在穿越处保存

\[
B_{3,{\rm ord,cross}},
\qquad
R_{\rm cross},
\qquad
V_{3,{\rm cross}},
\]

并在穿越后按体积和半径变化推进有序分量。数值验收应检查

\[
B_{3,{\rm ord}}(R)
=
B_{3,{\rm ord,cross}}
\frac{V_{3,{\rm cross}}}{V_3(R)}
\frac{R}{R_{\rm cross}},
\]

在穿越附近连续，且

\[
B_{3,{\rm ord}}\rightarrow0
\quad{\rm when}\quad
\sigma\rightarrow0.
\]

### 4.5 区域 2 与区域 3 的耦合

反向激波不是独立光源；它与正向激波共享爆波动力学背景。数值上需要同时追踪：

\[
M_2,\quad M_3,\quad \Gamma,\quad R,\quad U_3,\quad V_3.
\]

其中区域 2 影响整体运动演化，区域 3 决定反向激波电子注入、反向激波同步辐射、反向激波 SSC 和跨区逆康普顿种子光子场。跨区逆康普顿需要两个区域的光子场在同一观测者设置中组装：

\[
Q_{\rm IC}^{3\leftarrow2}
=
Q_{\rm IC}
\left[
N_{e,3},n_{\gamma,2}
\right],
\qquad
Q_{\rm IC}^{2\leftarrow3}
=
Q_{\rm IC}
\left[
N_{e,2},n_{\gamma,3}
\right].
\]

这要求正向/反向激波种子光子场的时间、半径和几何归一化一致。不能把反向激波光变作为一个后处理凸起叠加到正向激波光变上。

### 4.6 反向激波数值验收

反向激波路径需要比普通光变更细的验收：

| 检查项 | 目标 |
| --- | --- |
| \(m_{3,{\rm frac}}\) | 穿越前单调接近 1，穿越后固定为完整已激波加热抛射物 |
| \(\gamma_{34}(R)\) | 注入阶段平滑，穿越附近无孤立尖峰 |
| \(U_3/V_3\) | 穿越前后连续，穿越后热演化平滑 |
| \(B_{3,{\rm turb}}\) | 跟随 \(e_3\) 平滑变化 |
| \(B_{3,{\rm ord}}\) | \(\sigma>0\) 时穿越后连续衰减，\(\sigma\to0\) 时消失 |
| 反向激波同步辐射 | 峰时与穿越/冷却时间尺度一致 |
| 反向激波 SSC / 跨区逆康普顿 | 种子光子场和电子状态对齐 |
| \(\sigma\to0\) | 回到非磁化基线 |

这些检查失败时，应查反向激波动力学、跳跃条件、区域 3 热状态、磁场压缩或种子光子场构建。禁止在反向激波光变层加 smoothing 或经验时间因子。

## 5. 对数洛伦兹因子守恒变量

电子和质子谱都在对数能量网格上推进。令

\[
x
=
\log_{10}\gamma,
\qquad
\Delta x
=
x_{i+1}-x_i.
\]

若物理谱为

\[
N_\gamma
=
\frac{\mathrm{d}N}{\mathrm{d}\gamma},
\]

则数值守恒变量是

\[
N_x
=
\frac{\mathrm{d}N}{\mathrm{d}x}
=
\gamma\ln10\,N_\gamma.
\]

总粒子数计算为

\[
N_{\rm tot}
\simeq
\sum_i N_{x,i}\Delta x.
\]

这种变量选择让跨很多数量级的电子/质子谱仍能做守恒更新，并避免在低能端或高能截断附近丢失雅可比因子。

## 6. `fullhide_1d` 隐式电子输运

默认公开接口基线是 `fullhide_1d`，入口为 `src/Electron/electron_forward_fullhide_1d.f90`。它求解的 \(R\)-坐标输运方程为

\[
\frac{\partial N_x}{\partial R}
+
\frac{\partial}{\partial x}
(v_xN_x)
=
Q_x.
\]

其中

\[
v_x
=
\frac{1}{\gamma\ln10}
\frac{\mathrm{d}\gamma}{\mathrm{d}R}.
\]

有限体积形式为

\[
N_{x,i}^{n+1}
=
N_{x,i}^{n}
-
\frac{\Delta R}{\Delta x}
\left(
F_{i+1/2}^{n+1}
-
F_{i-1/2}^{n+1}
\right)
+
\Delta R\,Q_{x,i}.
\]

冷却通常让 \(v_x<0\)，即谱向低能方向平流。隐式迎风通量可概念化为

\[
F_{i+1/2}^{n+1}
=
\begin{cases}
v_{i+1/2}N_{x,i}^{n+1}, & v_{i+1/2}>0,\\
v_{i+1/2}N_{x,i+1}^{n+1}, & v_{i+1/2}<0.
\end{cases}
\]

于是每个 step 需要解线性系统

\[
\mathbf{A}\mathbf{N}^{n+1}
=
\mathbf{b}.
\]

`fullhide_1d` 的优势是在刚性冷却下稳定，适合 MCMC 和正式光变生成。它的代价是一阶隐式格式有数值扩散，因此高能截断和尖锐冷却断点的形状需要用 `charint_1d`、`weno5_1d` 或网格加密做交叉检查。

## 7. 子步压缩与自适应子步

`fullhide_1d` 中一个壳层间隔

\[
\Delta R_i
=
R_i-R_{i-1}
\]

会拆成多个 substeps。固定子步模式中

\[
\delta R
=
\frac{\Delta R_i}{L_i}.
\]

均匀介质、非热电子的快速路径会把一串固定子步的注入和界面耦合压缩到一个时空序列中。这是性能优化；它不是改变物理方程，而是在相同源项/冷却结构下减少重复构造三对角系统。

自适应子步比较一次 full step 与两次 half step：

\[
\epsilon
=
\frac{
\|\mathbf{N}_{1/2+1/2}-\mathbf{N}_{1}\|
}{
\|\mathbf{N}_{1/2+1/2}\|
}.
\]

若 \(\epsilon\) 超过 `electron_substep_rtol`，减小步长；若小得多，允许后续步长变大。该机制的验收目标是谱和特征频率平滑，而不是让所有壳层使用同样数量的子步。

## 8. 注入归一化

每个壳层的非热电子注入满足数目预算

\[
\sum_i Q_{x,i}\Delta x
=
\xi_N\Delta N_{\rm sw}.
\]

能量预算为

\[
\sum_i
(\gamma_i-1)m_ec^2
Q_{x,i}\Delta x
=
\epsilon_e\Delta E_{\rm sh}.
\]

截断幂律源项为

\[
Q_x(\gamma)
\propto
\gamma\ln10\,
\gamma^{-p}
\exp\!\left(-\frac{\gamma}{\gamma_{\max}}\right).
\]

`electron_gamma_m_exact` 负责处理 \(p\) 接近 2、高能截断不远离 \(\gamma_m\) 等情况。文档中的渐近公式不能替代这个精确归一化。

## 9. 冷却组装

冷却组装门面位于 `src/Electron/electron_cooling_kernel.f90`，SSA、IC 和 Compton-Y 实现分别位于 `electron_cooling_ssa_kernel.f90`、`electron_cooling_ic_kernel.f90` 和 `electron_cooling_y_kernel.f90`。门面把不同物理项组装为界面速度：

\[
v_x
=
v_{\rm ad}
+
v_{\rm syn}
+
v_{\rm IC}
+
v_{\rm SSA}.
\]

同步冷却满足

\[
\frac{\mathrm{d}\gamma}{\mathrm{d}R}
\bigg|_{\rm syn}
=
-
\frac{\sigma_TB'^2}{6\pi m_ec}
\frac{\gamma^2}{\beta\Gamma c}.
\]

逆康普顿冷却依赖种子光子场：

\[
v_{\rm IC}
=
v_{\rm IC}\left(\gamma,n_\gamma,\Gamma,R\right).
\]

ASGARD 把光子场准备和最终冷却组装分开：同一个种子光子场可以同时进入逆康普顿冷却、SSC 光子源、强子靶光子场和吸收靶光子场。这样可以检查能量预算，而不是让每个模块各自构造一份不一致的光子场。

## 10. 多求解器互证

电子求解器不是单一路径：

| 求解器 | 方法 | 作用 |
| --- | --- | --- |
| `fullhide_1d` | 隐式迎风 | 公开接口稳定基线 |
| `charint_1d` | 特征线重映射 | 检查冷却断点和截断 |
| `dg_1d` | 多域 LGL 谱元 DG | FS/RS opt-in 高阶 1D 路径 |
| `slc1_1d` | 半拉格朗日 | 方法比较 |
| `t2g1_1d` | 旧隐式格式 | 回归比较 |
| `weno5_1d` | WENO5 + RK | 高阶谱形诊断 |
| `fullhide_2d` | \(x,\chi\) 隐式输运 | 厚壳电子基线 |
| `charint_2d` | \(x\) 特征线 + \(\chi\) 隐式 | 2D 加速混合路径 |

方法互证的判据不是逐点完全一致，而是物理诊断一致：

\[
\nu_m(R),\quad
\nu_c(R),\quad
\nu_a(R),\quad
F_\nu(t),\quad
\int N_e\,\mathrm{d}\gamma.
\]

如果只有一个求解器出现非光滑谱形，优先检查该求解器的步长、重映射或边界处理。

## 11. 2D \(\chi\) 电子输运

2D 路径在能量坐标 \(x=\log_{10}\gamma_e\) 外增加有限主动壳层质量坐标 \(q\)，推进

\[
U(x,q,R)
=
\frac{\mathrm{d}N_e}{\mathrm{d}x\,\mathrm{d}q}.
\]

\(q=0\) 是激波前沿，\(q=q_{\rm active}=1-(3/4)^4\) 是当前主动下游壳层外边界。`chi_grid` 是 \(q\) cell 的 BM 等效 \(\chi_{\rm BM}=(1-q)^{-(4-k)/(3-k)}\) 诊断值；真正进入投影几何的是 `chi_radius_cm`、`chi_gamma_bulk` 和 `chi_dvolume_weight`。

完整连续方程

\[
\partial_R U+\partial_x(A_xU)+\partial_q(A_qU)
=
\partial_q(D_q\partial_q U)+S
\]

以及 \(A_q\)、有限 \(q\)-shell 下游速度、隐式三对角系统、characteristic remap 和 CFL 子步公式见 `doc/shock_shell_adaptive_algorithms.md`。该路径不表示强子、SSC 或对级联也变成 \(\chi\) 分辨；这些过程当前仍是壳层级契约。

## 12. 辐射积分与种子光子场复用

同步辐射积分为

\[
L'_{\nu_j}
\simeq
\sum_i
N_{x,i}
P'_{\nu_j}(\gamma_i,B')
\Delta x.
\]

固定网格同步辐射路径是公开接口快速路径。自适应同步辐射积分只作为诊断，因为正式拟合需要稳定、可复现且低开销的热路。

同步自吸收光深为

\[
\tau_{\nu}
=
\alpha_\nu\ell'.
\]

对于有限厚度网格单元，逃逸因子用

\[
S_{\nu,{\rm cell}}
=
\exp(-\tau_{\rm front})
\frac{1-\exp(-\tau_{\rm cell})}{\tau_{\rm cell}}.
\]

该公式避免只取网格单元前边界造成系统性透射偏差。它是 `chi_eats_2d` 投影相对薄壳投影的核心差异之一。

## 13. PhotonFieldState

`PhotonFieldState` 是 ASGARD 连接轻子和强子过程的数值接口。它不是一张观测通量表，而是壳层局域光子场契约。

关键变换是

\[
L_\nu
\rightarrow
u_\nu
\rightarrow
n_\epsilon,
\qquad
\epsilon=h\nu.
\]

若从每赫兹量转为每能量量，必须使用

\[
n_\epsilon
=
\frac{n_\nu}{h}.
\]

所有使用光子场的算子必须明确单位：

- 逆康普顿冷却使用光子能量密度。
- SSC 源项使用同一个种子光子场的散射发射率。
- p-gamma 和 BH 使用靶光子数密度。
- \(\gamma\gamma\) 使用吸收靶光子密度。

这使不同过程的源项/汇项可以在同一壳层上追踪。

## 14. 强子微物理算子化

正式强子路径把复杂强子过程拆成算子：

```text
质子输运
-> 质子同步辐射
-> p-gamma 响应
-> Bethe-Heitler 对产生
-> pp 源项/损失
-> 二级粒子输运
-> 二级粒子辐射
-> 中微子输出
-> 对产生分支
```

质子输运可写为

\[
\frac{\partial N_p}{\partial R}
+
\frac{\partial}{\partial E_p}
\left(
\dot{E}_{p,R}N_p
\right)
=
Q_{p,R}
+
Q_{p,{\rm reinj},R}
-
\frac{N_p}{R_{\rm esc}}.
\]

每个微物理算子返回的 \({\rm s}^{-1}\) 率要通过

\[
\frac{\mathrm{d}t'}{\mathrm{d}R}
=
\frac{1}{\beta\Gamma c}
\]

转换为每 \(R\) 源项或损失项。这样 BH、pp 和 \(\gamma\gamma\) 的二级粒子源可以进入电子方程，而 p-gamma/BH 的光子损失可以进入光子存活率。

## 15. 联合反馈的数值闭环

`electron_photon_coupling="joint"` 把电子、光子和强子反馈放在同一壳层序列中：

```text
主电子求解
-> 光子场构建
-> 正式强子算子
-> 二级 e± 源项 + 光子汇/源
-> 耦合电子求解
-> 光子场重建
```

核心方程是

\[
Q_{e,{\rm total},R}
=
Q_{e,{\rm shock},R}
+
\frac{1}{\beta\Gamma c}
\left(
Q_{e,{\rm BH}}
+
Q_{e,pp}
+
Q_{e,\gamma\gamma}
\right).
\]

光子存活率可抽象为

\[
n_{\gamma}^{\rm next}
=
n_{\gamma}^{\rm src}
\exp[-\tau_{\gamma\gamma}-\tau_{\rm BH}-\tau_{p\gamma}]
+
n_{\gamma}^{\rm add}.
\]

这里的表达只说明源-汇结构；实际实现使用正式算子输出的损失/存活率，不用经验衰减。验收极限包括弱反馈回到分离基线，以及强反馈时光子存活率、二级粒子源和电子谱都随 \(R\) 平滑。

## 16. 对级联的壳层序列

当前对级联是 \(\gamma\gamma\) 对产生/同步辐射级联：

\[
n_\gamma^{(m)}
\xrightarrow{\gamma\gamma}
Q_{e^\pm}^{(m)}
\xrightarrow{\rm synch}
n_{\gamma,{\rm pair}}^{(m)}
\xrightarrow{\rm add}
n_\gamma^{(m+1)}.
\]

数值上它沿壳层序列更新，而不是把所有壳层混成一个全局稳态。这样可以保留含时激波演化对对源项和对同步辐射种子光子场的影响。

未实现的逆康普顿介导电磁级联需要额外闭环：

\[
Q_{e^\pm}
\xrightarrow{\rm IC}
n_{\gamma,{\rm IC}}
\xrightarrow{\gamma\gamma}
Q_{e^\pm}^{\rm new}.
\]

该过程当前不属于实现边界，不能在文档或基准测试中当作已有功能解释。

## 17. 等到达时间面和多普勒插值

观测投影源文件包括：

- `src/Interpolation/SED_interpolation.f90`
- `src/Interpolation/SED_interpolation_structured.f90`
- `src/Interpolation/interpolation_common.f90`

角向到达时间为

\[
t_{\rm obs}(R,\mu)
=
t_{\rm obs,axis}(R)
+
(1+z)\frac{R(1-\mu)}{c}.
\]

对每个观测时间，算法在相邻壳层间做时间插值，再在对数频率空间做谱插值。实现中定义

\[
D
=
\Gamma(1-\beta\mu)
=
\delta^{-1},
\]

因此

\[
\nu'
=
(1+z)D\nu_{\rm obs},
\qquad
F_{\nu_{\rm obs}}
\propto
D^{-3}L'_{\nu'}.
\]

代码变量 `doppler` 对应的是 \(D\)，不是常见记号 \(\delta\)。这是阅读投影代码时最容易误解的点。

## 18. \(\chi\) 分辨等到达时间面投影

`sed_interpolation_chi` 在普通角向等到达时间面外加入下游厚度维度：

\[
(R_i,\theta_j,\phi_j)
\rightarrow
(R_i,\theta_j,\phi_j,\chi_k).
\]

每个 \(\chi\) cell 使用自己的 \(R_{\chi,k,i}\)、\(\Gamma_{\chi,k,i}\)、\(W_{\chi,k,i}\) 和 optical-depth survival。到达时间、\(D^{-3}\) Doppler 权重、SSA cell-average survival、transport-to-projection remap 和薄壳极限公式集中维护在 `doc/shock_shell_adaptive_algorithms.md`。这个薄壳极限仍是有限厚壳层投影最重要的回归测试。

## 19. 结构化喷流与面元重用

结构化喷流通过角向面元聚合：

\[
F_\nu^{\rm total}
=
\sum_j
F_{\nu,j}\Delta\Omega_j.
\]

每个面元有自己的

\[
E_{\rm iso}(\theta_j,\phi_j),
\qquad
\Gamma_0(\theta_j,\phi_j).
\]

当只扫描观测角且物理求解状态未改变时，可复用底层求解，只重跑投影。这是性能优化，不是物理近似。只有当介质、频率网格、时间网格、求解器和壳层状态完全一致时，这种重用才成立。

## 20. 缓存和运行时间报告

ASGARD 公开接口 `Model` 会缓存相同查询的原始求解。性能报告必须区分：

\[
t_{\rm cold}
=
t_{\rm dyn}
+t_{\rm electron}
+t_{\rm rad}
+t_{\rm had}
+t_{\rm proj},
\]

和

\[
t_{\rm warm}
\approx
t_{\rm cached\ query}.
\]

拟合或基准测试中如果把热缓存查询当作冷启动求解，会错误估计正式采样成本。文档和基准测试记录必须说明是否预热过缓存。

## 21. 编译与行截断检查方法

Fortran 重要改动后的最低验证包括：

```text
affected build_extensions.py --force
gfortran -Wall -Wline-truncation -Werror=line-truncation -fsyntax-only
minimal end-to-end check
```

行截断检查必须从干净 `/tmp` 工作目录执行，并显式指定模块输出目录：

```bash
rtk bash -lc 'source ~/.wsl_env && rm -rf /tmp/asgard_linecheck && mkdir -p /tmp/asgard_linecheck && cd /tmp && gfortran -cpp -fopenmp -Wall -Wline-truncation -Werror=line-truncation -fsyntax-only -J /tmp/asgard_linecheck -I /tmp/asgard_linecheck ...'
```

不要从仓库根目录跑该检查，因为旧 `.mod` 文件可能污染诊断。

## 22. 验证矩阵

| 数值处理 | 验证重点 |
| --- | --- |
| \(R\)-坐标输运 | 率是否用 \(\mathrm{d}t'/\mathrm{d}R\) 转换 |
| `fullhide_1d` 隐式输运 | 刚性冷却下 \(N_e\)、\(\nu_c\)、高能截断是否平滑 |
| `dg_1d` troubled positive-kernel | smooth-cell 空间阶数为 P12 的 \(O(\Delta y^{13})\)，端到端电子推进受 BE 限制为 \(O(\Delta R)\)；同时检查非负、活动支撑连续、无元素边界零洞、无多重 grid-scale sawtooth turns、辐射结果平滑；当前 RS DG sawtooth-turn 判据失败，作为待修诊断入口保留 |
| 子步压缩 | 与未压缩小步结果的粒子数和谱形一致 |
| 多求解器互证 | `fullhide_1d`、`charint_1d`、`weno5_1d` 的 break frequency 趋势一致 |
| `PhotonFieldState` | 每赫兹/每能量雅可比因子与壳层体积、逃逸时间一致 |
| 正式强子算子 | 源项、汇项、损失率来自同一微物理算子 |
| 联合反馈 | 弱反馈回归分离基线 |
| 对级联 | 对源项、对同步辐射种子光子场和光子存活率随壳层平滑 |
| \(\chi\)-EATS | 薄壳极限回归旧等到达时间面 |
| 结构化喷流投影重用 | 只在求解状态不变时复用投影 |

## 23. 不做的数值处理

ASGARD 的数值方法刻意不包含以下处理：

- 不在光变后处理层平滑非平滑结果。
- 不用经验时间平移修正峰时。
- 不给内部不可能状态添加防御性兜底。
- 不把观测者光度直接当局域光子数密度。
- 不把 \(\chi\) 分辨电子投影推广解释为 \(\chi\) 分辨强子输运。
- 不在逆康普顿介导级联尚未闭合前把它写成已实现功能。

这些限制不是保守措辞，而是为了保持物理方程和数值算法一一对应。真正的 bug 应在产生它的阶段暴露并修复。

## 24. 推荐阅读顺序

理解数值方法时建议按以下顺序阅读：

1. `doc/physical_processes.md`：先理解每个物理过程。
2. `doc/algorithm_workflow.md`：再看公开 API 到 `SolveState` 的路径。
3. `doc/electron_solver_algorithms.md`：深入电子输运器。
4. `doc/joint_secondary_feedback_algorithm.md`：理解联合反馈。
5. `doc/validation_and_benchmarks.md`：学习如何运行保留的 benchmark 入口和记录验收。

读代码时优先看这些入口：

- `src/Dynamics/Dynamics_forward.f90`
- `src/Electron/electron_forward_fullhide_1d.f90`
- `src/Electron/electron_transport_common.f90`
- `src/Electron/electron_cooling_kernel.f90`
- `src/Electron/electron_cooling_ssa_kernel.f90`
- `src/Electron/electron_cooling_ic_kernel.f90`
- `src/Electron/electron_cooling_y_kernel.f90`
- `src/Hadronic/hadronic_forward_1d.f90`
- `src/Interpolation/SED_interpolation.f90`
