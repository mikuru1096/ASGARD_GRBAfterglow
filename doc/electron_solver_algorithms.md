# 电子求解器算法说明

本文说明当前 ASGARD 正向/反向激波电子输运求解器的方程、离散方式和适用角色。相关 Fortran 源码位于 `src/Electron/`。代码标识保持英文，物理解释使用中文。

## 1. 统一变量

电子求解器对外输出仍使用电子 Lorentz 因子的对数网格

\[
x_\gamma=\log_{10}\gamma_e .
\]

辐射核需要的固定网格变量是

\[
dN_{x_\gamma}
=\frac{\mathrm{d}N}{\mathrm{d}\log_{10}\gamma_e}
=\gamma_e\ln 10\frac{\mathrm{d}N}{\mathrm{d}\gamma_e}.
\]

当前 1D transport 公共层允许内部能量坐标 \(y\) 不等于 \(x_\gamma\)。`fullhide_1d` 和 `dg_1d` 的正式 1D 路径使用均匀 log-four-velocity 坐标

\[
y
=
\log_{10}
\left[
1+\frac{\gamma_e^2-1}{\gamma_*^2-1}
\right],
\qquad
\gamma_*=2 .
\]

低能区接近四速度/动量空间，高能区渐近于 \(2x_\gamma-\log_{10}(\gamma_*^2-1)\)。任意内部坐标 \(y\) 上的守恒变量写为

\[
dN_y
=
\frac{\mathrm{d}N}{\mathrm{d}y}
=
dN_{x_\gamma}
\frac{\mathrm{d}x_\gamma}{\mathrm{d}y}.
\]

对 log-four-velocity 坐标，

\[
\gamma_e(y)
=
\left[
1+(\gamma_*^2-1)(10^y-1)
\right]^{1/2},
\]

\[
J_y
\equiv
\frac{\mathrm{d}x_\gamma}{\mathrm{d}y}
=
\frac{(\gamma_*^2-1)10^y}{2\gamma_e^2}.
\]

因此输出给辐射核时先做守恒投影，再除以 \(J_y\) 回到 \(dN_{x_\gamma}\) 或 \(\mathrm{d}N/\mathrm{d}\gamma_e\)。

2D 求解器额外展开下游有限壳层质量坐标

\[
q\in[0,q_{\rm active}],
\qquad
q_{\rm active}=1-\left(1-\frac14\right)^4 .
\]

\(q=0\) 是激波前沿，\(q=q_{\rm active}\) 是当前主动下游壳层的有限外边界。代码仍把投影输出字段命名为 `chi_*`，但 2D transport 主坐标是 `q_grid/q_face/dq`。投影和诊断中使用的 BM 等效坐标由

\[
q_{\rm tail}=1-q,
\qquad
\chi_{\rm BM}(q)=q_{\rm tail}^{-\alpha},
\qquad
\alpha=\frac{4-k}{3-k}
\]

给出，其中 \(k=0\) 是 ISM，\(k=2\) 是 wind。电子守恒变量为

\[
U(x,q,R)
=\frac{\mathrm{d}N}
{\mathrm{d}\log_{10}\gamma_e\,\mathrm{d}q}.
\]

这里 \(R\) 是激波半径，\(\Gamma_{\rm sh}\) 是激波 Lorentz 因子。

## 2. 共享微物理

所有电子求解器共享同一套壳层微物理：

- 外部介质密度 \(n_e\)，由 ISM、wind、密度跳变或注入历史给出。
- 局域磁场

\[
B=0.39\sqrt{\epsilon_B n_e\Gamma(\Gamma-1)}.
\]

- 最大电子 Lorentz 因子

\[
\gamma_{\max}
=\frac{3m_ec^2}{\sqrt{8Be^3}}.
\]

- \(\gamma_m\) 由 `electron_gamma_m_exact` 计算，覆盖 \(p>2\)、\(1<p<2\) 和 \(p=2\)。
- \(\gamma_c\) 用于诊断；2D 路径还可以从实际损失系数反推出 \(\gamma_c\)。
- 注入源项使用幂律加指数截断：

\[
Q_e(\gamma_e)\propto \gamma_e^{-p}
\exp\!\left(-\frac{\gamma_e}{\gamma_{\max}}\right).
\]

冷却系数由 `electron_cooling_kernel.f90` 组装。它包含同步辐射、同步自康普顿/逆康普顿和同步自吸收回热修正；绝热项通过输运系数中的几何项进入。`index_y` 的含义是：

- `0`：只使用同步辐射主项

\[
dE_\ell =
\left(f_r-\dot{\gamma}_{\rm SSA}C_{\rm scale}\right)\gamma_e.
\]

- `1`：数值逆康普顿冷却。
- `2`：Nakar 型 \(Y\) 参数。
- `3`：Fan 型 \(Y\) 参数。

## 3. 1D 连续方程

定义

\[
N(\gamma_e,R)=\frac{\mathrm{d}N}{\mathrm{d}\gamma_e},
\qquad
\gamma'_e=\frac{\mathrm{d}\gamma_e}{\mathrm{d}R}.
\]

连续方程写作

\[
\frac{\partial N}{\partial R}
+\frac{\partial}{\partial\gamma_e}
\left(\gamma'_eN\right)
=Q.
\]

变换到 \(x_\gamma=\log_{10}\gamma_e\) 后，

\[
\frac{\partial dN_{x_\gamma}}{\partial R}
+\frac{\partial}{\partial x_\gamma}
\left(A_{x_\gamma}dN_{x_\gamma}\right)
=S_x,
\]

其中

\[
A_{x_\gamma}=\frac{\gamma'_e}{\gamma_e\ln 10},
\qquad
S_x=\gamma_e\ln 10\,Q.
\]

若内部使用 \(y\) 坐标，方程乘以雅可比 \(J_y=\mathrm{d}x_\gamma/\mathrm{d}y\) 后仍保持守恒形式：

\[
\frac{\partial dN_y}{\partial R}
+\frac{\partial}{\partial y}
\left(A_y dN_y\right)
=S_y,
\]

\[
A_y=\frac{A_{x_\gamma}}{J_y},
\qquad
S_y=J_yS_x.
\]

实际离散中的 log-gamma 面速度近似为

\[
A_{x_\gamma,{\rm face}}
\simeq
\frac{dE_{\ell,{\rm mean}}+\alpha_{\rm ad}}{\ln 10},
\]

其中 \(dE_\ell\) 是 cooling kernel 汇总后的辐射、SSA 和 IC 损失系数，\(\alpha_{\rm ad}\) 是 wrapper 显式传入的绝热率。FS 使用 \(\alpha_{\rm ad}=1/R\)，RS 使用区域 3 体积演化给出的对应 branch rate；DG 核本身不重新推导 shock 物理。

## 4. 2D 连续方程

2D 路径推进 \(U(x,q,R)\)：

\[
\frac{\partial U}{\partial R}
+\frac{\partial(A_xU)}{\partial x}
\frac{\partial(A_q U)}{\partial q}
=
\frac{\partial}{\partial q}
\left(D_q\frac{\partial U}{\partial q}\right)
+S.
\]

\(x\) 方向是能量冷却平流，\(q\) 方向包含下游流体平流、扩散和局域源项。当前 `electron_transport_2d_kernel.f90` 中的几何关系为

\[
A_q(q,R)
=
\frac{3-k}{R}(1-q),
\]

\[
\log\frac{r(q)}{R}
=
w\left[-\frac{\chi_{\rm BM}(q)-1}{4(4-k)\Gamma_{\rm sh}^2}\right]
+
(1-w)\left[\frac{\log(1-q)}{4(3-k)}\right],
\]

\[
w=\frac{u_{\rm sh}^2}{1+u_{\rm sh}^2},
\qquad
u_{\rm sh}=\sqrt{\Gamma_{\rm sh}^2-1},
\qquad
\beta_f=\frac{u_{\rm sh}}{\Gamma_{\rm sh}}.
\]

局域下游 Lorentz 因子也用同一个 BM/Sedov-Taylor rapidity bridge：

\[
u_{\rm BM}(q)=\frac{u_{\rm sh}}{\sqrt{\chi_{\rm BM}(q)}},
\qquad
\beta_{\rm ST}(q)=\beta_f
\exp\left[\frac{\log(1-q)}{4(3-k)}\right],
\qquad
u_{\rm ST}(q)=\frac{\beta_{\rm ST}}{\sqrt{1-\beta_{\rm ST}^2}},
\]

\[
u(q)=\sinh\left[
w\sinh^{-1}u_{\rm BM}
+
(1-w)\sinh^{-1}u_{\rm ST}
\right].
\]

代码由 \(\beta_{\rm ST}\) 得到 \(u_{\rm ST}\)。这组公式的目标是让高 Lorentz 因子时回到 BM 下游，低 Lorentz 因子时保持有限正半径壳层。

`charint_2d` 中只有 shock-front 侧的 \(q\) cell 接收新注入源；已经 advect/diffuse 到下游的其它 \(q\) cell 在能量方向调用公共无源特征线冷却 primitive。`chi_radius_cm`、`chi_gamma_bulk` 和 `chi_dvolume_weight` 是从 \(q\) cell 映射到 observer projection 的有限厚壳层几何字段，不表示强子或 SSC 已经 \(\chi\) 局域。

## 5. 求解器角色

| 求解器 | 能量方向 | 厚壳方向 | 当前角色 |
| --- | --- | --- | --- |
| `fullhide_1d` | 一阶隐式迎风 | 无 | 默认稳定基线 |
| `slc1_1d` | 半拉格朗日 + Strang splitting | 无 | 方法比较 |
| `charint_1d` | 特征线保守重映射 | 无 | 1D 高保真路径 |
| `t2g1_1d` | BDF2 三层推进 | 无 | 时间推进比较 |
| `weno5_1d` | WENO5 + SSP RK3 | 无 | 高阶显式比较 |
| `dg_1d` | 多域 LGL 谱元 DG | 无 | FS/RS opt-in 高阶路径 |
| `fullhide_2d` | 一阶隐式迎风 | \(q\)-mass 隐式平流-扩散 | 2D 物理基线 |
| `charint_2d` | 特征线重映射 | \(q\)-mass 隐式平流-扩散 | 2D 加速混合版 |

`fullhide_1d` 最稳健，但数值扩散较强。`dg_1d` 使用 moving multi-domain LGL 谱元，在固定 1D 电子网格数下保留更尖锐的冷却断点和高能截断；它是 opt-in 路径，不替代默认 `fullhide_1d`。`dg_1d` 默认启用 troubled-cell positive-kernel 滤波，局部衰减高阶 Legendre 模态并保留 cell average，用于控制 Gibbs 振荡；尖锐曲率本身不作为失败条件。`charint_1d` 对高能截止和冷却断点也更锋利，但子步选择更敏感。`fullhide_2d` 是当前最完整的有限 \(q\)-shell 电子路径。`charint_2d` 只把能量方向换成特征线重映射，\(q\) 方向仍使用隐式平流-扩散，因为扩散项不能写成单纯特征线更新。

## 6. `dg_1d` 谱元路径

`dg_1d` 的核心守恒变量是 \(dN_y\)，当前 FS/RS 正式路径使用上节的 log-four-velocity 坐标。固定输出网格仍由 `Numerics.num_electron_gamma` 控制，当前基准使用 121 个电子网格；DG 内部使用 P12 LGL 谱元和移动分段网格，不要求输出格点与谱元节点一一对应。

设某个谱元 \(I_k=[y_L^k,y_R^k]\)，线性映射为

\[
y(r)
=
y_L^k+\frac{r+1}{2}\Delta y^k,
\qquad
\Delta y^k=y_R^k-y_L^k,
\qquad
r\in[-1,1].
\]

P12 表示每个谱元有 \(N+1=13\) 个 LGL 节点 \(r_i\) 和权重 \(w_i\)。节点插值写为

\[
U_h^k(r,R)
=
\sum_{i=0}^{N}U_i^k(R)\ell_i(r),
\qquad
U_i^k=dN_y(y(r_i)).
\]

当前网格把每个物理段再切成 6 个谱元。物理段边界来自当前活动的 \(\gamma_m\)、\(\gamma_c\)、高能活动边界和低能断点；落在计算域外或彼此太近的断点不会生成独立段。这个设计的目的不是减少到最小自由度，而是在 121 个正式输出格点下让冷却断点、高能指数尾和 RS 低能 kinetic bump 不被单个宽 cell 控制。

对守恒方程

\[
\frac{\partial U}{\partial R}
+\frac{\partial (aU)}{\partial y}
=S
\]

乘以 \(\ell_i\) 并在 LGL 节点上做质量矩阵求逆，得到半离散式

\[
\frac{\mathrm{d}U_i^k}{\mathrm{d}R}
=
\frac{2}{\Delta y^k w_i}
\sum_{j=0}^{N}
w_jD_{ji}\,a_jU_j^k
+{\cal F}_i^k
+S_i^k .
\]

这里 \(D_{ji}=\ell_i'(r_j)\)。边界项只作用在端点：

\[
{\cal F}_0^k
=
\frac{2}{\Delta y^k w_0}\hat f_{L}^{\,k},
\qquad
{\cal F}_N^k
=
-\frac{2}{\Delta y^k w_N}\hat f_{R}^{\,k}.
\]

冷却主导时 \(a<0\)，信息从高能流向低能，右端面取本谱元右端状态，左端面取右侧相邻谱元的状态。若局部 SSA 回热或其它项让 \(a>0\)，代码切到全局 dense assemble，按每个端面速度符号自动选择迎风状态。时间推进当前固定为后向 Euler：

\[
\left(I-\Delta R K\right)U^{n+1}
=
U^n+\Delta R S^{n+1}.
\]

FS 与 primary RS 共用同一个 DG transport kernel，但物理源项保持各自定义：

- FS 源项仍使用 \(\gamma^{-p}\)。
- RS kinetic 源项使用 \((\gamma-1)^{-p}\)。当前 primary RS DG 为避免 \(\gamma_m\simeq1\) 附近的 kinetic singular 被节点采样放大，先在正式 log-four-velocity cell 上构造与 fullhide 相同的 cell-averaged kinetic 源项，再投影到 DG 节点并按注入粒子数守恒归一化；不使用节点解析源项作为默认路径。
- SSA/IC/同步冷却仍由 cooling kernel 给出；DG transport 层只接收最终损失率、绝热率和源项。

FS DG 在壳层之间持续保存 DG state，并在 \(\gamma_m,\gamma_c,\gamma_{\max}\) 移动时重投影到新谱元。RS DG 同样跨 shell 持久保存 DG state；只在输出给固定电子网格和辐射核时使用守恒 cell 投影。RS DG 对断点附近的正谱多项式使用 cell-average preserving positivity limiter，保持元素平均粒子数，不作为后处理 smoothing。

默认 troubled positive-kernel 在每次 DG advance 和 remesh projection 后执行。每个谱元先投影到 Legendre 模态：

\[
U_h^k(r)
=
\sum_{m=0}^{N}\hat U_m^kP_m(r),
\qquad
\hat U_m^k
=
\frac{2m+1}{\Delta y^k}
\int_{y_L^k}^{y_R^k}
U_h^k(y)P_m(r(y))\,\mathrm{d}y .
\]

若节点值出现负值，或最高 6 个模态满足

\[
\frac{\sum_{m=m_{\rm hi}}^{N}|\hat U_m^k|}
{\sum_{m=0}^{N}|\hat U_m^k|}
>2\times10^{-2},
\qquad
m_{\rm hi}=\max(1,N-5),
\]

该谱元被标记为 troubled。默认模式只滤波 troubled 谱元及左右相邻谱元。滤波不改变 0 阶模态：

\[
\hat U_0^{k,{\rm new}}=\hat U_0^k,
\qquad
\hat U_m^{k,{\rm new}}
=
\sigma_m\hat U_m^k,\quad m\ge1 .
\]

默认 Jackson 因子为

\[
\sigma_m^J
=
\frac{
(N-m+2)\cos\theta_m
+\sin\theta_m/\tan\theta_1
}{N+2},
\qquad
\theta_m=\frac{\pi m}{N+2}.
\]

诊断模式 `ASGARD_DG1D_POSITIVE_KERNEL=fejer` 使用

\[
\sigma_m^F=1-\frac{m}{N+1}.
\]

环境变量未设置时使用 `troubled`。`0/off/false/none` 关闭滤波，`1/on/true/jackson` 使用全域 Jackson，`fejer` 使用全域 Fejer。默认 `troubled` 是局部高阶模态衰减，不是输出端 smoothing；它保留 cell average，因此不改变谱元平均粒子数。

当前时间推进基线不改变：DG 每个 shell 使用 10 个基础子步，密度跳变时 FS DG 仍由 jump-aware limiter 自适应缩步。RS fullhide 保持原 100--1000 冷却子步，RS DG 与 FS DG 对齐为 10 个基础子步。

### 6.1 当前收敛阶

`dg_1d` 的收敛阶必须按“能量空间离散”和“半径/时间推进”分开读。当前谱元阶数为 \(N=12\)。在光滑、未触发 positive-kernel 的谱元内，LGL-DG 对线性守恒平流的空间离散具有

\[
\|U-U_h\|_{L^2(I_k)}
=
{\cal O}\!\left((\Delta y_k)^{N+1}\right)
=
{\cal O}\!\left((\Delta y_k)^{13}\right),
\]

前提是速度、源项和映射在该谱元内足够光滑，且断点没有落在谱元内部。LGL 节点插值本身是 12 次多项式，单元平均由 \(m=0\) 模态给出；源项和投影在当前实现中按谱元守恒积分归一化，因此 smooth-cell 的主误差来自高阶截断而不是输出网格采样。

端到端电子推进的当前有效阶数不是 13 阶。原因是半径方向仍使用后向 Euler：

\[
(I-\Delta R K)U^{n+1}=U^n+\Delta R S^{n+1},
\]

所以在 cooling/source 系数随 \(R\) 平滑变化的普通 shell-to-shell 推进中，全局时间/半径阶数为

\[
U(R)-U_{\Delta R}(R)
=
{\cal O}(\Delta R).
\]

10 个 DG 基础子步只是降低 \(\Delta R\) 的误差系数，不改变这一阶数。密度跳变附近的 jump-aware 缩步同样服务于解析物理变化率，而不是把 BE 改成高阶时间格式。

在 troubled 谱元或真实断点附近，Jackson/Fejer kernel 会衰减高阶模态。此时不声明点值的 13 阶收敛；验收口径改为守恒平均不变、支撑连续、无负值和辐射结果平滑。断点被移动谱元边界捕捉时，断点两侧的 smooth-cell 仍按高阶空间收敛，断点本身按分段光滑问题处理。

Primary RS crossing 后的纯冷却段不是 BE 时间推进。若冷却系数在累计步内固定，解析特征线映射对冷却方程本身是精确的；剩余误差来自 crossing 谱的守恒重映射和冷却系数随壳层变化的取样，而不是 repeated upwind diffusion。

RS 高能尾的当前验收口径是分辨率收敛而不是强贴合 121 格 fullhide。121 格 fullhide 在 post-crossing 高能截断处有明显一阶隐式上风数值扩散；`num_gam_e=241` 时 fullhide 截断向 DG 收缩。因此 `dg_1d` 的较窄高能尾不视为单独错误，除非高分辨率 fullhide、解析特征线或守恒谱形诊断也显示同向偏差。

Primary RS 在完全 crossing 之后没有新的反向激波注入。`fullhide_1d` 和 `dg_1d` 都不再继续用各自的 shell-to-shell 数值输运推进这段纯冷却演化，而是在 crossing 后缓存固定网格上的 \(dN_y\)，随后只累计当前网格边界回溯到 crossing 网格的 characteristic edge map；输出给辐射核时从 crossing 谱一次做保守重映射。若冷却律写成

\[
\frac{\mathrm{d}\gamma}{\mathrm{d}\tau}
=
-a\gamma^2-b\gamma,
\]

则从 crossing 到当前时刻的解析映射为

\[
\gamma(\tau)
=
\frac{\gamma_0e^{-b\tau}}
{1+(a/b)\gamma_0(1-e^{-b\tau})},
\qquad b\ne0,
\]

\[
\gamma(\tau)
=
\frac{\gamma_0}{1+a\gamma_0\tau},
\qquad b=0.
\]

代码用这个映射把当前 cell edge 回溯到 crossing edge，再做守恒积分。低能端使用闭合边界，冷却到电子网格下界以下的数目保留在最低能单元。这样避免 repeated projection 或一阶隐式上风扩散生成的硬截断、平台、过宽高能尾和 late-time 清零。当前有效支撑图保留在 `output/asgard_doc/dg_radiation_stability_scan/postcross_direct_map_effective_support/`。

## 7. 子步与近似

2D 子步取

\[
\Delta R_{\rm try}
=\min(\Delta R_x,\Delta R_q,\Delta R_D).
\]

quick/formal 2D 路径使用隐式算子的精度步长：

\[
\Delta R_x \propto \frac{4\Delta x}{|A_x|},
\qquad
\Delta R_q \propto \frac{4\Delta q}{|A_q|}.
\]

2D 冷却使用 reduced frequency grid：

\[
N_{\nu,{\rm cool}}=\min(6,N_\nu).
\]

Python 调度层限制有效线程数：

\[
N_{\rm thread,eff}=\min(N_{\rm thread},N_\chi,4),
\]

避免小 \(\chi\) 或小频率网格上 OpenMP 开销支配运行时间。

## 8. 代码映射

| 方程项 | 主要实现 |
| --- | --- |
| \(\gamma_e\) 网格 | `electron_common.f90` |
| \(q/\chi_{\rm BM}\) 几何 | `electron_transport_2d_kernel.f90` |
| 注入源项 \(Q_e\) | `electron_injection_profiles.f90` |
| 冷却系数 \(A_x\) | `electron_cooling_kernel.f90` |
| FS/RS shared 1D transport wrapper | `electron_shell_transport_common.f90` |
| 1D 隐式迎风 | `electron_transport_common.f90` |
| 1D LGL-DG | `electron_transport_dg_1d_kernel.f90` |
| 1D 特征线 | `electron_transport_common.f90` |
| 2D \(q\) 输运 | `electron_transport_2d_kernel.f90` |
| 2D 能量输运 | `electron_transport_2d_kernel.f90` |
| 同步辐射谱 | `electron_radiation_kernel.f90` |
| 历史光子场 | `electron_seed_history_kernel.f90` |

## 9. 验收风险

需要重点检查：

- \(\nu_m\)、\(\nu_c\)、\(\nu_a\) 是否随时间平滑。
- 相邻观测时刻的总谱是否出现锯齿。
- 电子总数是否在无物理注入或逃逸事件时突变。
- 平滑参数扫描中高频端是否出现台阶式截断。
- `dg_1d` 输出给辐射核前的固定网格谱是否出现元素边界零洞或多重 grid-scale sawtooth turns。
- FS density-jump 场景使用默认 troubled positive-kernel 后，仍应检查支撑连续、非负和辐射结果平滑，不把真实尖断点的高曲率当作失败。当前 RS DG 谱形诊断暴露 sawtooth-turn 问题，保留为待修真实问题，不作为绿色门槛。
- `charint_2d` 的 \(x\) 子步是否过粗导致高能端提前消失。
- 2D reduced cooling bands 在窄频带问题中是否引入系统偏差。

如果这些诊断不平滑，应回到动力学、冷却源项、网格映射和投影检查，不能在后处理层 smoothing。
