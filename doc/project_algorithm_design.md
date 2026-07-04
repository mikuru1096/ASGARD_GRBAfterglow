# ASGARD 项目算法设计总纲

本文档给出 ASGARD 当前工作树的公式级算法设计。它从公开 API 的语义入口开始，说明这些参数如何被归一化为 runtime 配置，随后写清网格和 Jacobian、动力学离散、电子输运离散、光子场与强子算子、joint feedback、等到达时间面投影、缓存和验证。物理图像的完整说明见 `doc/project_physics_design.md`；这里侧重“代码如何把物理方程离散并组织成求解流程”。

ASGARD 的主变量不是观测光变，而是壳层半径 \(R\) 上的一组局域状态：

\[
\Gamma(R),\quad
M_{\rm sw}(R),\quad
B'(R),\quad
N_e(\gamma_e,R),\quad
N_p(\gamma_p,R),\quad
n_\gamma(\nu,R).
\]

观测通量

\[
F_\nu(t_{\rm obs})
\]

只在最后通过等到达时间面投影得到。这个分层很重要：输运、冷却、强子相互作用和二级反馈必须先在同一局域坐标中闭合，不能用观测者时间上的后处理去修补局域物理。

## 1. 公开 API 到 `RuntimeConfig`

典型入口为：

```python
model = Model(
    jet=top_hat_jet(...),
    medium=UniformMedium(...),
    observer=Observer(...),
    fwd_rad=Radiation(...),
    numerics=Numerics(...),
    observer_grid=ObserverGrid(...),
    solver_options=SolverOptions(...),
    reverse_shock=ReverseShock(...),
    hadronic=Hadronic(...),
)
flux = model.flux_density_grid(times_s, nu_hz)
```

内部对象链是：

```text
public dataclasses
-> RuntimeConfig
-> SimulationSetup
-> SolveState
-> ObsState / FluxComponents
```

`RuntimeConfig` 是唯一进入 runtime 和 Fortran wrapper 的归一化配置。公开 API 的字段映射遵循物理语义，而不是简单复制用户输入名：

| Public 字段 | Runtime 字段 | 算法含义 |
| --- | --- | --- |
| `jet.energy_iso_erg` | `e_iso` | 各向同性等效能量，进入动力学初始能量和结构化喷流 patch 权重。 |
| `jet.initial_lorentz_factor` | `eta_0` | 初始 bulk Lorentz factor，决定早期 coasting 和反向激波相对速度。 |
| `jet.opening_angle_rad` | `opening_angle_jet` | top-hat 边界或 structured patch 截断边界。 |
| `medium.number_density_cm3` | `d_ne` | ISM 数密度，进入 swept mass 和磁场归一化。 |
| `medium.a_star` | `a_star` | \(k=2\) wind 归一化；`a_star < 0` 表示非 wind baseline。 |
| `observer.z` | `z` | 频率、时间和 luminosity distance 的红移因子。 |
| `observer.viewing_angle_rad` | `theta_v` | Doppler factor 和 EATS 几何投影角。 |
| `radiation.epsilon_e` | `epsilon_e` | 电子注入能量分数。 |
| `radiation.epsilon_B` | `epsilon_b` | 湍流磁场能量分数。 |
| `radiation.accelerated_electron_fraction` | `f_e` | 参与非热加速的电子比例 \(\xi_N\)。 |
| `radiation.p` | `p` | 注入谱指数。 |
| `radiation.proton_energy_fraction` | `hadronic.epsilon_p` | 质子注入能量分数。 |
| `numerics.num_radius` | `num_r` | 动力学和壳层输运的半径/壳层数。 |
| `numerics.num_electron_gamma` | `num_gam_e` | 电子对数能量网格数。 |
| `numerics.num_photon_frequency` | `num_nu` | 光子频率/seed 网格数。 |
| `solver_options.electron_solver` | `electron_solver` | 选择电子输运离散。 |
| `solver_options.geometry_projection` | `geometry_kernel` | 选择 observer projection kernel。 |
| `solver_options.electron_photon_coupling` | `electron_photon_coupling` | 选择 separated 或 joint feedback 主链。 |

`SimulationSetup` 从 `RuntimeConfig` 派生四个基础数组：

\[
\texttt{boundary},\quad
\nu_j,\quad
t_{{\rm obs},q},\quad
d_L.
\]

其中 `boundary` 是 Fortran 动力学/电子核使用的紧凑参数向量；\(\nu_j\) 是本地辐射和 seed photon 频率网格；\(t_{{\rm obs},q}\) 是用户查询或内部投影的观测者时间网格；\(d_L\) 是由红移或用户覆盖值决定的 luminosity distance。

`SolveState` 是一次求解的不可混淆状态包：

```text
dynamics
electron
photon_field
hadronic
reverse_emission
observer
components
adapter_reports
timings
```

算法约束是：上游 stage 只通过 dataclass 字段向下游传值，observer projection 不回写 `electron`、`photon_field` 或 `hadronic`。

## 2. Runtime 主状态机

默认 separated 主链为：

```text
solve_dynamics
-> solve_electron
-> build_photon_field
-> solve_hadronic
-> solve_rsemission
-> assemble_observer_stage
-> SolveState
```

joint 主链为：

```text
solve_dynamics
-> primary electron solve
-> photon-field build
-> fixed internal feedback iterations:
      hadronic solve
      photon survival/source update
      secondary e± source assembly
      coupled electron solve
      photon-field rebuild
-> reverse shock stage if allowed by config
-> observer stage
```

每个 stage 的数学角色是：

\[
\mathcal{D}:\mathbf{p}\mapsto
\{R_i,\Gamma_i,t_{{\rm axis},i},M_{{\rm sw},i}\},
\]

\[
\mathcal{E}:
\{\mathcal{D},\nu_j,\mathbf{p}\}
\mapsto
\{N_e(\gamma_k,R_i),L_{\nu,{\rm syn}}(R_i),n_{\nu,{\rm syn}}(R_i)\},
\]

\[
\mathcal{H}:
\{N_p,n_\gamma,N_e\}
\mapsto
\{L_{\nu,{\rm had}},Q_{e^\pm},\tau_\nu,L_{\nu,\nu}\},
\]

\[
\mathcal{O}:
\{L'_{\nu'}(R_i),\Gamma_i,R_i\}
\mapsto
F_\nu(t_{{\rm obs},q}).
\]

不支持的组合在进入高代价 kernel 前由配置边界报错。这是系统边界校验，不是 fallback。

## 3. 网格、变量和 Jacobian

### 3.1 半径和时间网格

动力学输出离散壳层

\[
R_0<R_1<\cdots<R_{N_R-1}.
\]

壳层宽度使用

\[
\Delta R_i =
\begin{cases}
R_1-R_0, & i=0,\\
R_i-R_{i-1}, & i>0.
\end{cases}
\]

共动时间步不是观测者时间步，而是

\[
\Delta t'_i
=\frac{\Delta R_i}{\beta_i\Gamma_i c},
\qquad
\beta_i=\sqrt{1-\Gamma_i^{-2}}.
\]

因此任何共动率 \(\lambda'_i\,[{\rm s}^{-1}]\) 进入 \(R\) 坐标时都写成

\[
\lambda_{R,i}
=
\lambda'_i\frac{{\rm d}t'}{{\rm d}R}
=
\frac{\lambda'_i}{\beta_i\Gamma_i c}.
\]

这条公式同时用于电子冷却、BH/pγ/pp 损失、二级粒子注入和 photon survival。

### 3.2 对数能量网格和保守变量

电子、质子和中微子相关网格都使用对数坐标。以电子为例：

\[
x_k=\ln\gamma_{e,k},
\qquad
\gamma_{e,k}=\exp x_k.
\]

若物理分布定义为 \(N_e(\gamma)\)，对数网格上的守恒 cell content 是

\[
\mathcal{N}_{e,k}
=
\int_{x_{k-1/2}}^{x_{k+1/2}}
N_e(\gamma)\,\frac{{\rm d}\gamma}{{\rm d}x}\,{\rm d}x
=
\int_{x_{k-1/2}}^{x_{k+1/2}}
N_e(\gamma)\,\gamma\,{\rm d}x.
\]

因此连续变量转换必须带 Jacobian：

\[
N_e(\gamma)\,{\rm d}\gamma
=
\widetilde{N}_e(x)\,{\rm d}x,
\qquad
\widetilde{N}_e(x)=\gamma N_e(\gamma).
\]

光子频率同理：

\[
y_j=\ln\nu_j,\qquad
L_\nu\,{\rm d}\nu
=
\widetilde{L}_y\,{\rm d}y,
\qquad
\widetilde{L}_y=\nu L_\nu.
\]

当代码在 \(\nu L_\nu\)、\(L_\nu\)、photon number density 或 GeV energy density 之间转换时，必须显式使用这些 Jacobian。不能把每 Hz 的量直接当成每 log-bin 的量。

### 3.3 光子数密度单位转换

同步辐射或 SSC luminosity \(L_\nu\) 转为壳层局域 photon number density 时，需要壳层几何和逃逸时间。抽象形式为

\[
u_\nu'
\simeq
\frac{L_\nu' t'_{\rm esc}}{V'},
\qquad
n_\nu'
=
\frac{u_\nu'}{h\nu'}.
\]

其中 \(V'\) 是共动壳层体积，\(t'_{\rm esc}\) 是由 shell thickness/geometry 定义的逃逸时间。进入 pγ/BH 的 energy-space photon density 时还要做

\[
\epsilon_\gamma = h\nu_\gamma,
\qquad
n_{\epsilon}'\,{\rm d}\epsilon
=
n_\nu'\,{\rm d}\nu,
\qquad
n_{\epsilon}'
=
\frac{n_\nu'}{h}.
\]

这就是为什么 observer luminosity component 不能直接作为 local target photon density 使用；中间缺少体积、逃逸时间和变量 Jacobian。

## 4. 动力学离散

### 4.1 正向激波

正向激波求解器输出

\[
\{R_i,\Gamma_i,t_{{\rm axis},i},M_{{\rm sw},i}\}.
\]

外部介质质量增量为

\[
{\rm d}M_{\rm sw}
=
4\pi R^2 \rho_{\rm ext}(R)\,{\rm d}R
\]

对 top-hat 等效能量路径，动力学右端可以抽象写成

\[
\frac{{\rm d}\mathbf{y}}{{\rm d}s}
=
\mathbf{f}\!\left(\mathbf{y},R;\mathbf{p}_{\rm jet},\rho_{\rm ext},L_{\rm inj}\right),
\qquad
\mathbf{y}=(\Gamma,M_{\rm sw},R,t_{\rm axis},\ldots).
\]

这里 \(s\) 是 Fortran 内部的推进变量，可为 \(R\)、实验室时间或与观测轴时间相关的参数；对后续 stage 的正式合同只有输出数组，不依赖内部参数化。验收时使用解析标度作为物理直觉：

\[
\Gamma_{\rm BM,ISM}(R)\propto R^{-3/2},
\qquad
\Gamma_{\rm BM,wind}(R)\propto R^{-1/2}.
\]

能量注入以 \(L_{\rm inj}(t)\) 改变 blast-wave 能量：

\[
\frac{{\rm d}E}{{\rm d}t_{\rm lab}}=L_{\rm inj}(t_{\rm lab}),
\]

密度跳变以 \(\rho_{\rm ext}(R)\) 改变 \(M_{\rm sw}\) 的增长率。两者都应该在 \(\Gamma(R)\) 上产生连续响应；若光变有尖峰，先检查 \(\rho(R)\)、\(\Gamma(R)\) 和 \(B'(R)\) 是否平滑。

### 4.2 磁场和特征频率

正向激波下游能量密度的量级为

\[
e'\sim 4\Gamma^2 n_{\rm ext}m_pc^2.
\]

湍流磁场由

\[
\frac{B'^2}{8\pi}=\epsilon_B e'
\]

给出，即

\[
B'=\sqrt{8\pi\epsilon_B e'}.
\]

电子最小注入 Lorentz factor 的基线量级为

\[
\gamma_m
=
1+
\frac{\epsilon_e}{\xi_N}
\frac{p-2}{p-1}
\frac{m_p}{m_e}
(\Gamma_{\rm rel}-1),
\]

其中正向激波通常取 \(\Gamma_{\rm rel}\approx\Gamma\)，反向激波使用 \(\gamma_{34}\)。同步辐射特征频率为

\[
\nu_{\rm syn}(\gamma)
=
\frac{3q_e}{4\pi m_e c}
B'\gamma^2
\frac{\delta}{1+z},
\]

局域 stage 保存的是未投影或按内部约定归一化的壳层谱；observer 因子在投影 stage 统一处理。

### 4.3 反向激波和次级反向激波

反向激波用 shocked ejecta 区域 3 的热状态推进：

\[
e_3'=\frac{U_3}{V_3},
\qquad
B_{3,{\rm turb}}'=\sqrt{8\pi\epsilon_{B,3}e_3'}.
\]

磁化上游使用

\[
\sigma=\frac{B_4^2}{4\pi\Gamma_4^2\rho_4c^2},
\]

它同时改变 baryonic ejecta mass、pressure-balance 触发、MHD jump compression 和下游热比内能。当前有限强度 jump 给出

\[
C=\frac{u_{4s}}{u_{3s}},
\qquad
\epsilon_{\rm th,3}=\frac{h_3-1}{\hat\gamma},
\]

并用 \(\mathrm{d}U_{3,{\rm sh}}=\epsilon_{\rm th,3}\mathrm{d}M_3c^2\) 推进 region-3 shock heating。总磁场按湍流和有序分量合成：

\[
B_3'=\sqrt{B_{3,{\rm turb}}'^2+B_{3,{\rm ord}}'^2}.
\]

反向激波电子注入能标使用 shock-front 相对 Lorentz factor：

\[
\gamma_{34}
=
\Gamma_3\Gamma_4(1-\beta_3\beta_4),
\]

\[
\gamma_{m,3}
=
1+
\frac{\epsilon_{e,3}}{\xi_{N,3}}
\frac{p_3-2}{p_3-1}
\frac{m_p}{m_e}
\epsilon_{\rm th,3}.
\]

穿越前质量变量满足

\[
0<M_3/M_{\rm ej}<1,
\]

穿越后

\[
\frac{{\rm d}M_3}{{\rm d}R}=0.
\]

因此代码不能让一个积分步同时跨越“有新 ejecta 被 shock 加热”和“已无新 ejecta 注入”两个方程组。multi-density 次级反向激波在密度增强处新建局域 branch；该 branch 的耗散能、体积和磁场是独立诊断和辐射源，不覆盖原 forward-shock 电子状态。

## 5. 电子输运离散

### 5.1 连续方程

电子输运在 \(R\) 坐标上写成

\[
\frac{\partial N_e(\gamma,R)}{\partial R}
=
Q_{e,R}(\gamma,R)
-
\frac{\partial}{\partial\gamma}
\left[
\dot{\gamma}_{R}(\gamma,R)N_e(\gamma,R)
\right].
\]

其中

\[
\dot{\gamma}_{R}
=
\frac{\dot{\gamma}'}{\beta\Gamma c}
=
\dot{\gamma}_{{\rm ad},R}
\dot{\gamma}_{{\rm syn},R}
\dot{\gamma}_{{\rm IC},R}
\dot{\gamma}_{{\rm SSA},R}.
\]

shock 注入源项满足粒子数和能量归一化：

\[
\int Q_{e,R}(\gamma)\,{\rm d}\gamma
=
\frac{{\rm d}N_{e,{\rm acc}}}{{\rm d}R},
\]

\[
\int (\gamma-1)m_ec^2 Q_{e,R}(\gamma)\,{\rm d}\gamma
=
\epsilon_e\frac{{\rm d}E_{\rm diss}}{{\rm d}R}.
\]

非热形状为

\[
Q_{e,R}(\gamma)=A_e(R)\gamma^{-p}
\exp\!\left(-\frac{\gamma}{\gamma_{\max}}\right),
\qquad
\gamma\ge\gamma_m.
\]

`accelerated_electron_fraction` \(\xi_N\) 同时影响粒子数和 \(\gamma_m\)，不能只改谱归一化。

### 5.2 对数网格有限体积形式

令 \(x=\ln\gamma\)，\(\widetilde{N}=\gamma N\)。连续方程转为

\[
\frac{\partial \widetilde{N}}{\partial R}
=
\widetilde{Q}_{R}
-
\frac{\partial}{\partial x}
\left[
v_x\widetilde{N}
\right],
\]

其中

\[
v_x
=
\frac{{\rm d}x}{{\rm d}R}
=
\frac{1}{\gamma}\frac{{\rm d}\gamma}{{\rm d}R}
=
\frac{\dot{\gamma}_{R}}{\gamma},
\qquad
\widetilde{Q}_R=\gamma Q_R.
\]

第 \(k\) 个 cell 的保守更新为

\[
\mathcal{N}_{k}^{i+1}
=
\mathcal{N}_{k}^{i}
\Delta R_i\,\mathcal{Q}_{k}^{i}
-
\frac{\Delta R_i}{\Delta x_k}
\left(
\mathcal{F}_{k+1/2}^{i+\theta}
-
\mathcal{F}_{k-1/2}^{i+\theta}
\right).
\]

\(\theta=0\) 是显式更新，\(\theta=1\) 是隐式更新；`fullhide_1d` 的关键思想是把强冷却区域的通量项放入稳定的隐式/保守 stencil，而不是在谱后端裁剪负值或平滑尖峰。

面通量的迎风选择为

\[
\mathcal{F}_{k+1/2}
=
\begin{cases}
v_{k+1/2}\mathcal{N}_{k}, & v_{k+1/2}>0,\\
v_{k+1/2}\mathcal{N}_{k+1}, & v_{k+1/2}<0.
\end{cases}
\]

对冷却问题通常 \(v_x<0\)，粒子从高 \(\gamma\) 向低 \(\gamma\) 移动。若谱随 \(R\) 出现锯齿或非物理反弹，应检查 \(v_{k+1/2}\)、Jacobian 和边界通量，而不是给输出光变加 smoothing。

### 5.3 冷却项公式

绝热冷却按体积膨胀给出量级：

\[
\dot{\gamma}'_{\rm ad}
=
-\frac{\gamma}{3}
\frac{{\rm d}\ln V'}{{\rm d}t'}.
\]

同步冷却为

\[
\dot{\gamma}'_{\rm syn}
=
-\frac{4}{3}
\frac{\sigma_T c}{m_ec^2}
U_B'\gamma^2,
\qquad
U_B'=\frac{B'^2}{8\pi}.
\]

IC/SSC 冷却的 Thomson 量级为

\[
\dot{\gamma}'_{\rm IC,T}
=
-\frac{4}{3}
\frac{\sigma_T c}{m_ec^2}
U_{\rm ph}'\gamma^2.
\]

KN/Jones 路径使用核函数把电子和 photon seed 同时约束：

\[
\dot{\gamma}'_{\rm IC}
=
-c\int{\rm d}\epsilon\,n_\epsilon'
\int{\rm d}\epsilon_s\,
(\epsilon_s-\epsilon)\,
\frac{{\rm d}\sigma_{\rm IC}}{{\rm d}\epsilon_s}
(\gamma,\epsilon,\epsilon_s).
\]

对应 emissivity 为

\[
j_{\epsilon_s,{\rm IC}}'
=
c\epsilon_s
\int{\rm d}\gamma\,N_e(\gamma)
\int{\rm d}\epsilon\,n_\epsilon'
\frac{{\rm d}\sigma_{\rm IC}}{{\rm d}\epsilon_s}.
\]

joint 预算路径要求冷却和 photon production 使用同一个 \(N_e\)、\(n_\epsilon'\) 和 IC kernel；否则电子能量损失与 SSC luminosity 不能闭合。

### 5.4 电子 solver 的算法含义

`electron_solver` 的枚举值不是简单性能选项，而是不同离散假设：

| solver | 离散含义 | 适用判断 |
| --- | --- | --- |
| `fullhide_1d` | 1D \(\gamma\)-space 保守输运，固定或受控子步，默认生产路径。 | 光变、拟合、joint feedback 的基线。 |
| `slc1_1d` | 半拉格朗日特征移动，减少强冷却 CFL 压力。 | 用于交叉检查冷却轨迹。 |
| `charint_1d` | 沿特征线积分 \( {\rm d}\gamma/{\rm d}R=\dot{\gamma}_R \)。 | 用于检查谱峰移动和冷却 break。 |
| `dg_1d` | P12 LGL-DG，高阶谱元输运，默认 troubled positive-kernel。 | FS/RS opt-in 高阶路径；光滑谱元空间 \(O(\Delta y^{13})\)，端到端电子推进当前 \(O(\Delta R)\)。 |
| `t2g1_1d` | legacy 隐式输运。 | 回归比较路径。 |
| `weno5_1d` | 高阶重构的谱输运。 | 用于检查数值耗散和振荡。 |
| `fullhide_2d` / `charint_2d` | \((\gamma,\chi)\) 分辨电子输运。 | 服务有限厚度/chi 投影研究，不自动扩展强子反馈。 |

## 6. 同步辐射、SSA 和 photon seed

单电子同步辐射功率写为

\[
P_{\nu}'(\gamma)
=
\frac{\sqrt{3}q_e^3B'}{m_ec^2}
F\!\left(\frac{\nu'}{\nu_c'(\gamma)}\right),
\]

\[
\nu_c'(\gamma)
=
\frac{3q_eB'}{4\pi m_ec}\gamma^2.
\]

壳层同步 luminosity 为

\[
L_{\nu,{\rm syn}}'(R_i)
=
\int N_e(\gamma,R_i)P_{\nu}'(\gamma,R_i)\,{\rm d}\gamma.
\]

在对数网格上实际计算为

\[
L_{\nu_j}'(R_i)
\simeq
\sum_k
N_{e,k}(R_i)
P_{\nu_j}'(\gamma_k,R_i)\Delta\gamma_k,
\qquad
\Delta\gamma_k=\gamma_k\Delta x_k.
\]

SSA absorption coefficient 的形式为

\[
\alpha_\nu'
=
-\frac{1}{8\pi m_e\nu'^2}
\int P_\nu'(\gamma)\gamma^2
\frac{\partial}{\partial\gamma}
\left[
\frac{N_e(\gamma)}{\gamma^2}
\right]
{\rm d}\gamma.
\]

光学深度为

\[
\tau_{\nu,{\rm SSA}}'=\alpha_\nu'\ell',
\]

逃逸 luminosity 使用 transfer factor

\[
L_{\nu,{\rm esc}}'
=
L_{\nu,{\rm em}}'
\frac{1-\exp(-\tau_\nu')}{\tau_\nu'}.
\]

这里的重点是频率 cell 之间的连续性：\(\alpha_\nu'\)、\(\tau_\nu'\) 和 escape factor 都应随 \(\nu\) 平滑变化。`index_syn_integr=2` / `fixed_grid` 是 public 快速路径；adaptive 同步积分只应作为显式诊断使用。

`PhotonFieldState` 中的字段含义为：

| 字段 | 数学对象 | 下游用途 |
| --- | --- | --- |
| `forward_syn_seed` | \(n_{\nu,{\rm syn}}'\) | 构造 SSC、absorption 和 hadronic target 的基础 seed。 |
| `hadronic_forward_ssc_seed` | SSC 相关 seed contribution | hadronic 和 observer 组装。 |
| `hadronic_target_seed` | \(n_{\nu,{\rm target}}'\) | pγ、BH、pp 相关 photon target 和 joint IC cooling。 |
| `absorption_syn_seed` | absorption 用 synch seed | \(\gamma\gamma\)、SSA 或 survival 相关路径。 |
| `absorption_ssc_seed` | absorption 用 SSC seed | 高能 photon survival。 |

joint 可以更新 seed 数值，但不能改变字段语义。

## 7. 强子算子

### 7.1 质子注入和输运

质子输运形式为

\[
\frac{\partial N_p(E_p,R)}{\partial R}
=
Q_{p,R}(E_p,R)
-
\frac{\partial}{\partial E_p}
\left[
\dot{E}_{p,R}N_p
\right]
-
\lambda_{p,R}N_p
+Q_{p,{\rm reinj},R}.
\]

质子注入归一化满足

\[
\int Q_{p,R}(E_p)\,{\rm d}E_p
=
\frac{{\rm d}N_{p,{\rm acc}}}{{\rm d}R},
\]

\[
\int E_p Q_{p,R}(E_p)\,{\rm d}E_p
=
\epsilon_p\frac{{\rm d}E_{\rm diss}}{{\rm d}R}.
\]

加速上限由加速时间和损失/动力学时间比较确定，量级为

\[
t'_{\rm acc}
=
\eta_{\rm acc}\frac{E_p'}{e c B'},
\qquad
t'_{\rm acc}
=
\min(t'_{\rm dyn},t'_{\rm syn},t'_{p\gamma},t'_{\rm BH},t'_{pp}).
\]

### 7.2 Proton synchrotron

质子同步辐射使用同样的磁场和同步辐射结构，但质量替换为 \(m_p\)：

\[
\nu_{c,p}'(\gamma_p)
=
\frac{3q_eB'}{4\pi m_pc}\gamma_p^2,
\]

\[
P_{\nu,p}'(\gamma_p)
\propto
\frac{q_e^3B'}{m_pc^2}
F\!\left(\frac{\nu'}{\nu_{c,p}'}\right).
\]

由于 \(m_p\gg m_e\)，同一 \(\gamma\) 下质子特征频率和功率都与电子差异很大；不能把 electron synch kernel 的结果简单重标为 proton synch。

### 7.3 pγ 响应

pγ 相互作用率的基本结构为

\[
t_{p\gamma}^{\prime -1}(E_p)
=
c
\int{\rm d}\epsilon\,n_\epsilon'
\int{\rm d}\bar{\epsilon}\,
\sigma_{p\gamma}(\bar{\epsilon})
\kappa_{p\gamma}(\bar{\epsilon})
\mathcal{K}(E_p,\epsilon,\bar{\epsilon}).
\]

Hummer response path 把这个积分离散为响应矩阵：

\[
Q_a(E_a)
=
\sum_{p,j}
\mathcal{R}_{a,pj}
N_p(E_{p})
n_{\epsilon_j}'
\Delta E_p\Delta\epsilon_j,
\]

其中 \(a\) 可以是 gamma、\(\pi^\pm\)、\(\mu^\pm\)、\(e^\pm\) 或 neutrino 相关产物。矩阵响应的意义是：同一个 \(N_p\) 和 \(n_\epsilon'\) 同时给出 proton loss、secondary source、photon output 和 neutrino output，不能分别用经验比例拼接。

### 7.4 Bethe-Heitler 和 pp

BH 过程为

\[
p+\gamma\rightarrow p+e^++e^-.
\]

formal kernel 输出三类量：

\[
\lambda_{\rm BH}'(E_p),
\qquad
Q_{e^\pm,{\rm BH}}'(E_e),
\qquad
\lambda_{\gamma,{\rm BH}}'(\nu).
\]

进入 \(R\) 坐标后：

\[
Q_{e^\pm,{\rm BH},R}
=
\frac{Q_{e^\pm,{\rm BH}}'}{\beta\Gamma c},
\qquad
\tau_{\nu,{\rm BH},i}
=
\lambda_{\gamma,{\rm BH},i}'\Delta t'_i.
\]

photon survival factor 为

\[
S_{\nu,i}^{\rm BH}
=
\exp(-\tau_{\nu,{\rm BH},i}).
\]

pp 过程使用 baryon target density：

\[
t_{pp}^{\prime -1}
=
n_p'c\sigma_{pp}K_{pp},
\]

并输出 secondary pairs、gamma 和 neutrino 分量。只有 formal kernel 已经给出并完成单位归一化的 \(Q_{e^\pm,R}\) 可以进入电子方程。

## 8. Joint feedback 算法

`electron_photon_coupling="separated"` 的语义是：

\[
N_e \rightarrow n_\gamma \rightarrow \mathcal{H}
\]

强子和 pair 输出不改变同一轮电子输运。`joint` 的语义是求一个 shell-level 近似闭合：

\[
N_e
\leftrightarrow
n_\gamma
\leftrightarrow
\{N_p,Q_{e^\pm},S_\nu\}.
\]

内部迭代可写为：

\[
N_e^{(0)}
=
\mathcal{E}(Q_{e,{\rm shock}},n_\gamma^{\rm syn}),
\]

\[
n_\gamma^{(m)}
=
\mathcal{P}(N_e^{(m)}),
\]

\[
\mathcal{H}^{(m)}
=
\mathcal{H}(N_p,n_\gamma^{(m)}),
\]

\[
Q_{e^\pm,R}^{(m)}
=
Q_{e^\pm,{\rm BH},R}^{(m)}
+Q_{e^\pm,pp,R}^{(m)}
+Q_{e^\pm,\gamma\gamma,R}^{(m)}
+Q_{e^\pm,p\gamma,R}^{(m)},
\]

\[
n_{\gamma,{\rm surv}}^{(m)}
=
n_\gamma^{(m)}
\exp[-\tau_{\gamma\gamma}^{(m)}-\tau_{\rm BH}^{(m)}-\tau_{p\gamma}^{(m)}],
\]

\[
N_e^{(m+1)}
=
\mathcal{E}_{\rm coupled}
\left(
Q_{e,{\rm shock}}+Q_{e^\pm,R}^{(m)},
n_{\gamma,{\rm surv}}^{(m)}
\right).
\]

其中所有 rate 都用

\[
\frac{{\rm d}t'}{{\rm d}R}
=
\frac{1}{\beta\Gamma c}
\]

换算到 \(R\) 坐标。当前 fixed small iteration count 是内部实现细节，不是公开收敛参数。joint 模式的配置边界是算法合同的一部分：forward shock、`fullhide_1d`、formal `am3_1d`、BH enabled、Hummer pγ response、数值 IC cooling 和固定电子子步必须同时成立。

## 9. Pair cascade

\(\gamma\gamma\) 吸收阈值由

\[
\epsilon_1\epsilon_2(1-\cos\psi)
\ge
2(m_ec^2)^2
\]

给出。光学深度抽象为

\[
\tau_{\gamma\gamma}(E_\gamma)
=
\int{\rm d}\ell
\int{\rm d}\epsilon\,n_\epsilon'
\int{\rm d}\Omega\,
(1-\cos\psi)
\sigma_{\gamma\gamma}(s).
\]

pair injection 满足能量映射：

\[
E_{e^+}+E_{e^-}\simeq E_\gamma+\epsilon,
\]

在离散网格中写作

\[
Q_{e^\pm,k}
=
\sum_j
\mathcal{M}_{kj}^{\gamma\gamma}
L_{\gamma,j}
\left(1-e^{-\tau_{\gamma\gamma,j}}\right).
\]

`pair_cascade_iterations > 1` 时使用 shell-sequence time-dependent pair/synch cascade：

```text
absorbed high-energy photons
-> pair injection on gamma_e grid
-> pair synchrotron luminosity
-> pair synch seed added to photon field
-> next shell / next cascade iteration
```

该路径只闭合 pair synch cascade。IC-mediated electromagnetic cascade 不是当前 contract，因此算法说明和用户脚本都不能把 pair-synch seed 当成完整电磁级联。

## 10. Reverse shock、cross-zone IC 和结构化喷流

反向激波 emission stage 使用 region-3 电子谱和磁场：

\[
L_{\nu,{\rm RS}}'
=
\int N_{e,3}(\gamma)P_\nu'(\gamma,B_3')\,{\rm d}\gamma.
\]

RS SSC 的形式与 FS SSC 相同，但 seed 和 electron distribution 来自 region 3。cross-zone IC 是两个 zone 的交叉卷积：

\[
j_{\epsilon_s,{\rm FS}\leftarrow{\rm RS}}'
=
c\epsilon_s
\int{\rm d}\gamma\,N_{e,{\rm FS}}(\gamma)
\int{\rm d}\epsilon\,n_{\epsilon,{\rm RS}}'
\frac{{\rm d}\sigma_{\rm IC}}{{\rm d}\epsilon_s},
\]

\[
j_{\epsilon_s,{\rm RS}\leftarrow{\rm FS}}'
=
c\epsilon_s
\int{\rm d}\gamma\,N_{e,{\rm RS}}(\gamma)
\int{\rm d}\epsilon\,n_{\epsilon,{\rm FS}}'
\frac{{\rm d}\sigma_{\rm IC}}{{\rm d}\epsilon_s}.
\]

结构化喷流把角向 profile 离散成 patches：

\[
E_{\rm iso}=E_{\rm iso}(\theta,\phi),
\qquad
\Gamma_0=\Gamma_0(\theta,\phi).
\]

每个 patch 独立运行局域 1D solve：

\[
F_{\nu,q}
=
\sum_a
w_a
F_{\nu,q}^{(a)}
\left[
E_{\rm iso}(\theta_a,\phi_a),
\Gamma_0(\theta_a,\phi_a),
\theta_{\rm obs}
\right].
\]

patch 权重 \(w_a\) 来自角向面积或 dominant-region 采样。结构化喷流是角向求和算法，不表示 hadronic 或 pair feedback 已经拥有 \(\chi\)-local/patch-local transport。

## 11. EATS 观测者投影

### 11.1 Doppler 和红移

壳层速度与视线夹角为 \(\alpha\)，Doppler factor 为

\[
\delta
=
\frac{1}{\Gamma(1-\beta\cos\alpha)}.
\]

频率变换为

\[
\nu
=
\frac{\delta}{1+z}\nu',
\qquad
\nu'
=
\frac{1+z}{\delta}\nu.
\]

时间变换的轴向部分由动力学给出；几何延迟为

\[
t_{\rm obs}
=
(1+z)
\left[
t_{\rm lab}
-\frac{R\cos\alpha}{c}
\right].
\]

薄壳近似下，若动力学 stage 已经给出轴向观测时间 \(t_{\rm axis}\)，则 off-axis patch 使用

\[
t_{\rm obs}
\simeq
t_{\rm axis}
+(1+z)\frac{R}{c}(1-\cos\alpha)
\]

或等价的 Fortran interpolation 形式。

### 11.2 Luminosity 到 flux

各向同性共动谱 luminosity 的投影量级为

\[
F_\nu(t_{\rm obs})
=
\frac{1+z}{4\pi d_L^2}
\int_{\rm EATS}
\delta^3
L_{\nu'}'(R,\Omega)
\,{\rm d}W,
\]

其中 \({\rm d}W\) 包含角向面积、壳层权重和 EATS 插值权重。若输入是 emissivity 而不是 shell luminosity，Doppler 幂次和体积权重会按实现约定重排，但频率和时间变换必须保持上述关系。

离散形式为

\[
F_{\nu_q,t_q}
=
\sum_{i,a}
W_{i,a\rightarrow q}
\frac{1+z}{4\pi d_L^2}
\delta_{i,a}^{3}
L_{\nu'_{q,i,a}}'(R_i,a)
\exp[-\tau_{\nu'_{q,i,a}}].
\]

频率插值在 log-log 空间执行：

\[
\log L_{\nu'}(\nu_q')
=
\operatorname{interp}
\left(
\log\nu_j',
\log L_{\nu_j'}',
\log\nu_q'
\right),
\]

时间/EATS 插值同样不能改变本地输运状态，只能把已有 shell luminosity 映射到查询网格。

### 11.3 `lightcurve`、`sed` 和 `chi_eats_2d`

`projection_kind="lightcurve"` 面向固定频率、多时刻查询：

\[
\{\nu_q\}_{q=1}^{N_\nu}
\quad{\rm fixed},\qquad
\{t_q\}_{q=1}^{N_t}
\quad{\rm scanned}.
\]

`projection_kind="sed"` 面向固定时刻、扫频率查询：

\[
t_q\quad{\rm fixed\ or\ few},
\qquad
\nu_q\quad{\rm scanned}.
\]

二者使用同一 `SolveState`，但投影插值目标不同。不能用 SED 插值误差解释光变峰时，也不能用 lightcurve 路径替代宽频谱积分。

`chi_eats_2d` 在 FS synchrotron+SSA 分量中把 finite \(q\)-shell 的局域半径、bulk Lorentz factor 和体积权重纳入投影。`chi_grid` 只是 \(q\) cell 的 BM 等效诊断坐标；投影实际读取 `chi_radius_cm`、`chi_gamma_bulk` 和 `chi_dvolume_weight`：

\[
F_{\nu,t}^{\rm FS,syn}
=
\sum_{i,k_\chi,a}
W_{i,k_\chi,a\rightarrow t}
\frac{1+z}{4\pi d_L^2}
\delta_{i,k_\chi,a}^3
L_{\nu'}'(R_{i,k_\chi},q_{k_\chi},a)
S_{\nu'}(R_i,q_{k_\chi}).
\]

非 chi 分量仍遵守 shell-level projection。SSC、强子、pair cascade 和 cross-zone IC 不因为 `chi_eats_2d` 自动获得 \(q\)-local 反馈。

axisymmetric structured `chi_eats_2d` 将上述求和按 theta ring 拆开：每个 active ring 独立求解 2D electron state，Python 侧只准备该 ring 的 `F_ring = l_syn_spec_chi(1+z)/(4\pi d_L^2)` 和 `Tau_ring = tau_syn_chi`，观测者投影由 `sed_chi_ring` 累加 `theta_lo/theta_hi` 与 `structured_num_phi`。`structured_adaptive_rtol > 0` 只改变 theta-ring window 的投影采样；它不改变动力学、电子输运、强子输运或 photon-field closure。

## 12. 缓存、状态复用和性能边界

`Model` 查询可以复用已经构建的 `SolveState`。复用的充分条件是底层物理状态不变：

\[
\mathbf{p}_{\rm jet},\mathbf{p}_{\rm medium},\mathbf{p}_{\rm rad},
\mathbf{p}_{\rm solver},
\nu_{\rm seed},
R{\rm\ grid}
\]

均未改变。只改变查询频率、查询时间或 observer projection 目标时，可以复用 transport state 并重跑投影。

缓存边界如下：

| 变化 | 是否可复用 transport state | 原因 |
| --- | --- | --- |
| 只改变 `nu_hz` 查询点 | 可以 | 本地 solve 已经在 seed 网格上完成。 |
| 只改变 `times_s` 查询点且覆盖在已求解时间范围内 | 可以 | 只需 EATS/interpolation。 |
| 只改变 `viewing_angle_rad` | 取决于路径 | top-hat finite \(q\)-shell benchmark 可在同一 `SolveState` 上重跑 `project_flux`；若角向结构、patch 选择或 observer cache key 参与 transport solve，则需重新求解或重建 patch。 |
| 改变 `epsilon_e`、`epsilon_B`、`p` | 不可复用 | 电子谱和 photon field 改变。 |
| 改变 `include_pgamma`、BH、pp、pair cascade | 不可复用 | hadronic/source-sink 方程改变。 |
| 改变 `num_radius`、`num_electron_gamma`、`num_photon_frequency` | 不可复用 | 离散网格和 Jacobian 改变。 |

性能报告必须区分 cold solve 和 cached query：

\[
t_{\rm cold}
=
t_{\rm dynamics}+t_{\rm electron}+t_{\rm photon}+t_{\rm hadronic}+t_{\rm projection},
\]

\[
t_{\rm cached}
\simeq
t_{\rm projection}+t_{\rm interpolation}.
\]

把 cached query 当作完整求解时间会误导算法评估。

CPU wavefront 优化只允许重排满足因果依赖的 stencil：

\[
\mathcal{N}_{k}^{i+1}
\leftarrow
\{\mathcal{N}_{k-1}^{i+1},\mathcal{N}_{k}^{i},\mathcal{N}_{k+1}^{i}\}
\]

这类依赖决定了 tile/wavefront 的合法并行方向。不能打破依赖后再用谱平滑或数值保护掩盖误差。

## 13. 验证算法

文档和代码验证分为七层：

1. 文本编码和 Markdown 构建。
2. Public API 配置边界检查。
3. Fortran wrapper build。
4. `-Wline-truncation` source closure 检查。
5. 最小端到端验证。
6. 物理量连续性和预算闭合。
7. 文献或独立 backend benchmark。

物理连续性不是审美要求，而是方程对真实时空演化的约束。无注入、无密度跳变的基线中：

\[
\Gamma(R),\quad
B'(R),\quad
\nu_m(R),\quad
\nu_c(R),\quad
F_\nu(t)
\]

应当平滑。若有密度跳变，允许导数改变，但状态变量不应出现没有物理来源的离散断裂。能量预算检查包括：

\[
E_{e,{\rm inj}}
\simeq
E_{e,{\rm stored}}
+E_{\rm rad}
+E_{\rm ad},
\]

\[
E_{p,{\rm inj}}
\simeq
E_{p,{\rm stored}}
+E_{p{\rm syn}}
+E_{p\gamma}
+E_{\rm BH}
+E_{pp}
+E_{\nu}
+E_{\rm ad}.
\]

joint feedback 还必须检查同一 microphysics operator 的源汇一致性：

\[
\Delta E_{\gamma,{\rm absorbed}}
\simeq
\Delta E_{e^\pm,{\rm injected}}
+\Delta E_{\rm escaped/diagnostic},
\]

\[
\Delta E_{e,{\rm IC\ loss}}
\simeq
\Delta E_{\gamma,{\rm IC\ emitted}}
\]

在核函数允许的数值误差内成立。

## 14. 输出和诊断设计

`FluxComponents` 保存 observer-side 分量：

```text
total
fwd_sync
fwd_ssc
fwd_hadronic_gamma
fwd_hadronic_bethe_heitler
fwd_hadronic_inverse_compton
fwd_hadronic_pair_production
rev_sync
rev_ssc
cross_ic
```

总通量是已启用分量的和：

\[
F_{\nu}^{\rm total}
=
F_{\nu}^{\rm FS,syn}
+F_{\nu}^{\rm FS,SSC}
+F_{\nu}^{\rm had}
+F_{\nu}^{\rm RS}
+F_{\nu}^{\rm crossIC}.
\]

`details()` 和 diagnostics 应能追踪回上游物理量，例如：

\[
\nu_m,\nu_c,\nu_a,\Gamma,R,B',\tau_{\gamma\gamma},Q_{e^\pm},L_\nu.
\]

诊断字段只能解释主求解器状态，不能替代正式源项；输出层不能裁剪、平滑或重新归一化物理分量来制造连续光变。

## 15. 开发准则

算法改动必须满足：

- 先明确物理方程和离散变量，再改接口或实现。
- Python 负责 orchestration、配置归一化、缓存、benchmark 和 API glue；高代价微物理核优先在 Fortran。
- 新增源项必须说明单位、网格、Jacobian、进入哪个方程以及是否反馈。
- 不支持组合直接在配置边界报错。
- 不为内部不可能状态添加防御性 fallback。
- 非光滑物理输出先查动力学、Jacobian、输运通量、seed 归一化和 EATS 投影。
- 文档只描述当前能力、边界和验收口径，不分散维护未完成事项。
