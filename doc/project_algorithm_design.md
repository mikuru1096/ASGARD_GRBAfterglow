# ASGARD 项目算法设计总纲

本文档是公式级算法总纲，集中说明动力学、电子、辐射、强子和二级反馈。公开
配置、源码索引、observer 投影与验证流程分别见 `public_api.md`、
`code_overview.md`、`numerical_methods.md` 和 `developer_guide.md`。

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

只在最后通过等到达时间面投影得到。输运、冷却、强子相互作用和二级反馈先在
同一局域坐标闭合，不能用观测者时间后处理修补局域物理。

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

`electron_solver` 表示不同离散，不是失败时的性能 fallback。默认生产路径为
`fullhide_1d`；characteristic、SLC、DG、T2G 和 WENO 用于独立互证；二维路径
解析有限 q-shell。完整选择依据见 `electron_solver_algorithms.md`。

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
