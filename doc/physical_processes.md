# 物理过程详解

本文面向第一次认真读 ASGARD 物理链的新手。它从爆波、粒子输运和辐射转移出发，解释每个物理过程在代码中的位置、输入输出和验收方式。更高层的边界见 `doc/project_physics_design.md`，实现入口见 `doc/code_overview.md`。

## 1. 坐标、参考系和状态量

ASGARD 的输运主坐标是激波半径 \(R\)，不是观测者时间。每个壳层上的本地量先在激波/共动参考系中求解，再投影到观测者平面。常用时间关系为

\[
\mathrm{d}t'
=
\frac{\mathrm{d}R}{\beta\Gamma c},
\qquad
\beta
=
\sqrt{1-\Gamma^{-2}} ,
\]

其中 \(t'\) 是共动系时间，\(\Gamma\) 是 bulk Lorentz factor。观测者时间还包含几何延迟：

\[
t_{\rm obs}
=
(1+z)
\left[
t_{\rm lab}
-
\frac{R\mu}{c}
\right],
\qquad
\mu=\cos\alpha .
\]

对 on-axis 薄壳近似，常用量级关系是

\[
t_{\rm obs}
\sim
\frac{1+z}{2\Gamma^2}
\frac{R}{c}.
\]

代码中动力学、电子谱、强子谱和局域光子场都沿 \(R_i\) 网格组织；`flux_density_grid` 只在最后阶段把本地辐射投影到用户请求的 \((t_{\rm obs},\nu_{\rm obs})\)。

## 2. 外部介质与 swept-up mass

外部介质给出被正向激波扫过的重子/电子数密度。当前公开后端支持均匀星际介质和 \(k=2\) 恒星风介质：

\[
n_{\rm ISM}(R)=n_0 ,
\qquad
n_{\rm wind}(R)
=
\frac{3.0\times10^{35}A_\ast}{R^2}.
\]

若使用 density jump，密度是显式参数化的 \(n(R)\)。正向激波 swept-up mass 的壳层增量满足量级关系

\[
\mathrm{d}M_{\rm sw}
=
4\pi R^2 n(R)m_p\,\mathrm{d}R
\]

在代码中，实际右端项使用观测者时间/半径映射、介质剖面、`index_dyn` 分支和能量注入项共同推进。密度不连续会在物理上允许光变或特征断频出现响应；没有密度跳变、注入事件或激波穿越时，非连续结构优先视为 bug。

## 3. 正向激波动力学

正向激波动力学的核心是能量、动量和 swept-up mass 的耦合演化。ASGARD 的 `Dynamics_forward.f90` 推进的主变量包括 \(\Gamma\)、扫掠质量、内能相关变量和半径。常用绝热 blast-wave 图像可写为

\[
E_{\rm iso}
\sim
\Gamma^2 M_{\rm sw}c^2
\]

在均匀介质且绝热的 Blandford-McKee 标度下有

\[
\Gamma(R)
\propto R^{-3/2},
\qquad
\Gamma(t_{\rm obs})
\propto t_{\rm obs}^{-3/8}.
\]

在 wind 介质中

\[
\Gamma(R)
\propto R^{-1/2},
\qquad
\Gamma(t_{\rm obs})
\propto t_{\rm obs}^{-1/4}.
\]

这些标度不是替代代码右端项的公式，而是基准测试和物理直觉检查。若均匀星际介质基线在减速后不接近 \(\Gamma\propto R^{-3/2}\) 的趋势，应检查动力学输入、密度、能量注入开关和时间网格。

能量注入以源项形式进入动力学。若注入 luminosity 记作 \(L_{\rm inj}(t)\)，总能量变化可表示为

\[
\frac{\mathrm{d}E}{\mathrm{d}t_{\rm lab}}
=
L_{\rm inj}(t_{\rm lab})
-
L_{\rm rad}(t_{\rm lab}) .
\]

当前公开文档只把能量注入作为动力学参数说明；不要把拟合中需要的早期平台强行归因于注入，除非残差、物理时间尺度和参数后验共同支持。

## 4. 反向激波和磁化 ejecta

反向激波描述 ejecta 被反向激波加热后的 region 3。关键相对 Lorentz factor 是

\[
\gamma_{34}
=
\Gamma_3\Gamma_4
\left(
1-\beta_3\beta_4
\right),
\]

其中 3 表示已激波加热抛射物，4 表示未激波加热抛射物。ASGARD 当前反向激波电子注入能标来自激波前沿 \(\gamma_{34}\)，而不是从区域平均 \(\Gamma\) 经验替代。

区域 3 内能密度可抽象写成

\[
u_3
=
\frac{U_3}{V_3}.
\]

turbulent 磁场由

\[
B_{\rm turb}
=
\sqrt{8\pi\epsilon_{B,3}u_3}
\]

给出。若开启 upstream magnetization \(\sigma\)，有序磁场与 turbulent field 一起贡献总场：

\[
B_3
=
\sqrt{B_{\rm turb}^2+B_{\rm ord}^2}.
\]

这里的 \(\sigma\) 不是单独的辐射后处理因子。`ReverseShock.upstream_sigma` 先把总 ejecta 能量拆成 baryonic rest-mass 与 Poynting 分量：

\[
M_{\rm ej,b}
=
\frac{E_{\rm iso}}{(1+\sigma)\Gamma_0c^2}.
\]

随后 region 4 的 \(n_4\)、pressure-balance 条件、MHD jump compression、下游热比内能、有序场压缩和磁压焓惯性都随 \(\sigma\) 改变。当前代码使用有限强度 MHD jump：

\[
C=\frac{u_{4s}}{u_{3s}},
\qquad
\epsilon_{\rm th,3}
=
\frac{h_3-1}{\hat\gamma},
\]

并用

\[
\mathrm{d}U_{3,{\rm sh}}
=
\epsilon_{\rm th,3}\mathrm{d}M_3c^2
\]

代替非磁化极限的 \((\gamma_{34}-1)\mathrm{d}M_3c^2\)。因此 \(\sigma\to0\) 的验收同时约束 baryonic mass、shock triggering、压缩比、热能源项、ordered field 和动力学 inertia。完整公式和代码映射见 `doc/shock_shell_adaptive_algorithms.md`。

验收极限是

\[
\lim_{\sigma\rightarrow0}B_3
=
B_{\rm turb},
\qquad
\lim_{\sigma\rightarrow0}F_{\nu,{\rm RS}}(\sigma)
=
F_{\nu,{\rm RS}}(0).
\]

如果这个极限失败，问题在反向激波 jump、region-3 thermal state 或磁场归一化，而不是观测投影层。

## 5. 电子注入

正向激波把一部分激波能量放入非热电子。注入谱通常写成截断幂律：

\[
Q_e(\gamma_e)
=
Q_0
\gamma_e^{-p}
\exp\!\left(-\frac{\gamma_e}{\gamma_{e,\max}}\right)
H(\gamma_e-\gamma_m).
\]

电子数归一化满足

\[
\int_{\gamma_m}^{\gamma_{\max}}
Q_e(\gamma_e)\,\mathrm{d}\gamma_e
=
\xi_N\,\mathrm{d}N_{\rm sw},
\]

电子能量归一化满足

\[
m_ec^2
\int_{\gamma_m}^{\gamma_{\max}}
(\gamma_e-1)Q_e(\gamma_e)\,\mathrm{d}\gamma_e
=
\epsilon_e\,\mathrm{d}E_{\rm sh}.
\]

当 \(p>2\) 且 \(\gamma_{\max}\gg\gamma_m\) 时，常用估计为

\[
\gamma_m
\simeq
1+
\frac{\epsilon_e}{\xi_N}
\frac{p-2}{p-1}
\frac{m_p}{m_e}
(\Gamma_{\rm rel}-1).
\]

代码使用 `electron_gamma_m_exact` 处理截断和参数边界，因此文档中的闭式表达只作为物理解释。若 \(p\) 接近 2 或 cutoff 接近 \(\gamma_m\)，必须依赖精确归一化，而不是使用上述渐近式。

## 6. 电子输运方程

电子谱 \(N_e(\gamma_e,R)\) 的连续性方程是

\[
\frac{\partial N_e}{\partial R}
+
\frac{\partial}{\partial\gamma_e}
\left[
\dot{\gamma}_{e,R}N_e
\right]
=
Q_{e,R}.
\]

这里

\[
Q_{e,R}
=
Q_{e,{\rm shock},R}
+
Q_{e,{\rm secondary},R},
\]

而冷却速度由多项组成：

\[
\dot{\gamma}_{e,R}
=
\dot{\gamma}_{\rm ad,R}
+
\dot{\gamma}_{\rm syn,R}
+
\dot{\gamma}_{\rm IC,R}
+
\dot{\gamma}_{\rm SSA,R}.
\]

同步冷却在共动时间中可写为

\[
\left.
\frac{\mathrm{d}\gamma_e}{\mathrm{d}t'}
\right|_{\rm syn}
=
-
\frac{\sigma_T B'^2}{6\pi m_ec}
\gamma_e^2 .
\]

换成 \(R\) 坐标：

\[
\dot{\gamma}_{\rm syn,R}
=
-
\frac{\sigma_T B'^2}{6\pi m_ec}
\frac{\gamma_e^2}{\beta\Gamma c}.
\]

绝热冷却量级为

\[
\left.
\frac{\mathrm{d}\gamma_e}{\mathrm{d}t'}
\right|_{\rm ad}
=
-
\frac{\gamma_e}{3}
\nabla_\mu u^\mu .
\]

在壳层模型中该项由几何膨胀和 shell state 给出。物理验收看的是 \(N_e(\gamma_e,R)\) 随 \(R\) 平滑推进，以及总电子数与注入量的预算是否合理。

## 7. 磁场与特征频率

激波下游的磁能密度与内能密度关系为

\[
\frac{B'^2}{8\pi}
=
\epsilon_B u'.
\]

单个电子的同步特征频率是

\[
\nu'_{\rm syn}(\gamma_e)
=
\frac{3eB'}{4\pi m_ec}
\gamma_e^2 .
\]

观测者频率包含 Doppler factor \(\delta\) 和红移：

\[
\nu_{\rm obs}
=
\frac{\delta}{1+z}\nu',
\qquad
\delta
=
\frac{1}{\Gamma(1-\beta\mu)} .
\]

因此

\[
\nu_m
\sim
\frac{\delta}{1+z}
\frac{3eB'}{4\pi m_ec}
\gamma_m^2,
\]

\[
\nu_c
\sim
\frac{\delta}{1+z}
\frac{3eB'}{4\pi m_ec}
\gamma_c^2.
\]

冷却 Lorentz factor 由冷却时间和动力学时间相等给出：

\[
t'_{\rm cool}(\gamma_c)
=
t'_{\rm dyn}.
\]

同步自吸收频率由 optical depth 条件定义：

\[
\tau_\nu(\nu_a)
=
\alpha_\nu(\nu_a)\ell'
=
1.
\]

runtime 不再把 \(\nu_m\)、\(\nu_c\)、\(\nu_a\) 作为默认 details 字段输出；稳定性检查应回到电子谱、磁场、SSA kernel 和最终光变。

## 8. 同步辐射与 SSA

总同步谱来自电子谱对单粒子功率的卷积：

\[
L'_{\nu'}
=
\int
N_e(\gamma_e)
P'_{\nu'}(\gamma_e,B')
\,\mathrm{d}\gamma_e .
\]

单粒子同步谱可写为

\[
P'_{\nu'}
=
\frac{\sqrt{3}e^3B'}{m_ec^2}
F\!\left(\frac{\nu'}{\nu'_c}\right),
\]

\[
F(x)
=
x
\int_x^\infty
K_{5/3}(y)\,\mathrm{d}y .
\]

SSA 的吸收系数满足

\[
\alpha'_{\nu'}
=
-
\frac{1}{8\pi m_e\nu'^2}
\int
P'_{\nu'}(\gamma_e)
\gamma_e^2
\frac{\partial}{\partial\gamma_e}
\left[
\frac{N_e(\gamma_e)}{\gamma_e^2}
\right]
\mathrm{d}\gamma_e .
\]

观测投影前，本地谱会乘上 survival factor：

\[
S_\nu
=
\exp(-\tau_\nu)
\]

或在有限厚度网格单元中使用单元平均逃逸因子：

\[
S_{\nu,{\rm cell}}
=
\exp(-\tau_{\rm front})
\frac{1-\exp(-\tau_{\rm cell})}{\tau_{\rm cell}} .
\]

后一式是 `chi_eats_2d` 中有限厚度 SSA 投影的物理含义。

## 9. SSC 与 Compton 冷却

SSC 是同一电子群散射本地同步光子场。电子冷却中常用写法为

\[
\left.
\frac{\mathrm{d}\gamma_e}{\mathrm{d}t'}
\right|_{\rm IC}
=
-
\frac{4}{3}
\frac{\sigma_T c}{m_ec^2}
U'_{\rm ph}
\gamma_e^2
f_{\rm KN}(\gamma_e,\epsilon').
\]

Thomson 极限下 \(f_{\rm KN}\rightarrow1\)。Klein-Nishina 区域中，高能电子的散射截面和能量转移下降。ASGARD 的 IC/SSC 路径要求电子冷却种子光子场与光子源种子光子场一致，否则能量预算会错。

常用 Compton 参数写成

\[
Y
=
\frac{U'_{\rm ph}}{U'_B},
\qquad
U'_B=\frac{B'^2}{8\pi}.
\]

当 \(Y\) 大时，\(\gamma_c\) 会降低，进而改变 \(\nu_c\) 和整个 SED。拟合时不能只打开 SSC photon output 而忽略 `ssc_cooling` 对电子谱的反馈含义。

## 10. 局域光子场

ASGARD 文档中必须区分两类光子量：

- 局域靶光子场：参与 IC、p-gamma、BH、gamma-gamma。
- 观测者光度：投影给观测者的辐射输出。

若 \(L'_{\nu'}\) 要转为壳层内 photon density，基本关系是

\[
u'_{\nu'}
\sim
\frac{L'_{\nu'}t'_{\rm esc}}{V'},
\qquad
n'_{\epsilon'}
=
\frac{u'_{\epsilon'}}{\epsilon'}.
\]

其中 \(t'_{\rm esc}\)、\(V'\) 和路径长度必须由壳层几何定义。不能把观测者通量直接当作靶光子数密度。

## 11. Gamma-gamma 吸收与 pair 产生

两光子对产生阈值为

\[
\epsilon_1\epsilon_2(1-\cos\theta)
\ge
2(m_ec^2)^2 .
\]

光深可写为

\[
\tau_{\gamma\gamma}(\epsilon_1)
=
\int \mathrm{d}s
\int \mathrm{d}\epsilon_2
\int \mathrm{d}\Omega
n_{\epsilon_2}
(1-\cos\theta)
\sigma_{\gamma\gamma}
(\epsilon_1,\epsilon_2,\theta).
\]

吸收后的 photon survival 是

\[
f_{\rm surv}
=
\exp(-\tau_{\gamma\gamma}).
\]

被吸收的能量进入二级 \(e^\pm\)：

\[
Q_{e^\pm}(\gamma_e)
\propto
\int
N_\gamma(\epsilon_1)
n_\gamma(\epsilon_2)
c\sigma_{\gamma\gamma}
\mathcal{K}_{\gamma\gamma}(\epsilon_1,\epsilon_2,\gamma_e)
\mathrm{d}\epsilon_1
\mathrm{d}\epsilon_2 .
\]

当前 `pair_cascade_iterations > 1` 覆盖的是 \(\gamma\gamma\) 对产生/同步辐射级联；逆康普顿介导电磁级联尚未实现，因为它需要把二级 \(e^\pm\)、IC 光子、新的 \(\gamma\gamma\) 吸收和多代源-汇方程同时闭合。

## 12. Proton synchrotron

质子同步辐射形式与电子同步类似，但质量不同：

\[
\nu'_{p,{\rm syn}}(\gamma_p)
=
\frac{3eB'}{4\pi m_pc}
\gamma_p^2 .
\]

单粒子总功率按质量缩放：

\[
P'_{p,{\rm syn}}
=
\frac{4}{3}\sigma_T cU'_B\gamma_p^2
\left(\frac{m_e}{m_p}\right)^2 .
\]

因此质子同步辐射需要极高 \(\gamma_p\) 或强 \(B'\) 才能显著贡献高能辐射。它可以作为观测者分量，也可以作为强子能量损失通道，但不会自动成为局域靶光子场，除非有明确光子场回灌契约。

## 13. p-gamma photopion

p-gamma 反应阈值由质子静止系 photon energy 决定：

\[
\bar{\epsilon}
=
\gamma_p\epsilon(1-\beta_p\mu)
\gtrsim
\bar{\epsilon}_{\rm th}.
\]

反应率的一般形式为

\[
t_{p\gamma}^{-1}(\gamma_p)
=
c
\int \mathrm{d}\epsilon
\int \mathrm{d}\Omega
n_\gamma(\epsilon,\Omega)
(1-\beta_p\mu)
\sigma_{p\gamma}(\bar{\epsilon})
\kappa_{p\gamma}(\bar{\epsilon}).
\]

Hummer 型响应算子把一个质子能量格与光子场映射到多个产物：

\[
p+\gamma
\rightarrow
n+\pi^+,
\qquad
p+\gamma
\rightarrow
p+\pi^0,
\]

以及更多多 pion channel。中性 pion 给出 gamma-rays：

\[
\pi^0\rightarrow2\gamma .
\]

带电 pion/muon 链给出 neutrino 和二级电子：

\[
\pi^+
\rightarrow
\mu^+ + \nu_\mu,
\qquad
\mu^+
\rightarrow
e^+ + \nu_e + \bar{\nu}_\mu .
\]

正式路径的关键是同一个算子同时输出质子损失、二级粒子源、光子损失和中微子光度，不能用总能量守恒在后处理中猜一个谱形。

## 14. Bethe-Heitler pair production

BH 过程为

\[
p+\gamma
\rightarrow
p+e^+ + e^- .
\]

它的非弹性系数通常小于 photopion，但阈值更低，因此在某些光子场中可成为主要二级 \(e^\pm\) 来源。一般对源项写作

\[
Q_{e,{\rm BH}}(E_e)
=
\int
N_p(E_p)
n_\gamma(\epsilon)
c\sigma_{\rm BH}(E_p,\epsilon)
\mathcal{K}_{\rm BH}(E_p,\epsilon,E_e)
\mathrm{d}E_p
\mathrm{d}\epsilon .
\]

ASGARD 联合路径要求 BH 的质子损失、对源项和光子汇来自同一个微物理算子。若只补对源项而没有光子汇，会破坏光子能量预算。

## 15. pp 相互作用

pp 过程依赖已激波加热重子靶密度：

\[
p+p
\rightarrow
\pi^0+\pi^+ + \pi^- + X .
\]

反应率量级为

\[
t_{pp}^{-1}
=
n_p c\sigma_{pp}\kappa_{pp}.
\]

产物包含 gamma-rays、neutrinos 和二级 \(e^\pm\)。pp 强度主要受重子密度控制，因此在恒星风、致密壳层或特殊环境中更值得检查。若同一参数改变导致 pp 光度非平滑，应先检查靶密度和质子输运。

## 16. 二级粒子输运

中子、pion、muon 等二级粒子满足与电子/质子类似的源-损-衰变方程：

\[
\frac{\partial N_a}{\partial t'}
=
Q_a
-
\frac{\partial}{\partial E}
\left(
\dot{E}_a N_a
\right)
-
\frac{N_a}{t'_{\rm dec,a}}
-
\frac{N_a}{t'_{\rm esc,a}} .
\]

衰变时间为

\[
t'_{\rm dec,a}
=
\gamma_a\tau_{a,0}.
\]

当 \(t'_{\rm syn,a}<t'_{\rm dec,a}\) 时，pion/muon synchrotron 会改变 neutrino 与 photon 输出。对应的 synchrotron cooling 时间为

\[
t'_{\rm syn,a}
\propto
\frac{m_a^3c}{\sigma_Tm_e^2B'^2\gamma_a}.
\]

这也是强子二级辐射必须与粒子输运一起验收的原因。

## 17. Neutrino 输出

中微子由 p-gamma、pp 和衰变链产生。中微子不参与后续电磁反馈，因此它是输出诊断而非输运源项。观测能量由

\[
E_{\nu,{\rm obs}}
=
\frac{\delta}{1+z}E'_\nu
\]

给出。若开启 `neutrino=True`，应同时检查 proton loss、pion/muon chain 和 photon target field，而不是只看最终 neutrino spectrum。

## 18. Equal-arrival-time surface 投影

本地共动谱 \(L'_{\nu'}\) 投影到观测者时，频率和强度变换为

\[
\nu_{\rm obs}
=
\frac{\delta}{1+z}\nu',
\]

\[
F_{\nu_{\rm obs}}
=
\frac{1+z}{4\pi d_L^2}
\int
\delta^3
L'_{\nu'}
\frac{\mathrm{d}\Omega}{4\pi}.
\]

ASGARD 的 `SED_interpolation.f90` 中使用的变量 `doppler` 实际是

\[
D^{-1}
=
\Gamma(1-\beta\mu),
\]

因此实现中出现 \(D^{-3}\) 与 \(\nu'=\nu_{\rm obs}(1+z)D\) 的组合。阅读代码时要注意这个命名，不要把它误读为 \(\delta\)。

EATS 的角向延迟为

\[
\Delta t_{\rm obs}
=
(1+z)
\frac{R(1-\mu)}{c}.
\]

结构化喷流或 off-axis 观测时，不同 \((\theta,\phi)\) patch 的 \(\mu\)、\(\delta\)、arrival time 都不同，因此峰时、峰值和偏振角都可能显著变化。

## 19. \(\chi\) 分辨有限厚壳层

`chi_eats_2d` 只用于正向激波同步辐射+SSA 的有限厚度等到达时间面投影。它的物理对象是 \((R,\chi,\theta,\phi)\) 体元，其中 \(\chi=1\) 接近激波前沿，较大的 \(\chi\) 对应更深的下游。

当前强子、SSC 和对级联仍是壳层级契约。不能因为电子同步辐射的观测者投影有 \(\chi\) 维，就把强子输运解释成 \(\chi\) 局域。完整 \(\chi\) 几何、输运方程、SSA survival、投影权重和薄壳极限见 `doc/shock_shell_adaptive_algorithms.md`。

## 20. 多密度增强下的次级反向激波

多密度增强反向激波用于描述均匀 ISM 中多个平滑密度增强触发的次级反向激波电子同步辐射。它不是 wind bump、结构化喷流、强子、RS SSC 或 \(\chi\) 分辨输运的通用替代模型。

密度剖面使用数组合同：

```text
jump_r_cm
jump_factor
jump_width_log10
```

对应的光滑密度增强可写为

\[
n(R)
=
n_0
\left[
1+
\sum_j
(f_j-1)
\exp\left(
-
\frac{
(\log_{10}R-\log_{10}R_j)^2
}{
2w_j^2
}
\right)
\right].
\]

旧 `r_tr/f_jump/f_wide` 是单 bump 兼容入口；当多数组非空时，以数组合同为准。

次级反向激波采用四区图像。区域 1 是密度增强前的冷外介质，只提供 \(n_1(R)\)。区域 2 是透射正向激波下游。区域 4 是密度增强前已经被正向激波加热的旧 shocked shell。区域 3 是被次级反向激波再激波的旧 shocked shell，新增辐射只来自这部分额外耗散能，不能把旧 FS 电子辐射重复计入。

区域 4 是热上游，应使用其已有 \(p_4,e_4,n_4,\Gamma_4\)，不能套用抛射物反向激波的冷上游公式。局部 Riemann 问题、压缩比、\(u_{{\rm diss},3}=e_3-e_4C^{4/3}\)、\(\gamma_{m,3}\)、branch reservoir 和 `rev.sync` 合成公式见 `doc/shock_shell_adaptive_algorithms.md`。

物理验收应检查：

- 单 bump 数组与旧 `r_tr/f_jump/f_wide` 密度剖面等价。
- 无多 bump 时原有 RS benchmark 不变。
- 每个 bump 的 \(\Gamma_c,p_2,p_3,\gamma_{43},u_{{\rm diss},3},B_3\) 随 \(R\) 平滑。
- 区域 3 注入能量积分等于 RH 给出的次级 RS 额外耗散能乘 \(\epsilon_e\)。
- FS 不出现没有物理来源的尖锐重亮；RS 对 bump 更敏感，但必须随 bump 宽度和幅度连续变化。

该图像对应 Nakar & Granot 2007 的 density-jump 四区结构和平滑观测响应约束，也与 Uhm & Zhang 2014 中密度增强/空洞触发 long-lived reverse shock 的动力学机制相容。代码验收仍以当前实现的状态量连续性和能量预算为准。

## 21. 物理验收清单

每个正式结果至少检查：

- \(\Gamma(R)\)、\(B'(R)\)、\(\nu_m(R)\)、\(\nu_c(R)\)、\(\nu_a(R)\) 是否连续。
- \(N_e(\gamma,R)\) 和 \(N_p(\gamma,R)\) 是否随 \(R\) 平滑演化。
- \(F_\nu(t)\) 和 SED 是否在物理事件之外无孤立尖峰。
- SSC 冷却与 SSC 光子输出是否来自同一种子光子场。
- BH、pp、\(\gamma\gamma\) 的二级粒子源与光子/质子汇是否来自同一算子。
- `upstream_sigma -> 0` 是否回到非磁化反向激波基线；内部旧名 `reverse_sigma` 指同一量。
- 弱反馈联合算例是否回到分离基线。

这些检查失败时，不在投影层或后处理中加 smoothing、time shift 或经验 floor。应回到产生异常的物理状态量。
