# 2D 壳层、反向激波与自适应网格

本文集中说明 ASGARD 当前网页文档中最容易混用的三个模块：\(\chi\) 分辨有限厚壳层、抛射物反向激波与密度增强触发的次级反向激波、自适应网格算法。这里写的是当前代码契约；未实现的 \(\chi\)-local 强子输运、\(\chi\)-local SSC 和 IC-mediated electromagnetic cascade 不在本页扩展。

代码标识保持英文，例如 `chi_eats_2d`、`fullhide_2d`、`reverse_sigma`。正文术语按 `doc/terminology.md`：反向激波、次级反向激波、密度增强、壳层级、等到达时间面。

## 1. \(\chi\) 分辨有限厚壳层

### 1.1 物理坐标

正向激波半径为 \(R\)，激波 Lorentz 因子为 \(\Gamma_{\rm sh}\)。当前 2D 电子输运的主坐标不是无限延伸的 \(\eta=\log_{10}\chi\)，而是有限主动壳层质量坐标 \(q\)：

\[
q\in[0,q_{\rm active}],
\qquad
q_{\rm active}
=
1-\left(1-\frac{1}{4}\right)^4 .
\]

\(q=0\) 是激波前沿，\(q=q_{\rm active}\) 是强激波压缩给出的有限主动下游外边界量级。代码入口是 `compute_q_geometry`：

\[
q_{k+1/2}=k\Delta q,
\qquad
\Delta q=\frac{q_{\rm active}}{N_\chi},
\qquad
q_k=\left(k-\frac12\right)\Delta q .
\]

为了和 BM 解析图像、投影层字段名以及已有 `chi_eats_2d` 术语相连，代码为每个 \(q\) cell 定义 BM 等效坐标

\[
\chi_{\rm BM}(q)
=
(1-q)^{-\alpha},
\qquad
\alpha=\frac{4-k}{3-k},
\]

其中 \(k=0\) 表示 ISM，\(k=2\) 表示 wind。`chi_grid`、`chi_radius_cm`、`chi_gamma_bulk` 和 `chi_dvolume_weight` 保留历史命名，但当前语义是 finite \(q\)-shell cell 的投影几何。

cell 半径用 BM 极限和低 Lorentz 因子有限半径近似做 rapidity bridge。设

\[
u_f=\sqrt{\Gamma_{\rm sh}^2-1},
\qquad
\beta_f=\frac{u_f}{\Gamma_{\rm sh}},
\qquad
w=\frac{u_f^2}{1+u_f^2}.
\]

\[
\log x_{\rm BM}
=
-
\frac{\chi_{\rm BM}-1}
{4(4-k)\Gamma_{\rm sh}^2},
\qquad
\log x_{\rm ST}
=
\frac{\log(1-q)}{4(3-k)} .
\]

\[
\frac{r(q)}{R}
=
\exp\left[
w\log x_{\rm BM}+(1-w)\log x_{\rm ST}
\right].
\]

局域下游四速度用同一个权重桥接：

\[
\Gamma_{\rm BM}(q)
=
\sqrt{1+\frac{u_f^2}{\chi_{\rm BM}}},
\qquad
u_{\rm BM}(q)=\frac{u_f}{\sqrt{\chi_{\rm BM}}},
\qquad
\beta_{\rm ST}(q)
=
\beta_f\exp\left[\frac{\log(1-q)}{4(3-k)}\right].
\]

\[
u_{\rm ST}
=
\frac{\beta_{\rm ST}}{\sqrt{1-\beta_{\rm ST}^2}},
\qquad
u(q)
=
\sinh\left[
w\sinh^{-1}u_{\rm BM}
+
(1-w)\sinh^{-1}u_{\rm ST}
\right],
\]

\[
\Gamma(q)=\sqrt{1+u(q)^2},
\qquad
\beta(q)=\frac{u(q)}{\Gamma(q)}.
\]

共动厚度由相邻 \(q\) face 半径差得到：

\[
\Delta x'_k
=
\Gamma(q_k)
\left[
r(q_{k-1/2})-r(q_{k+1/2})
\right].
\]

这套坐标保证投影使用正半径、有限厚度和局域 \(\Gamma(q)\)，同时在 ultra-relativistic 极限保留 BM 下游结构。

### 1.2 电子输运方程

2D 电子状态使用

\[
U(x,q,R)
=
\frac{\mathrm{d}N_e}
{\mathrm{d}x\,\mathrm{d}q},
\qquad
x=\log_{10}\gamma_e .
\]

连续方程写成

\[
\frac{\partial U}{\partial R}
+
\frac{\partial(A_xU)}{\partial x}
+
\frac{\partial(A_q U)}{\partial q}
=
\frac{\partial}{\partial q}
\left(
D_q
\frac{\partial U}{\partial q}
\right)
+
S(x,q,R).
\]

能量方向系数来自冷却和绝热项：

\[
A_x
=
\frac{1}{\ln 10}
\frac{1}{\gamma_e}
\frac{\mathrm{d}\gamma_e}{\mathrm{d}R},
\]

\[
\frac{\mathrm{d}\gamma_e}{\mathrm{d}R}
=
\frac{\mathrm{d}t'}{\mathrm{d}R}
\left(
\dot{\gamma}'_{\rm syn}
+
\dot{\gamma}'_{\rm IC}
+
\dot{\gamma}'_{\rm ad}
+
\dot{\gamma}'_{\rm SSA}
\right).
\]

\(q\) 方向输运系数在 2D electron transport 的步长估计、隐式矩阵组装和特征线推进处内联计算：

\[
A_q(q,R)
=
\frac{3-k}{R}(q_{\rm active}-q),
\quad
q_{\rm active}
=
1-\left(1-\frac{1}{\sigma_{\rm strong}}\right)^{\sigma_{\rm strong}}.
\]

局域绝热冷却系数使用有限 \(q\)-shell 的下游速度场散度。代码在相邻 cell 中对 \(\beta(q)\) 与 \(r(q)\) 做有限差分：

\[
\left(\nabla\cdot \mathbf{v}\right)_q
=
c
\left[
\frac{2\beta(q)}{r(q)}
+
\frac{\partial \beta}{\partial r}
\right],
\]

\[
\left(\frac{\mathrm{d}x}{\mathrm{d}R}\right)_{\rm ad}
=
\frac{1}{3\beta_f c\ln 10}
\left(\nabla\cdot\mathbf{v}\right)_q .
\]

### 1.3 离散推进

`fullhide_2d` 对 \(q\) 方向使用隐式有限体积形式。设 cell 平均量为 \(U_k^n\)，步长为 \(\Delta R\)，通量为

\[
F_{k+1/2}
=
A^+_{k+1/2}U_k^{n+1}
+
A^-_{k+1/2}U_{k+1}^{n+1}
-
D_{k+1/2}
\frac{U_{k+1}^{n+1}-U_k^{n+1}}{\Delta q}.
\]

则每个能量 bin 上求解三对角系统

\[
U_k^{n+1}
+
\frac{\Delta R}{\Delta q}
\left(
F_{k+1/2}^{n+1}
-
F_{k-1/2}^{n+1}
\right)
=
U_k^n+\Delta R\,S_{k,{\rm shock}},
\]

其中 \(S_{k,{\rm shock}}\) 只注入到 shock-front 侧的 \(q\) cell。`charint_2d` 在 \(q\) 对流部分用特征线回溯：

\[
\frac{\mathrm{d}q}{\mathrm{d}R}=A_q(q,R),
\qquad
q_{\rm back}
=
q_{\rm face}(R_{n+1})
-
\int_{R_n}^{R_{n+1}}A_q\,\mathrm{d}R,
\]

然后用守恒 PPM 前缀积分重映射：

\[
U_k^{n+1}
=
\frac{
\int_{q_{{\rm back},k-1/2}}^{q_{{\rm back},k+1/2}}
U^n(q)\,\mathrm{d}q
}{\Delta q}.
\]

2D 子步限制来自 \(q\) 方向 CFL：

\[
\Delta R_q
=
C_{\rm CFL}
\frac{\Delta q}
{\max_k |A_{q,k+1/2}|}.
\]

### 1.4 \(\chi\) 分辨等到达时间面投影

`chi_eats_2d` 只替换正向激波同步辐射和同步自吸收的 lightcurve projection。对角向面元 \((\theta_j,\phi_j)\) 定义

\[
\mu_j
=
\cos\theta_{\rm obs}\cos\theta_j
+
\sin\theta_{\rm obs}\sin\theta_j\cos\phi_j .
\]

每个 \((R_i,\chi_k)\) 体元的观测者时间为

\[
t_{{\rm obs},i k j}
=
t_{{\rm axis},i}
+
(1+z)
\frac{
R_{{\rm front},i}
-
R_{\chi,ik}\mu_j
}{c}.
\]

代码变量 `doppler` 是

\[
D_{ikj}
=
\Gamma_{\chi,ik}(1-\beta_{\chi,ik}\mu_j)
=
\delta_{ikj}^{-1},
\]

所以

\[
\nu'
=
(1+z)D_{ikj}\nu_{\rm obs},
\qquad
L_{\nu_{\rm obs}}
\propto
D_{ikj}^{-3}
L'_{\nu'} .
\]

离散投影为

\[
F_{\nu_{\rm obs}}(t_m)
=
\sum_{i,k,j}
W_{\Omega,j}
W_{\chi,ik}
D_{ikj}^{-3}
L'_{\nu'}(R_i,\chi_k)
S_{\nu'}(R_i,\chi_k)
\mathcal{I}_{ikj\rightarrow m},
\]

其中 \(\mathcal{I}_{ikj\rightarrow m}\) 是沿 \(R_i\rightarrow R_{i+1}\) 段对目标观测时间 bin 的线性插值权重。

SSA survival 不取 cell 前边界单点值，而取 optical-depth coordinate 的 cell 平均：

\[
S_{\nu',k}
=
\exp(-\tau_{{\rm front},k})
\frac{1-\exp(-\tau_{{\rm cell},k})}
{\tau_{{\rm cell},k}},
\]

\[
\tau_{{\rm front},k}
=
\sum_{q<k}\tau_q .
\]

transport \(\chi\) 网格到 projection \(\chi\) 网格使用 cell-overlap 守恒重映射。若 transport cell 为 \(m\)，projection cell 为 \(k\)，重叠长度为

\[
\Delta\chi_{km}
=
\max\left[
0,\,
\min(\chi_{k+1/2}^{\rm p},\chi_{m+1/2}^{\rm t})
-
\max(\chi_{k-1/2}^{\rm p},\chi_{m-1/2}^{\rm t})
\right],
\]

则

\[
P_k^{\rm p}
=
\frac{1}{\Delta\chi_k^{\rm p}}
\sum_m P_m^{\rm t}\Delta\chi_{km},
\qquad
\tau_k^{\rm p}
=
\sum_m \tau_m^{\rm t}
\frac{\Delta\chi_{km}}{\Delta\chi_m^{\rm t}} .
\]

薄壳极限要求

\[
W_{\chi,k}\rightarrow \delta_{k k_0},
\qquad
R_{\chi,k}\rightarrow R,
\qquad
\Gamma_{\chi,k}\rightarrow\Gamma_{\rm shell},
\qquad
\tau_{\chi,k}\rightarrow0,
\]

此时 `sed_interpolation_chi` 回到壳层级等到达时间面投影。

## 2. 抛射物反向激波

### 2.1 四区图像

抛射物反向激波使用标准区域编号：

| 区域 | 含义 |
| --- | --- |
| 1 | 未激波外介质 |
| 2 | 正向激波下游 |
| 3 | 反向激波加热后的抛射物 |
| 4 | 未激波抛射物 |

区域 3 和区域 2 共享接触间断速度：

\[
\Gamma_2=\Gamma_3,
\qquad
\beta_2=\beta_3.
\]

反向激波电子注入能标必须使用 shock-front 相对 Lorentz 因子

\[
\gamma_{34}
=
\Gamma_3\Gamma_4(1-\beta_3\beta_4).
\]

当前动力学实现中 \(\Gamma_4=\eta_0\)，并写成等价形式

\[
\gamma_{34}
=
\frac{\Gamma_3^2+\eta_0^2-1}
{\eta_0\Gamma_3+u_3u_4},
\qquad
u=\sqrt{\Gamma^2-1}.
\]

### 2.2 区域 3 质量、能量和体积

区域 3 的状态量是显式 reservoir：

\[
M_3(R),\qquad U_3(R),\qquad V_3(R),
\qquad
u_3(R)=\frac{U_3}{V_3}.
\]

crossing 前，新扫入抛射物质量为

\[
\mathrm{d}M_3
=
4\pi R^2 m_p
(\beta_4-\beta_{\rm rs})
\eta_0 n_4\,\mathrm{d}R ,
\]

非磁化极限下，新增热能为

\[
\mathrm{d}U_{3,{\rm sh}}
=
(\gamma_{34}-1)\,\mathrm{d}M_3c^2 .
\]

有限 upstream magnetization 时，\((\gamma_{34}-1)\) 被 MHD jump 给出的下游热比内能 \(\epsilon_{\rm th,3}(\gamma_{34},\sigma)\) 替代，具体闭合见下一节。

区域 3 体积变化由新扫入体积和膨胀组成：

\[
\mathrm{d}V_3
=
\frac{\mathrm{d}M_3}{n_3m_p}
+
V_3
\left(
3\frac{\mathrm{d}R}{R}
-
\frac{\mathrm{d}\Gamma_3}{\Gamma_3}
\right).
\]

绝热项为

\[
\mathrm{d}U_{3,{\rm ad}}
=
-(\hat{\gamma}_3-1)U_3
\frac{\mathrm{d}V_{3,{\rm exp}}}{V_3},
\qquad
\hat{\gamma}_3
=
\frac{4\gamma_{\rm th,3}+1}
{3\gamma_{\rm th,3}}.
\]

crossing 后不再有新抛射物注入：

\[
\mathrm{d}M_3=0,
\qquad
\mathrm{d}U_{3,{\rm sh}}=0,
\]

但 \(U_3/V_3\) 仍随膨胀演化，不能用固定 crossing 值代替。

### 2.3 磁化反向激波

用户 API 字段是 `ReverseShock.upstream_sigma`；内部配置和 Fortran 参数中也会看到历史名 `reverse_sigma` 或 `sigma_r`。三者描述同一个未激波抛射物上游磁化：

\[
\sigma
=
\frac{B_4^2}{4\pi\Gamma_4^2\rho_4c^2}.
\]

`E_iso` 在 ASGARD 中是总 ejecta 能量；有限 \(\sigma\) 时只有 \(1/(1+\sigma)\) 是 baryonic rest-mass 分量。因此

\[
M_{\rm ej,b}
=
\frac{E_{\rm iso}}
{(1+\sigma)\eta_0c^2}.
\]

这一步很重要：增加 \(\sigma\) 不只是增强磁场，也降低同一总能量下的 baryonic ejecta mass，并改变 region 4 数密度 \(n_4\)。

反向激波是否能形成由两条条件共同决定。第一条是上游相对于 shock 的四速度必须超过 fast magnetosonic 条件：

\[
u_{4s}^2>\sigma .
\]

第二条是 region 2 的 shocked external pressure 必须压过 ejecta 有序磁压。代码使用的归一化压力为

\[
p_2^{\rm norm}
=
\frac{(4\Gamma_{\rm cd}+3)(\Gamma_{\rm cd}-1)n_1}{3},
\qquad
p_{B,4}^{\rm norm}
=
\frac12\sigma n_4 ,
\]

所以触发条件是

\[
\frac{p_2^{\rm norm}}{p_{B,4}^{\rm norm}}\ge1.
\]

等价的临界磁化可写为

\[
\sigma_{\rm cr}
=
\frac{2(4\Gamma_{\rm cd}+3)(\Gamma_{\rm cd}-1)n_1}{3n_4}.
\]

通过条件后，压缩比不再使用 ultra-relativistic \(4\gamma_{34}+3\) 近似。当前实现使用有限强度 MHD jump 的下游四速度 \(u_{3s}\) 和上游 shock-frame 四速度 \(u_{4s}\)：

\[
C(\gamma_{34},\sigma)
=
\frac{u_{4s}}{u_{3s}} .
\]

\(\sigma\to0\) 时该 jump 回到当前有限强度 hydrodynamic baseline；这比把 \(4\gamma_{34}+3\) 当成所有强度下的基线更严格。下游热比内能由同一个 MHD jump 给出。设

\[
\hat\gamma
=
\frac{4}{3}+\frac{1}{3\gamma_{34}},
\qquad
\gamma_{is}=\sqrt{1+u_{is}^2},
\]

\[
h_3
=
\frac{(1+\sigma)\gamma_{4s}}{\gamma_{3s}}
-C\sigma ,
\]

则代码推进的 shock-heating 源项使用

\[
\epsilon_{\rm th,3}
=
\frac{h_3-1}{\hat\gamma},
\qquad
\mathrm{d}U_{3,{\rm sh}}
=
\epsilon_{\rm th,3}\,\mathrm{d}M_3c^2 .
\]

\(\sigma=0\) 时 \(\epsilon_{\rm th,3}=\gamma_{34}-1\)。

上游有序磁场定义为

\[
B_{4,{\rm ord}}
=
\sqrt{4\pi c^2\sigma\rho_4},
\qquad
B_{3,{\rm ord}}
=
C(\gamma_{34},\sigma)B_{4,{\rm ord}}.
\]

总区域 3 磁场为 turbulent 分量和 ordered 分量的平方和：

\[
B_3
=
\sqrt{
8\pi\epsilon_{B,3}\frac{U_3}{V_3}
+
B_{3,{\rm ord}}^2
}.
\]

因此

\[
\lim_{\sigma\to0}B_3
=
\sqrt{8\pi\epsilon_{B,3}U_3/V_3}.
\]

穿越后不再有新的 region 4 物质进入。代码在 crossing 处保存 \(B_{3,{\rm ord,cross}}\)、\(R_{\rm cross}\) 和 \(V_{3,{\rm cross}}\)，随后按

\[
B_{3,{\rm ord}}(R)
=
B_{3,{\rm ord,cross}}
\frac{V_{3,{\rm cross}}}{V_3(R)}
\frac{R}{R_{\rm cross}}
\]

推进有序场。这个项还进入动力学惯性。若

\[
p_{B,3}=\frac{B_{3,{\rm ord}}^2}{8\pi},
\qquad
E_{B,3}=p_{B,3}V_3,
\]

则有序场给 bulk 方程的有效惯性质量为

\[
M_{B,{\rm eff}}
=
\frac{E_{B,3}+p_{B,3}V_3}{c^2}
+
\frac{p_{B,3}V_3}{\Gamma_3^2c^2}.
\]

因此当前磁化 RS 的完整链条是：

1. `upstream_sigma` 改变 \(M_{\rm ej,b}\) 和 \(n_4\)。
2. fast-wave 与 pressure-balance 条件决定 RS 何时真正形成。
3. MHD jump 同时给出 \(C\) 和 \(\epsilon_{\rm th,3}\)。
4. \(C\) 压缩上游 ordered field，\(\epsilon_{\rm th,3}\) 决定 \(\gamma_{m,3}\) 和区域 3 thermal state。
5. crossing 后 \(U_3/V_3\)、\(V_3\) 和 \(B_{3,{\rm ord}}\) 继续演化，并通过磁压焓进入 inertia。

VegasAfterglow 在这里只作为 jump-condition comparison backend，不是光变目标或全局动力学基准。

磁化反向激波动力学在每个平滑分支内由 `advance_reverse_logtime` 系列 RK4 步推进，右端项平滑且没有事件端点时全局为四阶：

\[
Y(R)-Y_{\Delta R}(R)
=
{\cal O}(\Delta R^4).
\]

waiting-to-shock、pressure-balance 触发和 crossing 是物理分支切换，代码用端点捕获而不是跨步后修补。事件附近不对单个跨事件步声明四阶；验收口径是 \(M_3,U_3,V_3,\gamma_{34},B_{3,{\rm ord}}\) 和光变分量在物理事件两侧连续。耦合到 `dg_1d` 电子输运后，端到端电子谱普通注入/冷却段仍受后向 Euler 限制为一阶半径收敛。

### 2.4 反向激波电子和辐射

电子注入谱为

\[
Q_{e,3}(\gamma)
=
Q_{0,3}(\gamma-1)^{-p_3}
H(\gamma-\gamma_{m,3})
\exp\!\left(-\frac{\gamma}{\gamma_{M,3}}\right),
\]

数和能量归一化满足

\[
\int Q_{e,3}\,\mathrm{d}\gamma
=
\xi_{N,3}\frac{\mathrm{d}M_3}{m_p},
\]

\[
m_ec^2
\int(\gamma-1)Q_{e,3}(\gamma)\,\mathrm{d}\gamma
=
\epsilon_{e,3}\,\mathrm{d}U_{3,{\rm sh}}.
\]

在 \(p_3>2\) 且高能截断足够远时，

\[
\gamma_{m,3}
\simeq
1+
\frac{\epsilon_{e,3}}{\xi_{N,3}}
\frac{p_3-2}{p_3-1}
\frac{m_p}{m_e}
\epsilon_{\rm th,3}.
\]

非磁化极限下 \(\epsilon_{\rm th,3}=\gamma_{34}-1\)。

反向激波同步辐射、RS SSC 和 FS/RS 跨区逆康普顿都使用同一套区域 3 电子谱、\(B_3\) 和 seed photon field。RS 强子 light path 只表示壳层级质子注入/输运和质子同步辐射；full-chain RS path 复用 formal 1D 强子核，不是 \(\chi\) 分辨 RS 强子输运。

## 3. 密度增强触发的次级反向激波

### 3.1 密度增强数组

多密度增强使用数组合同：

```text
jump_r_cm
jump_factor
jump_width
```

外介质密度写成

\[
n(R)
=
n_{\rm base}(R)
\left[
1+
\sum_j
(f_j-1)
\exp\left(
-
\frac{(R-R_j)^2}{2(w_j R_j)^2}
\right)
\right].
\]

当前次级反向激波分支只使用每个密度增强的上升段作为局部 source：

\[
n_{{\rm exc},j}(R)>0
\quad
{\rm for}
\quad
-4w_jR_j\le R-R_j<0.
\]

这使每个 branch 有明确的能量来源，避免在下降段重复注入。

### 3.2 热上游四区 Riemann 问题

密度增强触发的次级反向激波不是抛射物冷上游问题。它的区域 4 是已经被正向激波加热的旧 shocked shell，必须使用

\[
n_4,\qquad e_4,\qquad p_4,\qquad \Gamma_4
\]

作为热上游。对第 \(j\) 个增强，若没有父 branch，则当前实现取

\[
n_4=4\Gamma_4 n_{\rm pre},
\qquad
e_4=4\Gamma_4(\Gamma_4-1)n_{\rm pre}m_pc^2,
\qquad
p_4=\frac{e_4}{3}.
\]

若前一个增强的 branch 已存在，则第 \(j\) 个增强可使用父 branch 的 reservoir：

\[
n_4=\frac{M_{3,j-1}}{m_pV_{3,j-1}},
\qquad
e_4=\frac{U_{3,j-1}}{V_{3,j-1}},
\qquad
p_4=\frac{e_4}{3}.
\]

接触间断 Lorentz 因子 \(\Gamma_c\) 由压力平衡确定：

\[
p_2(\Gamma_c)=p_3(\Gamma_c).
\]

透射正向激波侧压力为

\[
p_2(\Gamma_c)
=
\frac{4}{3}
(\Gamma_c^2-1)n_1m_pc^2.
\]

次级反向激波侧相对 Lorentz 因子为

\[
\gamma_{43}
=
\Gamma_c\Gamma_4(1-\beta_c\beta_4),
\]

压力为

\[
p_3(\Gamma_c)
=
p_4+
\frac{4}{3}
(\gamma_{43}^2-1)(e_4+p_4).
\]

数值上用二分法在 \(1<\Gamma_c<\Gamma_4\) 中求根：

\[
f(\Gamma_c)=p_2(\Gamma_c)-p_3(\Gamma_c)=0.
\]

shock-frame 压缩比 \(C=n_3/n_4\) 和次级反向激波速度 \(\beta_s\) 来自动量通量平衡：

\[
C=\frac{u_{4s}}{u_{3s}},
\]

\[
w_4u_{4s}^2+p_4
=
w_3u_{3s}^2+p_3,
\]

\[
w_4=n_4m_pc^2+e_4+p_4,
\qquad
w_3=Cn_4m_pc^2+e_3+p_3,
\qquad
e_3=3p_3 .
\]

### 3.3 新耗散能和电子注入

次级反向激波只允许把新增耗散能注入新电子。先扣除热上游的绝热压缩贡献：

\[
e_{\rm ad}
=
e_4 C^{4/3},
\qquad
u_{{\rm diss},3}
=
e_3-e_{\rm ad}.
\]

只有 \(u_{{\rm diss},3}>0\) 时 branch 才有新注入。区域 3 数密度为

\[
n_3=Cn_4.
\]

最小电子 Lorentz 因子为

\[
\gamma_{m,3}^{\rm sec}
=
1+
\frac{\epsilon_{e,3}}{\xi_{N,3}}
\frac{p_3-2}{p_3-1}
\frac{u_{{\rm diss},3}}
{n_3m_ec^2}.
\]

branch 的扫入质量、体积和热能增量为

\[
\frac{\mathrm{d}M_{3,j}}{\mathrm{d}R}
=
4\pi R^2m_p n_4\Gamma_4
\frac{\beta_2-\beta_s}{\beta_c},
\]

\[
\mathrm{d}V_{3,j}
=
\frac{\mathrm{d}M_{3,j}}{n_3m_p}
+
V_{3,j}
\left(
3\frac{\mathrm{d}R}{R}
-
\frac{\mathrm{d}\Gamma_c}{\Gamma_c}
\right),
\]

\[
\mathrm{d}U_{3,j}
=
u_{{\rm diss},3}\,\mathrm{d}V_{{\rm shock},j}
-
\frac{1}{3}U_{3,j}
\frac{\mathrm{d}V_{{\rm exp},j}}{V_{3,j}}.
\]

磁场为

\[
B_{3,j}
=
\sqrt{8\pi\epsilon_{B,3}e_3}.
\]

输出的 `rev.sync` 是抛射物反向激波同步辐射和所有 active 次级反向激波 branch 同步辐射之和：

\[
L_{\nu}^{\rm rev}
=
L_{\nu}^{\rm ejecta\,RS}
+
\sum_j L_{\nu,j}^{\rm secondary\,RS}.
\]

### 3.4 验收量

次级反向激波不是拟合残差的经验 bump。每个 branch 必须检查

\[
\Gamma_c(R),\quad p_3(R),\quad \gamma_{43}(R),
\quad u_{{\rm diss},3}(R),\quad B_{3,j}(R)
\]

随 \(R\) 连续平滑。能量预算应满足

\[
E_{e,j}^{\rm inj}
=
\epsilon_{e,3}
\int u_{{\rm diss},3}\,\mathrm{d}V_{{\rm shock},j}.
\]

若 FS 光变出现无物理来源的尖峰，优先检查密度增强宽度、event root、branch 权重和 reservoir 演化，不在观测投影层做 smoothing。

## 4. 自适应网格算法

ASGARD 中“自适应网格”的物理含义不同，不能互相替代。

### 4.1 电子辐射积分的曲率重采样

旧的 `adaptive_resampling_log` 曲率重采样没有生产调用者，已从 Electron 构建闭包移除。当前同步辐射积分直接使用正式电子谱网格；若以后需要重新引入辐射采样压缩，必须先说明误差预算、能量/粒子矩保持方式和触达入口，而不是作为后处理降采样。

### 4.2 自适应子步误差

电子输运可用 full step 与 two half steps 比较误差。设

\[
\mathbf{N}_1
=
\Phi_{\Delta R}(\mathbf{N}^n),
\qquad
\mathbf{N}_{1/2,1/2}
=
\Phi_{\Delta R/2}
\left[
\Phi_{\Delta R/2}(\mathbf{N}^n)
\right],
\]

误差估计为

\[
\epsilon
=
\frac{
\left\|\mathbf{N}_{1/2,1/2}-\mathbf{N}_1\right\|_w
}{
\left\|\mathbf{N}_{1/2,1/2}\right\|_w
}.
\]

若 \(\epsilon\) 超过目标阈值，减小 \(\Delta R\) 后重算；若 \(\epsilon\) 很小，可增大下一步。这个误差控制服务输运方程本身，不允许作为后处理去修光变尖峰。

### 4.3 自适应观测时间网格

`flux_density_grid_adaptive` 的观测时间网格不改变求解状态，只改变查询时间集合。先保留用户输入时间

\[
\mathcal{T}_{\rm user}=\{t_m\},
\]

再加入 baseline log grid

\[
\mathcal{T}_{\rm base}
=
\left\{
10^{\log t_{\min}+q\Delta\log t}
\right\},
\]

并从每条发射轨道加入等到达时间 knot：

\[
t_{{\rm knot},ij}
=
t_{{\rm axis},i}
+
(1+z)
\frac{R_i(1-\mu_j)}{c}.
\]

合并、排序、去重后，再加入相邻点的 log midpoint：

\[
t_{{\rm mid},q}
=
\sqrt{t_qt_{q+1}}.
\]

最终查询集合为

\[
\mathcal{T}_{\rm adaptive}
=
\operatorname{unique}
\left(
\mathcal{T}_{\rm user}
\cup
\mathcal{T}_{\rm base}
\cup
\mathcal{T}_{\rm knot}
\cup
\mathcal{T}_{\rm mid}
\right).
\]

该算法的目标是解析 EATS 和 RS/次级 RS crossing 的时间结构；返回给用户的原始时间轴语义不变。

### 4.4 曝光时间平均

曝光平均使用 Gauss-Legendre 节点。对第 \(m\) 个观测点，曝光区间为

\[
[t_m-\Delta t_m/2,\;t_m+\Delta t_m/2].
\]

节点和权重由

\[
t_{mq}
=
\frac{t_+ + t_-}{2}
+
\frac{t_+ - t_-}{2}\xi_q,
\qquad
w_{mq}
=
\frac{t_+ - t_-}{2}\omega_q
\]

给出，平均通量为

\[
\bar{F}_{\nu,m}
=
\frac{1}{\Delta t_m}
\sum_q w_{mq}F_\nu(t_{mq}).
\]

这也是观测查询层算法，不改变动力学、电子输运或反向激波物理状态。

## 5. 不允许的解释

- \(\chi\) 分辨同步投影不等于 \(\chi\)-local 强子输运。
- 次级反向激波不是光变后处理 bump；它必须来自热上游 Riemann 问题和新增耗散能。
- 自适应网格不是 smoothing；若真实物理参数随时间或半径不平滑，应回查源项、事件定位、输运矩阵或投影权重。
- `doppler` 代码变量在插值核中表示 \(D=\delta^{-1}\)，物理 Doppler factor 仍是 \(\delta=[\Gamma(1-\beta\mu)]^{-1}\)。
