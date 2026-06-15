# 2D 壳层、反向激波与自适应网格

本文集中说明 ASGARD 当前网页文档中最容易混用的三个模块：\(\chi\) 分辨有限厚壳层、抛射物反向激波与密度增强触发的次级反向激波、自适应网格算法。这里写的是当前代码契约；未实现的 \(\chi\)-local 强子输运、\(\chi\)-local SSC 和 IC-mediated electromagnetic cascade 不在本页扩展。

代码标识保持英文，例如 `chi_eats_2d`、`fullhide_2d`、`reverse_sigma`。正文术语按 `doc/terminology.md`：反向激波、次级反向激波、密度增强、壳层级、等到达时间面。

## 1. \(\chi\) 分辨有限厚壳层

### 1.1 物理坐标

正向激波半径为 \(R\)，激波 Lorentz 因子为 \(\Gamma_{\rm sh}\)，当前实现的 ISM/BM 厚壳几何使用

\[
\chi
=
1+
8\Gamma_{\rm sh}^2
\left(1-\frac{r}{R}\right),
\qquad
1\le \chi \le 1+8\Gamma_{\rm sh}^2 .
\]

\(\chi=1\) 是激波前沿，较大的 \(\chi\) 是更深的下游。数值网格使用

\[
\eta=\log_{10}\chi,
\qquad
\eta_{k+1/2}=k\Delta\eta,
\qquad
\chi_{k+1/2}=10^{\eta_{k+1/2}},
\]

\[
\eta_k=\left(k-\frac12\right)\Delta\eta,
\qquad
\chi_k=10^{\eta_k}.
\]

代码中的几何系数为

\[
a(R)=\frac{8\Gamma_{\rm sh}^2}{R},
\qquad
\chi=1+a(R)(R-r).
\]

激波系下游距离和共动距离为

\[
x_{{\rm sh},k+1/2}
=
\frac{(\chi_{k+1/2}-1)R}{8\Gamma_{\rm sh}^2},
\]

\[
\Delta x'_k
=
\gamma_{{\rm rel},k}
\left(x_{{\rm sh},k+1/2}-x_{{\rm sh},k-1/2}\right),
\qquad
\gamma_{{\rm rel},k}
=
\left(1-\beta_{2,{\rm sh}}^2(\chi_k)\right)^{-1/2}.
\]

BM 下游速度近似为

\[
\Gamma_2(\chi)
=
\frac{\Gamma_{\rm sh}}{\sqrt{2\chi}},
\qquad
\beta_2(\chi)
=
\sqrt{1-\Gamma_2^{-2}(\chi)} ,
\]

\[
\beta_{2,{\rm sh}}(\chi)
=
\frac{\beta_{\rm sh}-\beta_2(\chi)}
{1-\beta_{\rm sh}\beta_2(\chi)} .
\]

### 1.2 电子输运方程

2D 电子状态使用

\[
U(x,\eta,R)
=
\frac{\mathrm{d}N_e}
{\mathrm{d}x\,\mathrm{d}\eta},
\qquad
x=\log_{10}\gamma_e .
\]

连续方程写成

\[
\frac{\partial U}{\partial R}
+
\frac{\partial(A_xU)}{\partial x}
+
\frac{\partial(A_\eta U)}{\partial \eta}
=
\frac{\partial}{\partial\eta}
\left(
D_\eta
\frac{\partial U}{\partial\eta}
\right)
+
S(x,\eta,R).
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

\(\eta\) 方向输运系数与代码 `eta_face_transport_coeff` 对应：

\[
A_\eta(\chi,R)
=
\frac{1}{\ln 10}
\left[
\frac{a(R)\beta_{2,{\rm sh}}(\chi)}
{\chi\beta_{\rm sh}}
+
\frac{\chi-1}{\chi}
\frac{\mathrm{d}\ln a}{\mathrm{d}R}
\right].
\]

局域绝热冷却系数使用 BM 下游速度场散度。当前实现的有效形式是

\[
\left(\nabla\cdot \mathbf{v}\right)_\chi
=
c
\left[
\frac{2\beta_2(\chi)}
{R_\chi/R}
+
8\frac{1}{\beta_2(\chi)}
\right]\frac{1}{R},
\qquad
R_\chi=R\left[1-\frac{\chi-1}{8\Gamma_{\rm sh}^2}\right],
\]

\[
\left(\frac{\mathrm{d}x}{\mathrm{d}R}\right)_{\rm ad}
=
\frac{1}{3\beta_{\rm sh}c\ln 10}
\left(\nabla\cdot\mathbf{v}\right)_\chi .
\]

### 1.3 离散推进

`fullhide_2d` 对 \(\eta\) 方向使用隐式有限体积形式。设 cell 平均量为 \(U_k^n\)，步长为 \(\Delta R\)，通量为

\[
F_{k+1/2}
=
A^+_{k+1/2}U_k^{n+1}
+
A^-_{k+1/2}U_{k+1}^{n+1}
-
D_{k+1/2}
\frac{U_{k+1}^{n+1}-U_k^{n+1}}{\Delta\eta}.
\]

则每个能量 bin 上求解三对角系统

\[
U_k^{n+1}
+
\frac{\Delta R}{\Delta\eta}
\left(
F_{k+1/2}^{n+1}
-
F_{k-1/2}^{n+1}
\right)
=
U_k^n+\Delta R\,S_{k,{\rm shock}},
\]

其中 \(S_{k,{\rm shock}}\) 只注入到 shock-front 侧的 \(\eta\) cell。`charint_2d` 在 \(\eta\) 对流部分用特征线回溯：

\[
\frac{\mathrm{d}\eta}{\mathrm{d}R}=A_\eta(\eta,R),
\qquad
\eta_{\rm back}
=
\eta_{\rm face}(R_{n+1})
-
\int_{R_n}^{R_{n+1}}A_\eta\,\mathrm{d}R,
\]

然后用守恒 PPM 前缀积分重映射：

\[
U_k^{n+1}
=
\frac{
\int_{\eta_{{\rm back},k-1/2}}^{\eta_{{\rm back},k+1/2}}
U^n(\eta)\,\mathrm{d}\eta
}{\Delta\eta}.
\]

2D 子步限制来自 \(\eta\) 方向 CFL：

\[
\Delta R_\eta
=
C_{\rm CFL}
\frac{\Delta\eta}
{\max_k |A_{\eta,k+1/2}|}.
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

新增热能为

\[
\mathrm{d}U_{3,{\rm sh}}
=
(\gamma_{34}-1)\,\mathrm{d}M_3c^2 .
\]

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

非磁化压缩比固定为 ASGARD 基线：

\[
C_0=4\gamma_{34}+3.
\]

若 `reverse_sigma=\sigma>0`，代码用 VegasAfterglow 的 \(\sigma\) 依赖相对修正，但保持 \(\sigma\to0\) 极限：

\[
C(\gamma_{34},\sigma)
=
C_0(\gamma_{34})
\frac{C_{\rm Vegas}(\gamma_{34},\sigma)}
{C_{\rm Vegas}(\gamma_{34},0)}.
\]

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

总区域 3 磁场为

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

### 2.4 反向激波电子和辐射

电子注入谱为

\[
Q_{e,3}(\gamma)
=
Q_{0,3}\gamma^{-p_3}
H(\gamma-\gamma_{m,3})
H(\gamma_{M,3}-\gamma),
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
(\gamma_{34}-1).
\]

反向激波同步辐射、RS SSC 和 FS/RS 跨区逆康普顿都使用同一套区域 3 电子谱、\(B_3\) 和 seed photon field。RS 强子 light path 只表示壳层级质子注入/输运和质子同步辐射；full-chain RS path 复用 formal 1D 强子核，不是 \(\chi\) 分辨 RS 强子输运。

## 3. 密度增强触发的次级反向激波

### 3.1 密度增强数组

多密度增强使用数组合同：

```text
jump_r_cm
jump_factor
jump_width_log10
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
\frac{(\log_{10}R-\log_{10}R_j)^2}{2w_j^2}
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

ASGARD 中“自适应网格”有三类，物理含义不同，不能互相替代。

### 4.1 电子辐射积分的曲率重采样

`adaptive_resampling_log` 用于压缩活跃电子谱辐射积分网格。给定 \(x_i=\log_{10}\gamma_i\) 和 \(f_i=\log_{10}N_i\)，先估计

\[
f'_i
\simeq
\frac{f_{i+1}-f_{i-1}}
{x_{i+1}-x_{i-1}},
\]

\[
f''_i
\simeq
\frac{f'_{i+1}-f'_{i-1}}
{x_{i+1}-x_{i-1}}.
\]

曲率权重为

\[
\kappa_i=|f''_i|,
\]

并按 cell 宽度修正和平滑：

\[
\tilde{\kappa}_i
=
\mathcal{M}_5
\left[
\kappa_i\Delta x_i
\right],
\]

其中 \(\mathcal{M}_5\) 是五点滑动平均。为避免只保留强曲率点而丢失宽区间背景，代码设置最小权重

\[
w_i
=
\max\left(
\tilde{\kappa}_i,\,
\frac{\sum_l\tilde{\kappa}_l}{m g}
\right),
\]

其中 \(m\) 是目标点数，\(g\) 是平滑/保底参数。累积权重为

\[
C_i=\sum_{q\le i}w_q\Delta x_q.
\]

目标点从等间隔累计权重中选出：

\[
C_{i_k}
\approx
\frac{k-1}{m-1}C_N,
\qquad
k=1,\ldots,m.
\]

重复索引会被去重，不足的点从缺失索引中均匀补足，然后排序。这个算法只改变辐射积分采样点，不改变电子输运状态本身。

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
