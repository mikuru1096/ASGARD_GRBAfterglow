# 磁化反向激波与 DG1D 教程

本文面向已经跑通 `quickstart.md` 的用户，说明如何启用修复后的磁化反向激波求解器和 `dg_1d` 高阶电子输运。它不是新手最小路径；默认公开基线仍是 `electron_solver="fullhide_1d"`，`dg_1d` 是需要显式启用的高阶路径。

## 1. 最小开关

完整 `Model(...)` 构造器见 `quickstart.md`。启用本专题路径时只需要改两个对象：

```python
solver_options = SolverOptions(
    electron_solver="dg_1d",
    dynamics_solver="forward_legacy",
    geometry_projection="sed_legacy",
    electron_photon_coupling="separated",
    ssc_cooling_mode="nakar_y_thomson",
    synchrotron_integration="fixed_grid",
    cooling_kernel="legacy",
    radiation_kernel="legacy",
    structured_backend="fortran_1d",
    patch_sampling="uniform",
    patch_projection="auto",
    patch_sampling_pilot_theta=0,
    patch_sampling_num_times=12,
    patch_sampling_beaming_factor=3.0,
    patch_sampling_beaming_resolution=8.0,
    structured_parallel_mode="outer",
    structured_outer_threads=None,
    structured_inner_threads=None,
    fullhide2d_transport_model="legacy",
    fullhide2d_stochastic_accel_norm=0.0,
    fullhide2d_escape_mode="closed",
)

reverse_shock = ReverseShock(
    enabled=True,
    shell_duration_s=30.0,
    upstream_sigma=0.1,
    include_cross_zone_ic=False,
    include_ssc=True,
)
```

需要对照 `fullhide_1d` 时，只把 `electron_solver` 改回 `"fullhide_1d"`，其余物理参数不变。`ReverseShock.upstream_sigma` 是公开 API 字段；内部配置和 Fortran 中的 `reverse_sigma`、`sigma_r` 指同一物理量。

受影响扩展的构建命令是：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --module electron_forward_dg_1d --module electron_reverse_kernel --module structured_jet_1d --force'
```

## 2. 磁化反向激波的物理闭合

反向激波使用四区图像。区域 3 是已激波加热的抛射物，区域 4 是未激波抛射物。反向激波电子注入能标来自激波前沿相对 Lorentz 因子：

\[
\gamma_{34}
=
\Gamma_3\Gamma_4(1-\beta_3\beta_4)
=
\frac{\Gamma_3^2+\eta_0^2-1}
{\eta_0\Gamma_3+u_3u_4}.
\]

这里 \(\Gamma_4=\eta_0\)，\(u_i=\sqrt{\Gamma_i^2-1}\)。不要用区域平均 \(\Gamma\) 经验替代 \(\gamma_{34}\)。

上游磁化定义为

\[
\sigma
=
\frac{B_4^2}{4\pi\Gamma_4^2\rho_4c^2}.
\]

ASGARD 中 `E_iso` 是总抛射物能量；有限 \(\sigma\) 时重子静质量是

\[
M_{\rm ej,b}
=
\frac{E_{\rm iso}}{(1+\sigma)\eta_0c^2}.
\]

所以增大 `upstream_sigma` 会同时降低 \(n_4\)，而不是只给辐射端增加一个有序磁场。

磁化反向激波真正形成前要通过两个条件。第一，反向波必须是激波分支：

\[
u_{4s}^2>\sigma .
\]

第二，已激波外介质压力要压过抛射物磁压：

\[
p_2^{\rm norm}
=
\frac{(4\Gamma_{\rm cd}+3)(\Gamma_{\rm cd}-1)n_1}{3},
\qquad
p_{B,4}^{\rm norm}
=
\frac12\sigma n_4,
\]

\[
\frac{p_2^{\rm norm}}{p_{B,4}^{\rm norm}}\ge1.
\]

等价地，可定义

\[
\sigma_{\rm cr}
=
\frac{2(4\Gamma_{\rm cd}+3)(\Gamma_{\rm cd}-1)n_1}{3n_4}.
\]

通过触发条件后，压缩比使用有限强度 MHD jump：

\[
C(\gamma_{34},\sigma)
=
\frac{u_{4s}}{u_{3s}}.
\]

这条式子替代旧文档里的 ultra-relativistic \(4\gamma_{34}+3\) 近似；\(\sigma\to0\) 的极限回到当前有限强度 hydrodynamic jump。

同一个 MHD jump 给出下游热比内能。设

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
\frac{(1+\sigma)\gamma_{4s}}{\gamma_{3s}}-C\sigma .
\]

则

\[
\epsilon_{\rm th,3}
=
\frac{h_3-1}{\hat\gamma},
\qquad
\mathrm{d}U_{3,{\rm sh}}
=
\epsilon_{\rm th,3}\mathrm{d}M_3c^2.
\]

\(\sigma=0\) 时 \(\epsilon_{\rm th,3}=\gamma_{34}-1\)。反向激波最小电子 Lorentz factor 因而使用

\[
\gamma_{m,3}
=
1+
\frac{\epsilon_{e,3}}{\xi_{N,3}}
\frac{p_3-2}{p_3-1}
\frac{m_p}{m_e}
\epsilon_{\rm th,3}.
\]

反向激波注入源项使用动能幂律：

\[
Q_{e,3}(\gamma)
\propto
(\gamma-1)^{-p_3}
\exp\!\left(-\frac{\gamma}{\gamma_{M,3}}\right),
\]

而 FS 注入仍使用 \(\gamma^{-p}\)。

有序场由上游场经 MHD 压缩给出：

\[
B_{4,{\rm ord}}
=
\sqrt{4\pi c^2\sigma\rho_4},
\qquad
B_{3,{\rm ord}}
=
C B_{4,{\rm ord}}.
\]

总磁场是

\[
B_3
=
\sqrt{
8\pi\epsilon_{B,3}\frac{U_3}{V_3}
+B_{3,{\rm ord}}^2
}.
\]

穿越后不再有新的区域 4 物质进入，代码保存 \(B_{3,{\rm ord,cross}}\)、\(R_{\rm cross}\)、\(V_{3,{\rm cross}}\)，并用

\[
B_{3,{\rm ord}}(R)
=
B_{3,{\rm ord,cross}}
\frac{V_{3,{\rm cross}}}{V_3(R)}
\frac{R}{R_{\rm cross}}
\]

推进有序场。有序场还给 bulk 方程贡献磁压焓惯性：

\[
M_{B,{\rm eff}}
=
\frac{E_{B,3}+p_{B,3}V_3}{c^2}
+
\frac{p_{B,3}V_3}{\Gamma_3^2c^2},
\qquad
p_{B,3}=\frac{B_{3,{\rm ord}}^2}{8\pi}.
\]

因此 `upstream_sigma -> 0` 的验收不是只看 \(B_3\)，而是同时检查质量、触发条件、压缩比、热状态、有序场和惯性项。

## 3. DG1D 高阶输运机制

`dg_1d` 对外仍输出固定电子网格上的 \(dN/\mathrm{d}\log\gamma_e\)，但内部使用 log-four-velocity 坐标

\[
y=\log_{10}\left[1+\frac{\gamma_e^2-1}{\gamma_*^2-1}\right],
\qquad \gamma_*=2.
\]

守恒变量是

\[
dN_y=dN_{x_\gamma}\frac{\mathrm{d}x_\gamma}{\mathrm{d}y},
\qquad
x_\gamma=\log_{10}\gamma_e.
\]

对 log-four-velocity 坐标，

\[
\frac{\mathrm{d}x_\gamma}{\mathrm{d}y}
=
\frac{(\gamma_*^2-1)10^y}{2\gamma_e^2}.
\]

输运方程写成

\[
\frac{\partial dN_y}{\partial R}
+\frac{\partial}{\partial y}(a_y dN_y)
=S_y.
\]

DG 在每个谱元上使用 P12 LGL 节点：

\[
dN_y^k(r)
=
\sum_{i=0}^{12}U_i^k\ell_i(r).
\]

弱式离散为

\[
\frac{\mathrm{d}U_i^k}{\mathrm{d}R}
=
\frac{2}{\Delta y^k w_i}
\sum_{j=0}^{12}w_jD_{ji}a_jU_j^k
+{\cal F}_i^k
+S_i^k .
\]

冷却主导时 \(a_y<0\)，信息从高能流向低能；端面通量按迎风方向取右侧状态。若局部回热使速度反号，代码改用稠密输运矩阵并按端面速度符号选择迎风态。时间推进为后向 Euler：

\[
(I-\Delta R K)U^{n+1}=U^n+\Delta R S^{n+1}.
\]

当前基线为 P12、每个物理段 6 个谱元、每个 shell 10 个基础子步；密度跳变时 FS DG 由 jump-aware limiter 缩小步长。RS DG 子步数与 FS 保持一致。

## 4. 当前收敛阶

当前 DG1D 的阶数不是一个单一数字。能量空间和半径方向需要分开判断。

在光滑且未触发正性核的谱元内，P12 LGL-DG 的空间误差满足

\[
\|dN_y-dN_{y,h}\|_{L^2}
=
{\cal O}\!\left((\Delta y)^{13}\right).
\]

这是 smooth-cell 空间阶数，来自 12 次 Legendre/Lagrange 谱元展开。它要求 \(\gamma_m\)、\(\gamma_c\)、高能 cutoff 或 kinetic bump 不在谱元内部造成不可解析断点；当前移动物理段边界就是为了让这些断点尽量落在谱元边界。

半径/时间方向当前仍是后向 Euler，所以普通注入和冷却推进的全局时间阶数是

\[
{\cal O}(\Delta R).
\]

每个 shell 的 10 个 DG 基础子步降低误差幅度，但不改变阶数。密度跳变自适应缩步的作用是解析跳变区宽度和变化率，不是把时间积分提升到二阶或四阶。

磁化 RS 动力学本身在平滑分支上由 RK4 推进，光滑段全局四阶；crossing、waiting-to-shock 和 pressure-balance 触发是物理事件分裂，事件附近用端点捕获保持分支状态连续。耦合到电子输运和光变后，当前端到端 DG 电子结果的限制阶通常仍来自 BE 电子推进，即一阶半径收敛。

troubled 正性核触发的谱元不声明 13 阶点值收敛。该区域的验收目标是保守平均、非负、支撑连续和辐射平滑；断点两侧未滤波的光滑谱元仍保留高阶空间收敛。

主反向激波穿越后纯冷却段使用解析特征线映射。若冷却系数在累计区间内固定，该段对冷却方程本身没有 BE 时间截断误差；实际误差来自穿越时刻电子谱重映射和冷却系数随壳层演化的取样。

## 5. 问题单元正性核

DG 高阶多项式在注入断点和冷却断点附近会出现 Gibbs 振荡。当前默认不是强行截断输出，而是在输运层做局部模态核。

每个谱元先投影到 Legendre 模态：

\[
dN_y^k(r)
=
\sum_{m=0}^{12}\hat U_m^kP_m(r).
\]

若节点值有负值，或最高 6 个模态占比

\[
\frac{\sum_{m=m_{\rm hi}}^{12}|\hat U_m^k|}
{\sum_{m=0}^{12}|\hat U_m^k|}
>2\times10^{-2}
\]

超过阈值，该谱元标记为问题单元。默认模式只滤波该谱元和相邻谱元，并保持 0 阶模态：

\[
\hat U_0^{\rm new}=\hat U_0,
\qquad
\hat U_m^{\rm new}=\sigma_m^J\hat U_m,\quad m\ge1.
\]

Jackson 因子为

\[
\sigma_m^J
=
\frac{(N-m+2)\cos\theta_m+\sin\theta_m/\tan\theta_1}{N+2},
\qquad
\theta_m=\frac{\pi m}{N+2}.
\]

环境变量控制诊断模式：

| 环境变量 | 行为 |
| --- | --- |
| 未设置或 `troubled` | 默认局部 troubled Jackson |
| `0`, `off`, `false`, `none` | 关闭 kernel |
| `1`, `on`, `true`, `jackson` | 全域 Jackson |
| `fejer` | 全域 Fejer |

这个核保留单元平均值，因此不是输出端平滑。验收重点是非负、活动支撑连续、无元素边界零洞、无多重网格尺度锯齿振荡、最终光变平滑。

## 6. 穿越后纯冷却

主反向激波完全穿越后没有新注入。`fullhide_1d` 和 `dg_1d` 都使用解析特征线更新，避免反复数值输运造成晚期平台、硬截断或过宽高能尾。

若冷却律为

\[
\frac{\mathrm{d}\gamma}{\mathrm{d}\tau}
=
-a\gamma^2-b\gamma,
\]

则

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

代码把当前单元边界回溯到穿越时刻的单元边界，再从穿越时刻的电子谱做一次守恒重映射。低能端为闭合边界，冷却到网格下界以下的电子数保留在最低能单元。

## 7. 结果检查

推荐先跑三组对照：

1. `electron_solver="fullhide_1d"`, `upstream_sigma=0`。
2. `electron_solver="dg_1d"`, `upstream_sigma=0`。
3. `electron_solver="dg_1d"`, `upstream_sigma>0`。

需要检查：

- \(\Gamma(R)\)、\(\gamma_{34}(R)\)、\(U_3/V_3\)、\(B_{3,{\rm ord}}\) 在穿越前后连续。
- `upstream_sigma -> 0` 回到非磁化基线。
- DG 电子谱活动支撑连续；低能峰附近的小振荡若不影响辐射且无零洞，可以接受。
- 多波段光变从 \(0.1\,{\rm s}\) 或科学问题要求的最早时间开始检查，Y 轴不要压缩到看不见 RS 峰。
- 结构化喷流 `fortran_1d` 后端可以使用 `fullhide_1d` 或 `dg_1d` 的同步路径；`dg_1d` 不支持热电子分支。

对应冒烟测试和基准测试入口见 `validation_and_benchmarks.md`。正式刷新基准测试前记录 git HEAD、命令、构建状态、输出路径和物理验收口径。
