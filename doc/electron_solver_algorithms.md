# 电子求解器算法说明

本文说明当前 ASGARD 正向/反向激波电子输运求解器的方程、离散方式和适用角色。相关 Fortran 源码位于 `src/Electron/`。代码标识保持英文，物理解释使用中文。

## 1. 统一变量

1D 求解器在能量对数坐标

\[
x=\log_{10}\gamma_e
\]

上推进电子谱。主守恒变量是

\[
dN_x
=\frac{\mathrm{d}N}{\mathrm{d}\log_{10}\gamma_e}
=\gamma_e\ln 10\frac{\mathrm{d}N}{\mathrm{d}\gamma_e}.
\]

2D 求解器额外展开下游厚壳坐标

\[
\eta=\log_{10}\chi,
\qquad
\chi=1+\frac{8\Gamma_{\rm sh}^2 x_{\rm downstream}}{R},
\]

并推进

\[
U(x,\eta,R)
=\frac{\mathrm{d}N}
{\mathrm{d}\log_{10}\gamma_e\,\mathrm{d}\log_{10}\chi}.
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

变换到 \(x=\log_{10}\gamma_e\) 后，

\[
\frac{\partial dN_x}{\partial R}
+\frac{\partial}{\partial x}
\left(A_xdN_x\right)
=S_x,
\]

其中

\[
A_x=\frac{\gamma'_e}{\gamma_e\ln 10},
\qquad
S_x=\gamma_e\ln 10\,Q.
\]

实际离散中的面速度近似为

\[
A_{\rm face}\simeq dE_{\ell,{\rm mean}}+\frac{1}{R\ln 10},
\]

第二项是球对称膨胀给出的绝热项。

## 4. 2D 连续方程

2D 路径推进 \(U(x,\eta,R)\)：

\[
\frac{\partial U}{\partial R}
+\frac{\partial(A_xU)}{\partial x}
+\frac{\partial(A_\eta U)}{\partial\eta}
=
\frac{\partial}{\partial\eta}
\left(D_\eta\frac{\partial U}{\partial\eta}\right)
+S.
\]

\(x\) 方向是能量冷却平流，\(\eta\) 方向包含下游流体平流、扩散和局域源项。`charint_2d` 中只有 \(\eta=1\) 的 shock-front 列接收新注入源；已经 advect/diffuse 到下游的其它 \(\chi\) 列在能量方向调用公共无源特征线冷却 primitive。

## 5. 求解器角色

| 求解器 | 能量方向 | \(\chi\) 方向 | 当前角色 |
| --- | --- | --- | --- |
| `fullhide_1d` | 一阶隐式迎风 | 无 | 默认稳定基线 |
| `slc1_1d` | 半拉格朗日 + Strang splitting | 无 | 方法比较 |
| `charint_1d` | 特征线保守重映射 | 无 | 1D 高保真路径 |
| `t2g1_1d` | BDF2 三层推进 | 无 | 时间推进比较 |
| `weno5_1d` | WENO5 + SSP RK3 | 无 | 高阶显式比较 |
| `dg_1d` | 多域 LGL 谱元 DG | 无 | FS/RS opt-in 高阶路径 |
| `fullhide_2d` | 一阶隐式迎风 | 隐式平流-扩散 | 2D 物理基线 |
| `charint_2d` | 特征线重映射 | 隐式平流-扩散 | 2D 加速混合版 |

`fullhide_1d` 最稳健，但数值扩散较强。`dg_1d` 使用 moving multi-domain LGL 谱元，在固定 1D 电子网格数下保留更尖锐的冷却断点和高能截断；它是 opt-in 路径，不替代默认 `fullhide_1d`。`charint_1d` 对高能截止和冷却断点也更锋利，但子步选择更敏感。`fullhide_2d` 是当前最完整的厚壳电子路径。`charint_2d` 只把能量方向换成特征线重映射，\(\eta\) 方向仍使用隐式平流-扩散，因为扩散项不能写成单纯特征线更新。

## 6. `dg_1d` 谱元路径

`dg_1d` 的核心空间变量仍是 \(dN_x\)，内部在多个 LGL 谱元上表示多项式状态。FS 与 primary RS 共用同一个 DG transport kernel，但物理源项保持各自定义：

- FS 源项仍使用 \(\gamma^{-p}\)。
- RS kinetic 源项使用 \((\gamma-1)^{-p}\)。当前 primary RS DG 为避免 \(\gamma_m\simeq1\) 附近的 kinetic singular 被节点采样放大，先在正式 \(\log\gamma\) 网格上构造与 fullhide 相同的 cell-averaged kinetic 源项，再插值到 DG 节点并按注入粒子数守恒归一化；不使用节点解析源项作为默认路径。
- SSA/IC/同步冷却仍由 cooling kernel 给出；DG transport 层只接收最终损失率、绝热率和源项。

FS DG 在壳层之间持续保存 DG state，并在 \(\gamma_m,\gamma_c,\gamma_{\max}\) 移动时重投影到新谱元。RS DG 同样跨 shell 持久保存 DG state；只在输出给固定 \(\log\gamma\) 网格和辐射核时使用守恒 cell 投影。RS DG 对断点附近的正谱多项式使用 cell-average preserving positivity limiter，保持元素平均粒子数，不作为后处理 smoothing。

当前时间推进基线不改变：DG 每个 shell 使用 10 个基础子步，密度跳变时 FS DG 仍由 jump-aware limiter 自适应缩步。RS fullhide 保持原 100--1000 冷却子步，RS DG 与 FS DG 对齐为 10 个基础子步。

RS 高能尾的当前验收口径是分辨率收敛而不是强贴合 121 格 fullhide。121 格 fullhide 在 post-crossing 高能截断处有明显一阶隐式上风数值扩散；`num_gam_e=241` 时 fullhide 截断向 DG 收缩。因此 `dg_1d` 的较窄高能尾不视为单独错误，除非高分辨率 fullhide、解析特征线或守恒谱形诊断也显示同向偏差。

Primary RS 在完全 crossing 之后没有新的反向激波注入。`fullhide_1d` 和 `dg_1d` 都不再继续用各自的 shell-to-shell 数值输运推进这段纯冷却演化，而是在 crossing 后缓存固定网格上的 \(dN_x\)，随后只累计当前网格边界回溯到 crossing 网格的 characteristic edge map；输出给辐射核时从 crossing 谱一次做保守重映射。低能端使用闭合边界，冷却到电子网格下界以下的数目保留在最低能单元。这样避免 repeated projection 或一阶隐式上风扩散生成的硬截断、平台、过宽高能尾和 late-time 清零。当前有效支撑图保留在 `output/asgard_doc/dg_radiation_stability_scan/postcross_direct_map_effective_support/`。

## 7. 子步与近似

2D 子步取

\[
\Delta R_{\rm try}
=\min(\Delta R_\xi,\Delta R_\eta,\Delta R_D).
\]

quick/formal 2D 路径使用隐式算子的精度步长：

\[
\Delta R_\xi \propto \frac{4\Delta x}{|A_\xi|},
\qquad
\Delta R_\eta \propto \frac{4\Delta\eta}{|A_\eta|}.
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
| \(\chi/\eta\) 几何 | `electron_transport_2d_kernel.f90` |
| 注入源项 \(Q_e\) | `electron_injection_profiles.f90` |
| 冷却系数 \(A_x\) | `electron_cooling_kernel.f90` |
| FS/RS shared 1D transport wrapper | `electron_shell_transport_common.f90` |
| 1D 隐式迎风 | `electron_transport_common.f90` |
| 1D LGL-DG | `electron_transport_dg_1d_kernel.f90` |
| 1D 特征线 | `electron_transport_common.f90` |
| 2D \(\eta\) 输运 | `electron_transport_2d_kernel.f90` |
| 2D 能量输运 | `electron_transport_2d_kernel.f90` |
| 同步辐射谱 | `electron_radiation_kernel.f90` |
| 历史光子场 | `electron_seed_history_kernel.f90` |

## 9. 验收风险

需要重点检查：

- \(\nu_m\)、\(\nu_c\)、\(\nu_a\) 是否随时间平滑。
- 相邻观测时刻的总谱是否出现锯齿。
- 电子总数是否在无物理注入或逃逸事件时突变。
- 平滑参数扫描中高频端是否出现台阶式截断。
- `dg_1d` 输出给辐射核前的固定网格谱是否出现元素边界零洞或断点 Gibbs 过冲。
- FS density-jump 场景当前仍会触发严格 `dg_1d_smoke.py` 的 3 个 sawtooth turns；这是未完成数值边界，不通过后处理 smoothing 处理。
- `charint_2d` 的 \(\xi\) 子步是否过粗导致高能端提前消失。
- 2D reduced cooling bands 在窄频带问题中是否引入系统偏差。

如果这些诊断不平滑，应回到动力学、冷却源项、网格映射和投影检查，不能在后处理层 smoothing。
