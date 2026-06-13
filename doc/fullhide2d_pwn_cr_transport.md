# `fullhide_2d` 的 PWN/CR 输运扩展

本文记录可选路径 `fullhide2d_transport_model="pwn_cr_v1"` 的物理契约。默认路径仍是 `"legacy"`，沿用现有 27 元素 `Boundary` 布局和 `fullhide_2d` 数值算子；PIC 后端不属于本实现。

## 1. 边界数组

现有公共解包器固定 `Boundary(27)=R_0`。因此扩展字段只能追加在 `R_0` 之后，不能占用第 27 位：

- `Boundary(28)`：`transport_model_selector`，`0=legacy`，`1=pwn_cr_v1`。
- `Boundary(29)`：`stochastic_accel_norm`，无量纲动量扩散强度，默认 `0`。
- `Boundary(30)`：`escape_mode_selector`，`0=closed`，`1=free_outer`。

这个编号故意不同于早期草稿，因为覆盖 `Boundary(27)` 会破坏外部密度半径尺度。

## 2. 输运方程基础

PWN 多区模型和宇宙线输运模型通常求解包含空间输运、辐射损失、绝热损失、注入和可选再加速的 Fokker-Planck 型方程。`pwn_cr_v1` 不把 ASGARD 改造成完整 PWN/CR 代码，而是在现有正向激波 2D 电子求解器上加入当前有明确物理意义的项：

- \(\chi\) 空间扩散边界控制。
- \(\chi\) 局域绝热损失。
- \(\chi\) 局域磁场冷却和 Bohm 型扩散系数。
- 可选的 \(\log_{10}\gamma_e\) 空间扩散。

没有加入 PIC 命中概率、双磁场 patch 发射或 PIC scattering 算子。

## 3. BM 局域绝热损失

旧能量算子使用均匀膨胀近似：

\[
\frac{\mathrm{d}\log_{10}\gamma_e}{\mathrm{d}R}
= \frac{1}{R\ln 10}.
\]

`pwn_cr_v1` 改为从 BM 下游速度剖面计算局域径向散度。设激波半径为 \(R\)，激波 Lorentz 因子为 \(\Gamma_{\rm sh}\)，BM 坐标为 \(\chi\)。在超相对论有效区域，

\[
\Gamma_2(\chi)=\frac{\Gamma_{\rm sh}}{\sqrt{2\chi}},
\qquad
\beta_2(\chi)=\sqrt{1-\Gamma_2^{-2}}.
\]

局域半径比写为

\[
\frac{r}{R}
=1-\frac{\chi-1}{8\Gamma_{\rm sh}^2}.
\]

径向速度散度按

\[
\frac{\nabla\cdot v}{c}
=\frac{2\beta_2}{(r/R)R}
 + \frac{\mathrm{d}\beta_2}{\mathrm{d}r},
\qquad
\frac{\mathrm{d}\beta_2}{\mathrm{d}r}
=\frac{8}{\beta_2 R}.
\]

因此每单位激波半径推进的绝热能量损失系数为

\[
A_{\rm ad}
=\frac{\nabla\cdot v}{3\beta_{\rm sh}c\ln 10}.
\]

BM 表达式只在 \(\Gamma_2\ge 2\) 的区域使用。更远下游不满足超相对论速度梯度假设，且 \(\Gamma_2\to 1\) 时公式会发散；这些单元回到均匀各向同性膨胀闭合：

\[
A_{\rm ad}=\frac{1}{R\ln 10}.
\]

这里的回退不是数值兜底，而是 BM 近似适用域之外的显式物理边界。

## 4. \(\eta\) 空间扩散

物理扩散通量为

\[
F_x=-\kappa\frac{\partial N}{\partial x}.
\]

现有 `fullhide_2d` 几何把下游物理距离映射到

\[
\eta=\log_{10}\chi.
\]

`pwn_cr_v1` 继续在 \(\eta\) 空间使用三对角隐式求解，系数来自 `advance_eta_logchi_implicit` 已使用的几何因子，避免在每个子步中切换物理网格。

边界条件：

- `closed`：外侧 \(\eta\) 面扩散通量为零。
- `free_outer`：外侧 ghost density 为零，最外层单元出现向外逃逸 sink。

## 5. 微湍动磁场闭合

`pwn_cr_v1` 复用现有 `epsilon_b_floor`、`magnetic_decay_alpha_t` 和 `magnetic_decay_t0_s`。同一个局域 \(B_\chi\) 用于同步辐射/SSA 诊断、辐射冷却和 Bohm 型扩散系数估计。这个选择对应 GRB 余辉中激波后微湍动可衰减的物理动机，同时保持实现绑定到现有 `fullhide_2d` 字段。

## 6. 随机加速

`stochastic_accel_norm=0` 时关闭该算子。若其为正数，代码使用 Strang splitting：

1. \(\log_{10}\gamma_e\) 空间半步守恒扩散。
2. 辐射冷却和绝热冷却全步推进。
3. \(\log_{10}\gamma_e\) 空间半步守恒扩散。

能量网格两端使用零通量边界。隐式扩散系数为

\[
C_{\rm stoch}
= \frac{\Delta R\;a_{\rm stoch}}
        {R\left(\Delta\log_{10}\gamma_e\right)^2},
\]

其中 \(a_{\rm stoch}\) 是 `stochastic_accel_norm`。它只是每个动力学半径上的研究型无量纲归一化，不是独立物理量 \(D_{pp}\)。当 \(a_{\rm stoch}=0\) 时，该分支不被调用，结果应与无随机加速的 `pwn_cr_v1` 路径完全一致。

## 7. 参考依据

- Van Rensburg, Krüger & Venter 2018, MNRAS 477, 3853：空间依赖 PWN 输运模型。
- 3D PWN/CR 输运模型：MNRAS 528, 2749。
- Lemoine 2013, MNRAS 428, 845：GRB 余辉微湍动衰减。
- GALPROP theory page：宇宙线扩散、再加速和逃逸边界分类。
- BM/\(\chi\) 余辉剖面：MNRAS 442, 3495。
