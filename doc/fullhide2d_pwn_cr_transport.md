# `fullhide_2d` PWN/CR 输运扩展

本文固定可选 `fullhide2d_transport_model="pwn_cr_v1"` 的算法边界。默认
`legacy` 不变；该扩展不是 PIC 或完整 PWN 模型。

## 1. ABI 字段

`Boundary(27)=R_0`，density-jump 槽位固定到 `Boundary(52)`；
`Boundary(53)=N` 后依次存放 N 个 profile radius 和 N 个 density。扩展只追加到动态数组末尾：

- `Boundary(n-2)`：transport，`0=legacy`、`1=pwn_cr_v1`；
- `Boundary(n-1)`：无量纲 stochastic acceleration 强度；
- `Boundary(n)`：escape，`0=closed`、`1=free_outer`。

不得占用旧半径或 density-jump 槽位。

## 2. 方程

`pwn_cr_v1` 在 finite-q electron transport 上加入：

- q 空间扩散及外边界；
- q-local 绝热和同步冷却；
- 同一局域磁场下的 Bohm 型扩散估计；
- 可选 `log10(gamma)` 随机加速。

q 是有限主动壳层质量坐标，不是无限 BM χ 网格。

## 3. 局域绝热损失

在超相对论 BM 有效域，

\[
\Gamma_2(\chi)=\Gamma_{\rm sh}/\sqrt{2\chi},\qquad
r/R=1-(\chi-1)/(8\Gamma_{\rm sh}^2),
\]

并由径向速度散度得到

\[
A_{\rm ad}=\frac{\nabla\!\cdot v}
{3\beta_{\rm sh}c\ln 10}.
\]

只在 `Gamma2 >= 2` 使用 BM 梯度。更远下游采用均匀各向同性膨胀
`A_ad=1/(R ln 10)`；这是两种物理近似的显式适用域，不是失败 fallback。

## 4. q 扩散和逃逸

主动域为

\[
q\in[0,q_{\rm active}],\qquad q_{\rm active}=1-(3/4)^4.
\]

`advance_q_pwncr_implicit` 使用三对角隐式通量。`closed` 令外侧面通量为零；
`free_outer` 令外侧 ghost density 为零，从最外 cell 产生逃逸 sink。

## 5. 磁场与随机加速

局域 `B_chi` 同时用于同步/SSA、辐射冷却和 Bohm 扩散。它复用
`epsilon_b_floor`、`magnetic_decay_alpha_t` 和 `magnetic_decay_t0_s`。

随机加速为正时采用 Strang splitting：能量扩散半步、cooling 全步、能量扩散
半步。能量边界零通量，离散强度

\[
C_{\rm stoch}=\frac{\Delta R\,a_{\rm stoch}}
{R(\Delta\log_{10}\gamma)^2}.
\]

`a_stoch` 是研究型无量纲归一化，不等同于独立物理量 `D_pp`。为零时不调用该
算子，必须回到同一 `pwn_cr_v1` 无随机加速结果。

## 6. 验收

- legacy 与关闭扩展项的结果不变；
- closed 边界守恒粒子数，free_outer 的损失等于边界通量；
- diffusion-only 展宽平滑且网格收敛；
- magnetic decay 同时一致地影响 cooling、emission 和 diffusion；
- 2D state 与 χ-resolved observer 使用同一 cell geometry；
- 强制构建 `electron_forward_transport_2d` 相关扩展并检查 line truncation。

参考物理包括 Blandford--McKee 下游剖面、Lemoine 2013 微湍动衰减，以及
Van Rensburg et al. 2018 的空间依赖 PWN 输运。实现只采用上面明确列出的算子。
