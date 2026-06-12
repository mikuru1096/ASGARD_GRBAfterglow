# 多密度增强下次级反向激波三区模型

## 目标

本计划固定 ASGARD v1 的多密度增强反向激波模型合同。v1 只处理 ISM 中多个平滑 density bump 触发的次级反向激波电子同步辐射，不处理强子过程、RS SSC、FS/RS cross-zone IC、wind、structured jet 或 2D/chi resolved transport。

## 区域定义

区域 1 是未激波外部介质，只作为密度边界 \(n_1(R)\) 输入，不作为演化变量。区域 2 是透射 forward shock 下游，和区域 3 共用接触间断速度与压强。区域 4 是 density bump 前已经被 forward shock 加热的旧 shocked shell，状态来自当前 ASGARD 动力学局部量。区域 3 是被次级 reverse shock 再激波的旧 shocked shell；新增辐射只来自区域 3 的额外耗散能，不能把旧 forward-shock 电子辐射重复计入。

接触间断条件为

\[
\Gamma_2=\Gamma_3=\Gamma_c,\qquad p_2=p_3 .
\]

区域 1 采用冷上游 \(p_1=0,\Gamma_1=1,n_1=n(R)\)。区域 4 采用热上游 \(p_4,e_4,n_4,\Gamma_4\)，不能套用原始 ejecta reverse shock 的冷上游公式。

## 密度合同

多 density bump 使用数组接口：

- `jump_r_cm`
- `jump_factor`
- `jump_width_log10`

ISM 密度为

\[
n(R)=n_0\left[1+\sum_j(f_j-1)\exp\left(-\frac{(\log_{10}R-\log_{10}R_j)^2}{2w_j^2}\right)\right].
\]

旧 `r_tr/f_jump/f_wide` 是单 bump 兼容入口；多数组非空时使用数组合同。v1 不从 callable medium 拟合多 bump，也不支持 wind bump。

## 跳变条件

每个 density bump 的上升段按局部 Riemann 问题求解。跨 shock 使用相对论 Rankine-Hugoniot 守恒：

\[
[\rho u^\mu]n_\mu=0,\qquad [T^{\mu\nu}]n_\nu=0 .
\]

数值实现以 \(\Gamma_c\) 为未知量求

\[
p_2(\Gamma_c)-p_3(\Gamma_c)=0 .
\]

无物理解时必须直接报错。禁止加入经验 smoothing、fallback 或事后光变修正。

## 电子性质

区域 3 新注入电子最小洛伦兹因子使用次级 shock 额外耗散能：

\[
\gamma_{m,3}=1+
\frac{\epsilon_{e,3}}{\xi_{N,3}}
\frac{p_3-2}{p_3-1}
\frac{u_{\rm diss,3}}{n_3 m_e c^2}.
\]

\(u_{\rm diss,3}\) 定义为 Rankine-Hugoniot 下游热能减去区域 4 绝热压缩贡献，避免把旧 shocked shell 的已有热能当作新注入能量。

区域 3 磁场为

\[
B_3=\sqrt{8\pi\epsilon_{B,3} e_{\rm int,3}} .
\]

区域 3 冷却、SSA 和同步谱复用现有 reverse electron transport/radiation kernel 的 shell-level 合同，但输入必须是每个次级 RS reservoir 的 \(M_3,U_3,V_3,B_3,\Gamma_3,\gamma_{43}\)。

## 输出与验收

`rev.sync` 等于原始 ejecta RS synch 加所有次级 RS synch。details 需暴露每个 bump 的 \(p_3,\Gamma_c,\gamma_{43},u_{\rm diss,3},B_3,\nu_m,\nu_c,\nu_a\)。

验收包括：

- 单 bump 数组与旧 `r_tr/f_jump/f_wide` density profile 等价。
- 无多 bump 时现有 RS benchmark 不变。
- 每个 bump 的 \(\Gamma_c,p_2,p_3,\gamma_{43},u_{\rm diss,3},B_3\) 连续平滑。
- 区域 3 注入能量积分等于 RH 给出的次级 RS 额外耗散能乘 \(\epsilon_e\)。
- FS 不出现非物理尖锐重亮；RS 对 bump 更敏感，但随 bump 宽度和幅度连续变化。
- 受影响 Fortran build、`-Wline-truncation` 和最小 multi-bump RS smoke 必须通过。

## 文献来源

- Nakar & Granot 2007, density jump 四区结构和观测平滑性约束。
- Uhm & Zhang 2014, density bump/void 对 long-lived reverse shock 响应的动力学机制。
