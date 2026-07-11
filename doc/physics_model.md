# 物理模型

本文是 ASGARD 物理定义的唯一总源。它说明代码实际求解的量、近似与支持边界；离散方法和调用关系分别见[数值方法](numerical_methods.md)与[代码概览](code_overview.md)。

## 1. 计算链

一次完整计算依次经过：

1. 外部介质与喷流角结构；
2. 正向激波（FS）及可选反向激波（RS）动力学；
3. 非热电子注入、冷却与输运；
4. 同步辐射、自吸收（SSA）与可选 SSC；
5. 可选质子输运、强子辐射和二级粒子；
6. 可选含时 electron–photon–hadron 联合反馈；
7. 等到达时面（EATS）投影。

内部量使用 cgs。观测频率以 Hz 输入，时间以 s 输入，距离以 cm 输入。

## 2. 坐标与参考系

实验室系半径、时间和激波速度记为

\[
R,\qquad t,\qquad \beta c=\frac{dR}{dt},\qquad
\Gamma=(1-\beta^2)^{-1/2}.
\]

撇号量位于局域共动系。对与视线夹角为 \(\alpha\) 的流体元，

\[
\delta_D=[\Gamma(1-\beta\cos\alpha)]^{-1},
\qquad
\nu'=\frac{(1+z)\nu_{\rm obs}}{\delta_D}.
\]

观测者时间满足

\[
t_{\rm obs}=(1+z)\left(t-\frac{R\cos\alpha}{c}\right).
\]

代码从动力学轨迹插值满足该式的发射事件，不以单一的 \(R/\Gamma^2c\) 近似替代 EATS。

## 3. 外部介质

外部质子数密度写为

\[
n_{\rm ext}(R)=n_0\left(\frac{R}{R_0}\right)^{-k}.
\]

当前介质类型包括均匀介质、风介质以及由公开参数定义的分段密度增强。扫掠静质量为

\[
dm=4\pi R^2 m_p n_{\rm ext}(R)\,dR
\]

的喷流立体角份额。结构化喷流逐角元独立演化；角元间没有横向能量交换，除非启用已有的 spreading 动力学处方。

密度跳变会改变局域扫掠率并可激发次级 RS。它不是简单的事后光变缩放。

## 4. 正向激波动力学

动力学状态包含 \(R,\Gamma\)、扫掠质量、内能和喷流开角等量。绝热情形满足总能量在扫掠、膨胀和几何演化间闭合；辐射效率进入相应能量损失项。

激波后数密度和内能密度由相对论 Rankine–Hugoniot 关系给出。常用超相对论极限

\[
n'_2\simeq4\Gamma n_{\rm ext},\qquad
e'_2\simeq4\Gamma^2 n_{\rm ext}m_pc^2
\]

只用于理解；内核保留跨相对论区间的表达式。

磁场由

\[
\frac{B'^2}{8\pi}=\epsilon_B e'
\]

闭合。该 \(B'\) 是辐射区的无序场强，不等同于 RS 上游的有序磁化参数。

## 5. 电子注入

非热电子获得耗散内能的份额 \(\epsilon_e\)，参与加速的电子份额为 \(\xi_e\)。注入谱通常为

\[
Q_e(\gamma)=Q_0\gamma^{-p},
\qquad \gamma_m\le\gamma\le\gamma_M.
\]

归一化同时满足注入粒子数和能量预算。对 \(p\) 接近 2 的情况使用相应积分极限，不用任意小偏移替代对数极限。

最小洛伦兹因子由每个被加速电子可得能量确定；最大洛伦兹因子由加速时间与辐射、动力学时间的竞争确定。

热电子开关增加独立的热分量；它不改变非热谱的定义。

## 6. 电子输运

壳层共动系电子谱满足

\[
\frac{\partial N_e}{\partial t'}+
\frac{\partial}{\partial\gamma}(\dot\gamma N_e)
=Q_e-\frac{N_e}{t'_{\rm esc}}+Q_{e,\rm sec}.
\]

能损率包括同步、逆康普顿和绝热项：

\[
\dot\gamma=\dot\gamma_{\rm syn}+
\dot\gamma_{\rm IC}+\dot\gamma_{\rm ad}.
\]

同步损失为

\[
\dot\gamma_{\rm syn}=-\frac{\sigma_TB'^2}{6\pi m_ec}\gamma^2.
\]

逆康普顿核在光子能量接近 Klein–Nishina 区域时直接降低截面；不存在可关闭的额外“KN 修正”后处理开关。

公开 1D solver 在单壳能量网格上演化 \(N_e(\gamma)\)。`fullhide_2d` 使用有限 q-mass 壳层坐标。

## 7. 有限厚壳层

2D 输运的第二坐标是有限壳层质量坐标 \(q\)，不是无限的 \(\eta=\log_{10}\chi\) 网格。每个 q 单元保存其粒子谱、半径、体积权重与流体状态。

`chi_grid` 是每个 q 单元的 Blandford–McKee 等价诊断坐标。观测投影实际使用

- `chi_radius_cm`；
- `chi_gamma_bulk`；
- `chi_dvolume_weight`。

`chi_eats_2d` 当前只替换 FS 同步辐射与 SSA 的 light-curve 投影。SSC、强子过程、pair cascade 和二级反馈仍是壳层级契约，不能被称为 chi-resolved 输运。

## 8. 同步辐射

单电子特征频率为

\[
\nu'_c=\frac{3eB'}{4\pi m_ec}\gamma^2\sin\alpha_p.
\]

各向同性 pitch-angle 分布的体发射率为

\[
j'_{\nu'}=\frac{1}{4\pi}\int N_e(\gamma)
P'_{\nu'}(\gamma)\,d\gamma.
\]

默认 `fixed_grid` 核直接积分同步函数。`cyclotron` 路径为深牛顿区 \(\gamma<2\) 提供连续的回旋–同步处理。

SSA 吸收系数由同一电子谱计算，均匀壳层传输因子为

\[
T_{\rm ssa}(\tau)=\frac{1-e^{-\tau}}{\tau}.
\]

`fullhide` 路径复用同一同步光深网格确定 \(\nu_a\)，避免另一次物理上重复的根搜索。

## 9. SSC 与外部种子场

SSC 使用壳层内同步光子作为种子，计算

\[
j'_{\epsilon_s}=c\int d\gamma\,N_e(\gamma)
\int d\epsilon\,n'_\gamma(\epsilon)
\frac{d\sigma_{\rm IC}}{d\epsilon_s}.
\]

核保留完整的 Klein–Nishina 能量依赖。`include_ssc` 控制是否输出 SSC 分量；冷却模型选择决定 IC 损失如何进入电子输运，两者不是同一个键。

FS–RS cross-IC 使用另一激波区的光子场，但仍在接收区电子分布上散射。几何稀释和参考系变换属于显式模型输入。

## 10. Gamma-gamma 吸收

高能光子与目标光子满足

\[
\epsilon_1\epsilon_2(1-\cos\theta)\ge2(m_ec^2)^2
\]

时可产生电子–正电子对。光深由角平均截面与目标光子谱积分得到：

\[
\tau_{\gamma\gamma}(E_\gamma)=
\ell'\int n'_\gamma(\epsilon)\,
\bar\sigma_{\gamma\gamma}(E_\gamma,\epsilon)\,d\epsilon.
\]

transfer 由几何所有者决定。共空间均匀发射与吸收使用

\[
\psi(\tau)=\frac{1-e^{-\tau}}{\tau},
\]

前景屏使用 \(e^{-\tau}\)，有序多区则沿传播方向逐区积分。生产路径中仍有从 \(R,\Gamma\) 推测 \(R/(12\Gamma)\) 的历史实现；该长度只对应均匀外介质的 one-zone 极限，wind、density jump、prompt 和多区路径尚未完成几何所有权闭合，见根目录 `BUG.md`。

吸收能量只有在启用对应二级注入时才进入 pair 源项；不能同时由 observer 层和 joint solver 重复拥有。

## 11. Pair cascade

`separated` 模式按已完成的一次辐射解计算 pair 注入和后续辐射；它不回写原电子–光子演化。

`joint` 的目标合同是让 pair 注入只属于含时 electron–photon 系统：光子 sink、pair 源和二级辐射共享同一能量预算，observer 只投影最终分量。当前 joint secondary pair 的坐标 Jacobian、累计状态与辐射所有权仍有已登记缺陷，因此相关组合尚不能视为完成验收。

当前级联是壳层平均的一维能量输运，不是空间分辨 Monte Carlo cascade，也不是 chi-resolved cascade。

## 12. 质子注入与输运

非热质子谱写为

\[
Q_p(E_p)=Q_{p,0}E_p^{-s_p},
\qquad E_{p,\min}\le E_p\le E_{p,\max}.
\]

质子输运满足

\[
\frac{\partial N_p}{\partial t'}+
\frac{\partial}{\partial E_p}(\dot E_pN_p)
=Q_p-\frac{N_p}{t'_{p,\rm esc}}.
\]

损失可包括绝热、质子同步、pγ、Bethe–Heitler（BH）和 pp。过程开关只控制已实现链路，不隐含空间扩散。

## 13. Proton synchrotron

质子同步公式由电子同步表达式作质量替换得到。固定电荷下特征频率随 \(m^{-1}\) 缩放，辐射功率随 \(m^{-2}\) 缩放。该过程同时进入质子损失与光子源，二者必须共享同一能量预算。

## 14. Photopion pγ

pγ 相互作用由共动目标光子谱驱动。当前正式响应模型为 `hummer_2010_response`；`disabled` 完全关闭 pγ 源与损失。

相互作用率具有形式

\[
t_{p\gamma}^{-1}(\gamma_p)=
\frac{c}{2\gamma_p^2}
\int_{\bar\epsilon_{\rm th}}^\infty d\bar\epsilon\,
\sigma_{p\gamma}\kappa\bar\epsilon
\int_{\bar\epsilon/(2\gamma_p)}^\infty
\frac{n'_\gamma(\epsilon)}{\epsilon^2}\,d\epsilon.
\]

响应表同时产生光子、电子/正电子和 neutrino。所有产物必须来自同一质子能量损失，不能分别归一化。

## 15. Bethe–Heitler

BH 过程

\[
p+\gamma\rightarrow p+e^-+e^+
\]

主要把质子能量注入二级 pair。它对目标光子谱的阈值积分同时决定质子损失和 pair 源；`joint` 模式将 pair 注入电子方程，`separated` 模式在主输运后组装。

## 16. Proton-proton

pp 总链始终计算现有的质子损失、neutrino 和二级 pair。π0 衰变 gamma 谱由

`Hadronic.pp_gamma_model` 选择：

- `delta`：默认，保持历史正式结果；
- `geant4`、`sibyll`、`qgsjet`、`pythia8`：Kafexhiu/AM3 详细谱，显式 opt-in。

模型 ID 只在 Python 边界映射一次。详细模型仅替换 π0 gamma 分量，不改变 pp 的 `ploss/qnu/qpair`。

详细参数化在若干生成器接合能量处继承上游分段跳变；代码不添加平滑、clamp 或事后归一化掩盖这一性质。

## 17. Neutrino

neutrino 来自 pγ 和 pp 的带电介子链。源端谱与电磁分量共享同一强子能量预算；传播振荡若由调用层处理，也不得回写源端输运。

neutrino 不受 SSA 或 gamma-gamma 吸收。红移和几何稀释仍按观测者投影处理。

## 18. 含时联合反馈

`electron_photon_coupling="joint"` 的目标是在同一壳层时间步内联合演化

\[
\{N_p(E_p),N_e(\gamma),n_\gamma(\epsilon)\}.
\]

方程通过以下源和 sink 闭合：

- 电子同步与 IC：电子损失对应光子源；
- gamma-gamma：光子 sink 对应 pair 源；
- BH：质子损失对应 pair 源；
- pγ/pp：质子损失对应 gamma、pair 与 neutrino 产物；
- 光子逃逸：只从局域光子场移除能量。

`separated` 和 `joint` 是互斥算法契约，不应在一个结果中叠加同一二级分量。完整定义见 [含时二级反馈](joint_secondary_feedback_physics.md)。

## 19. 反向激波

RS 由 ejecta 与接触面之间的相对运动驱动。上游磁化定义为

\[
\sigma=\frac{B_u'^2}{4\pi\rho_u'c^2}.
\]

主 RS 使用冷上游、横向有序场的一维 MHD jump，压力平衡决定接触面状态。\(\sigma=0\) 恢复流体极限；有限 \(\sigma\) 同时改变压缩、热化和下游磁能。

RS crossing 前持续注入新电子，crossing 后进入纯冷却/绝热演化。次级密度跳变 RS 当前是流体诊断，不继承主 RS 的有序场压缩与磁焓反馈。

## 20. Prompt 内部激波

`prompt/` 是两壳内部激波的快照与诊断工作流，不是稳定 `asgard_core` API 或拟合分支。快壳追上慢壳后，接触面两侧形成 FS/RS；接触面洛伦兹因子由两侧压力平衡求得。

当前入口为

- `prompt.internal_shock`；
- `prompt.radiation`；
- `prompt.eats`。

详细假设与例子见 [Prompt 内部激波教程](prompt_internal_shock_tutorial.md)。

## 21. 偏振

局域同步偏振由磁场几何、电子谱和视线方向决定。总 Stokes 参数在 EATS 上线性叠加：

\[
I=\sum I_i,\qquad Q=\sum Q_i,\qquad U=\sum U_i.
\]

偏振度和位置角为

\[
\Pi=\frac{\sqrt{Q^2+U^2}}{I},\qquad
\chi_P=\frac12\operatorname{atan2}(U,Q).
\]

不能先平均局域偏振度再乘总光度。

## 22. 观测者投影

观测通量由 EATS 上各角元和壳层贡献积分：

\[
F_{\nu_{\rm obs}}=
\frac{1+z}{4\pi d_L^2}
\int_{\rm EATS}\delta_D^3 L'_{\nu'}\,d\Omega,
\qquad
\nu'=\frac{(1+z)\nu_{\rm obs}}{\delta_D}.
\]

具体 Doppler 次数随被积量是体发射率、谱光度还是已含体积权重的数组而定；实现只在一个层级拥有该变换。

`total_only` 返回总量，`full_components` 保留过程分量。总量按固定顺序求和以保持可复现浮点结果。

## 23. 模型绑定关系

若两个键共同定义一个物理近似，应在公开边界绑定，而不是允许无意义笛卡尔积。当前重要关系是：

- `chi_eats_2d` 只与有限 q-mass 的 `fullhide_2d` FS 同步/SSA 路径配套；
- `joint` 拥有二级 pair 反馈，observer 不重复注入；
- pp 详细模型只控制 π0 gamma；
- Fan–Y 路径尚未公开：必须先按原始文献修正有限 \(\gamma_M\)、分支连续性和自洽 \(\gamma_c\)；
- RS 与 prompt 的磁化 jump 是不同计算对象，不可互换 selector。

## 24. 支持边界

当前代码不声称支持：

- 全三维 MHD 或角元间磁流体耦合；
- 空间分辨的 hadronic/SSC/pair chi 输运；
- Monte Carlo 电磁级联；
- prompt 与 afterglow 的统一含时拟合；
- 未经原始文献推导验证的 Fan–Y 公开 selector；
- 通过平滑或经验重归一化修复微物理参数化。

## 25. 物理验收

每次物理修改至少检查：

1. 粒子数、能量和单位闭合；
2. 阈值以下源项为零，允许区有限且非负；
3. 关闭新过程时旧路径不变；
4. 共享损失的多种产物使用同一预算；
5. 真实时空演化随参数连续且平滑；
6. EATS、Doppler、红移和距离因子只应用一次；
7. 数值差异若来自运算顺序，应同时给出性能收益与误差尺度；
8. 上游公式本身的不连续必须明确记录，不能用后处理掩盖。

物理问题先记录到 `BUG.md`。修复必须同时通过真实入口、编译和文献核对。
