# 含时二级反馈的物理契约

本文说明 `electron_photon_coupling="joint"` 的物理含义。它固定电子、光子、强子和二级对在同一半径坐标 \(R\) 上闭合的要求；算法执行顺序见 `doc/joint_secondary_feedback_algorithm.md`。

## 1. 为什么使用半径 \(R\)

ASGARD 的动力学主状态沿激波半径 \(R\) 存储。若电子、质子、二级 \(e^\pm\) 和光子场在同一个壳层内相互反馈，它们必须使用同一个自变量。

设壳层 bulk Lorentz 因子为 \(\Gamma\)，速度为

\[
\beta=\sqrt{1-\Gamma^{-2}}.
\]

激波半径满足

\[
\frac{\mathrm{d}R}{\mathrm{d}t_{\rm lab}}=\beta c.
\]

共动时间与实验室系时间的关系为

\[
\mathrm{d}t'=\frac{\mathrm{d}t_{\rm lab}}{\Gamma}.
\]

因此

\[
\frac{\mathrm{d}t'}{\mathrm{d}R}
=\frac{1}{\beta\Gamma c}.
\]

所有自然单位为 \({\rm s}^{-1}\) 的微物理率进入 \(R\) 坐标输运时必须转换为

\[
\Lambda_R=\frac{\Lambda_t}{\beta\Gamma c},
\qquad
Q_R=\frac{Q_t}{\beta\Gamma c},
\]

\[
\dot{\gamma}_R=\frac{\dot{\gamma}_t}{\beta\Gamma c},
\qquad
\dot{E}_R=\frac{\dot{E}_t}{\beta\Gamma c}.
\]

这里 \(\Lambda\) 是损失或吸收率，\(Q\) 是注入率，\(\dot{\gamma}\) 和 \(\dot{E}\) 是连续冷却漂移速度。

## 2. 变量归一化

公式使用连续谱记号：

- \(N_p(E_p,R)\)：质子谱。
- \(N_e(\gamma_e,R)\)：电子/正电子谱。
- \(n_\gamma(\epsilon_\gamma,R)\)：目标光子谱。
- \(Q_{x,R}\)：已经换算到 per-\(R\) 的源项。
- \(\Lambda_{\gamma,R}\)：已经换算到 per-\(R\) 的光子吸收率。

实现中同时存在 shell-integrated spectra 和 local density spectra。任一源项进入某个输运方程前，必须与该方程推进的谱变量使用同一归一化。

## 3. 质子方程

formal 1D 强子路径在每个壳层推进

\[
\frac{\partial N_p}{\partial R}
=Q_{p,{\rm shock},R}
-\frac{\partial}{\partial E_p}
\left(\dot{E}_{p,{\rm loss},R}N_p\right)
+Q_{p,{\rm reinj},R}.
\]

其中

\[
\dot{E}_{p,{\rm loss},R}
=\dot{E}_{{\rm ad},R}
+\dot{E}_{p{\rm syn},R}
+\dot{E}_{{\rm BH},R}
+\dot{E}_{p\gamma,R}
+\dot{E}_{pp,R}.
\]

若某个过程只输出 photon/neutrino luminosity，而没有给出可回灌的输运源项，则它只能作为辐射输出或诊断，不能伪造成反馈源。

## 4. 电子/正电子方程

`joint` 模式下，电子方程直接接收强子和 \(\gamma\gamma\) 过程产生的二级源项：

\[
\frac{\partial N_e}{\partial R}
=Q_{e,{\rm shock},R}
+Q_{e,{\rm secondary},R}
-\frac{\partial}{\partial \gamma_e}
\left[
\left(
\dot{\gamma}_{{\rm ad},R}
+\dot{\gamma}_{{\rm syn},R}
+\dot{\gamma}_{{\rm IC},R}[n_\gamma]
\right)N_e
\right].
\]

当前实现的二级源项为

\[
Q_{e,{\rm secondary},R}
=Q_{e,{\rm BH},R}
+Q_{e,pp,R}
+Q_{e,\gamma\gamma,R}.
\]

未接入的项不能用总能量外推构造谱形。尤其是 \(p\gamma/\pi/\mu\) decay chain 的 \(e^\pm\) 注入，必须等待 formal kernel 输出谱形和归一化后才能进入 \(Q_{e,{\rm secondary},R}\)。

## 5. 光子场方程

joint 光子场同时服务电子 IC 冷却、强子 target field 和吸收。连续形式可写为

\[
\frac{\partial n_\gamma}{\partial R}
=Q_{\gamma,e{\rm syn},R}
+Q_{\gamma,e{\rm IC},R}
+Q_{\gamma,{\rm pair\,syn},R}
+Q_{\gamma,{\rm had,out},R}
-\Lambda_{\gamma,R}n_\gamma.
\]

必须区分两类光子：

- 参与 joint seed 闭合的目标光子场：正向激波电子同步辐射 seed、SSC/IC seed、\(\gamma\gamma\) pair synch seed，以及 \(p\gamma\)、BH、\(\gamma\gamma\) survival 作用后的 field。
- 作为 observer component 输出的强子辐射 luminosity：质子同步辐射、\(p\gamma\) gamma、BH pair radiation、hadronic IC、pair-production radiation 等。

并非所有 observer-side hadronic photon luminosity 都已经回灌成 target photon density。后续若要加入某个强子 photon source，必须先定义 luminosity 到壳层 photon density 的几何、逃逸时间和吸收归一化。

## 6. 光子 sink

光子吸收率在半径坐标中的形式为

\[
\Lambda_{\gamma,R}
=\frac{
\alpha_{p\gamma}
+\alpha_{\rm BH}
+\alpha_{\gamma\gamma}
+t_{\rm esc}^{-1}
}{\beta\Gamma c}.
\]

当前代码中：

- \(p\gamma\) photon depletion 使用 formal \(p\gamma\) kernel 输出的 photon survival。
- BH photon loss 使用 BH kernel 输出的 photon loss rate，并与 proton loss 和 pair rate 同一归一化。
- \(\gamma\gamma\) absorption 使用 pair-production branch 的 photon loss rate 或 cascade optical depth。
- 若某路径没有显式建模 \(t_{\rm esc}^{-1}\)，不能用经验 sink 代替。

## 7. 能量预算

### IC 预算

同一组 \(N_e\) 和 \(n_\gamma\) 必须同时决定

\[
P_{e,{\rm IC\,loss}}
\quad{\rm and}\quad
P_{\gamma,{\rm IC\,source}}.
\]

如果只改变 electron cooling 而没有同步改变 IC photon source，电子能量会从系统中消失。当前 joint 路径通过同一 seed 传入 coupled electron pass；`tests/electron_photon_joint_secondary_feedback_smoke.py` 是对应诊断入口，但当前工作树会在 formal hadronic electron-energy grid contract 处失败，需先修复该契约后再把它作为绿色验收。

### BH 预算

Bethe-Heitler kernel 同时给出

\[
\dot{N}_{p,{\rm loss}}(E_p),
\qquad
\dot{N}_{e^\pm}(E_e),
\qquad
\dot{n}_{\gamma,{\rm loss}}(\epsilon_\gamma).
\]

三者必须来自同一个微物理算子。允许存在由网格截断、离散积分和未观测逃逸项造成的可解释差异；不允许用经验比例补齐 photon sink 或 pair source。

### \(\gamma\gamma\) 预算

\(\gamma\gamma\) absorption 的当前物理链条是

\[
\gamma+\gamma
\rightarrow e^+ + e^-
\rightarrow \text{pair synchrotron photons}.
\]

当前接入的是 pair/synch cascade。IC-mediated electromagnetic cascade 需要额外的 photon/\(e^\pm\) source-sink 方程和 IC kernel 预算，尚未实现。

## 8. `separated` 与 `joint`

`separated` 是默认模式：

```text
primary electron solve
-> photon field
-> hadronic solve
-> BH/secondary post-merge or diagnostic output
```

它适合 weak-feedback 情况，因为二级对不显著改变电子冷却和目标光子场。

`joint` 是 opt-in 模式：电子、光子场、强子输运和二级对在同一 \(R\) 网格内迭代闭合。它的目标不是改变强子微物理，而是让二级 \(e^\pm\) 与 photon sink/source 在同一壳层坐标下反馈到 electron/photon state。

## 9. 支持边界

当前 `joint` 支持：

- 正向激波。
- `electron_solver="fullhide_1d"`。
- `hadronic_solver="am3_1d"`。
- `bethe_heitler=True`。
- `ssc=True` 且 `index_y=1`。
- fixed electron substeps。

当前 `joint` 不支持：

- 反向激波 full-chain feedback。
- \(\chi\)-resolved / 2D 强子输运。
- 结构化喷流后端。
- `fullhide_1d` 之外的 electron solver。
- 非 `hummer_2010_response` 的 \(p\gamma\) scheme。
- IC-mediated electromagnetic cascade。
- formal \(p\gamma/\pi/\mu\) \(e^\pm\) feedback。

这些边界必须显式报错，不能 fallback 到 `separated` 或静默忽略。

## 10. 物理验收

合格结果应满足：

- 默认 `separated` 回归不变。
- weak-feedback 极限下 `joint` 接近 `separated`。
- strong-wind / strong-BH 场景中 \(\tau_{\rm BH}\)、二级 pair 谱、photon source 和 light curve 随 \(R\) 或 observer time 平滑演化。
- \(\gamma\gamma\) 打开时，\(\tau_{\gamma\gamma}\)、pair luminosity 和 cascade synchrotron luminosity 连续。
- IC electron loss 与 IC photon source 的壳层积分闭合。
- BH proton loss、BH pair injection 和 BH photon loss 同量级、同谱形来源可追踪。

若图中出现孤立尖峰、锯齿或不可解释断崖，应先查坐标换算、源项归一化、网格映射和 photon survival，而不是后处理 smoothing。
