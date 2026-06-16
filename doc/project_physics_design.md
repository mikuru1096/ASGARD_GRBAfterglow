# ASGARD 项目物理设计总纲

本文档给出 ASGARD 当前工作树的项目级物理设计。它回答三个问题：

- 这个项目在物理上求解什么系统。
- 每个重要辐射/粒子过程在总方程中扮演什么角色。
- 哪些过程已经自洽反馈，哪些只作为输出或诊断，哪些明确不支持。

更细的专题说明见 `doc/physical_processes.md`、`doc/physics_model.md`、`doc/joint_secondary_feedback_physics.md`、`doc/hadronic_chi_transport_decision.md`、`doc/pair_cascade_extension_boundary.md` 和 `doc/public_backend_limits.md`。

## 1. 总体物理对象

ASGARD 建模的是 GRB afterglow 中相对论外流与外部介质相互作用后形成的辐射系统：

```text
relativistic ejecta
-> forward shock / reverse shock
-> non-thermal particles
-> local radiation and photon fields
-> absorption, pair production, cascade, hadronic channels
-> equal-arrival-time observer projection
```

主状态沿 shock radius `R` 演化。观测者时间 `t_obs` 是输出投影坐标，不是所有输运方程的共同自变量。任何需要在同一 shell 中闭合的粒子和光子反馈，都必须先投影到 `R` 坐标。

## 2. 几何与动力学

### Forward shock

正向激波是默认主线。动力学求解 blast-wave radius、bulk Lorentz factor、swept-up mass 和 observer-time mapping。支持的外部介质是：

- ISM：常数密度 `n_ism`。
- Wind：`rho ∝ R^-2`，当前 backend 支持 `k=2`。
- density jump / transition：作为显式介质参数进入动力学。
- energy injection / magnetar-like injection：通过 jet/dynamics 参数进入。

正向激波的物理输出为后续电子、强子和观测投影提供：

```text
R_i, Gamma_i, t_obs,i, swept mass, external density
```

### Reverse shock

反向激波描述 ejecta 内部 shocked region。当前物理基线包括：

- shock-front `gamma34` 作为电子注入能标来源。
- region-3 thermal state `U3/V3`。
- turbulent magnetic field。
- optional upstream magnetization `sigma`，并加入 ordered magnetic component。
- crossing 前后分段动力学，避免 RK step 跨越物理分支。
- 多密度增强触发的次级反向激波电子同步分支；该分支使用热上游区域 4 和区域 3 新耗散能，不重复计算旧 FS 电子辐射。

`sigma -> 0` 必须回到 unmagnetized baseline。VegasAfterglow 是 comparison backend，不是目标真值。

多密度增强反向激波的四区物理图像、数组密度合同和验收条件见 `doc/physical_processes.md` 的“多密度增强下的次级反向激波”。

### 喷流几何

Public jet constructors cover tophat、Gaussian 和 power-law profiles。结构化喷流通过 patch/structured backend 聚合每个角向单元的局域动力学与辐射，再做 observer projection。

当前 hadronic 和 pair feedback 的 formal transport contract 是壳层级 1D。结构化喷流可以复用壳层级 radiation/hadronic outputs，但不等于已经有 \(\chi\)-local 或 patch-local hadronic feedback。

## 3. 电子物理

电子谱是正向激波和反向激波辐射的核心状态。典型电子源项是 shock acceleration 注入的非热幂律：

\[
Q_{e,{\rm shock}}(\gamma_e)
\propto
\gamma_e^{-p}
\exp\!\left(-\frac{\gamma_e}{\gamma_{\max}}\right).
\]

归一化由 swept-up electron number、`epsilon_e`、`accelerated_electron_fraction` 和局域 shock energy 决定。可选 thermal electron 分支使用同一物理 shell 的 electron count 和 shock-front thermal scale。

电子连续方程的物理形式是

\[
\frac{\partial N_e}{\partial R}
=Q_{e,{\rm shock},R}
+Q_{e,{\rm secondary},R}
-\frac{\partial}{\partial\gamma_e}
\left(\dot{\gamma}_{e,R}N_e\right).
\]

其中

\[
\dot{\gamma}_{e,R}
=\dot{\gamma}_{{\rm ad},R}
+\dot{\gamma}_{{\rm syn},R}
+\dot{\gamma}_{{\rm IC},R}
+\dot{\gamma}_{{\rm SSA},R}.
\]

当前 `Q_e,secondary,R` 只在 joint feedback 路径中作为电子方程源项进入；默认 separated 路径保留后处理合并语义。

## 4. 磁场与同步辐射

局域磁场主要由 shock internal energy 和 `epsilon_B` 确定。对正向激波，磁场是壳层级 state；对 2D/\(\chi\) observer projection，当前磁场在一个 shell 的投影网格内保持当前实现约定，不引入未定义的 \(\chi\)-local 磁场演化。

同步辐射输出包括：

- electron synchrotron luminosity。
- synchrotron seed photon field。
- characteristic frequencies `nu_m`, `nu_c`, `nu_a`。
- synchrotron self-absorption transfer。
- polarization Stokes emissivity。

同步辐射是多个后续过程的 seed：

```text
electron synch seed
-> SSC / IC cooling
-> hadronic pγ / BH target photons
-> gamma-gamma absorption target
```

因此 seed field 的单位、几何归一化和 shell 对齐必须保持稳定。

## 5. IC / SSC 物理

SSC 是电子与本地 synchrotron photon field 的 inverse-Compton scattering。它在项目中有两个角色：

- 电子冷却项 `dotgamma_IC`。
- observer-side IC photon source。

joint feedback 要求同一个 `N_e` 和同一个 photon seed 同时决定这两个量。只修改 electron cooling 而不生成对应 photon source 会破坏能量预算。

当前 IC kernel 包含 KN/Jones 类型截面约束。`index_y=1` 是 joint 预算路径要求的数值 IC cooling 模式。

## 6. 光子场物理

ASGARD 中 photon field 有两种不同用途，必须区分：

- local target photon field：参与 IC cooling、pγ、BH、gamma-gamma absorption 和 joint feedback。
- observer luminosity component：投影给观测者的辐射输出。

一个 observer luminosity 不能自动当作 local target photon density。若要把某个 luminosity 回灌进 photon equation，必须定义：

```text
shell volume
escape time
angular / path geometry
absorption and survival
frequency or energy variable conversion
```

当前 `PhotonFieldState` 是壳层级 photon-field contract。joint 模式只更新其 target seed、absorption seed 和 diagnostics，不改变字段语义。

## 7. Gamma-gamma absorption 与 pair cascade

\(\gamma\gamma\) absorption 的基本链条是

\[
\gamma+\gamma
\rightarrow e^+ + e^-
\rightarrow \text{pair synchrotron photons}.
\]

当前实现包含：

- observer-side gamma-gamma attenuation。
- pair-production branch。
- `pair_cascade_iterations > 1` 时的 shell-sequence time-dependent gamma-gamma pair/synch cascade。

当前没有实现 IC-mediated electromagnetic cascade。原因是该过程需要同时闭合 secondary e±、IC photons、new gamma-gamma absorption 和多代 cascade 的 source-sink 方程，不能用经验补项替代。

## 8. 强子物理

### 质子注入与输运

质子由公开 API 字段 `proton_energy_fraction` 控制注入能量预算；下式中的物理符号记为 \(\epsilon_p\)。formal hadronic path 在 log-gamma grid 上推进

\[
\frac{\partial N_p}{\partial R}
=Q_{p,R}
-\frac{\partial}{\partial E_p}
\left(\dot{E}_{p,R}N_p\right)
+Q_{p,{\rm reinj},R}.
\]

损失项包括 adiabatic、proton synchrotron、pγ、BH、pp 等已启用过程。

### Proton synchrotron

高能质子在局域磁场中产生 synchrotron radiation。该分支可作为 observer component 输出，也可作为 hadronic diagnostics。它不自动变成 local target photon field，除非有明确回灌契约。

### Photopion pγ

pγ 过程使用 photon target field 和 proton distribution 计算：

- gamma-ray photons。
- charged pion / muon chain products。
- neutrino luminosity。
- proton energy loss。
- photon depletion / survival。

当前 formal pγ path 支持 Hummer 2010 response 作为 transport feedback path。若 pγ/π/μ 链要把 e± 注入反馈到 electron equation，必须先由 formal kernel 输出谱形和能量预算。

### Bethe-Heitler

BH 过程是 proton-photon pair production：

```text
p + gamma -> p + e+ + e-
```

当前 BH kernel 同时输出：

```text
proton_loss_rate(E_p)
pair_rate(E_e)
photon_loss_rate(epsilon_gamma)
```

joint 模式使用同一个 microphysics operator 的 pair source 和 photon sink，避免凭经验构造 photon depletion。

### pp

pp 过程由 baryon target density、proton spectrum 和 pp kernel 决定。当前可输出 pp secondary pairs 和相关辐射/损失诊断。joint 中只有 formal kernel 已明确输出并归一化的 e± source 进入 `Q_e,secondary,R`。

### Neutrino

`neutrino=True` 输出 neutrino luminosity。中微子逃逸，不反馈到 electron/photon equations。

## 9. 含时二级反馈

`SolverOptions(electron_photon_coupling="joint")` 是 opt-in 壳层级 feedback 模式。其物理目标是把以下对象放在同一 \(R\) 坐标中闭合：

```text
primary electrons
protons
BH / pp / gamma-gamma secondary e±
photon target field
photon survival / absorption
```

所有自然单位为 \({\rm s}^{-1}\) 的 rate 进入求解器前必须换算：

\[
\frac{\mathrm{d}t'}{\mathrm{d}R}
=\frac{1}{\beta\Gamma c}.
\]

joint 当前只支持 forward-shock、1D fullhide electron solver、formal `am3_1d` hadronic path、BH enabled、SSC cooling enabled。详细契约见 `doc/joint_secondary_feedback_physics.md`。

## 10. Reverse-shock 辐射与 cross-zone IC

反向激波辐射包括：

- RS electron synchrotron。
- RS SSC。
- FS/RS cross-zone IC。
- optional RS hadronic light path。
- optional full-chain RS hadronic dispatch through formal 1D kernels。

Cross-zone IC 的物理含义是一个 zone 的 photons 被另一个 zone 的 electrons 散射。它要求两个 zone 的 photon geometry 和 shock state 在同一 observer setup 中组装。

## 11. 偏振

偏振只对同步辐射类分支有物理意义。当前 Stokes projection 覆盖：

- FS electron synchrotron。
- RS electron synchrotron。
- FS/RS hadronic synchrotron。

非同步辐射分支不混入 Stokes。Lan 2023 overlay 中峰值幅度和峰时是两个不同诊断；禁止用 projection-layer time shift 或 smoothing 修正动力学层问题。

## 12. 观测者投影

本地 shell radiation 最终通过 equal-arrival-time surface 投影到观测者：

```text
local luminosity / seed / absorption
-> Doppler transform
-> redshift
-> luminosity distance dilution
-> observer time-frequency grid interpolation
```

`projection_kind="lightcurve"` 是光变和拟合默认路径。`projection_kind="sed"` 是固定时刻扫频率和频段积分默认路径。

`solver_options.geometry_projection="chi_eats_2d"` 只让 FS synchrotron+SSA 的 lightcurve projection 使用 \(\chi\) 分辨有限厚壳层；该 public API 字段在 runtime 配置中对应 `geometry_kernel`。SSC、hadronic、pair cascade 仍是壳层级契约。

## 13. 拟合中的物理含义

Fitting 不改变物理模型。`Fitter` 只是把参数路径映射到 `Model` / `RuntimeConfig`，执行同一个 solve-state 和 projection chain，再用观测数据计算 likelihood 或 redchi。

因此 fit 中发现的非平滑 light curve、非连续 break frequency 或异常 hadronic diagnostics，应回到动力学、输运、源项归一化或投影网格检查，而不是在 likelihood 层修补。

## 14. 不支持边界

当前明确不支持或只部分支持的物理边界包括：

- Jet spreading backend dynamics。
- 用户自定义 `Medium` 的 Fortran kernel dispatch。
- Wind `k != 2`。
- `fullhide_1d` 之外的 thermal electron branch。
- 非轴对称喷流上的 toroidal polarization。
- \(\chi\) 分辨 hadronic transport。
- IC-mediated electromagnetic cascade。
- formal pγ/π/μ e± feedback。

完整准入条件见 `TODO.md` 和 `doc/public_backend_limits.md`。

## 15. 物理验收原则

任何真实物理轨迹若不连续或不光滑，优先视为 bug。例外必须来自明确物理事件，例如 density jump、injection event 或 shock crossing。

验收时重点检查：

- `Gamma(R)`、`B(R)`、`nu_m/nu_c/nu_a` 连续。
- electron/proton spectra 随 shell 平滑演化。
- light curves 和 SED 不出现孤立尖峰。
- photon survival 和 optical depths 单调/连续地响应源项。
- IC、BH、gamma-gamma 的能量预算在可解释误差内闭合。
- weak-feedback 极限回到 separated 或 legacy baseline。

禁止使用经验 smoothing、fallback、后处理时间平移或无物理来源的 floor 来掩盖失败。
