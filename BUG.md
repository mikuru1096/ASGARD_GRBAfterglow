# 未修复缺陷

## γγ/电磁级联的壳层几何与传输所有权未闭合

### 第一性原理

当前 shell-level `annihilation` 与 pair cascade 固定采用共动路径
`R/(12 Gamma)`。该式不是一般 Blandford--McKee 壳宽，只能由 `k=0`
homogeneous slab 的质量守恒推出。若瞬时激波后 proper density 取
`n2'=4 Gamma n1(R)`，则

```text
Msw = 4 pi R^2 Delta' mp n2',
Delta' = Msw/[16 pi R^2 Gamma n1(R) mp].
```

对从原点延伸到当前半径的幂律外介质 `n1=A R^-k`，

```text
Msw = 4 pi mp n1(R) R^3/(3-k),
Delta' = R/[4(3-k) Gamma].
```

所以仅 `k=0` 给出 `Delta'=R/(12 Gamma)`；纯风 `k=2` 给出
`Delta'=R/(4 Gamma)`。对 density jump、tabulated profile 和
transrelativistic 流，`Msw+n1(R)` 至多定义把全部旧物质瞬时压到当前
激波后密度的 one-zone 等效宽度，不能替代演化后的多区体积与空间柱密度。

光深与 transfer 也必须分层。局域 γγ opacity 由目标光子柱密度
`Nnu=nnu Delta'` 决定；共空间、均匀发射与吸收使用
`psi(tau)=(1-exp(-tau))/tau`，前景屏使用 `exp(-tau)`，多区问题需要按
空间顺序积分。不能用一个固定路径或额外常数同时代表这三种几何。

### 受影响的真实调用

- `src/Radiation/pair_absorption.f90::annihilation`：从 `R,Gamma` 猜测
  `R/(12 Gamma)`，并在同一过程内同时组装 opacity 与有限单元 transfer。
- `asgard_core/asgard_state.py`：forward、reverse 和 hadronic target 可被合并后
  共用一个 absorption；pair cascade 开启时局域 target 置零，实际结果完全由
  `tau_extra` 中的 cascade 几何控制。
- `src/Structured/structured_jet_1d.f90::apply_absorption`：已有 `R_mass` 与局域
  density 状态，但仍只向辐射核传 `R,Gamma`；reverse/hadronic 输出复用同一个
  transfer。
- `prompt/radiation.py::compute_branch_radiation`：internal-shock branch 已有
  `comoving_volume_cm3`，仍错误套用外激波的 `R/(12 Gamma)` 路径。
- `src/Hadronic/hadronic_cascade.f90::cascade_seq`：当前存在三套没有共同所有者的
  时钟。pair electron 状态以 `Delta tobs/nsub` 推进，`electron_loss` 使用
  `R/(Gamma c)`，而 `tesc=R/(12 Gamma c)` 控制 `tau_pair`、光子源归一化和
  光子逃逸。

`src/Electron/electron_forward_transport_2d.f90` 通过显式 `dx_comov_hist`
调用 `pair_tau`，不属于这个固定路径缺陷；它说明路径应由几何所有者提供。

### 修复时序

该缺陷必须在强子阶段最后做一次原子修复。若只修改 shell annihilation，
pair-cascade 开启的公开路径仍由旧 `tau_extra` 控制；若只修改 cascade，普通
leptonic、structured 和 prompt 路径仍不一致。禁止在中间阶段加入 `k` switch、
经验系数、clamp、平滑或用 `R,Gamma` 猜测历史相关壳宽。

### 验收条件

- 建立纯函数式 opacity 核：输入目标光子柱密度或显式 density/path，输出
  `tau_gamma_gamma`；核内不读取 density profile，不决定 shell 几何，也不应用
  transfer。
- 几何所有者分别构造 FS、RS、prompt、structured 与 cascade 的光子柱密度；
  支持的 jump/多区路径必须使用演化后的分区体积或 `dx_comov`，没有物理状态的
  组合在公开边界明确拒绝，不能退回 `R/(12 Gamma)`。
- transfer 所有者明确区分共空间 `psi(tau)`、前景 `exp(-tau)` 与有序多区积分；
  FS+RS target 不再无条件合并后套一个共同路径。
- 保留现有 `Radiation.annihilation` 的 public f2py 参数顺序与数组 shape；通过薄包装
  明确保留其 `k=0 homogeneous slab` 语义，真实生产调用使用新的显式几何核心。
- `hadronic_cascade` 的 pair opacity、光子源与逃逸、pair electron 演化和
  `electron_loss` 时钟统一由显式几何/时钟所有者定义，不再各自使用互不相干的
  `Delta tobs`、`R/(Gamma c)` 或固定 `R/(12 Gamma c)`。
- 通过真实构建入口 `pair_absorption`、`structured_jet_1d`、
  `hadronic_forward_1d` 及受影响 source closure 的严格
  `-Wline-truncation -Werror=line-truncation` 检查。
- 直接运行并比较至少 `k=0`、wind、prompt、FS+RS、density jump、pair-cascade
  开/关路径；验证 `k=0` homogeneous 极限、wind 解析比例、prompt 体积路径和
  cascade 光深闭合。光深、频谱与时序曲线不得出现人为跳变，性能比较使用至少
  3 次 median。

## 强子总通量与公开子分量在观测端被重复求和

### 第一性原理

`_hadronlum` 返回的 `fwd_hadronic_gamma` 已经是质子同步、光介子以及所有启用的
Bethe--Heitler、强子逆康普顿和次级粒子辐射之和。公开的
`fwd_hadronic_bethe_heitler`、`fwd_hadronic_inverse_compton` 与
`fwd_hadronic_pair_production` 是该总量的诊断分解，不是额外辐射源。因此总通量只能
包含 `fwd_hadronic_gamma` 一次。

### 当前缺陷

`asgard_core/asgard_state.py::_sumobs` 遍历 `_OBSERVED_COMPONENT_ATTRS`，同时累加
`fwd_hadronic` 总量及其公开子分量。于是 `full_components["total"]` 会重复计入已包含在
强子总量中的通道；`total_only` 的逐分量 batch 只累加 `fwd_hadronic` 总量一次，两种公开模式
在强子开启后仍存在一项独立于投影算法的系统差异。pair cascade 开启时
`pair_lum_total` 还会同时写入总强子量和 pair-production 诊断量，重复同样发生。

### 修复时序与验收

- 与强子/电磁级联最终阶段原子修复，不在当前非强子批量 EATS 改动中重定义公开字段。
- 明确区分“互斥总通量成员”和“总量内部诊断分解”；`total_only` 与
  `full_components["total"]` 只累加前者，仍原样返回所有诊断分量。
- 用至少一个同时启用 Bethe--Heitler、强子逆康普顿和 pair cascade 的真实算例逐项验证
  代数闭合，并检查频谱/光变连续、非负；不得通过删除诊断输出或事后相减掩盖重复计数。

## 强子壳层时间与扫掠质量所有权错误

### 当前缺陷

- `src/Dynamics/Dynamics_forward.f90:55` 输出的 `R_Tobs` 已含红移因子 `(1+z)`，但
  `src/Hadronic/hadronic_base.f90:96-108::shell_dt` 直接取其相邻差；默认 legacy FS 在
  `src/Hadronic/hadronic_forward_1d.f90:249-250`、RS 在
  `src/Hadronic/hadronic_reverse_1d.f90:59-66` 又把该 observer-time 差用于共动粒子冷却。
  对相同的 `R,Gamma` 轨迹，共动步长必须由 `dR/(beta Gamma c)` 决定，不能随红移改变。
- formal 路径虽在 `src/Hadronic/hadronic_shell.f90:724-742::shell_geom` 使用
  `dR/(beta Gamma c)`，但第一个输出点借用 `R(2)-R(1)`，第二个输出点再次使用同一径向间隔，
  因而首个相邻区间被推进两次；初始状态没有独立的时间所有权。
- `asgard_core/asgard_runtime.py:1017-1028::solve_hadronic` 以“几何球壳体积乘当前局域密度”
  重新构造 proton 注入质量，而 dynamics 已提供累计 `swept_mass_g`。纯 wind 的首壳会少算三倍，
  density jump/tabulated profile 还会把区间积分错误替换为端点密度。
  `src/Structured/structured_jet_1d.f90:477-485` 已持有累计 `R_mass`，却重复了同一重积分。

### 修复边界与验收

- 粒子输运统一消费由 `R,Gamma` 构造的共动 proper-time 历史；每个相邻径向区间恰好推进一次，
  第一个输出点只由明确的初始条件定义，不借用未来区间。删除强子输运对 `shell_dt(R_Tobs)` 的依赖，
  保持现有 FS/RS f2py 参数顺序与数组 shape。
- FS 注入质量使用相邻累计扫掠质量差
  `Delta Msw(i)=Msw(i)-Msw(i-1)`；structured 使用 `R_mass` 的同一差分，不再从局域密度重积分。
- 直接比较同一物理 `R,Gamma,B` 轨迹在不同 `z` 下的局域强子解，结果必须不变；比较至少三组
  `Num_R` 验证粒子数、注入能量、冷却谱和强子光度收敛。分别运行 `k=0`、wind、density jump
  真实入口，逐壳验证 `sum(Delta Msw)=Msw` 与 proton 注入能量闭合，并完成受影响模块的
  `-Wline-truncation -Werror=line-truncation` 检查。

## joint secondary pair 的坐标单位与状态/辐射重复所有权

### 当前缺陷

- fullhide 的输运坐标是 `x=ln(gamma)`：`src/Electron/electron_injection_profiles.f90:168-182`
  构造 `ln(gamma)` 边界，`src/Electron/electron_transport_common.f90:969-985` 以
  `dR*dF/dln(gamma)` 更新状态。因此 `dN/dgamma -> dN/dln(gamma)` 只应乘 `gamma`。
  但 `src/Hadronic/hadronic_shell.f90:372::electron_seq` 与
  `asgard_core/asgard_state.py:1189-1196::_sourcer` 都额外乘了 `ln(10)`，使 joint BH/pp/cascade
  secondary source 固定放大 `ln(10)`。
- `src/Hadronic/hadronic_cascade.f90:182-184::cascade_seq` 输出 `pden=cur` 后以 `prev=cur`
  把它作为下一壳累计 pair 状态；`_sourcer` 却把完整 `pden` 除以当前 `Delta R` 当作本壳新注入，
  使旧 pair 在以后每个半径再次进入主电子输运。
- `asgard_core/asgard_state.py:413-430::_jointfeedback` 同时让 cascade 推进 pair、产生 synch seed/lum，
  又把同一 pair 状态注入主电子 synch。它还在 `:426` 写入 `hadronic.l_had_pair_production`；
  `_hadronlum` 已将该量计入强子 aggregate 后，observer 在 `:1008-1020` 重跑 cascade，并在 `:1057`
  再加一次 `pair_lum_total`。固定 `JOINT_ITERS` 的中间迭代态由此被当作额外物理辐射，而不是被最终态替代。

### 修复所有权与验收

- gamma-gamma secondary pair 由 cascade 独占状态推进和 synch 辐射；主电子输运仍只接收 BH/pp
  的真实增量 source。cascade 的累计 `pden` 只作为状态/诊断，不再经过 `_sourcer` 注入主电子。
  最终收敛的 `PairCascade` 随 solve state 保存，observer 直接复用其 luminosity、seed 与 tau，
  不重算 cascade、不累加 Picard 中间态，也不再通过 `hadronic.l_had_pair_production` 复制所有权。
- BH/pp source 按自然对数坐标传递，验证
  `integral dR dln(gamma) Q_R` 等于各壳实际注入 pair 数；禁止用重标定、clamp 或后处理相减修正。
- 对 gamma-gamma 链逐壳验证 absorbed photon power、pair injected power、pair 储能、synch 与绝热损失闭合；
  同一真实算例下 cascade synch 在 `fwd_sync`、强子 aggregate 和 pair diagnostic 中只能有一个物理所有者。
  比较 pair 开/关、joint 迭代收敛、至少三组径向网格与 cascade substep，要求粒子数、频谱、tau
  和时序曲线收敛且连续非负；构建 `hadronic_forward_1d` 与受影响 electron source closure，并执行严格
  line-truncation 检查。
