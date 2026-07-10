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

## joint secondary pair 的坐标单位与状态/辐射重复所有权

### 当前缺陷

- fullhide 的输运坐标是 `x=ln(gamma)`：`src/Electron/electron_injection_profiles.f90:168-182`
  构造 `ln(gamma)` 边界，`src/Electron/electron_transport_common.f90:969-985` 以
  `dR*dF/dln(gamma)` 更新状态。因此 `dN/dgamma -> dN/dln(gamma)` 只应乘 `gamma`。
  `src/Hadronic/hadronic_shell.f90::electron_seq` 已按该 Jacobian 输出 BH/pp 径向源；但
  `asgard_core/asgard_state.py:1189-1196::_sourcer` 仍额外乘 `ln(10)`，使 joint cascade
  secondary source 固定放大 `ln(10)`。
- `src/Hadronic/hadronic_cascade.f90:182-184::cascade_seq` 输出 `pden=cur` 后以 `prev=cur`
  把它作为下一壳累计 pair 状态；`_sourcer` 却把完整 `pden` 除以当前 `Delta R` 当作本壳新注入，
  使旧 pair 在以后每个半径再次进入主电子输运。
- `asgard_core/asgard_state.py:413-430::_jointfeedback` 同时让 cascade 推进 pair、产生 synch seed/lum，
  又把同一 pair 状态注入主电子 synch。它还在 `:426` 写入 `hadronic.l_had_pair_production`；
  `_hadronlum` 已将该量计入强子 aggregate 后，observer 在 `:1008-1020` 重跑 cascade，并在 `:1057`
  再加一次 `pair_lum_total`。固定 `JOINT_ITERS` 的中间迭代态由此被当作额外物理辐射，而不是被最终态替代。
- `src/Hadronic/hadronic_pair.f90` 的 absorbed/injected power 离散积分都漏乘能格宽度
  `Delta ln E`；两侧同时漏因子让内部相对闭合检查仍通过，但绝对功率被放大 `1/Delta ln E`。
  一个直接算例中当前值为 `2.7982`，真实离散积分为 `0.64431`，比值恰为 `4.34294=1/dln`。
- `cascade_seq` 推进后的唯一 photon 状态是 `phden`，当前却把 `cphden` 输出写成尚未传播的
  `pseed/h`，并在最后一次 pair operator 之前保存 power 诊断；Python 随后丢弃该输出，再用
  `syn_seed + survival(tau)` 近似重建，重复且不等价于 Fortran 已完成的 transfer。
- joint 路径已把 BH/pp 增量 source 注入主 electron，故其 synch 已包含在
  `ElectronSolution.l_syn_spec/seed_syn`；但 `l_had_bethe_heitler/seed_had_bethe_heitler` 仍再次进入
  `_hadronlum/_hadronseed`。separated 路径已清空这两个字段，joint 尚未执行相同的唯一所有权边界。

### 修复所有权与验收

- gamma-gamma secondary pair 由 cascade 独占状态推进和 synch 辐射；主电子输运仍只接收 BH/pp
  的真实增量 source。cascade 的累计 `pden` 只作为状态/诊断，不再经过 `_sourcer` 注入主电子。
  最终收敛的 `PairCascade` 随 solve state 保存，observer 直接复用其 luminosity、seed 与 tau，
  不重算 cascade、不累加 Picard 中间态，也不再通过 `hadronic.l_had_pair_production` 复制所有权。
- Fortran 直接输出最终传播后的 cascade photon density、最终 power 与 native pair state；下一轮 Picard
  直接消费这一 photon state，不再以 seed 相加和平均 survival 重建 transfer。删除无调用者的第二套
  `cascade_step/cool_deposit`，能量推进复用唯一保守 remap kernel。
- BH/pp source 按自然对数坐标传递，验证
  `integral dR dln(gamma) Q_R` 等于各壳实际注入 pair 数；禁止用重标定、clamp 或后处理相减修正。
- 对 gamma-gamma 链逐壳验证 absorbed photon power、pair injected power、pair 储能、synch 与绝热损失闭合；
  同一真实算例下 cascade synch 在 `fwd_sync`、强子 aggregate 和 pair diagnostic 中只能有一个物理所有者。
  比较 pair 开/关、joint 迭代收敛、至少三组径向网格与 cascade substep，要求粒子数、频谱、tau
  和时序曲线收敛且连续非负；构建 `hadronic_forward_1d` 与受影响 electron source closure，并执行严格
  line-truncation 检查。

## Bethe--Heitler pair kernel 与 proton-loss kernel 能量不闭合

### 当前缺陷

formal 已把 fractional proton loss 正确转换为 `|dγ/dt|=-γ_p ploss`，并把单电荷 `qbh` 以 `2 qbh`
注入代表 electron 与 positron 的主电子谱。直接离散审计中，单电荷 pair 数率与 absorbed photon 数率闭合到机器精度；
但总 pair 功率与 proton loss 功率之比在能格 `32, 48, 64, 96` 上分别为
`0.857143, 0.855777, 0.855861, 0.855711`，不向 1 收敛。

把临时副本的 BH 内外积分阶数从 `5` 提到 `10` 和 `20` 后，该比值仅从 `0.855861` 变为
`0.855589` 和 `0.855507`。因此约 14.4% 的差异不是固定积分阶数误差，而是 differential pair kernel 与
proton-loss `phi` 近似或其积分边界语义不一致。

### 修复边界与验收

- 从同一 BH 微分截面推导 proton energy moment，或证明现有 `phi` 近似与 pair kernel 的适用区间和截断完全一致；
  禁止用全局乘子强制闭合。
- 同时加密 proton、photon、electron 能格并扩大物理积分边界；总 pair 功率与 proton loss 功率误差必须单调下降，
  单电荷数率闭合仍保持机器精度，再运行 formal proton/electron 连续非负检查。
