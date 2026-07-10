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

## prompt internal-shock 入口导入已删除的反向激波私有函数

`prompt/internal_shock.py` 在模块导入时仍从 `asgard_core.asgard_runtime` 请求
`_rs_fast_shock_allowed`、`_rs_shock_upstream_u` 和 `_rs_vegas_downstream_u`，但当前
运行时已不再定义或导出这三个名字；仓库内也没有其他实现。因而任何 prompt 直接入口都会在
物理计算前以 `ImportError` 终止，当前无法执行 prompt fixed/adaptive EATS 回归。

修复不能重新复制一套私有 Python jump 公式或添加兼容兜底。应先确认当前反向激波跳跃条件的
唯一物理所有者，把 prompt 两壳碰撞显式接到同一纯函数核，并用磁化与非磁化 FS/RS 两支真实
算例验证压缩比、下游四速度、过激波判据和 EATS 光变连续性；完成后删除本条。
