# 对产生级联扩展边界

本文记录当前是否应把 pair cascade 从 shell-sequence γγ pair/synch 扩展到 inverse-Compton-mediated electromagnetic cascade。结论是：**暂不扩展**。当前代码已经有一个可检验的 γγ pair/synch 合同；IC-mediated cascade 需要新的 photon/e± source-sink 系统，不能作为现有 synch branch 的后处理项添加。

## 当前契约

- `asgard_core/hadronic_cascade.py` 提供两条路径：
  - `compute_time_dependent_pair_cascade_sequence(...)`：主链使用的 shell-sequence time-dependent γγ pair/synch cascade。
  - `compute_iterative_pair_cascade(...)`：legacy single-shell iterative synch branch，仅保留作诊断。
- `_compute_pair_production_branch(...)` 在 `pair_cascade_iterations > 1` 时调用 shell-sequence path，返回 `pair_syn_luminosity_hz`、`pair_syn_seed_per_hz` 和 `tau_pair_path`。
- Fortran `hadronic_pair_cascade_kernel.f90` 的单步物理是 `γγ -> e± -> synchrotron photons`。它不演化 IC photons，也不把 IC photons 再反馈进 γγ target field。
- 当前 public contract 中，`Radiation.pair_production=True` 加 `Setups.pair_cascade_iterations>1` 的含义是 γγ pair/synch cascade substeps，不是完整 electromagnetic cascade。

## 物理问题

IC-mediated electromagnetic cascade 的闭合变量不同于当前 synch branch。需要同时演化 photon density、pair density、synch source、IC source、γγ absorption、escape 和 Klein-Nishina 截面。若只把 IC emissivity 加到当前 pair synch luminosity 上，会出现两个问题：

1. IC photons 的 source-sink 没有回写到 target photon field，不能决定下一步 γγ pair injection。
2. pair energy loss 被 synch 和 IC 共同分配时，必须在同一时间步内守恒；不能先按 synch-only 冷却算电子谱，再额外叠加 IC luminosity。

因此 IC-mediated cascade 必须是新的耦合方程，而不是当前分支的输出修饰。

## 当前决策

当前主链保持：

- γγ absorption / pair injection 来自正式 pair-production operator。
- `pair_cascade_iterations > 1` 使用 shell-sequence time-dependent γγ pair/synch cascade。
- single-shell iterative cascade 只作低层诊断。
- IC-mediated electromagnetic cascade 继续列为未实现边界。

进入实现前必须先完成：

1. 写出 photon/e± source-sink 方程，包含 synch、IC、γγ、escape 和 pair cooling。
2. 明确 IC kernel：Thomson/KN 处理、target photon field、输出 photon grid 与 `PhotonFieldState` 的耦合。
3. 定义能量预算：每个 shell 的 absorbed photon power、pair internal energy、synch luminosity、IC luminosity、escaping photon power 必须在数值误差内闭合。
4. 给出最小 benchmark：synch-dominated limit 回到当前 γγ pair/synch path；IC-dominated optically thin limit 的谱和能量分配随 radius 平滑。
5. 明确 observer projection：IC cascade photons 是 hadronic component、absorption seed 还是新的 emission component，不能在 flux assembly 里隐式混合。

这些条件满足前，扩展 IC-mediated cascade 不会提高物理可信度。
