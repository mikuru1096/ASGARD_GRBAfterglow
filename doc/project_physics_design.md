# 项目物理设计导航

本页只回答“物理定义在哪里”。完整公式与支持边界统一放在 [物理模型](physics_model.md)，避免多份总纲漂移。

## 主计算链

| 阶段 | 物理内容 | 进一步阅读 |
|---|---|---|
| 环境与动力学 | 外部介质、结构化喷流、FS/RS、密度跳变 | [物理模型](physics_model.md#3) |
| 轻子输运 | 注入、同步/IC/绝热冷却、有限壳层 | [电子求解器](electron_solver_algorithms.md) |
| 辐射 | 同步、SSA、SSC、gamma-gamma | [物理模型](physics_model.md#8) |
| 强子 | proton synchrotron、pγ、BH、pp、neutrino | [物理模型](physics_model.md#12) |
| 二级反馈 | electron–photon–hadron 含时闭合 | [二级反馈物理](joint_secondary_feedback_physics.md) |
| 投影 | EATS、Doppler、红移、偏振 | [物理模型](physics_model.md#22) |

## 专题页

- [磁化 RS 与 DG1D](magnetized_rs_dg1d_tutorial.md)：MHD jump、crossing 和高阶电子输运。
- [Prompt 内部激波](prompt_internal_shock_tutorial.md)：两壳快照诊断；不是稳定 public API。
- [有限厚壳层输运](fullhide2d_pwn_cr_transport.md)：q-mass 网格与 chi 诊断坐标。
- [代码概览](code_overview.md)：物理阶段对应的 Python/Fortran 模块。
- [后端能力边界](public_backend_limits.md)：公开组合与未实现组合。

## 关键所有权

- 动力学产生壳层状态，不产生观测通量。
- electron kernel 拥有电子谱及同步/SSA 局域量。
- hadronic kernel 拥有质子损失及强子产物。
- `joint` solver 拥有含时 pair 反馈；observer 不重复注入。
- observer 只负责分量组装与 EATS 投影。

## 绑定而非独立的选择

- `chi_eats_2d` 只投影 `fullhide_2d` 的 FS 同步与 SSA。
- pp gamma model 只替换 π0 gamma；默认 `delta`，详细 AM3 模型显式启用。
- `include_ssc` 与 IC 冷却近似含义不同，不能相互代替。
- Fan–Y 尚未公开，需先完成原文推导与内核修正。

## 验收原则

物理改动必须满足单位、粒子数和能量闭合；过程关闭时旧路径不变；阈值与分支符合原始文献；真实时空量连续平滑。禁止 fallback、clamp、经验平滑或事后重归一化。
