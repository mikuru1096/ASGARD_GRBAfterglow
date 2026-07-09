# ASGARD TODO

本文档是当前工作树唯一的 TODO / 未完成项入口。新任务必须写清动机、物理假设、验收口径和受影响 build/test 范围；已完成任务不保留在这里，只在提交历史中追踪。

## 原则

- 当前活动任务必须从物理动机、验收口径和受影响 build/test 范围出发定义。
- 不新建 OOP manager、context controller、registry、compose/pipe/Option 层，不为“统一风格”增加薄 wrapper。
- 不新增 fallback、clamp、smoothing、heuristic post-processing 或非边界输入验证。
- public f2py ABI、模块名、入口名、参数顺序和数组 shape 默认不变；若确需收窄 public 面，必须先确认 Python/外部调用。
- 物理公式和数值离散行为保持不变；重构只移动、内联、重排或重命名低价值抽象。
- Fortran 重要改动后必须运行受影响 `build_extensions.py --force`、干净 `/tmp` source closure 的 `-Wline-truncation` 检查和最小相关直接运行或正式 benchmark 入口。
- 不新增测试文件；通过真实构建、直接相关运行和物理连续性检查暴露问题。

## 当前活动任务

- 重新审计 Fortran 主计算链并按 `Dynamics -> Electron -> Radiation/Interpolation -> Hadronic` 顺序优化。动机是删除已证明的重复计算和低价值临时状态，同时提高真实求解速度、数值精度与物理连续性；只实施能由调用图、基准或编译诊断证明收益的改动。
- 验收口径：public f2py ABI 与输出 shape 不变；受影响扩展强制构建和干净 `/tmp` `-Wline-truncation` 闭包通过；直接模型输出有限、非负且动力学/谱随时间连续；每个性能改动用同一输入至少 3 次 median 比较；总代码行数不因薄 wrapper 或无效抽象增加。
- 受影响范围按阶段确定：`Dynamics_forward/Dynamics_reverse`，默认 1D/2D electron targets，`ssc_spectrum/pair_absorption/SED_interpolation`，最后 `hadronic_forward_1d`；不新增测试文件。

## 长期未完成物理边界

这些不是当前重构任务，除非目标观测或物理问题明确需要，否则不进入实现：

- 2D / chi-resolved hadronic transport：formal hadronic path 仍保持 1D shell 契约。
- IC-mediated electromagnetic cascade：当前 pair cascade 只覆盖 gamma-gamma pair/synch contract。
- Formal pγ / π / μ 二级电子谱输出：joint 电子方程只接入 formal kernel 已直接输出且归一化明确的二级 e± 源项。
- Jet spreading、自定义 `Medium` Fortran dispatch、wind `k != 2`、`fullhide_1d`/`fullhide_1d_hz`/`dg_1d` 之外 thermal electrons、非轴对称 toroidal polarization 仍是 unsupported boundary。
- Lan 2023 polarization 峰时偏早问题指向 dynamics/jet-evolution benchmark，禁止 projection 层经验修正。
- FS formal hadronic benchmark refresh 只在假设和决策价值明确时做。
- RS full-chain hadronic 的 `electron_energy_gev must be logarithmically uniform` blocker 需要从输入网格合同证明或修正 kernel 假设；禁止重采样 fallback。

## 不做

- 不删除 `tests/*.npz` baseline、`output/asgard_doc/**` benchmark artifacts 或物理验收图，除非先按 benchmark refresh protocol 证明可复现且无记录价值。
- 不清理 `.venv/`、`.vscode/`、`.codex-remote-attachments/` 等本地目录到 git diff。
- 不把短小纯函数改成类层级。
- 不为两个物理合同不同的路径强行抽象统一。
- 不做无目标驱动的 `RuntimeConfig -> SimulationConfig` 迁移。
