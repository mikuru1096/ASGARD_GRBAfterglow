# ASGARD TODO

本文档是当前工作树唯一的 TODO / 未完成项入口。其他文档只描述当前能力、历史决策或 backend 边界；新任务进入前先在这里补充动机、物理假设、验收口径和受影响 build/test 范围。

已完成或已过期的旧 public alias、旧观测 builder、重复 enum member、冗余 property、旧 kernel alias 表、零调用 hadronic alias、physics helper 薄 wrapper、旧 AM3 process-label alias、pgamma 两阶段查找和偏振模块阶段计划不再作为待办。

## 原则

- 不新建无必要工具模块，不引入 compose / pipe / registry / Option 等抽象层。
- 不恢复旧 public alias，不新增 wrapper、fallback、compat shim。
- 物理学内核代码不因清理任务改动；如后续触及 Fortran，必须做受影响编译和 `-Wline-truncation` 检查。
- `build_extensions.py` 的 `convert_utf8_to_ascii` 源码回写按当前项目决策保留，清理任务不得修改。

## 当前活动重构计划

### Fortran kernel readability refactor

目标：按 Fortran 技能的顶层路线收窄声明块和主流程复杂度，只做行为保持的作用域、命名 helper、contained stage 或 workspace 重构。不得改 public/f2py 入口、数组形状、物理公式、数值策略或 wrapper 调用契约。

当前未提交批次已覆盖：

- f2py 临时签名源处理：`build_extensions.py` 仅在签名扫描副本中剥离 `block/associate` 结构行，真实 gfortran 构建仍使用原始 Fortran。
- 2D electron / PIC / transport common：谱峰、active support、substep 系数扫描、η 推进矩阵阶段已命名化。
- 1D electron drivers：`fullhide_1d`、`fullhide_1d_hybrid`、`charint_1d` 的壳层准备、hybrid/source 构造和末端诊断已收成 contained helper；输运主体仍留在父过程。
- electron radiation/cooling/hybrid：同步辐射单频点、SSA/IC/Y 积分和 hybrid spectrum point 已收成局部阶段 helper。
- structured / interpolation / radiation：structured patch solve、χ-resolved 和 structured EATS segment projection、SSC、γγ absorption、reverse seed、radiation common seed core 已完成局部阶段命名。
- dynamics：reverse RK event trial-state、reverse RHS 的状态解码/区域 2/区域 3 场与辐射效率阶段已命名化。
- hadronic：pγ Hummer、interaction、forward driver、pp delta、pp models、acceleration、Bethe-Heitler、decay、pair cascade/production、secondary radiation、species transport、hadronic IC、radiation、transport/remap 已完成 contained helper 或小型共享 helper 重构。

剩余不动边界：

- `electron_forward_weno5_1d.f90` 含 NaN/非负裁剪等数值策略，需单独物理/数值审计，不在本轮“只读性重构”中用 helper 包起来。
- `electron_forward_slc1_1d.f90`、`electron_forward_t2g1_1d.f90`、`Dynamics_forward.f90`、`hadronic_common.f90`、`interpolation_common.f90`、`quantum_synchrotron_kernel.f90`、`synchrotron_polarization_kernel.f90` 等文件短且职责单一，继续拆分的边际信息增益低。
- `adaptive_resampling_mod.f90` 已有内部 helper；若未来要改采样策略或模块私有过程边界，需另开数值审计，不能夹在行为保持重构中。

验收口径：每个批次必须通过受影响 `build_extensions.py --force`、干净 `/tmp` module 目录的 `gfortran -Wall -Werror=line-truncation -Wline-truncation -cpp -fopenmp -fsyntax-only` source closure 检查、最小相关 smoke test；每累计约 5 个代码文件或高风险批次后运行 CCreview。最终完成标记只记录提交哈希。

## 当前未完成边界

这些条目不是可随手实现的 backlog。只有在目标观测或物理问题明确需要，并且先写清契约与验收口径后，才允许进入实现。

### 1. 2D / chi-resolved hadronic transport

当前 formal hadronic path 保持 1D shell 契约。FS synchrotron + SSA 已有 `geometry_kernel="chi_eats_2d"` 的 chi-resolved observer projection；hadronic 仍暂不实现 2D / chi-resolved transport，直到 chi-local photon field、hadron density、secondary feedback 和 hadronic observer projection 的物理契约完成。决策记录见 `doc/hadronic_chi_transport_decision.md`。

### 2. IC-mediated electromagnetic cascade

当前 pair branch 停在 shell-sequence time-dependent gamma-gamma pair/synch cascade。超出当前 gamma-gamma pair/synch contract 的 inverse-Compton-mediated electromagnetic cascade 暂不实现，直到 photon/e± source-sink 方程、IC kernel 契约和 energy-budget benchmark 完成。边界见 `doc/pair_cascade_extension_boundary.md`。

### 3. Public/backend unsupported boundaries

以下 public API 或配置入口已暴露但 backend 明确不支持或只部分支持，不能静默 fallback：

- Jet spreading backend dynamics。
- 用户自定义 `Medium` 的 Fortran kernel dispatch。
- Wind `k != 2`。
- `fullhide_1d` 之外的 thermal electron branch。
- 非轴对称喷流上的 toroidal polarization。

完整边界和实现准入条件见 `doc/public_backend_limits.md`。

### 4. Polarization timing diagnostic

Lan 2023 overlay 的峰值幅度已匹配，峰时仍偏早。当前证据指向 dynamics/jet-evolution benchmark，而不是 surface-element EATS 或 patch solid-angle 权重。禁止在 polarization projection 层使用经验 time shift、smoothing 或投影层补丁修正。诊断记录见 `doc/polarization_timing_diagnostic.md`。

### 5. FS formal hadronic benchmark refresh

baseline Vegas comparison 已按 `doc/benchmark_refresh_protocol.md` 全量刷新；含 AM3 对照或 hadronic-dominated scenario 的 FS formal hadronic benchmark figures 仍需在目标明确时单独刷新。刷新前后必须记录 HEAD、tracked diff、完整命令、受影响 Fortran build 状态、输出路径和物理验收口径。

## 不做

- 不删除 `tests/*.npz` baseline、`output/asgard_doc/**` benchmark artifacts 或文献/物理验收图，除非先按 benchmark refresh protocol 证明可复现且无记录价值。
- 不清理 `.venv/`、`.vscode/`、`.codex-remote-attachments/` 等本地目录到 git diff。
- 不把短小的纯函数改成类层级，也不为两个不同物理契约强行抽象统一。
- 不做无目标驱动的 `FitConfig -> SimulationConfig` 主链迁移；`FitConfig` 仍是 runtime、state、postprocess、tests 和 scripts 的主输入类型。
- 不做 `ISM`、`Wind`、`TophatJet` 等 public constructor alias 的破坏性移除；这些别名是当前文档化公开入口。
