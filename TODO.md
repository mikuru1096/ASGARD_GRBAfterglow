# ASGARD TODO

本文档是当前工作树唯一的 TODO / 未完成项入口。其他文档只描述当前能力、历史决策或 backend 边界；新任务进入前先在这里补充动机、物理假设、验收口径和受影响 build/test 范围。

## 原则

- 不新建无必要工具模块，不引入 compose / pipe / registry / Option 等抽象层。
- 不恢复旧 public alias，不新增 wrapper、fallback、compat shim。
- 物理学内核代码不因清理任务改动；如后续触及 Fortran，必须做受影响编译和 `-Wline-truncation` 检查。
- `build_extensions.py` 的 `convert_utf8_to_ascii` 源码回写按当前项目决策保留，清理任务不得修改。

## 当前活动任务

当前没有活动任务。

## 当前未完成边界

这些条目不是可随手实现的 backlog。只有在目标观测或物理问题明确需要，并且先写清契约与验收口径后，才允许进入实现。

### 1. 2D / chi-resolved hadronic transport

当前 formal hadronic path 保持 1D shell 契约；hadronic 暂不实现 2D / chi-resolved transport，直到 chi-local photon field、hadron density、secondary feedback 和 hadronic observer projection 的物理契约完成。决策记录见 `doc/hadronic_chi_transport_decision.md`。

### 2. IC-mediated electromagnetic cascade

超出当前 gamma-gamma pair/synch contract 的 inverse-Compton-mediated electromagnetic cascade 暂不实现，直到 photon/e± source-sink 方程、IC kernel 契约和 energy-budget benchmark 完成。边界见 `doc/pair_cascade_extension_boundary.md`。

### 3. Formal pγ / π / μ 二级电子谱输出

joint 电子方程只接入 formal kernel 已直接输出且归一化明确的二级 e± 源项。若 pγ/π/μ 链需要反馈到电子方程，必须先在 formal hadronic kernel 中提供 e± 注入谱及其能量预算 benchmark；禁止用总能量守恒外推临时构造谱形。

### 4. Public/backend unsupported boundaries

以下 public API 或配置入口已暴露但 backend 明确不支持或只部分支持，不能静默 fallback：

- Jet spreading backend dynamics。
- 用户自定义 `Medium` 的 Fortran kernel dispatch。
- Wind `k != 2`。
- `fullhide_1d` 之外的 thermal electron branch。
- 非轴对称喷流上的 toroidal polarization。

完整边界和实现准入条件见 `doc/public_backend_limits.md`。

### 5. Polarization timing diagnostic

Lan 2023 overlay 的峰时偏早问题应进入 dynamics/jet-evolution benchmark，而不是 surface-element EATS 或 patch solid-angle 权重。禁止在 polarization projection 层使用经验 time shift、smoothing 或投影层补丁修正。

### 6. FS formal hadronic benchmark refresh

含 AM3 对照或 hadronic-dominated scenario 的 FS formal hadronic benchmark figures 只有在目标明确时才单独刷新。新增正式 benchmark 必须先说明假设、决策价值、受影响 Fortran build 状态、输出路径和物理验收口径。

### 7. Formal hadronic validation blockers

当前未完成项不是迁移到 Fortran，而是修复或重新验证 formal hadronic validation blocker：

- full-chain RS hadronic 分支报 `electron_energy_gev must be logarithmically uniform`。

修复契约：

- 不添加 fallback、网格重采样补丁或忽略错误；必须从输入 `electron_gamma` 到 `electron_energy_gev` 的构造路径证明 Fortran kernel 看到的是严格等比能量网格，或者修正 kernel 对当前契约的错误假设。
- RS full-chain 必须继续使用 RS seed photons、RS magnetic field、RS shell energy 和 RS baryon target density。
- 若再触及 joint feedback，必须保持 shell-level electron/photon/hadronic 闭合，不恢复 Python substep 推进。
- 修复后 FS 1D hadronic、RS full-chain hadronic 和相关 joint feedback 端到端验证必须全部通过。
- 与迁移前同一输入的 proton synch、pγ gamma、neutrino、BH/pp secondary luminosity 和 energy budget 在数值误差内一致，且随半径连续平滑。
- 受影响 `hadronic_forward_1d` / `structured_jet_1d` build、干净 `/tmp` source closure `-Wline-truncation` 检查通过。

## 不做

- 不删除 `tests/*.npz` baseline、`output/asgard_doc/**` benchmark artifacts 或文献/物理验收图，除非先按 benchmark refresh protocol 证明可复现且无记录价值。
- 不清理 `.venv/`、`.vscode/`、`.codex-remote-attachments/` 等本地目录到 git diff。
- 不把短小的纯函数改成类层级，也不为两个不同物理契约强行抽象统一。
- 不做无目标驱动的 `RuntimeConfig -> SimulationConfig` 主链迁移；`RuntimeConfig` 仍是 runtime、state、postprocess、tests 和 scripts 的主输入类型。
