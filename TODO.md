# ASGARD 待优化项

本文档只记录当前工作树仍成立的清理候选。已完成或已过期的旧 public alias、旧观测 builder、重复 enum member、冗余 property、旧 kernel alias 表、零调用 hadronic alias、physics helper 薄 wrapper、旧 AM3 process-label alias 和 pgamma 两阶段查找不再作为待办。

## 原则

- 不新建无必要工具模块，不引入 compose / pipe / registry / Option 等抽象层。
- 不恢复旧 public alias，不新增 wrapper、fallback、compat shim。
- 物理学内核代码不因清理任务改动；如后续触及 Fortran，必须做受影响编译和 `-Wline-truncation` 检查。
- `build_extensions.py` 的 `convert_utf8_to_ascii` 源码回写按当前项目决策保留，清理任务不得修改。

## 已知候选

### 1. 2D electron Boundary 解包重复

`src/Electron/electron_transport_2d_kernel.f90` 与 `src/Electron/electron_seed_history_kernel.f90` 可能仍有可合并的 Boundary 解包逻辑。执行前必须先证明重复代码完全等价，并确认不会改变 `fullhide_2d` 与 `charint_2d` 的物理路径。

### 2. 2D electron extension 文档持续同步

当前源码没有独立的 charint 2D Fortran 源文件。`electron_forward_charint_2d` 是 f2py extension 名称，由 `src/Electron/electron_forward_transport_2d.f90` 中的 `fs_electron_transport_2d_core` 构建，并通过 `use_charint_transport` 启用 charint 2D path。后续文档、测试和构建说明必须保持这个事实。

### 3. `FitConfig -> SimulationConfig` 主链迁移

`FitConfig` 仍是当前 runtime、state、postprocess、tests 和 scripts 的主输入类型。不要只因为注释里写了 legacy 就删除它；只有当 `SimulationConfig` 贯穿运行主链并完成等价 smoke/benchmark 后，才能进入破坏性迁移。

### 4. Public constructor alias 破坏性移除

`ISM`、`Wind`、`TophatJet` 等 constructor aliases 是当前文档化公开入口。若后续要移除，必须作为 public API breaking change 单独处理，并同步 README、public API 文档、示例和 tests。

## 不做

- 不删除 `tests/*.npz` baseline、`output/asgard_doc/**` benchmark artifacts 或文献/物理验收图，除非先按 benchmark refresh protocol 证明可复现且无记录价值。
- 不清理 `.venv/`、`.vscode/`、`.codex-remote-attachments/` 等本地目录到 git diff。
- 不把短小的纯函数改成类层级，也不为两个不同物理契约强行抽象统一。
