# ASGARD TODO

本文档是当前工作树唯一的 TODO / 未完成项入口。新任务必须写清动机、物理假设、验收口径和受影响 build/test 范围；已完成任务不保留在这里，只在提交历史中追踪。

## 原则

- 当前主任务是 Fortran 可读性反抽象重构：恢复初始提交 `04f2bf4` 那种从入口顺读物理链的体验，但不回滚当前功能。
- 不新建 OOP manager、context controller、registry、compose/pipe/Option 层，不为“统一风格”增加薄 wrapper。
- 不新增 fallback、clamp、smoothing、heuristic post-processing 或非边界输入验证。
- public f2py ABI、模块名、入口名、参数顺序和数组 shape 默认不变；若确需收窄 public 面，必须先确认 Python/外部调用。
- 物理公式和数值离散行为保持不变；重构只移动、内联、重排或重命名低价值抽象。
- Fortran 重要改动后必须运行受影响 `build_extensions.py --force`、干净 `/tmp` source closure 的 `-Wline-truncation` 检查和最小相关 smoke。
- 不新增测试文件；通过已有 smoke、真实构建和物理连续性检查暴露问题。

## 当前活动任务：Fortran 主链恢复可读性

### 动机

初始提交的 Fortran 主链可读性好，是因为物理阶段线性集中：`Dynamics -> Electron -> Radiation -> Absorption -> Interpolation`。当前 HEAD 功能更完整，但同一条物理链被拆成大量跨文件 helper、thin wrapper、diagnostic wrapper 和局部 stage，读者需要跨文件跳转才能恢复公式和数组演化。

目标不是继续抽象，而是反抽象：让 f2py 主入口按真实物理阶段顺读；低价值薄层内联；模块边界按物理职责，不按抽象类型或历史 wrapper 堆叠。

### 验收口径

- 每个主入口顶部保留短注释，列出真实执行顺序。
- 主入口读下来必须能看到主要物理阶段和状态数组如何演化；最多一次跳转能到核心公式或数值 primitive。
- 只调用一次、只转发参数、名字比被调过程更抽象的 wrapper 必须删除或内联，除非它承担 ABI、单位转换、数组 shape 转换或明确物理边界。
- 反抽象后代码总跳转数减少，正式 driver 不被 diagnostic wrappers、unit helpers 或历史 public 面打断。
- 不改变物理量随时间/半径的连续平滑性；若 smoke 暴露不连续，先查 solver 或投影公式，不做后处理补丁。

## 实施顺序

### 1. Electron 1D：先恢复 `fullhide_1d` 主链

目标文件：
- `src/Electron/electron_forward_fullhide_1d.f90`
- 关联审计：`electron_common.f90`、`electron_injection_profiles.f90`、`electron_shell_transport_common.f90`、`electron_transport_common.f90`

目标阅读结构：

```text
unpack boundary/config
-> build energy coordinate
-> initialize electron distribution
-> loop radius shells
   -> local density / B / gamma_m / gamma_c
   -> injection / cooling coefficients
   -> electron transport
   -> synchrotron / SSA diagnostics
-> return gam_e, dN, luminosity, seed, nu_m/c/a
```

具体 TODO：
- 审计只服务 `fullhide_1d` 的 helper；若 helper 只转发参数或隐藏一两行公式，内联回 `electron_forward_fullhide_1d.f90`。
- `electron_shell_transport_common` 只保留真正可复用的数值离散 primitive；不要把 shell 主推进逻辑藏进去。
- `fullhide_1d` 与 `fullhide_1d_coupled` 共享阶段只用本文件 `contains` 小过程表达，不新建 manager/context。
- 保留 fast uniform-density fixed-substep 路径；不得把已合并的 shell substep 退回逐 substep 慢路径。

验证：
- `build_extensions.py --module electron_forward_fullhide_1d --force`
- `ELECTRON_HISTORY_SOURCES + electron_forward_fullhide_1d.f90` 干净 `/tmp` `gfortran -Wall -Wextra -Wimplicit-interface -ffree-line-length-none -Wline-truncation -fsyntax-only`
- `tests/dg_1d_smoke.py`
- fullhide 相关最小 smoke 或现有 public solve smoke

### 2. Electron 2D：入口保留完整有限厚壳层物理推进

目标文件：
- `src/Electron/electron_forward_transport_2d.f90`
- `src/Electron/electron_transport_2d_kernel.f90`
- `src/Electron/electron_seed_history_kernel.f90`

目标阅读结构：

```text
unpack config
-> construct finite-q shell geometry
-> initialize shock-front electron state
-> loop shells/substeps
   -> update shock state
   -> assemble chi-local cooling/radiation/history fields
   -> advance energy coordinate
   -> advance q advection/diffusion
   -> accumulate shell and chi-resolved spectra
-> reduce to shell-level outputs
-> return chi projection geometry
```

具体 TODO：
- `electron_forward_transport_2d.f90` 保留主流程，不再把 shell/substep 物理推进继续下沉到泛化 helper。
- `electron_transport_2d_kernel.f90` 只保留 q geometry、q advection/diffusion matrix、energy advance primitive。
- 继续清理只为统一命名存在的 q/chi thin wrappers；保留 public procedure 名称和数组 shape。
- `chi_*` 输出继续标注为 observer projection geometry，不是 hadronic chi-local 合同。

验证：
- `build_extensions.py --module electron_forward_transport_2d --force`
- 2D electron ordered source closure 的 `/tmp` line-truncation 检查
- `tests/fullhide_2d_smoke_bench.py`
- `tests/structured_chi_2d_smoke.py`

### 3. Hadronic：隔离 diagnostic ABI，正式 shell sequence 顺读

目标文件：
- `src/Hadronic/hadronic_forward_1d.f90`
- `src/Hadronic/hadronic_forward_formal_1d.f90`
- `build_extensions.py`

目标阅读结构：

```text
diagnostic f2py wrappers
-> fs_hadronic_1d light driver
-> fs_hadronic_formal_transport_1d ABI driver
-> formal shell sequence implementation
   -> proton injection/transport
   -> pγ / BH / pp
   -> secondary species
   -> radiation
   -> photon survival
   -> secondary e± source
```

具体 TODO：
- 将 diagnostic wrappers 与正式 driver 分离；可拆成单独 wrapper 源文件，但必须保持同一个 f2py module 的 public ABI，除非先确认收窄 public 面。
- `hadronic_forward_1d.f90` 不再被 unit/grid helper 和 diagnostic wrapper 打断正式主链。
- `hadronic_forward_formal_1d.f90` 保持 formal shell sequence 内部可顺读；不把过程再拆成无物理含义的 operator 层。
- 不删除底层物理 kernel；只整理入口和 wrapper 边界。

验证：
- `build_extensions.py --module hadronic_forward_1d --force`
- hadronic ordered source closure `/tmp` line-truncation 检查
- `tests/hadronic_1d_smoke.py`
- `tests/electron_photon_joint_secondary_feedback_smoke.py`

### 4. Dynamics_reverse：按物理块切分长文件

目标文件：
- `src/Dynamics/Dynamics_reverse.f90`
- 可新增物理块文件：`reverse_secondary_events.f90`、`reverse_jump_conditions.f90`、`reverse_rhs.f90`

目标阅读结构：

```text
resolve RS start/crossing
-> integrate main RS dynamics
-> update secondary reverse events
-> compute hydro/MHD jump state
-> output region-3 thermal/magnetic records
```

具体 TODO：
- 主入口只保留 reverse shock 时间推进和输出组装。
- secondary reverse event start/end/source 放入独立物理块。
- RH/MHD jump、`sigma -> 0` hydrodynamic limit 放入独立物理块。
- RK RHS 与状态导数放入独立物理块。
- 不改变 `U3/V3` thermal state、`gamma34` 注入合同、有限 upstream `sigma` 惯性和 ordered+turbulent `B3` 合同。

验证：
- `build_extensions.py --module Dynamics_forward --module Dynamics_reverse --force`
- `Constants.f90 + dynamics_common.f90 + Dynamics_reverse.f90 + 新增物理块` `/tmp` line-truncation 检查
- `tests/reverse_shared_solver_smoke.py`
- `tests/reverse_shock_smoke.py`
- `tests/hadronic_reverse_shock_smoke.py`
- 检查 radius、Gamma、M3、B3、U3/V3 随半径有限、连续、平滑

### 5. SED / Projection：projection 类型分文件，投影文件不放 radiation 公式

目标文件：
- `src/Interpolation/SED_interpolation.f90`
- 可拆分：`SED_interpolation_shell.f90`、`SED_interpolation_adaptive_theta.f90`、`SED_interpolation_chi.f90`、`SED_interpolation_structured_chi.f90`

目标阅读结构：

```text
shell-level projection
adaptive theta projection
top-hat chi projection
direct-electron chi projection
structured precomputed-ring chi projection
```

具体 TODO：
- 每个 projection entry 内部顺序固定为 EATS geometry -> Doppler -> frequency interpolation -> SSA survival -> accumulation。
- `chi_synch_point` 这类 radiation 公式移回 radiation kernel 或删除重复路径；projection 文件只调用 radiation kernel。
- 不把 geometry segment 计算继续泛化成难读 helper；重复少量几何公式可以接受。
- `radiation_syn_kernel_value` 继续作为同步核唯一中心。

验证：
- `build_extensions.py --module SED_interpolation --force`
- `Constants.f90 + radiation_common.f90 + interpolation_common.f90 + SED_interpolation*.f90` `/tmp` line-truncation 检查
- `tests/fullhide_2d_smoke_bench.py`
- `tests/structured_chi_2d_smoke.py`

### 6. Radiation / Electron kernel：按公式族重排，不继续泛化

目标文件：
- `src/Electron/electron_radiation_kernel.f90`
- `src/Radiation/radiation_common.f90`
- 相关 cooling/history kernel

具体 TODO：
- 将同步 emissivity / SSA tau / break frequencies / chi batch radiation 按公式族分段或拆文件。
- 不复制同步核公式，不在 SED、Electron、Structured 中各自实现局部版本。
- 只保留真实复用的公式入口；单调用薄 wrapper 内联。
- 清理历史遗留 `implicit real` 时只做局部必要声明，不做大规模 style migration。

验证：
- 受影响 electron/radiation build
- 对应 ordered source closure `/tmp` line-truncation 检查
- `tests/fullhide_2d_smoke_bench.py`
- `tests/dg_1d_smoke.py`

## 长期未完成物理边界

这些不是当前重构任务，除非目标观测或物理问题明确需要，否则不进入实现：

- 2D / chi-resolved hadronic transport：formal hadronic path 仍保持 1D shell 契约。
- IC-mediated electromagnetic cascade：当前 pair cascade 只覆盖 gamma-gamma pair/synch contract。
- Formal pγ / π / μ 二级电子谱输出：joint 电子方程只接入 formal kernel 已直接输出且归一化明确的二级 e± 源项。
- Jet spreading、自定义 `Medium` Fortran dispatch、wind `k != 2`、`fullhide_1d` 之外 thermal electrons、非轴对称 toroidal polarization 仍是 unsupported boundary。
- Lan 2023 polarization 峰时偏早问题指向 dynamics/jet-evolution benchmark，禁止 projection 层经验修正。
- FS formal hadronic benchmark refresh 只在假设和决策价值明确时做。
- RS full-chain hadronic 的 `electron_energy_gev must be logarithmically uniform` blocker 需要从输入网格合同证明或修正 kernel 假设；禁止重采样 fallback。

## 不做

- 不删除 `tests/*.npz` baseline、`output/asgard_doc/**` benchmark artifacts 或物理验收图，除非先按 benchmark refresh protocol 证明可复现且无记录价值。
- 不清理 `.venv/`、`.vscode/`、`.codex-remote-attachments/` 等本地目录到 git diff。
- 不把短小纯函数改成类层级。
- 不为两个物理合同不同的路径强行抽象统一。
- 不做无目标驱动的 `RuntimeConfig -> SimulationConfig` 迁移。
