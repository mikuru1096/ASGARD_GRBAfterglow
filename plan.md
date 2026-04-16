# ASGARD 推断链与新用户入口清单

## 1. 当前已经固定的事实

- 生产推断快路径已经接入，不再依赖旧的逐数据集重复求值。
- 当前默认推断链：
  - `ASGARD/api.py::observe`
  - `asgard_fit.py`
  - `asgard_state.py` 的 `mode="total_only"`
- `Fitter.loglike()`、`mcmc_fit.py`、`multinest_fit.py` 都已经切到 compiled fast-path。
- 新用户入口已经收束到少数脚本：
  - `lc_spec_demo.py`
  - `tests/readme_smoke_bench.py`
  - `tests/slc1_vs_fullhide_bench.py`
  - `tests/wind_vs_fullhide_bench.py`
  - `tests/doc_figures.py`
  - `tests/ic_doc_plots.py`
  - `tests/sed_electron_compare.py`

## 2. 已经做完，后续不再重复做

- 不再新增单独的 inference profile / check 脚本。
- 不再恢复旧的深层验证测试树。
- 不再为生产推断保留逐分量 observer 插值。
- 不再为生产推断保留 `deepcopy(model)` 的逐点评估路径。
- 不再引入新的中间包装层，除非它能替代至少两处现有重复代码。
- 不再把诊断图、验证图默认挂到新用户主入口里。

## 3. 当前代码结构

1. 配置层：`asgard_models.py`
2. 运行时绑定层：`asgard_runtime.py`
3. 观测与后处理层：`ASGARD/api.py`、`asgard_postprocess.py`
4. 采样与推断层：`asgard_fit.py`、`mcmc_fit.py`、`multinest_fit.py`
5. 最小示例层：`lc_spec_demo.py`
6. benchmark 与绘图层：`tests/`

## 4. 后续只保留的工作

### 4.1 推断链

- 继续维持 compiled fast-path 的简单结构。
- 如果要改推断逻辑，优先改 `asgard_fit.py`，不要再分裂出新入口。
- 如果观测块的复用策略要调整，只允许在当前结构内改，不回到旧的多层调用链。

### 4.2 新用户入口

- `lc_spec_demo.py` 保持为最短可运行 demo。
- `tests/` 只保留两类入口：
  - benchmark
  - plot
- 新增脚本前先问自己一件事：
  - 它是否直接帮助新用户理解“如何跑”和“如何看结果”。
  - 如果不是，就不要加。

### 4.3 文档

- `README.md` 只保留最短上手路径。
- `AGENTS.md` 记录当前事实，不记录已经删除的旧入口为主流程。
- 若代码结构变化，先更新 `AGENTS.md`，再更新 `plan.md`。

## 5. 执行顺序

1. 优先保持当前主入口可运行。
2. 只做和 benchmark / plot / 新用户入口直接相关的改动。
3. 修改推断链时，先看是否能复用现有 compiled fast-path。
4. 任何会增加结构深度的改动，默认先拒绝。

## 6. 明确不做

- 不恢复旧的验证脚本树。
- 不恢复旧的 profile/check 脚本。
- 不把内部诊断工具包装成默认新用户入口。
- 不为未确认的收益新增抽象层。
- 不把 `tests/` 再扩成“所有诊断都放进去”的集合。

## 7. 已记录的待修复物理问题

- **高优先：charint 高频电子端异常抬起（fast 测试集中）**
  - 表现：在快网格测试下，`electron_solver="charint"` 在电子高能端或等效高频辐射端出现非物理抬起。
  - 现判：更像 Fortran 核层问题，优先排除接口/外层采样误差。
  - 暂行要求（本周内）：
    - 定位 `src/Electron/FS_electron_charint.f90` + `src/Electron/electron_common.f90` 中高能前沿/重映射路径（含 `electron_find_high_energy_front`、`electron_anchor_high_energy_edges`、可变步长步进逻辑）。
    - 在不改物理公式前提下，只做最小核逻辑修复。
    - 加入一个回归检查（只在快测试或专用回归中触发），确保高能尾不再出现单调性反转。
