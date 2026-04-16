# ASGARD 当前主线计划

## 0. 使用规则

- 这份文件只保留当前仍有效的主线任务、约束和验收口径。
- 已完成的离线诊断、已放弃路线、已降级事项不再长期堆在这里，改由 `AGENTS.md` 留痕。
- 后续只要主线目标、验收标准或默认工作流发生变化，都必须同步更新本文件。

## 1. 当前稳定基线

- 默认稳定电子求解器仍是 `fullhide`。
- 当前还不能把实验链更新为新的全项目主基准。
- 稳定链的基础回归门槛仍是：
  - `python tests/regression_check.py`
  - `python tests/ASGARD_comprehensive_validation.py`
  - `python tests/bug_audit_check.py`
  - `python tests/physical_closure_check.py`

## 2. 当前主线任务

### 2.1 `charint` 继续作为实验主线

- `charint` 已经接入为新 solver，公开名称固定为 `electron_solver="charint"`。
- 当前主线目标不再是继续追 `peak-gamma`、低能前沿结构或 `g_lo` 跳变。
- 当前主线目标改为：
  - 在不改坏正式算法口径的前提下，继续降低 `charint` 的单线程计算开销；
  - 尤其减少 `L1` 子步内重复预处理、重复 tracing、重复 remap、以及 `exp/log/pow` 开销；
  - 保持粗阶数精度不退化。

### 2.2 `slc1`、`fullhide`、`weno5` 当前不作为主攻

- `slc1` 当前不再继续追 `peak-gamma` / 低能前沿问题。
- `fullhide` 和 `weno5` 当前不作为主线重构目标。
- 这些链路的历史问题和诊断结论保留在 `AGENTS.md`，不再占用当前计划主区。

## 3. 当前有效验收口径

### 3.1 `charint` 的正式 done-when

- 速度目标：
  - 正式 benchmark 上比 `fullhide` 快 `10x`
  - 正式 benchmark 上比 `slc1` 快 `5x`
- 精度目标：
  - 至少与当前基线相同，或更高
  - 不接受靠牺牲连续性、平滑性或整体物理一致性换速度

### 3.2 当前不再作为主线验收项

- `peak-gamma`
- 低能前沿结构
- `g_lo` 跳变
- 细碎的峰位锁格诊断

这些量可以保留历史记录，但不再驱动当前开发优先级。

### 3.3 当前保留的诊断口径

- `tests/order_convergence_check.py` 只保留粗阶数检查：
  - `electron-spectrum`
  - `radiation-bands-flux`
  - `dynamics-forward-gamma/radius`
- 理论口径：
  - `charint` 按整体二阶算法看待；
  - 平滑区 remap 可表现出更高阶风格，但不作为全局 formal order 声称。

## 4. 当前 benchmark 口径

- benchmark 默认只比较：
  - `fullhide`
  - `slc1`
  - `charint`
- 不再测 `weno5`。
- 当前电子网格口径固定为：
  - `fullhide.num_gam_e = 121`
  - `slc1.num_gam_e = 41`
  - `charint.num_gam_e = 41`
- benchmark 很耗时，默认尽量少做：
  - 平时优先做轻量 one-shot timing；
  - 只有在阶段性优化完成后才跑正式 benchmark。

## 5. 当前实现约束

- `charint` 的正式算法口径保持为：
  - 固定子步主链
  - `4` 点 Gauss-Legendre 的 Duhamel 源项积分
- `adaptive_substeps` 仍只视为实验入口，不默认混入正式 benchmark、谱图或精度口径。
- 所有 Fortran 改动都必须做：
  - 编译检查
  - `line truncation` 检查

## 6. 当前默认工作流

### 6.1 每轮 `charint` 优化的最低检查

- `python build_extensions.py --module FS_electron_charint --force`
- `python tests/order_convergence_check.py`

### 6.2 只有在必要时才追加

- 轻量 one-shot timing
- `python tests/electron_solver_comparison.py`
- 光变图 / 谱图输出

### 6.3 不再默认每轮都做

- 重 benchmark
- 大量峰位、前沿、`g_lo` 诊断

## 7. 当前最值得做的优化方向

1. 继续降低 `L1` 子步循环单位成本
   - 减少子步内重复构造的几何量、source 预处理和 remap 预处理
   - 尽量把壳层不变量、介质不变量搬到子步循环外

2. 继续减少高代价标量函数
   - 少做 `exp/log/pow`
   - 减少反复除法
   - 减少小函数热点调用造成的开销

3. 继续减少 tracing / remap 的重复工作
   - 复用 cell hint、sweep、prepared remap
   - 避免对固定 source shape 或固定几何做重复预处理

4. 优先找结构性降本，而不是直接粗暴减小物理步长精度
   - 不优先靠减少 `L1` 或降低积分阶数换速度
   - 先找等价重用、缓存、批量化和循环外提

## 8. 当前明确不做的事

- 不把 `charint` 当前实验结果直接升级成全项目默认基线。
- 不继续围绕 `peak-gamma`、低能前沿、`g_lo` 跳变做主线排障。
- 不在没有必要时重跑重 benchmark。
- 不把历史阶段性的离线研究结果继续堆积在 `plan.md`。

## 9. 下一步直接行动项

1. 继续压 `charint` 主核的单线程算量。
2. 重点审查 `L1` 子步循环内还能外提或复用的固定对象。
3. 优先寻找新的结构性优化点，不只盯现有两三个热点。
4. 每轮只保留必要的编译和粗阶数检查，benchmark 按需触发。
