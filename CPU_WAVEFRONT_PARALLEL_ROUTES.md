# CPU fullhide 波前并行路线

## 目标

CPU 版本应对齐 GPU 版本的物理离散和代数依赖，而不是照搬 GPU 的执行形状。GPU 端适合 cell-level anti-diagonal wavefront，因为设备有大量线程且 kernel launch/设备驻留是主要约束。CPU 端若直接对单个 shell 的 `(gamma, substep)` cell 做每条反对角 OpenMP 同步，`num_gam_e ~= 96` 时每个 wave 的并行宽度太小，同步成本会吞掉收益。

因此，坚持 CPU 波前并行时，需要把并行粒度放大。下面三条路线都保留 fullhide spacetime stencil 的因果依赖，但采用 CPU 可承受的任务粒度。

## 路线 1：多 shell / 多 patch 合并波前

把单个 shell 的 wavefront 扩展为批量 wavefront，在同一个反对角层上并行多个独立 patch、theta ring、shell 或 sample：

```text
for wave:
    parallel over active (sample, patch, shell, gamma) cells on this wave
```

有效并行宽度从

```text
min(num_gam_e, num_substeps)
```

提升为

```text
num_active_items * min(num_gam_e, num_substeps)
```

其中 `num_active_items` 可以来自 theta patches、多个独立样本、或一组已完成前置辐射/冷却计算的 shell。

优点：

- 最符合 CPU 多核：每个 wave 的任务数足够多，线程不会只抢 96 个 gamma lane。
- 保留 GPU 同源 wavefront 依赖图。
- 适合生产网格和拟合批量评估。

风险：

- 需要重排 runtime 调度，让多个独立 transport work item 一起进入同一个批处理入口。
- 必须保证不同 work item 的冷却、注入、边界状态已经预先完成，不能跨物理依赖乱序。

验收：

- 生产网格 `num_r=300, num_theta=300, num_phi=1, num_gam_e=96, num_nu=128, num_tobs=72`。
- 比较 `fullhide_1d` baseline 与批量 wavefront 的 electron stage 和 total wall time。
- ISM/Wind flux 必须 finite、non-negative，抽样光变曲率不能出现非物理尖峰。

## 路线 2：tile-wavefront，块内串行、块间波前并行

不要把单个 `(gamma, substep)` cell 作为 OpenMP 任务。把 spacetime 平面切成 tile，例如：

```text
tile_gamma = 16 或 32
tile_step  = 16 或 32
```

每个 tile 内按 CPU cache 友好的 step-major 顺序串行回代；tile 之间按反对角依赖并行：

```text
for tile_wave:
    parallel over active tiles:
        solve_tile_step_major()
```

这样同步次数从

```text
num_gam_e + num_substeps
```

降到

```text
ceil(num_gam_e / tile_gamma) + ceil(num_substeps / tile_step)
```

每个并行任务也从一个 cell 变成一个 tile，算术量和内存局部性都更适合 CPU。

优点：

- 仍然是严格波前并行，但同步粒度变粗。
- tile 内数据连续，cache 命中优于 cell-level wavefront。
- 可以和路线 1 组合，形成 `(batch item, tile)` 级并行。

风险：

- tile 边界需要明确传递左/上依赖状态，不能用经验补丁修边界。
- tile 太小同步仍然贵；tile 太大并行宽度不足，需要生产网格调参。

验收：

- 先用单 shell 代数等价测试证明 tile 结果与 step-major 标量结果一致。
- 再用生产网格测 `tile_gamma/tile_step` 的少量候选，不做低信息增益穷举。
- 接受标准是 electron stage 明确快于 baseline，且 ISM/Wind 都不退化。

## 路线 3：OpenMP task depend 波前流水

用 OpenMP task dependency 表达 tile 之间的因果图，而不是手写：

```text
for wave:
    parallel do
    barrier
```

任务形式：

```fortran
!$omp task depend(in: previous_step_tile, upper_gamma_tile) depend(out: current_tile)
call solve_tile(...)
```

runtime 可以在某个 tile 的依赖满足后立即执行下游 tile，不必等待同一 wave 的所有 tile 全部完成。这个方法仍是波前并行，但调度是流水式的。

优点：

- 避免显式每条 wave 的全局 barrier。
- 对不规则 `num_substeps`、Wind/jump 介质更自然。
- 与 tile-wavefront 兼容，适合作为路线 2 的调度实现。

风险：

- 单 shell tile 数不足时 task runtime 开销仍然过大。
- Fortran OpenMP task depend 的数组 section 写法和编译器支持需要单独验证。
- 依赖对象必须真实表示物理/代数边界，不能用粗糙 token 隐藏数据竞争。

验收：

- 先做独立 Fortran probe：相同随机正系数系统，task-depend tile wavefront 与标量 step-major 在 roundoff 内一致。
- 再接入 `electron_transport_common.f90`，只替换执行调度，不改物理系数公式。
- 生产网格下报告 1/8 线程 electron stage、total wall time、flux 平滑性。

## 推荐顺序

1. 先实现路线 2 的 tile-wavefront 标量等价 probe。
2. 若 tile 级同步成本仍偏高，把路线 1 的 batch item 维度合并进 tile-wavefront。
3. 最后用路线 3 替代显式 wave barrier，前提是 task-depend probe 和生产网格基准都证明收益。

不应接受的路线：

- 直接照搬 GPU cell-level anti-diagonal OpenMP `parallel do`。
- 只在 smoke 网格上证明加速。
- 用 smoothing、截断、经验因子修正不平滑输出。
- 为了避免退化添加隐藏 fallback。若某介质或线程数不快，必须直接报告并重新设计并行粒度。

## 2026-06-09 调度 probe 结果

当前生产 `fullhide_1d/fullhide_2d` 尚未接入路线 3/2/1 的 runtime selector；本次只测试独立 Fortran 调度 probe
`scripts/benchmarks/cpu_wavefront_routes_probe.f90`。probe 使用同一个正系数 spacetime fullhide 回代系统，把当前 CPU step-major
序列推进作为基线，然后按 3-2-1 顺序测试：

```bash
rtk wsl -e bash -c 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && gfortran -O3 -march=native -fopenmp -Wall -Wline-truncation scripts/benchmarks/cpu_wavefront_routes_probe.f90 -o /tmp/cpu_wavefront_routes_probe && /tmp/cpu_wavefront_routes_probe 8 300 96 300 32 32 3 | tee /tmp/cpu_wavefront_routes_probe_8t.csv'
rtk wsl -e bash -c 'source ~/.wsl_env && /tmp/cpu_wavefront_routes_probe 1 300 96 300 32 32 3 | tee /tmp/cpu_wavefront_routes_probe_1t.csv'
```

基线提交：`e6952163cf1255716f005439d13c78de0fca872a`。网格形状为 `items=300, gamma=96, steps=300`，
tile 为 `32x32`，每个 timing 为 3 次重复均值。相对当前 step-major 基线的加速比如下：

| route | threads | seconds | baseline seconds | speedup | max abs err | max rel err |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 3_task_depend | 8 | 7.989372e-02 | 9.883206e-02 | 1.2370 | 3.330669e-16 | 5.870908e-16 |
| 2_tile_wave | 8 | 2.214191e-02 | 9.883206e-02 | 4.4636 | 3.330669e-16 | 5.870908e-16 |
| 1_item_batch | 8 | 1.340803e-02 | 9.883206e-02 | 7.3711 | 0.0 | 0.0 |
| 3_task_depend | 1 | 1.460849e-01 | 1.016890e-01 | 0.6961 | 3.330669e-16 | 5.870908e-16 |
| 2_tile_wave | 1 | 1.257861e-01 | 1.016890e-01 | 0.8084 | 3.330669e-16 | 5.870908e-16 |
| 1_item_batch | 1 | 1.022198e-01 | 1.016890e-01 | 0.9948 | 0.0 | 0.0 |

结论：路线 1 的批量 item 粗粒度并行在 CPU 上最强，路线 2 的 tile-wavefront 次之，路线 3 的 `omp task depend`
调度开销较大，只在 8 线程下略快于标量基线。三条路线在 probe 上与 step-major 结果在 roundoff 内一致。
