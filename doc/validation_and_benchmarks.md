# 验证与基准

本文档固定 ASGARD 的验证层级。核心原则：验证应回答明确假设或支持决策，不做低信息增益的穷举。

## 验证层级

### 文档或 Python-only 轻量改动

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git diff --check'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run --with "mkdocs<2" --with "mkdocs-material>=9.5" mkdocs build --strict --site-dir /tmp/asgard_mkdocs_site'
```

如果改动涉及示例命令，还应运行对应 smoke 或最小命令片段。

### Fortran 语法和行截断

Fortran 重要改动必须跑：

- 受影响 `build_extensions.py --module ... --force`
- `gfortran -Wall -Wline-truncation -fsyntax-only`
- 最小相关 smoke test

Line check 必须从 `/tmp` 执行，并指定临时 `-J` / `-I` 目录。

### 冒烟测试

Forward/electron：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/readme_smoke_bench.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/fitter_public_api_smoke.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/eats_adaptive_projection_smoke.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/fullhide_2d_smoke_bench.py'
```

Reverse shock：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/reverse_shock_smoke.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/reverse_shared_solver_smoke.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/structured_shared_solver_smoke.py'
```

Hadronic：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/hadronic_1d_smoke.py'
```

Prompt snapshot：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/internal_shock_prompt_smoke.py'
```

`tests/dg_1d_smoke.py` 是 DG 谱形诊断入口。它检查有限值、非负、活动支撑无零洞、无多重 grid-scale sawtooth turns、粒子数和同步光度量级；尖锐曲率本身不判失败，因为 troubled positive-kernel 的目标是消除 Gibbs 振荡，同时保留真实冷却断点和高能 cutoff。当前工作树中该脚本在 RS DG spectrum sawtooth-turn 判据处失败，因此不属于绿色基线。这个 smoke 不是阶数测试。当前阶数口径是：未滤波光滑谱元的能量空间形式阶数为 P12 的 \(O(\Delta y^{13})\)，普通 shell-to-shell 电子推进受后向 Euler 限制为 \(O(\Delta R)\)，主 RS crossing 后纯冷却解析映射对固定冷却系数无 BE 时间截断误差。需要验证阶数时，应分开做能量空间光滑谱测试和半径步长减半测试，不用真实断点附近的局部 Gibbs 区域拟合单一阶数。

### 当前已知失败诊断

以下入口当前不属于“绿色基线”，但必须保留为真实问题入口：

- `tests/hadronic_reverse_shock_smoke.py`：base 和 RS light proton-synch 分支可执行；RS full-chain hadronic 分支报 `electron_energy_gev must be logarithmically uniform`。
- `tests/electron_photon_joint_secondary_feedback_smoke.py`：同样在 formal hadronic electron-energy grid contract 处失败。
- `tests/dg_1d_smoke.py`：当前在 RS DG spectrum sawtooth-turn 判据处失败。

修复这些失败时不得删除断言、放宽阈值或添加 fallback；必须回到对应 hadronic grid contract 或 DG RS 电子谱演化路径。

## 按区域划分的构建门槛

Electron 1D：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_dg_1d --force'
```

Electron 2D：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_charint_2d --force'
```

反向激波电子：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_reverse_kernel --force'
```

Hadronic：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_forward_1d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_reverse_1d --force'
```

Dynamics：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module Dynamics_forward --module Dynamics_reverse --force'
```

## 基准刷新

最少元数据：

- 刷新前后 Git HEAD。
- `git status --short --branch`。
- `git diff --stat`。
- 完整命令文本。
- 是否重建 Fortran extension。
- 输出路径。
- 物理验收口径。

旧 `scripts/benchmarks/` 入口已删除。新增 formal benchmark 必须先进入 `tests/`，并在脚本内固定物理参数、网格、输出路径和验收口径。

当前保留的 q-shell diagnostic benchmark：

| 入口 | 决策问题 | 输出默认目录 | 注意事项 |
| --- | --- | --- | --- |
| `tests/benchmark_theta_j_multiples_magnetic_decay.py` | 有限 q-shell 与磁场衰减在离轴角度扫描中如何改变 lightcurve/spectrum。 | `output/benchmark_1d_vs_qshell_theta_j_multiples_bdecay_alpha04/` | 2D 分支复用同一 solve state 只重跑 projection；只适用于底层物理状态不变的 top-hat 诊断。 |
| `tests/benchmark_skymap_centroid_motion.py` | 1D thin shell 与 2D q-shell 的 sky map、flux centroid 和 apparent speed 差异。 | `output/benchmark_1d_vs_qshell_skymap_motion_bdecay_alpha04/` | 是诊断 benchmark，不是已验证天图科学产品。 |

## 产物策略

可以提交：

- 作为文档资产的可复现 benchmark figures。
- 由已提交脚本生成且物理验收通过的 CSV summaries。
- 复现图像所需的脚本和文档。

不要提交：

- `.buildcache/`
- `__pycache__/`
- 临时 debug 脚本
- 失败占位图
- 不完整 CSV
- 一次性 notebook export
- 本地编译 `.so`, `.o`, `.mod`

## 物理验收

Forward-shock：

- Light curves 应平滑，除非物理 density jump 或 injection event 产生已记录特征。
- 电子谱、磁场和最终光变应连续演化；需要断频诊断时通过 `nu_callback` 临时收集。
- SSA breaks 不应出现 grid-cell discontinuity。
- `solver_options.geometry_projection="chi_eats_2d"` 只验收 FS synchrotron+SSA；图中 forward SSC 仍是 shell-level 总通量贡献。Finite q-shell 投影必须使用正半径 `chi_radius_cm`、局域 `chi_gamma_bulk` 和非零 `chi_dvolume_weight`，SSA survival 必须按 emitting cell 的 optical-depth coordinate 平均。图中不得出现由负半径、负通量、孤立尖峰、全部 `chi_dvolume_weight` 同时归零或源项截断造成的光变断崖。2D/1D SED 与 top-hat 角度扫描允许离轴情况下出现 order-unity 以上差异，但光变和频谱方向应保持连续。
- ISM q-shell `num_chi` convergence scan 中 `num_chi=512` 是参考曲线；`F_chi/F_512` 的时间演化应连续，不允许用 smoothing 或显示裁剪掩盖孤立尖峰。
- q-shell magnetic-decay angle benchmark 允许 1D thin shell 与 2D q-shell 在离轴峰时/峰值上出现 order-unity 差异，但每个角度和频段的光变、谱形、centroid offset 和 \(\beta_{\rm app}\) 应随时间平滑。

Reverse-shock：

- Pre-crossing 的 `M3` crossing 端点应由 `m3_frac=1` 给出，不允许 RK step 跨越 pre/post 方程分支。
- `sigma -> 0` 必须恢复 unmagnetized baseline。
- `B3`, `gamma34`, `U3/V3` 应平滑；断频只作为可选 `nu_callback` 诊断。
- `dg_1d` 的 RS primary electron 路径必须通过 `reverse_shared_solver_smoke.py` 和 `structured_shared_solver_smoke.py`。RS crossing 后的纯冷却段在 `fullhide_1d` 和 `dg_1d` 下都使用 post-crossing direct characteristic map，并采用闭合低能边界保持冷却到网格下界以下的电子数；121 格 fullhide 的旧 post-crossing 宽尾不作为 DG 真值，需用 direct-map 有效支撑图、高分辨率对照或守恒谱形诊断判断。
- VegasAfterglow 是 comparison backend，不是目标真值。

Hadronic：

- 启用过程必须明确列出。
- Formal path 在定义 \(\chi\)-local contract 前保持壳层级 1D。
- Pair cascade 必须保持 gamma-gamma pair/synch branch 的解释，不扩展成 IC-mediated cascade。

Polarization：

- Peak amplitude 和 peak time 是分开的诊断。
- 不用经验 time shift 或 smoothing 修正 timing mismatch。

Prompt snapshot：

- `prompt/` 内部激波只按 snapshot 诊断验收：\(\sigma\rightarrow0\) 接近 hydrodynamic baseline，磁化 case 的 baryonic mass 降低且 ordered field 上升，fast shock 不允许时对应分量为零。
- `fs_sync`、`fs_ssc`、`rs_sync`、`rs_ssc` 和 `total` 必须有限、非负，活动光变段分段平滑。
- 不把 prompt snapshot smoke 当作观测 GRB prompt light curve 拟合验收。

## 失败处理

验证失败时：

1. 记录准确命令和 artifact。
2. 判断失败属于 compile、API、numerical 还是 physical。
3. 检查最小 source closure。
4. 除非失败发生在真实系统边界，不添加 fallback 或 guard code。
5. 只重跑能覆盖修复的最窄检查。
