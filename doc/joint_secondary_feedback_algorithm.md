# 含时二级反馈的算法契约

本文档描述 `electron_photon_coupling="joint"` 的实现流程、数据结构、函数入口和验证要求。物理方程与源汇预算见 `doc/joint_secondary_feedback_physics.md`。

## 1. public 接口

joint 模式只新增一个 public switch：

```python
SolverOptions(electron_photon_coupling="separated" | "joint")
```

默认值是 `separated`。`hadronic_solver="am3_1d"` 仍只表示 formal hadronic 微物理核；`electron_photon_coupling` 只控制 electron / photon / hadronic 阶段的耦合方式。

不新增 `joint_solver`、`bh_time_solver` 或类似 public 名称。后续扩展仍应复用现有 `Radiation` / `Hadronic` flags：

```text
proton_synch
include_pgamma
pp
bethe_heitler
hadronic_inverse_compton
pair_production
pair_cascade_iterations
neutrino
```

## 2. separated 主链

默认主链保持原语义：

```text
solve_dynamics
-> solve_electron
-> build_photon_field
-> solve_hadronic
-> separated BH merge / recompute seed
-> reverse shock
-> observer assembly
```

其中 `_mergebh` 只属于 separated 路径。该路径不能被 joint 复用为真实反馈，因为它发生在 electron solve 之后，无法改变本 shell 中电子冷却与 IC photon source。

## 3. joint 主链

joint 主链位于 `asgard_core/asgard_state.py::_jointstage`：

```text
validate joint config
primary_electron = solve_electron()
photon_field = build_photon_field(primary_electron)

repeat fixed small iteration count:
    hadronic = solve_hadronic(primary_electron, photon_field)
    photon_field = rebuild electron photon field
    apply pγ/BH photon survival
    secondary_source_R = assemble secondary e± source
    add gamma-gamma pair/synch feedback when enabled
    primary_electron = solve_coolingseed(
        cooling_seed = photon_field.hadronic_target_seed,
        secondary_source_R = secondary_source_R,
    )
    photon_field = rebuild from updated electron
    apply secondary photon feedback for final observer state

observer assembly
```

当前 fixed iteration count 是内部常量，不是 public convergence knob。第一版目标是 shell-level closure，而不是暴露一个可调非线性求解器接口。

## 4. 配置校验

`joint` 在进入主链前校验：

- `electron_solver == "fullhide_1d"`。
- 无 reverse shock。
- 非 chi-resolved / 2D transport。
- `structured_backend == "fortran_1d"`。
- hadronic enabled 且 `epsilon_p > 0`。
- `bethe_heitler=True`。
- `hadronic_solver == "am3_1d"`。
- `pgamma_scheme == "hummer_2010_response"`。
- `index_y == 1`。
- `electron_adaptive_substeps=False`。

不满足这些条件时直接报错。这里的报错是系统边界校验，不是 fallback。

## 5. `PhotonFieldState` 语义

`PhotonFieldState` public 字段语义保持不变：

- `forward_syn_seed`：forward electron synchrotron seed。
- `hadronic_forward_ssc_seed`：由 forward synch seed 计算的 SSC/IC seed contribution。
- `hadronic_target_seed`：hadronic 和 joint electron cooling 使用的目标 photon field。
- `absorption_syn_seed`：synchrotron absorption 使用的 seed。
- `absorption_ssc_seed`：SSC absorption 使用的 seed。

joint 内部允许更新：

```text
hadronic_target_seed
absorption_syn_seed
absorption_ssc_seed
hadronic_forward_ssc_seed
```

但不改变字段含义。pγ/BH/gamma-gamma photon survival 作用在这些 seed field 上；pair-synch cascade source 以 seed contribution 加回。

## 6. radius 坐标步长

hadronic Python path 不再用 observer-time shell spacing 推进 formal BH/proton/secondary steps。radius shell 的 comoving 时间步是：

```text
Delta t'_i = Delta R_i / (beta_i Gamma_i c)
```

其中

```text
Delta R_i =
  R_1 - R_0,          i = 0
  R_i - R_{i-1},      i > 0
```

对应实现：

- `asgard_core/asgard_runtime.py::_hadronic_shell_comoving_dt_from_radius`

shell-local 推进使用 shell 间隔对应的 comoving time。

## 7. hadronic 输出契约

`HadronicSolution` 在 formal path 中增加内部字段：

```text
secondary_electron_source_r(gamma_e, R)
tau_bh(nu, R)
bh_photon_loss_rate(nu, R)
```

含义：

- `secondary_electron_source_r` 是已经换算为 electron equation 可直接使用的 `Q_e,secondary,R`。
- `tau_bh` 是 BH photon loss 在当前 shell path 上的 optical depth。
- `bh_photon_loss_rate` 保留 comoving microphysics rate，用于诊断和预算。

这些字段不是新的 public solver 接口；它们是 `joint` state 内部闭合所需的正式 dataclass contract。

## 8. 电子 coupled pass

joint 电子入口：

```text
solve_coolingseed(...)
  -> fs_fullhide_coupled(...)
```

Fortran entry 位于：

```text
src/Electron/electron_forward_fullhide_1d.f90
```

新增输入：

```text
Seed_cooling(num_nu, num_R)
sec_source(num_gam_e, num_R)
```

`Seed_cooling` 用于 IC cooling。`sec_source` 以 `Q_e,secondary,R` 加到每个 shell substep 的电子源项中：

```text
dF1 = Q_e,shock + sec_source(:, I_tobs)
```

`fs_fullhide_coupled` 只接受 fixed substeps 和 `index_y=1`，因为当前预算测试依赖同一个数值 IC kernel。

## 9. IC 预算 kernel

IC 冷却预算 kernel 位于：

```text
src/Electron/electron_cooling_ic_kernel.f90
electron_cooling_ic_loss_emissivity_budget
```

目的：用与 SSC emissivity 一致的 Jones/KN kernel 计算 electron IC loss。这样同一个 `N_e,n_gamma` 能同时约束：

```text
electron cooling
IC photon production
```

该局部一致性不再保留独立脚本，joint feedback 的可执行验收收敛到端到端入口。

## 10. gamma-gamma pair branch

joint secondary feedback 入口：

```text
asgard_core/asgard_state.py::_jointfeedback
```

当 `pair_production=True`：

```text
solve_branch(photon field)
-> pair injection / photon loss
-> advance pair spectrum in shell sequence
-> compute pair synchrotron seed
-> add pair source to Q_e,secondary,R
-> apply photon survival to photon field
```

如果 `pair_cascade_iterations > 1`，使用 shell-sequence time-dependent pair/synch cascade path；否则使用 single-pass pair production branch。两条路径都只覆盖 gamma-gamma pair/synch feedback，不覆盖 IC-mediated electromagnetic cascade。

## 11. observer assembly

observer stage 仍在 joint forward stage 之后统一执行：

```text
_observerstage
```

hadronic luminosity components、absorption factors、reverse shock 和 final projection 的 public 输出语义不变。joint 改变的是进入 observer stage 前的 electron/photon/hadronic state，而不是 observer projection API。

## 12. benchmark 入口

旧含时 BH / joint Python benchmark 脚本已删除。重新建立 formal benchmark 前必须先说明要回答的物理假设，并把可复用入口放入 `tests/`。

## 13. 最小验证集合

Fortran 改动后：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force'
```

必须额外跑 line-truncation source closure 检查，使用干净 `/tmp` module 目录。

当前工作树中 joint 端到端验证会在 formal hadronic electron-energy grid contract 处失败，错误为 `electron_energy_gev must be logarithmically uniform`。修复时必须回到 hadronic input grid 契约和 joint shell state 构造，不能删除断言、跳过 formal hadronic 分支或添加 fallback。

文档/格式：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && git diff --check'
```

## 14. 开发边界

- 不新增 public solver name 表示 joint；继续使用 `electron_photon_coupling`。
- 不把 separated 的 BH post-merge 当成 joint feedback。
- 不用 smoothing 修复 light curve 或 shell diagnostics。
- 不用经验 photon sink/source 补齐 formal kernel 没有输出的项。
- 不把 observer luminosity 直接当作 local photon density；必须先定义逃逸、体积和吸收归一化。
- 不在没有 \(\chi\)-local photon/hadron contract 前实现 \(\chi\) 分辨 hadronic transport。
