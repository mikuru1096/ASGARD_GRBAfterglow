# 含时 BH / 二级反馈与光子谱联立方案

## 当前契约

- 所有输运自变量固定为半径 `R`，不以时间 `t` 推进。
- 微物理冷却、反应率、源项若自然单位是 `s^-1` 或 `cm^-3 GeV^-1 s^-1`，进入求解器前统一换算到 `R` 坐标：

```text
d t' / dR = 1 / (beta Gamma c)
Lambda_R = Lambda_t / (beta Gamma c)
Q_R = Q_t / (beta Gamma c)
dotgamma_R = dotgamma_t / (beta Gamma c)
```

- public 接口保持一个最小 switch：

```python
Setups(electron_photon_coupling="separated" | "joint")
```

- `electron_photon_coupling="separated"` 是默认值，保留原流程：

```text
electron solver -> photon_field_stage -> hadronic -> separated BH merge/recompute seed
```

- `electron_photon_coupling="joint"` 在同一 `R` 网格上迭代闭合：

```text
electron state -> photon field -> hadronic transport
-> secondary electron source + photon source/sink
-> electron solve with external Qe and IC seed
-> photon field rebuild
```

## joint 方程含义

在每个 shell 的 `R` 区间内推进：

```text
dN_p/dR =
  Q_p,R
  - d/dE_p[(dotE_ad,R + dotE_syn,R + dotE_BH,R + ...) N_p]
  + Q_p,reinj,R

dN_e/dR =
  Q_e,shock,R + Q_e,secondary,R
  - d/dgamma_e[(dotgamma_ad,R + dotgamma_syn,R + dotgamma_IC,R[n_gamma]) N_e]

dn_gamma/dR =
  Q_syn,R[N_e,B] + Q_IC,R[N_e,n_gamma]
  + Q_hadronic_gamma,R
  - Lambda_gamma,R n_gamma
```

其中 `Q_e,secondary,R` 由已启用且 kernel 已正式输出的二级电子源项构成：BH pairs、pp pairs、gamma-gamma pairs。pγ/π/μ 链若当前 formal kernel 没有给出可归一化的 e± 注入谱，则不以能量守恒外推补项。

光子源项在 joint 内使用同一组 `N_e,n_gamma` 输入生成 IC 冷却和 IC photon source，避免只改变电子冷却而丢失光子源。光子 sink 使用 formal kernel 输出的 pγ photon survival、BH photon loss 和 gamma-gamma absorption；不引入经验 sink。

## 支持边界

`joint` 第一版只支持：

- forward shock。
- `electron_solver="fullhide_1d"`。
- `hadronic_solver="am3_1d"`。
- `bethe_heitler=True`。
- `ssc=True` 且 `index_y=1`。
- 非 adaptive electron substeps。

直接报错的组合：

- reverse shock full-chain。
- chi-resolved / 2D hadronic transport。
- structured backend。
- `fullhide_1d` 之外的 electron solver。
- 非 `hummer_2010_response` pγ scheme。
- IC-mediated electromagnetic cascade。

## 实现位置

- `ASGARD/api_model.py` / `asgard_core/asgard_config.py`：新增 `electron_photon_coupling` public switch。
- `asgard_core/asgard_state.py`：新增 joint stage、secondary feedback 汇总、photon survival/sink 应用、gamma-gamma pair/synch cascade 接入。
- `asgard_core/asgard_runtime.py`：formal hadronic path 输出 `secondary_electron_source_r`，并把 BH/pp source 从 `s^-1` 换算到 `R^-1`。
- `asgard_core/asgard_types.py`：`HadronicSolution.secondary_electron_source_r`。
- `src/Electron/electron_forward_fullhide_1d.f90`：电子求解入口接受外部 `Secondary_source(gamma_e,R)`，作为 `Q_e,secondary,R` 加入电子方程。
- `src/Hadronic/hadronic_bethe_heitler_kernel.f90`：BH kernel 输出与 pair/proton loss 同归一化的 photon loss rate。
- `scripts/benchmarks/time_dependent_bh_photon_benchmark.py`：生成 weak-feedback、BH-active、strong-wind-BH 三组 separated/joint 对比图、CSV 和 metadata。

## 验收记录

已执行的验证：

- `build_extensions.py --module electron_forward_fullhide_1d --force`
- `gfortran -Wall -Wline-truncation -Werror=line-truncation -fsyntax-only` source closure 检查。
- `tests/electron_photon_coupling_config_smoke.py`
- `tests/electron_photon_joint_smoke.py`
- `tests/electron_photon_joint_secondary_feedback_smoke.py`
- `tests/electron_photon_separated_regression_smoke.py`
- `tests/electron_photon_ic_consistency_smoke.py`
- `tests/hadronic_r_coordinate_smoke.py`
- `tests/hadronic_bethe_heitler_smoke.py`
- `tests/hadronic_pair_cascade_smoke.py`
- `tools/check_text_encoding.py`
- `git diff --check`

formal benchmark 命令：

```bash
uv run python scripts/benchmarks/time_dependent_bh_photon_benchmark.py --mode formal
```

formal benchmark 摘要：

```text
scenario        tau_BH,max joint   median LC joint/separated   max LC joint/separated
weak_feedback   3.60e-13           1.045                       1.324
bh_active       9.08e-11           1.108                       2.671
strong_wind_bh  1.31e-10           1.023                       1.917
```

生成的 PNG/PDF/CSV/metadata 位于 `output/asgard_doc/time_dependent_bh_photon_benchmark/`，属于 benchmark artifact，可由脚本复现；清理提交时不纳入版本控制。
