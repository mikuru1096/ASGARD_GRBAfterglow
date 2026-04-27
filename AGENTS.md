# ASGARD Agent Notes

## Working Rules

- 用中文回复。
- Shell 命令前缀 `rtk`，必须在 WSL Ubuntu 内执行。
- 默认环境 `WSL + uv`。不要用 `wsl.exe`、PowerShell、Windows CMD 或 Windows 侧 Python。
- 非模拟代码禁止加数值保护。
- 物理量随时间演化若不连续不平滑，优先怀疑数值或物理 bug。
- Fortran 代码保持 idiomatic Fortran 风格。
- Fortran 重要改动后必须跑：受影响的 `build_extensions.py --force`、`-Wline-truncation` 检查、最小相关 smoke test。
- 不要提交 `.buildcache/`、临时 debug 脚本、失败占位图。

## Build Commands

```bash
rtk uv run python build_extensions.py --module FS_electron_fullhide_2d --force
rtk uv run python build_extensions.py --module FS_electron_charint_2d --force
rtk /usr/bin/gfortran --version
rtk uv run python tests/readme_smoke_bench.py
rtk uv run python tests/reverse_shock_smoke.py
rtk uv run python tests/hadronic_1d_smoke.py
rtk uv run python tests/hadronic_species_transport_smoke.py
```

## Declaration Policy

**模块级**: 统一 `implicit none`。

**子程序级**: 按性质分两类。

**A 类 — 保留 `implicit REAL(8)(A-H,O-Z)`**（物理/数值密集子程序）:
- `src/Electron/electron_cooling_kernel.f90`（全部）、`electron_radiation_kernel.f90`（大部分）
- `src/Radiation/radiation_common.f90`（辐射核）、`SSC_spec.f90`、`Annihilation.f90`、`Seed_reverse.f90`
- `src/Hadronic/hadronic_radiation_kernel.f90`、`hadronic_transport_kernel.f90`、`hadronic_common.f90`（物理子程序）、`hadronic_pair_production_kernel.f90`
- `FS_electron_*.f90`（驱动子程序）
- `src/Interpolation/SED_interpolation*.f90`

**B 类 — 使用 `implicit none`**（验证/编排子程序，整数/逻辑比例高）:
- `src/Electron/electron_transport_common.f90`、`electron_seed_history_kernel.f90`、`electron_transport_2d_kernel.f90`、`electron_injection_profiles.f90`、`electron_reverse_kernel.f90`、`adaptive_resampling_mod.f90`
- `src/Hadronic/hadronic_acceleration_kernel.f90`、`hadronic_bethe_heitler_kernel.f90`、`hadronic_decay_kernel.f90`、`hadronic_hadronic_ic_kernel.f90`、`hadronic_interaction_kernel.f90`、`hadronic_pp_kernel.f90`、`hadronic_secondary_radiation_kernel.f90`、`hadronic_species_transport_kernel.f90`
- `src/Interpolation/interpolation_common.f90`

**声明块压缩规则**: 同类型声明合并到一行，语义相近的量分组。B 类子程序声明块 ≤15 行。禁止每行只声明一个变量。

## Public Runtime Status

- **Public API**: `observe(model, config=...)`, `run_fit(config)`
- **Electron solvers**: `fullhide_1d`, `slc1_1d`, `charint_1d`, `charint_2d`, `t2g1_1d`, `weno5_1d`, `fullhide_2d`
- **Aliases**: `fullhide→fullhide_1d`, `slc1→slc1_1d`, `charint→charint_1d`, `t2g1→t2g1_1d`, `weno5→weno5_1d`
- `*_1d → num_chi=1`, `*_2d → num_chi=64` (default)
- **`pgamma_scheme`**: `hummer_2010_response`, `ka2008_reference`, `disabled`
- **Hadronic coupling** (all active in formal 1D runtime):
  - proton injection, adiabatic/synchrotron cooling, proton synchrotron
  - `hummer_2010_response`: `alpha_p` + `Q_p^reinj` feedback to proton transport; `alpha_gamma^{pγ}` as local shell photon-loss closure (`tau_pg = alpha * R/(12 Γ c)`), survival factor `(1-e^{-τ})/τ` applied before observer projection
  - Bethe-Heitler: proton continuous cooling + secondary e± → merge into forward electron distribution → recompute seed_syn
  - Pair production: observer-side `tau_pair` attenuation + pair synchrotron branch
  - Hadronic IC: proton channel active; pion/muon IC via explicit secondary species
  - pp: gamma, neutrino, pair source, proton loss
  - Secondary species transport: neutron, π±, μ± (left/right)
  - Secondary radiation: pion/muon synchrotron + IC
  - Hadronic acceleration/injection operators
- **Not yet implemented**: reverse-shock hadronic, 2D/chi-resolved hadronic transport, full time-dependent pair cascade PDE

## Electron Kernel Layout

- `src/Electron/electron_common.f90`
- `src/Electron/electron_radiation_kernel.f90`
- `src/Electron/electron_cooling_kernel.f90`
- `src/Electron/electron_seed_history_kernel.f90`
- `src/Electron/electron_transport_2d_kernel.f90`
- `src/Electron/electron_injection_profiles.f90`
- `src/Electron/adaptive_resampling_mod.f90`
- `src/Electron/electron_reverse_kernel.f90`
- `src/Electron/electron_transport_common.f90`

## 2D Solver State

- `fullhide_2d`: χ-dependent downstream velocity, implicit η/log-χ transport, χ-resolved historical photon fields, shell-reduced public `P_syn/Seed_syn`, shell-level diagnostic `ν_m/ν_c/ν_a`, reduced 6-band cooling grid, shell cooling assembly once per shell
- `charint_2d`: shared `fullhide_2d` outer physics + history + shell diagnostics; characteristic η/log-χ advection + implicit η diffusion + characteristic ξ/log-γ advance

## Current Test Entrypoints

- `tests/readme_smoke_bench.py`
- `tests/reverse_shock_smoke.py`
- `tests/fullhide_2d_smoke_bench.py`
- `tests/fullhide_2d_medium_diag.py`
- `tests/vegas_afterglow_comparison.py`
- Hadronic: `hadronic_1d_smoke.py`, `hadronic_proton_synch_1d_diag.py`, `hadronic_pg_neutrino_1d_diag.py`, `hadronic_species_transport_smoke.py`, `hadronic_secondary_radiation_smoke.py`, `hadronic_acceleration_smoke.py`, `hadronic_bethe_heitler_smoke.py`, `hadronic_hadronic_ic_smoke.py`, `hadronic_pair_production_smoke.py`, `hadronic_pp_smoke.py`, `hadronic_pgamma_benchmark_report.py`

## Documentation Entrypoints

- `doc/code_overview.md` — 代码总览、API、编排层、Fortran 层、hadronic 边界、构建/测试入口
- `doc/source_tree.md` — 源码树索引
- `doc/call_chain.md` — Python/Fortran 调用链
- `doc/electron_solver_algorithms.md` — 电子求解器算法与离散推导
- `doc/hadronic_pgamma_notes.md` — hadronic/pγ 状态与边界
- `doc/am3_migration_plan.md` — AM3 迁移计划（历史参考）

## AM3 / ASGARD Coexistence

- AM3 本地参考: `/mnt/c/Users/jia/Documents/New project/_external/am3_reference` (HEAD: `7aba970b`)
- AM3 仅作微观核和过程拆分参考，不替换 ASGARD 的 dynamics/electron/observer 主链
- ASGARD 主链: `dynamics → electron solver → photon target field → hadronic add-on → observer projection`
- 所有 AM3-derived 最终微物理实现必须在 `src/Hadronic/*.f90`
- Python 只做 orchestration、wrapping、benchmark、API glue
