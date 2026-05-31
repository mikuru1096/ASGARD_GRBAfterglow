# Physics Model

本文档描述当前 ASGARD 物理模型和实现边界。它不是教科书推导，而是代码中已实现契约的索引。

## 总体图像

ASGARD 的主线是：

```text
relativistic blast wave
  -> electron / proton distribution evolution
  -> local radiation and photon fields
  -> absorption / cascade / cross-zone coupling
  -> equal-arrival-time observer projection
```

Python 层组织状态机、配置和观测投影；Fortran 层求解电子、辐射、动力学和强相互作用微物理核。

## Forward Shock

Forward-shock branch 包含：

- blast-wave dynamics
- electron injection and transport
- synchrotron radiation
- synchrotron self-absorption
- SSC and Compton cooling
- gamma-gamma absorption
- optional hadronic emission and feedback
- observer-frame projection

重要微物理参数：

- `eps_e`: 电子能量分数。
- `eps_B`: 磁场能量分数。
- `p`: 注入电子谱指数。
- `xi_N`: 非热电子数分数。
- `ssc`: 是否输出 SSC。
- `ssc_cooling` / `index_y`: 是否在冷却中包含 Compton 项以及采用的近似。

电子谱演化是核心物理状态，不能用后处理 smoothing 修复不连续结果。

## Electron Transport

当前登记求解器：

- `fullhide_1d`: 默认稳定基线，适合一般 public runtime 和拟合。
- `slc1_1d`: semi-Lagrangian / characteristic family path。
- `charint_1d`: characteristic integration path。
- `t2g1_1d`: legacy implicit transport path。
- `weno5_1d`: 高阶电子谱解析路径。
- `fullhide_2d`: energy + chi resolved electron transport。
- `charint_2d`: 2D characteristic path。

电子输运输出：

- `gam_e`
- `d_n_gam_e`
- `l_syn_spec`
- `seed_syn`
- `nu_m`
- `nu_c`
- `nu_a`

2D path 额外输出：

- `d_n_gam_e_chi`
- `chi_grid`

## Synchrotron / SSA / SSC

Synchrotron kernels:

- `src/Electron/electron_radiation_kernel.f90`
- `src/Radiation/radiation_common.f90`

SSA 冷却和 transfer：

- `src/Electron/electron_cooling_kernel.f90`
- `src/Radiation/radiation_common.f90`

SSC:

- `src/Radiation/SSC_spec.f90`
- `asgard_core/asgard_ssc.py`

当前同步辐射积分选择中，`index_syn_integr=1/2` 是固定网格快速路径；adaptive path 只作为显式诊断路径使用，不作为 public 默认。

## Reverse Shock

Reverse shock 当前基线：

- electron synchrotron
- RS SSC
- FS/RS cross-zone IC
- optional RS hadronic light path
- optional full-chain RS hadronic dispatch through formal 1D hadronic kernels

关键物理约束：

- 注入能标使用 shock-front `gamma34`。
- 区域 3 turbulent field 和 crossing 后热演化使用显式 `U3/V3` thermal state。
- `reverse_sigma` 引入 upstream magnetization；磁化 jump 使用 VegasAfterglow 的 jump-condition 形式作为来源和 comparison backend。
- `B3` 是 turbulent + ordered total field。
- `sigma -> 0` 必须回到当前非磁化 baseline。

VegasAfterglow 在当前项目中是 comparison backend，不是光变目标或 RS 物理基准。

## Hadronic

Hadronic forward-shock solver:

- `legacy_1d`: proton transport + proton synchrotron baseline。
- `am3_1d`: formal research path，覆盖 p-gamma、BH、pp、hadronic IC、secondary species transport、secondary radiation、pair production branch、neutrino。

Hadronic process switches:

- `epsilon_p`
- `proton_synch`
- `pg`
- `bethe_heitler`
- `hadronic_inverse_compton`
- `pp`
- `neutrino`
- `pair_production`
- `pgamma_scheme`
- `eta_acc`

RS hadronic:

- `FS_hadronic_reverse_1d`: light proton injection/transport + proton synchrotron。
- Full-chain RS path reuses `FS_hadronic_1d` formal kernels with RS seed photons, RS `B3`, shell energy and baryon target density.

Current hadronic boundary:

- Hadronic transport remains 1D shell-level.
- 2D / chi-resolved hadronic transport is intentionally not implemented.
- Required future contract: chi-local photon field, hadron density, secondary feedback and observer projection.

## Pair Production and Cascade

Current implemented path:

- observer-side gamma-gamma attenuation
- pair-production branch
- shell-sequence time-dependent gamma-gamma pair/synch cascade when `pair_cascade_iterations > 1`

Not implemented:

- IC-mediated electromagnetic cascade.

Reason: that extension requires a photon/e± source-sink equation, IC kernel contract and energy-budget benchmark. See `doc/pair_cascade_extension_boundary.md`.

## Polarization

Polarization path:

- synchrotron Stokes projection
- FS/RS electron synch
- FS/RS hadronic synch

Non-synchrotron branches are not mixed into polarization Stokes. Current Lan 2023 comparison records that amplitude is matched while peak time remains early; the current evidence points to dynamics/jet-evolution benchmark before projection-layer changes.

## Observer Projection

Observer stage combines:

- EATS/Doppler interpolation
- redshift and luminosity distance
- synchrotron/SSC/hadronic/RS/cross-zone components
- absorption factors
- optional patch integration for structured jets and sky image

Main Fortran interpolation:

- `src/Interpolation/SED_interpolation.f90`
- `src/Interpolation/SED_interpolation_structured.f90`

## Physical Acceptance Rules

- Time/space evolution of physical rates should be continuous and smooth.
- Non-smooth final physical parameter tracks are treated as likely bugs until proven otherwise.
- Do not fix physical discontinuities with empirical smoothing or projection-layer time shifts.
- Do not add numerical guards outside true system boundaries.
- Python should not replace final numerical microphysics kernels; Python orchestrates, wraps and benchmarks.
