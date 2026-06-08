# ASGARD 代码总览

本文档按当前工作树整理代码结构、运行主链和关键边界。电子算法细节见 `doc/electron_solver_algorithms.md`；当前唯一 TODO / 未完成项入口见根目录 `TODO.md`。

## 1. 公开 API

- `ASGARD/api_model.py`：`Model`, `Medium`, `JetProfile`, `ISM`, `Wind`, `TophatJet`, `GaussianJet`, `PowerLawJet`, `TwoComponentJet`, `StepPowerLawJet`, `Ejecta`, `Observer`, `Radiation`, `Setups`。`ISM/Wind` 与 named jet constructors 是 factory-style public constructors，返回带 `kind` 标记的 `Medium` / `JetProfile`。
- `Model.flux_density_grid(times_s, nu_hz)`, `flux_density(times_s, nu_hz)`, `spectrum(time_s, nu_hz)`, `flux(time_s, nu_min, nu_max)`, `sky_image(t_obs, nu_obs, fov)`, `details()`。
- `Model.polarization(times_s, nu_hz, magnetic_geometry=..., local_emissivity=...)`。
- Hadronic public switches：`Radiation.pair_production`, `Radiation.pg`, `Radiation.bethe_heitler`, `Radiation.pp`, `Radiation.neutrino`, `Radiation.reverse_epsilon_p`；cascade substeps 使用 `Setups.pair_cascade_iterations`。
- Reverse-shock magnetization switch：`ReverseShockConfig.sigma` / `Setups.reverse_sigma`。
- `ASGARD/api_observe.py`：`observe(model, config)`, `run_fit(config)`。
- `ASGARD/api_fit.py`：`Fitter`, `Param`, `FitResult`。
- Electron solver names：`fullhide_1d`, `slc1_1d`, `charint_1d`, `charint_2d`, `t2g1_1d`, `weno5_1d`, `fullhide_2d`。
- Aliases：`fullhide`, `slc1`, `charint`, `t2g1`, `weno5` 映射到对应 `*_1d`。

## 2. 运行时主链

```text
Model / observe / run_fit
  -> FitConfig -> SimulationSetup
  -> solve_state_from_setup
  -> solve_dynamics -> solve_electron -> photon_field_stage
  -> solve_hadronic -> solve_reverse_shock_emission
  -> observer assembly -> Radiation.annihilation
  -> Interpolation.sed_interpolation -> API result
```

核心状态对象位于 `asgard_core/asgard_types.py`：

- `DynamicsSolution`：`r_tobs`, `r_gamma`, `radius`, `swept_mass_g`。
- `ReverseShockDynamics`：`M3`, total `B3`, ordered crossing field `B3_ordered_cross`, `U3/V3`, `gamma34` 和 crossing thermal records；`B3` 是 turbulent `sqrt(8 pi epsilon_B,r U3/V3)` 与可选 ordered upstream field 的总和。
- `ElectronSolution`：`gam_e`, `d_n_gam_e`, `l_syn_spec`, `seed_syn`, `nu_m`, `nu_c`, `nu_a`；2D 额外包含 `d_n_gam_e_chi`, `chi_grid`；BH 额外包含 `d_n_gam_e_bh`。
- `PhotonFieldState`：forward synch seed、hadronic target field、absorption seed field。
- `HadronicSolution`：1D hadronic proton/secondary/radiation results。
- `ObserverState`：absorption factors、`tau_pair`、flux components。
- `FluxComponents`：`total`、FS synch/SSC、hadronic、RS synch/SSC、cross-zone IC。

状态机位于 `asgard_core/asgard_state.py`：

- `solve_state_from_setup`：dynamics -> electron -> photon field -> hadronic -> reverse shock -> observer。
- `_build_photon_field_stage`：复制 electron `seed_syn`；hadronic SSC seed 写入 target field。
- `_solve_hadronic_stage`：调用 `solve_hadronic`；BH 次级 e± 并入 forward electron 后重算 `l_syn_spec/seed_syn`；pγ photon survival 写回 photon field。
- `_assemble_observer_stage`：组装 FS synch/SSC、RS synch/SSC、cross-zone IC 和 hadronic components；hadronic photons 使用 electron Fortran kernel 的 SSA transfer。

拟合最短路径：

```text
Fitter.loglike -> compile_problem -> eval_loglike -> solve_state_from_setup
  -> project_flux_grid -> combine_multiband_flux -> compute_light_curve_redchi
```

## 3. Python 编排层

- `asgard_setup.py`：`FitConfig -> SimulationSetup`。
- `asgard_config.py`：`FitConfig`, `SimulationSetup`, `PhysicalSolution`, `FitResult` 和 runtime config dataclasses；旧 `asgard_models.py` compatibility shim 已移除。
- `asgard_runtime.py`：backend selection、Fortran extension dispatch、array wrapping。
- `asgard_state.py`：主状态机和跨阶段耦合。
- `asgard_ssc.py`：forward SSC auxiliary grid 与 seed。
- `asgard_coupling.py`：FS/RS cross-zone IC geometry 与 seed field coupling。
- `asgard_postprocess.py`：observer projection、band aggregation、fit postprocessing。
- `asgard_fit.py`：fit problem compilation 和 likelihood path。
- `asgard_types.py`：runtime dataclass contracts。
- `structured_jet_kernel.py`：结构化喷流 Fortran backend 的薄中间层，负责采样结构化参数、选择轴对称/非轴对称分支、调用 `structured_jet_1d` 并组装 API 结果。

Hadronic Python 模块只做 orchestration、wrapping 和正式 reference backend：

- Fortran wrappers：`hadronic_hummer.py`, `hadronic_bethe_heitler.py`, `hadronic_hadronic_ic.py`, `hadronic_pp.py`, `hadronic_pair_production.py`, `hadronic_species_transport.py`, `hadronic_secondary_radiation.py`, `hadronic_acceleration.py`。
- Reverse shock wrapper：`hadronic_reverse.py`；开启 RS full-chain flags 时，runtime 通过 formal 1D hadronic kernels 处理 RS seed photons、RS `B3`、shell energy 和 baryon target density。
- Reference/backend：`hadronic_pgamma.py`, `hadronic_am3_solver.py`, `hadronic_cascade.py`。

最终 AM3-derived microphysics 位于 `src/Hadronic/*.f90`。

## 4. Fortran 数值核

### 动力学

- `src/Dynamics/Dynamics_forward.f90`：正激波动力学、ISM/wind、density jumps、energy injection。
- `src/Dynamics/Dynamics_reverse.f90`：反激波动力学，含显式 region-3 `U3/V3/gamma34`；注入使用 shock-front `gamma34`，turbulent field 和 post-crossing thermal evolution 使用 `U3/V3`，可选 upstream `sigma_r` 添加 MHD-jump compression 与 ordered magnetic component。
- `src/Dynamics/dynamics_common.f90`：共享动力学辅助函数。

### 电子

- Main entries：`electron_forward_fullhide_1d.f90`, `electron_forward_transport_2d.f90`, `electron_forward_charint_1d.f90`, `electron_forward_slc1_1d.f90`, `electron_forward_t2g1_1d.f90`, `electron_forward_weno5_1d.f90`。`electron_forward_charint_2d` extension 复用 `electron_forward_transport_2d.f90` 中的 `fs_electron_transport_2d_core`，通过 `use_charint_transport` 选择 charint 2D path。
- Shared kernels：`electron_common.f90`, `electron_cooling_kernel.f90`, `electron_radiation_kernel.f90`, `electron_seed_history_kernel.f90`, `electron_transport_2d_kernel.f90`, `electron_injection_profiles.f90`, `electron_transport_common.f90`, `electron_reverse_kernel.f90`, `adaptive_resampling_mod.f90`。

### 辐射与插值

- `src/Radiation/radiation_ssc_spectrum.f90`：SSC spectrum 和 seed。
- `src/Radiation/radiation_reverse_seed.f90`：反激波同步辐射和 seed。
- `src/Radiation/radiation_gamma_gamma_absorption.f90`：gamma-gamma absorption。
- `src/Radiation/radiation_common.f90`：Simpson weights、power-law interpolation、pair cross-section、synchrotron seed core、transfer factor。
- `src/Radiation/synchrotron_polarization_kernel.f90`：频率相关同步辐射偏振 emissivity。
- `src/Radiation/quantum_synchrotron_kernel.f90`：quantum synchrotron helper。
- `src/Interpolation/SED_interpolation.f90`：observer-frame EATS/Doppler interpolation。
- `src/Interpolation/SED_interpolation_structured.f90`：structured jet interpolation。
- `src/Structured/structured_jet_1d.f90`：结构化喷流 Fortran 聚合入口，调度 theta 或 theta-phi 网格并复用现有动力学、电子、辐射、强子和 SED 插值核。

### 强子

`src/Hadronic/hadronic_forward_1d.f90` 是 forward-shock f2py entry point，调度：

- `hadronic_transport_kernel.f90`：proton injection、adiabatic/synchrotron loss、log-gamma energy advance。
- `hadronic_transport_remap_kernel.f90`：强子 transport 网格 remap helper。
- `hadronic_radiation_kernel.f90`：proton synchrotron。
- `hadronic_interaction_kernel.f90`：Hummer 2010 photopion operator。
- `hadronic_pgamma_hummer_1d.f90`：Hummer pγ 1D aggregate helper，供 formal hadronic 和 structured jet path 复用。
- `hadronic_decay_kernel.f90`：pi0 -> gamma、pi/mu decay、neutrino emissivity。
- `hadronic_pair_production_kernel.f90`：gamma-gamma pair production。
- `hadronic_pair_cascade_kernel.f90`：pair-cascade synchrotron kernel。
- `hadronic_pp_kernel.f90`、`hadronic_pp_models_kernel.f90`：pp source/loss helpers。
- `hadronic_bethe_heitler_kernel.f90`：Bethe-Heitler pair source 与 proton loss。
- `hadronic_hadronic_ic_kernel.f90`：hadronic inverse Compton。
- `hadronic_species_transport_kernel.f90`：neutron、pi±、mu± explicit transport。
- `hadronic_acceleration_kernel.f90`：acceleration timescale、injection operator、gamma_max estimate。
- `hadronic_secondary_radiation_kernel.f90`：pion/muon synchrotron 与 IC。
- `hadronic_common.f90`：共享常量、grid builders、validation。

Reverse-shock hadronic light entry 是 `src/Hadronic/hadronic_reverse_1d.f90`。Full-chain RS hadronic dispatch 通过 Python runtime wrapper 复用 `hadronic_forward_1d` formal 1D kernels，使用 RS magnetic field、RS seed photons、RS shell energy 和 RS baryon target density。

2D / chi-resolved hadronic transport 有意不实现。当前 `chi_grid` 属于 2D electron transport contract，而 `PhotonFieldState` 与 `HadronicSolution` 是 shell-level contracts。边界见 `doc/hadronic_chi_transport_decision.md`。

## 5. 强子当前状态

- 配置：`Radiation.epsilon_p`, `.proton_synch`, `.pg`, `.bethe_heitler`, `.hadronic_inverse_compton`, `.pp`, `.neutrino`, `.eta_acc`, `.pgamma_scheme`；`Setups.hadronic_enabled`, `.hadronic_solver`, `.num_gam_p`, `.num_nu_nu`。
- Solver names：`legacy_1d` 只覆盖 proton transport + proton synchrotron；`am3_1d` 是当前 formal hadronic main path。
- `pgamma_scheme`：`hummer_2010_response` 含 transport feedback；`ka2008_reference` 仅 emission benchmark；`disabled` 禁用。
- Pair cascade：`pair_cascade_iterations > 1` 选择 shell-sequence time-dependent gamma-gamma pair/synch cascade path；IC-mediated electromagnetic cascade 不属于当前契约。

## 6. 构建和测试

默认：WSL Ubuntu + uv，命令使用 `rtk` 前缀。

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_charint_2d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_forward_1d --force'
```

Smoke tests：`tests/readme_smoke_bench.py`, `tests/fullhide_2d_smoke_bench.py`, `tests/polarization_smoke.py`。

Hadronic regressions：`tests/hadronic_1d_smoke.py`, `tests/hadronic_species_transport_smoke.py`, `tests/hadronic_secondary_radiation_smoke.py`, `tests/hadronic_acceleration_smoke.py`, `tests/hadronic_bethe_heitler_smoke.py`, `tests/hadronic_hadronic_ic_smoke.py`, `tests/hadronic_pair_production_smoke.py`, `tests/hadronic_pp_smoke.py`。

Benchmark：`tests/vegas_afterglow_comparison.py`, `tests/sed_electron_compare.py`。Benchmark refresh protocol 见 `doc/benchmark_refresh_protocol.md`。

## 7. 已知边界

当前未完成项和 public/backend unsupported boundaries 集中维护在根目录 `TODO.md` 与 `doc/public_backend_limits.md`。本节只保留架构边界。

架构边界：

- ASGARD = shell-evolving blast-wave + observer projection。
- AM3 = microphysics/kernel reference，不替代 ASGARD dynamics/electron/observer chain。
- 最终 AM3-derived microphysics 写入 `src/Hadronic/*.f90`；Python 只做 orchestration。
- 非光滑物理时间演化优先视为 bug。
- Public/backend unsupported boundaries 固定在 `doc/public_backend_limits.md`。
