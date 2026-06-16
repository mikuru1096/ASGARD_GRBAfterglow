# ASGARD 代码总览

本文档按当前工作树整理代码结构、运行主链和关键边界。电子算法细节见 `doc/electron_solver_algorithms.md`；当前唯一 TODO / 未完成项入口见根目录 `TODO.md`。

## 1. 公开 API

- `asgard_core/api_model.py`：`Model`, `Medium`, `JetProfile`, `UniformMedium`, `WindMedium`, `TabulatedMedium`, `top_hat_jet`, `gaussian_jet`, `power_law_jet`, `Observer`, `Radiation`, `Numerics`, `ObserverGrid`, `SolverOptions`, `ReverseShock`, `Hadronic`。介质和喷流 public constructors 返回带 `kind` 标记的 `Medium` / `JetProfile`。`Model` 查询路径在本文件内完成 direct top-hat、structured Fortran backend 和 Python patch backend 调度，并直接构造内部 `RuntimeConfig`。
- `Model.flux_density_grid(times_s, nu_hz, projection_kind="lightcurve")`, `flux_density(times_s, nu_hz, projection_kind="lightcurve")`, `spectrum(time_s, nu_hz, projection_kind="sed")`, `flux(time_s, nu_min, nu_max, projection_kind="sed")`, `sky_image(t_obs, nu_obs, fov)`, `details()`。
- `Model.polarization(times_s, nu_hz, magnetic_geometry=..., local_emissivity=...)`。
- Hadronic public switches：`Radiation.pair_production`, `Radiation.include_pgamma`, `Radiation.bethe_heitler`, `Radiation.pp`, `Radiation.neutrino`, `Radiation.reverse_proton_energy_fraction`；cascade substeps 使用 `Hadronic.pair_cascade_iterations`。
- Electron-photon coupling switch：`SolverOptions.electron_photon_coupling="separated" | "joint"`；`joint` 是正向激波 1D formal 强子壳层级反馈路径，物理契约见 `doc/joint_secondary_feedback_physics.md`。
- Reverse-shock magnetization switch：`ReverseShock.upstream_sigma`，控制反向激波 upstream magnetization。
- `asgard_core/api_observe.py`：内部/旧配置观测工具，以及 `Model.sky_image(...)` / `Model.polarization(...)` 复用的实现函数；`observe(model, config)` 和 `run_fit(config)` 不从 `asgard_core` 顶层导出，不作为新教程的公开入口。
- `asgard_core/api_fit.py`：`Fitter`, `Param`, `FitResult`。
- Electron solver names：`fullhide_1d`, `fullhide_1d_hz`, `slc1_1d`, `charint_1d`, `charint_2d`, `t2g1_1d`, `weno5_1d`, `fullhide_2d`, `fullhide_2d_pic`。public API 只使用这些完整名称。

## 2. 运行时主链

```text
Model.flux_density_grid / flux_density / spectrum / flux
  -> RuntimeConfig -> SimulationSetup
  -> solve_state_from_setup
  -> solve_dynamics -> solve_electron / joint electron-photon-hadronic stage
  -> solve_reverse_shock_emission
  -> observer assembly -> Radiation.annihilation
  -> project_flux_grid -> Interpolation.sed_interpolation[_chi] -> API result
```

`Fitter` 是当前公开拟合入口；低层 `api_observe.run_fit(config)` 仅服务旧 `RuntimeConfig` 测试和内部工具。二者最终进入同一 `RuntimeConfig -> SimulationSetup -> solve_state_from_setup -> project_flux_grid` 主链。

核心状态对象位于 `asgard_core/asgard_types.py`：

- `DynamicsSolution`：`r_tobs`, `r_gamma`, `radius`, `swept_mass_g`。
- `ReverseShockDynamics`：`M3`, total `B3`, ordered crossing field `B3_ordered_cross`, `U3/V3`, `gamma34` 和 crossing thermal records；`B3` 是 turbulent `sqrt(8 pi epsilon_B,r U3/V3)` 与可选 ordered upstream field 的总和。
- `ElectronSolution`：`gam_e`, `d_n_gam_e`, `l_syn_spec`, `seed_syn`, `nu_m`, `nu_c`, `nu_a`；2D 额外包含 `d_n_gam_e_chi`, `chi_grid`, `l_syn_spec_chi`, `seed_syn_chi`, `tau_syn_chi`, `chi_radius_cm`, `chi_gamma_bulk`, `chi_dvolume_weight`；BH 额外包含 `d_n_gam_e_bh`。
- `PhotonFieldState`：forward synch seed、hadronic target field、absorption seed field。
- `HadronicSolution`：1D hadronic proton/secondary/radiation results；joint path 额外使用 `secondary_electron_source_r`、`tau_bh` 和 `bh_photon_loss_rate` 做 shell-level feedback。
- `ObserverState`：absorption factors、`tau_pair`、flux components。
- `FluxComponents`：`total`、FS synch/SSC、hadronic、RS synch/SSC、cross-zone IC。

状态机位于 `asgard_core/asgard_state.py`：

- `solve_state_from_setup`：dynamics -> separated 或 joint forward stage -> reverse shock -> observer。
- `_build_photon_field_stage`：复制 electron `seed_syn`；hadronic SSC seed 写入 target field。
- `_solve_hadronic_stage`：调用 `solve_hadronic`；BH 次级 e± 并入 forward electron 后重算 `l_syn_spec/seed_syn`；pγ photon survival 写回 photon field。
- `_solve_joint_forward_stage`：在同一 `R` 网格上迭代 electron、photon field、formal hadronic transport、BH/pp/gamma-gamma 二级 e± 源项和 photon survival；不使用 separated BH post-merge。
- `_assemble_observer_stage`：组装 FS synch/SSC、RS synch/SSC、cross-zone IC 和 hadronic components；hadronic photons 使用 electron Fortran kernel 的 SSA transfer。
- `project_flux_grid`：按 `projection_kind` 选择观测投影。`lightcurve` 是光变/拟合默认路径；当 public API 设置 `solver_options.geometry_projection="chi_eats_2d"` 时，底层 `geometry_kernel` 对 FS synchrotron+SSA 使用 χ 分辨 `sed_interpolation_chi`，并将非 χ 分量保持 shell-level projection。`sed` 是 `spectrum()` / `flux()` 默认路径，使用通用 shell SED 插值器。

拟合最短路径：

```text
Fitter.loglike -> compile_problem -> eval_loglike -> solve_state_from_setup
  -> project_flux_grid -> combine_multiband_flux -> compute_light_curve_redchi
```

## 3. Python 编排层

- `asgard_core/api_model.py`：public model objects、`Model` 查询缓存、direct/patch solve 调度、`Model -> RuntimeConfig` 适配和 details 打包。
- `asgard_core/api_observe.py`：`observe`, `run_fit` 兼容入口，以及 sky image / polarization / observation dataset helpers。
- `asgard_setup.py`：`RuntimeConfig -> SimulationSetup`。
- `asgard_config.py`：`RuntimeConfig`, `SimulationSetup`, `FitResult` 和 runtime config dataclasses；旧 compatibility shim 和配置 preset 中转层均已移除。
- `asgard_runtime.py`：backend selection、Fortran extension dispatch、array wrapping。
- `asgard_state.py`：主状态机和跨阶段耦合。
- `asgard_ssc.py`：forward SSC auxiliary grid 与 seed。
- `asgard_coupling.py`：FS/RS cross-zone IC geometry 与 seed field coupling。
- `asgard_postprocess.py`：observer projection、band aggregation、fit postprocessing。
- `asgard_fit.py`：fit problem compilation 和 likelihood path。
- `asgard_types.py`：runtime dataclass contracts。
- `structured_jet_kernel.py`：结构化喷流 Fortran backend 的薄中间层，负责采样结构化参数、选择轴对称/非轴对称分支、调用 `structured_jet_1d` 并组装 API 结果。

强子 Python 模块只做编排、包装和正式参考后端：

- Fortran wrappers：`hadronic_hummer.py`, `hadronic_bethe_heitler.py`, `hadronic_hadronic_ic.py`, `hadronic_pp.py`, `hadronic_pair_production.py`, `hadronic_species_transport.py`, `hadronic_secondary_radiation.py`, `hadronic_acceleration.py`。
- Reverse shock wrapper：`hadronic_reverse.py`；开启 RS full-chain flags 时，runtime 通过 formal 1D 强子核处理 RS seed photons、RS `B3`、shell energy 和 baryon target density。
- Reference/backend：`hadronic_pgamma.py`, `hadronic_am3_solver.py`, `hadronic_cascade.py`。

最终 AM3-derived microphysics 位于 `src/Hadronic/*.f90`。

## 4. Fortran 数值核

### 动力学

- `src/Dynamics/Dynamics_forward.f90`：正向激波动力学、ISM/wind、density jumps、energy injection。
- `src/Dynamics/Dynamics_reverse.f90`：反向激波动力学，含显式 region-3 `U3/V3/gamma34`；注入使用 shock-front `gamma34`，turbulent field 和 post-crossing thermal evolution 使用 `U3/V3`，可选 upstream `sigma_r` 添加 MHD-jump compression 与 ordered magnetic component。
- `src/Dynamics/dynamics_common.f90`：共享动力学辅助函数。

### 电子

- Main entries：`electron_forward_fullhide_1d.f90`, `electron_forward_transport_2d.f90`, `electron_forward_charint_1d.f90`, `electron_forward_slc1_1d.f90`, `electron_forward_t2g1_1d.f90`, `electron_forward_weno5_1d.f90`。`electron_forward_charint_2d` extension 复用 `electron_forward_transport_2d.f90` 中的 `fs_electron_transport_2d_core`，通过 `use_charint_transport` 选择 charint 2D path。
- Shared kernels：`electron_common.f90`, `electron_cooling_kernel.f90`, `electron_radiation_kernel.f90`, `electron_seed_history_kernel.f90`, `electron_transport_2d_kernel.f90`, `electron_injection_profiles.f90`, `electron_transport_common.f90`, `electron_reverse_kernel.f90`, `adaptive_resampling_mod.f90`。

### 辐射与插值

- `src/Radiation/radiation_ssc_spectrum.f90`：SSC spectrum 和 seed。
- `src/Radiation/radiation_reverse_seed.f90`：反向激波同步辐射和 seed。
- `src/Radiation/radiation_gamma_gamma_absorption.f90`：gamma-gamma absorption。
- `src/Radiation/radiation_common.f90`：Simpson weights、power-law interpolation、pair cross-section、synchrotron seed core、transfer factor。
- `src/Radiation/synchrotron_polarization_kernel.f90`：频率相关同步辐射偏振 emissivity。
- `src/Radiation/quantum_synchrotron_kernel.f90`：quantum synchrotron helper。
- `src/Interpolation/SED_interpolation.f90`：observer-frame EATS/Doppler interpolation。
- `src/Interpolation/SED_interpolation_structured.f90`：structured jet interpolation。
- `src/Structured/structured_jet_1d.f90`：结构化喷流 Fortran 聚合入口，调度 theta 或 theta-phi 网格并复用现有动力学、电子、辐射、强子和 SED 插值核。

### 强子

`src/Hadronic/hadronic_forward_1d.f90` 是正向激波强子 f2py 入口，调度：

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

反向激波强子 light entry 是 `src/Hadronic/hadronic_reverse_1d.f90`。Full-chain RS hadronic dispatch 通过 Python runtime wrapper 复用 `hadronic_forward_1d` formal 1D kernels，使用 RS magnetic field、RS seed photons、RS shell energy 和 RS baryon target density。

2D / \(\chi\) 分辨 hadronic transport 有意不实现。当前 `chi_grid` 属于 2D electron transport contract，而 `PhotonFieldState` 与 `HadronicSolution` 是壳层级 contracts。边界见 `doc/hadronic_chi_transport_decision.md`。

## 5. 强子当前状态

- 配置：`Radiation.proton_energy_fraction`, `.proton_synch`, `.include_pgamma`, `.bethe_heitler`, `.hadronic_inverse_compton`, `.pp`, `.neutrino`, `.acceleration_efficiency`, `.pgamma_scheme`；`Hadronic.enabled`, `.solver`, `.num_proton_gamma`, `.num_neutrino_frequency`。
- Solver names：`legacy_1d` 只覆盖 proton transport + proton synchrotron；`am3_1d` 是当前 formal hadronic main path。
- `pgamma_scheme`：`hummer_2010_response` 含 transport feedback；`ka2008_reference` 仅 emission benchmark；`disabled` 禁用。
- Pair cascade：`pair_cascade_iterations > 1` 选择 shell-sequence time-dependent \(\gamma\gamma\) pair/synch cascade path；IC-mediated electromagnetic cascade 不属于当前契约。
- Joint secondary feedback：`electron_photon_coupling="joint"` 在正向激波 1D 上把 BH/pp/\(\gamma\gamma\) 二级 \(e^\pm\) 作为外部 `Q_e,secondary,R` 输入电子方程，并把 \(p\gamma\)/BH/\(\gamma\gamma\) photon survival 作用到 joint photon field。详细算法见 `doc/joint_secondary_feedback_algorithm.md`。

## 6. 构建和测试

默认：WSL Ubuntu + uv，命令使用 `rtk` 前缀。

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_charint_2d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_forward_1d --force'
```

Smoke tests：`tests/readme_smoke_bench.py`, `tests/fullhide_2d_smoke_bench.py`。

Hadronic regressions：`tests/hadronic_1d_smoke.py`, `tests/hadronic_species_transport_smoke.py`, `tests/hadronic_secondary_radiation_smoke.py`, `tests/hadronic_acceleration_smoke.py`, `tests/hadronic_bethe_heitler_smoke.py`, `tests/hadronic_hadronic_ic_smoke.py`, `tests/hadronic_pair_production_smoke.py`, `tests/hadronic_pp_smoke.py`。

新增正式 benchmark 入口前必须先明确假设、决策价值和物理验收口径。

## 7. 已知边界

当前未完成项和 public/backend unsupported boundaries 集中维护在根目录 `TODO.md` 与 `doc/public_backend_limits.md`。本节只保留架构边界。

架构边界：

- ASGARD = 壳层演化爆波 + 观测者投影。
- AM3 = 微物理/数值核参考，不替代 ASGARD dynamics/electron/observer chain。
- 最终 AM3-derived microphysics 写入 `src/Hadronic/*.f90`；Python 只做编排。
- 非光滑物理时间演化优先视为 bug。
- Public/backend unsupported boundaries 固定在 `doc/public_backend_limits.md`。
