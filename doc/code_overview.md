# ASGARD 代码总览

本文档按当前工作树整理代码结构、运行主链和关键边界。电子算法细节见 `doc/electron_solver_algorithms.md`；物理过程到 Fortran 程序单元的交叉索引见 `doc/physics_algorithm_crosswalk.md`；逐文件程序单元索引见 `doc/fortran_kernel_index.md`；当前唯一 TODO / 未完成项入口见根目录 `TODO.md`。

## 1. 公开 API

- `asgard_core/api_model.py`：`Model`, `Medium`, `JetProfile`, `UniformMedium`, `WindMedium`, `TabulatedMedium`, `top_hat_jet`, `gaussian_jet`, `power_law_jet`, `Observer`, `Radiation`, `Numerics`, `ObserverGrid`, `SolverOptions`, `ReverseShock`, `Hadronic`。介质和喷流 public constructors 返回带 `kind` 标记的 `Medium` / `JetProfile`。`Model` 查询路径在本文件内完成 direct top-hat、structured Fortran backend 和 Python patch backend 调度，并直接构造内部 `RuntimeConfig`。
- `Model.flux_density_grid(times_s, nu_hz, projection_kind="lightcurve")`, `flux_density_grid_adaptive(...)`, `flux_density(...)`, `flux_density_exposures(...)`, `spectrum(time_s, nu_hz, projection_kind="sed")`, `flux(time_s, nu_min, nu_max, projection_kind="sed")`, `sky_image(t_obs, nu_obs, fov)`, `details()`。
- `Model.polarization(times_s, nu_hz, magnetic_geometry=..., local_emissivity=...)`。
- Hadronic public switches：`Radiation.pair_production`, `Radiation.include_pgamma`, `Radiation.bethe_heitler`, `Radiation.pp`, `Radiation.neutrino`, `Radiation.reverse_proton_energy_fraction`；cascade substeps 使用 `Hadronic.pair_cascade_iterations`。
- Electron-photon coupling switch：`SolverOptions.electron_photon_coupling="separated" | "joint"`；`joint` 是正向激波 1D formal 强子壳层级反馈路径，物理契约见 `doc/joint_secondary_feedback_physics.md`。
- Reverse-shock magnetization switch：`ReverseShock.upstream_sigma`，控制反向激波 upstream magnetization。
- `asgard_core/api_observe.py`：内部/旧配置观测工具，以及 `Model.sky_image(...)` / `Model.polarization(...)` 复用的实现函数；`observe(model, config)` 和 `run_fit(config)` 不从 `asgard_core` 顶层导出，不作为新教程的公开入口。
- `asgard_core/api_fit.py`：`Fitter`, `Param`, `FitResult`。
- Electron solver names：`fullhide_1d`, `fullhide_1d_hz`, `slc1_1d`, `charint_1d`, `dg_1d`, `charint_2d`, `t2g1_1d`, `weno5_1d`, `fullhide_2d`。public API 只使用这些完整名称。`fullhide_2d_pic` 没有跟踪源码和 `build_extensions.py` 构建登记，运行时映射已删除。
- `prompt/`：内部激波 snapshot 研究入口，不从 `asgard_core` 顶层导出。当前对象包括 `InternalShockShell`, `simulate_internal_shock`, `compute_prompt_observed_flux`，用于两壳碰撞、磁化 jump、FS/RS sync/SSC 和 prompt EATS 诊断。

## 2. 运行时主链

```text
Model.flux_density_grid / flux_density / spectrum / flux
  -> RuntimeConfig -> SimulationSetup
  -> solve_setup
  -> solve_dynamics -> solve_electron / joint electron-photon-hadronic stage
  -> solve_rsemission
  -> observer assembly -> Radiation.annihilation
  -> project_flux -> Interpolation.sed_interpolation[_chi] / structured chi ring projection -> API result
```

`Fitter` 是当前公开拟合入口；低层 `api_observe.run_fit(config)` 仅服务旧 `RuntimeConfig` 测试和内部工具。二者最终进入同一 `RuntimeConfig -> SimulationSetup -> solve_setup -> project_flux` 主链。

核心状态对象位于 `asgard_core/asgard_types.py`：

- `DynamicsSolution`：`r_tobs`, `r_gamma`, `radius`, `swept_mass_g`。
- `ReverseShockDynamics`：`M3`, total `B3`, ordered crossing field `B3_ordered_cross`, `U3/V3`, `gamma34` 和 crossing thermal records；`B3` 是 turbulent `sqrt(8 pi epsilon_B,r U3/V3)` 与可选 ordered upstream field 的总和。
- `ElectronSolution`：`gam_e`, `d_n_gam_e`, `l_syn_spec`, `seed_syn`；2D finite \(q\)-shell 额外包含 `d_n_gam_e_chi`, `chi_grid`, `l_syn_spec_chi`, `seed_syn_chi`, `tau_syn_chi`, `chi_radius_cm`, `chi_gamma_bulk`, `chi_dvolume_weight`。其中 `chi_grid` 是 \(q\) cell 的 BM 等效诊断坐标，observer projection 以 `chi_radius_cm`、`chi_gamma_bulk` 和 `chi_dvolume_weight` 为准；BH 额外包含 `d_n_gam_e_bh`。
- `PhotonFieldState`：forward synch seed、hadronic target field、absorption seed field。
- `HadronicSolution`：1D hadronic proton/secondary/radiation results；joint path 额外使用 `secondary_electron_source_r`、`tau_bh` 和 `bh_phloss` 做 shell-level feedback。
- `ObserverState`：absorption factors、`tau_pair`、flux components。
- `FluxComponents`：`total`、FS synch/SSC、hadronic、RS synch/SSC、cross-zone IC。

状态机位于 `asgard_core/asgard_state.py`：

- `solve_setup`：dynamics -> separated 或 joint forward stage -> reverse shock -> observer。
- `_photonfield`：复制 electron `seed_syn`；hadronic SSC seed 写入 target field。
- `_hadronstage`：调用 `solve_hadronic`；BH 次级 e± 并入 forward electron 后重算 `l_syn_spec/seed_syn`；pγ photon survival 写回 photon field。
- `_jointstage`：在同一 `R` 网格上迭代 electron、photon field、formal hadronic transport、BH/pp/gamma-gamma 二级 e± 源项和 photon survival；不使用 separated BH post-merge。
- `_observerstage`：组装 FS synch/SSC、RS synch/SSC、cross-zone IC 和 hadronic components；hadronic photons 使用 electron Fortran kernel 的 SSA transfer。
- `project_flux`：按 `projection_kind` 选择观测投影。`lightcurve` 是光变/拟合默认路径；`solver_options.geometry_projection="sed_adaptive_theta"` 对壳层级 EATS 使用 θ 自适应积分；`geometry_projection="chi_eats_2d"` 对 FS synchrotron+SSA 使用 χ 分辨 `sed_interpolation_chi`，并将非 χ 分量保持 shell-level projection。`sed` 是 `spectrum()` / `flux()` 默认路径，使用通用 shell SED 插值器或显式选择的 shell-level adaptive kernel。

拟合最短路径：

```text
Fitter.loglike -> compile_problem -> eval_loglike -> solve_setup
  -> project_flux -> combine_flux -> light_chi
```

## 3. Python 编排层

- `asgard_core/api_model.py`：public model objects、`Model` 查询缓存、direct/patch solve 调度、`Model -> RuntimeConfig` 适配和 details 打包。
- `asgard_core/api_observe.py`：`observe`, `run_fit` 兼容入口，以及 sky image / polarization / observation dataset helpers。
- `asgard_setup.py`：`RuntimeConfig -> SimulationSetup`。
- `asgard_config.py`：`RuntimeConfig`, `SimulationSetup`, `FitResult` 和 runtime config dataclasses；旧 compatibility shim 和配置 preset 中转层均已移除。
- `asgard_runtime.py`：backend selection、Fortran extension dispatch、array wrapping。
- `asgard_state.py`：主状态机、跨阶段耦合和 FS/RS cross-zone IC geometry。
- `asgard_postprocess.py`：observer projection、band aggregation、fit postprocessing 和观测数据 χ² helpers。
- `api_fit.py`：public `Fitter`、fit problem compilation 和 likelihood path。
- `asgard_types.py`：runtime dataclass contracts。
- `structured_jet_kernel.py`：结构化喷流 Fortran backend 的薄中间层，负责采样结构化参数、选择轴对称/非轴对称分支、调用 `structured_jet_1d` 并组装 API 结果。`fullhide_2d + chi_eats_2d` 的 axisymmetric structured 路径逐 theta ring 求解并预计算 chi-local synchrotron/SSA spectra，再调用 `sed_chi_ring` 做观测者投影；外层并行使用 POSIX `fork` 进程。
- `prompt/internal_shock.py`、`prompt/radiation.py`、`prompt/eats.py`：prompt internal-shock snapshot 的 Python orchestration。它复用现有 Fortran jump/electron/radiation/interpolation 核，不是 afterglow `Model` 主链的一部分。

强子 Python 模块只做编排、包装和轻量 helper：

- Fortran wrappers：`hadronic_processes.py`。
- Formal FS/RS shell-sequence transport 由 `src/Hadronic/hadronic_forward_1d.f90::formal_transport_1d` 推进；Python 只展开配置、传入数组并组装 `HadronicSolution`。
- Reverse shock light wrapper 已并入 `asgard_runtime.py`；开启 RS full-chain flags 时，runtime 通过 formal 1D 强子核处理 RS seed photons、RS `B3`、shell energy 和 baryon target density。
- Process/backend glue：`hadronic_am3_solver.py`, `hadronic_cascade.py`；pγ 单位转换和共享 wrapper 校验位于 `hadronic_processes.py`。

最终 AM3-derived microphysics 位于 `src/Hadronic/*.f90`。

## 4. Fortran 数值核

本节只列主要文件和阶段责任。需要逐个查 `module`、`subroutine`、`function` 时，使用 `doc/fortran_kernel_index.md`；需要从物理过程反查实现入口和验收指纹时，使用 `doc/physics_algorithm_crosswalk.md`。

### 动力学

- `src/Dynamics/Dynamics_forward.f90`：正向激波动力学、ISM/wind、density jumps、energy injection。
- `src/Dynamics/reverse_shock.f90`：反向激波 `Dynamics_reverse.dynamics_reverse` f2py 入口；first RS 推进、crossing event split 和 density-jump multiple RS 分支在同一输出循环内同步记录。
- `src/Dynamics/dynamics_density_profile.f90`：共享外介质密度、density jump/profile 状态和插值。
- 正向 RK4 推进留在 `Dynamics_forward.f90`；反向 event-split RK 推进位于 `reverse_shock.f90`。
- `src/Dynamics/reverse_shock_state.f90`、`src/Dynamics/reverse_shock_mhd_jump.f90`：反向激波 phase/state 下标和有限强度 MHD jump 公式。

### 电子

- Main entries：`electron_forward_fullhide_1d.f90`, `electron_forward_dg_1d.f90`, `electron_forward_transport_2d.f90`, `electron_forward_charint_1d.f90`, `electron_forward_slc1_1d.f90`, `electron_forward_t2g1_1d.f90`, `electron_forward_weno5_1d.f90`。`electron_forward_charint_2d` extension 复用 `electron_forward_transport_2d.f90` 中的 `fs_transport_2d`，通过 `use_charint_transport` 选择 charint 2D path。
- Shared kernels：`electron_common.f90`, `electron_cooling_kernel.f90`（门面/组装）, `electron_cooling_ssa_kernel.f90`, `electron_cooling_ic_kernel.f90`, `electron_cooling_y_kernel.f90`, `electron_radiation_kernel.f90`, `electron_seed_history_kernel.f90`, `electron_transport_2d_kernel.f90`, `electron_injection_profiles.f90`, `electron_shell_transport_common.f90`, `electron_transport_common.f90`, `electron_dg_transport.f90`, `electron_reverse_kernel.f90`, `adaptive_resampling_mod.f90`。

### 辐射与插值

- `src/Radiation/ssc_spectrum.f90`：SSC spectrum 和 seed。
- `src/Radiation/pair_absorption.f90`：gamma-gamma absorption。
- `src/Radiation/rad_common.f90`：Simpson weights、power-law interpolation、pair cross-section、synchrotron seed core、transfer factor。
- `src/Radiation/syn_polarization.f90`：频率相关同步辐射偏振 emissivity。
- `src/Radiation/quantum_synch.f90`：quantum synchrotron helper。
- `src/Interpolation/SED_interpolation.f90`：observer-frame EATS/Doppler interpolation，包含 shell-level、adaptive theta、top-hat chi、direct-electron chi 和 structured ring-precomputed chi 投影入口。
- `src/Interpolation/SED_interpolation_structured.f90`：`structured_jet_1d` 内部使用的 shell-level structured jet interpolation；不再从 `src.Interpolation` 暴露旧 Python 绑定。
- `src/Structured/structured_jet_1d.f90`：结构化喷流 Fortran 聚合入口，调度 theta 或 theta-phi 网格并复用现有动力学、电子、辐射、强子和 SED 插值核。

### 强子

`src/Hadronic/hadronic_forward_1d.f90` 是正向激波强子 f2py 入口。公开 f2py 面只保留 Python 正式路径调用的 drivers、process wrappers 和少量单位/插值 helper；其余 shell operators 保留为 Fortran 内部实现，不再从 `build_extensions.py` 导出。

- `formal_transport_1d`：formal 1D shell sequence driver，闭合 proton transport、pγ photon survival、BH、pp、secondary species、hadronic IC 和 BH/pp 二级电子序列。
- `hadronic_transport_kernel.f90`：proton injection、adiabatic/synchrotron loss、log-gamma energy advance。
- `hadronic_transport_remap_kernel.f90`：强子 transport 网格 remap helper。
- `hadronic_rad.f90`：proton synchrotron。
- `hadronic_pg.f90`：Hummer 2010 photopion operator。
- `hadronic_hummer.f90`：旧 Hummer pγ 1D aggregate helper；formal path 不使用其 photon escape 时间尺度。
- `hadronic_decay.f90`：pi0 -> gamma、pi/mu decay、neutrino emissivity。
- `hadronic_pair.f90`：gamma-gamma pair production。
- `hadronic_cascade.f90`：pair-cascade synchrotron kernel。
- `hadronic_pp.f90`、`pp_models.f90`：pp source/loss helpers。
- `hadronic_bh.f90`：Bethe-Heitler pair source 与 proton loss。
- `hadronic_ic.f90`：hadronic inverse Compton。
- `hadronic_species.f90`：neutron、pi±、mu± explicit transport。
- `hadronic_accel.f90`：acceleration timescale、injection rate、gamma limit estimate。
- `hadronic_secondary.f90`：pion/muon synchrotron 与 IC。
- `hadronic_base.f90`：共享常量、grid builders、validation。

反向激波强子 light entry 是 `src/Hadronic/hadronic_reverse_1d.f90`。Full-chain RS hadronic dispatch 通过 Python runtime wrapper 复用 `hadronic_forward_1d` formal 1D kernels，使用 RS magnetic field、RS seed photons、RS shell energy 和 RS baryon target density。

2D / \(\chi\) 分辨 hadronic transport 有意不实现。当前 `chi_grid` 属于 2D electron transport contract，而 `PhotonFieldState` 与 `HadronicSolution` 是壳层级 contracts。边界见 `doc/hadronic_chi_transport_decision.md`。

## 5. 强子当前状态

- 配置：`Radiation.proton_energy_fraction`, `.proton_synch`, `.include_pgamma`, `.bethe_heitler`, `.hadronic_inverse_compton`, `.pp`, `.neutrino`, `.acceleration_efficiency`, `.pgamma_scheme`；`Hadronic.enabled`, `.solver`, `.num_proton_gamma`, `.num_neutrino_frequency`。
- Solver names：`legacy_1d` 只覆盖 proton transport + proton synchrotron；`am3_1d` 是当前 formal hadronic main path。
- `pgamma_scheme`：`hummer_2010_response` 含 transport feedback；`disabled` 禁用。
- Pair cascade：`pair_cascade_iterations > 1` 选择 shell-sequence time-dependent \(\gamma\gamma\) pair/synch cascade path；IC-mediated electromagnetic cascade 不属于当前契约。
- Joint secondary feedback：`electron_photon_coupling="joint"` 在正向激波 1D 上把 BH/pp/\(\gamma\gamma\) 二级 \(e^\pm\) 作为外部 `Q_e,secondary,R` 输入电子方程，并把 \(p\gamma\)/BH/\(\gamma\gamma\) photon survival 作用到 joint photon field。详细算法见 `doc/joint_secondary_feedback_algorithm.md`。

## 6. 构建和测试

默认：WSL Ubuntu + uv，命令使用 `rtk` 前缀。

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_charint_2d --force'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && TMPDIR=/tmp uv run python build_extensions.py --module hadronic_forward_1d --force'
```

端到端验证入口不在代码概览页逐项列出；构建门槛、benchmark refresh 和当前已知验证阻塞集中见 `doc/validation_and_benchmarks.md`。

新增正式基准测试入口前必须先明确假设、决策价值和物理验收口径。

当前 finite q-shell diagnostic benchmark 入口：

- `tests/benchmark_theta_j_multiples_magnetic_decay.py`：1D thin shell 与 2D q-shell magnetic-decay 多观测角 lightcurve/spectrum 对照，2D 分支复用同一 solve state 重跑 projection。
- `tests/benchmark_skymap_centroid_motion.py`：1D thin shell 与 2D q-shell magnetic-decay 天图、flux centroid 和 apparent speed 诊断。

## 7. 已知边界

当前未完成项和公开 API/后端不支持边界集中维护在根目录 `TODO.md` 与 `doc/public_backend_limits.md`。本节只保留架构边界。

架构边界：

- ASGARD = 壳层演化爆波 + 观测者投影。
- AM3 = 微物理/数值核参考，不替代 ASGARD 动力学/电子/观测者主链。
- 最终来自 AM3 的微物理写入 `src/Hadronic/*.f90`；Python 只做编排。
- 非光滑物理时间演化优先视为 bug。
- 公开 API/后端不支持边界固定在 `doc/public_backend_limits.md`。
