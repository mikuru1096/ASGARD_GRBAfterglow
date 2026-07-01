# Fortran 程序单元索引

本文是当前工作树的 Fortran kernel 索引。它面向需要逐个进入子程序读算法的人：先看 f2py 入口和物理阶段，再进入文件内的 `module`、`subroutine`、`function`。更高层的物理到算法映射见 `doc/physics_algorithm_crosswalk.md`，运行主链见 `doc/call_chain.md`。

当前索引按 ASGARD 自有 Fortran 数值核抽取，排除第三方固定格式特殊函数依赖，共 808 个程序单元：35 个 module、565 个 subroutine、208 个 function。行号是生成本页时的源文件位置。

## 读源码顺序

1. 从 public Python API 进入 `asgard_core/asgard_state.py`，确认阶段顺序。
2. 在本页查 f2py module 和主 entry。
3. 进入对应 Fortran 文件，先读入口参数和数组形状，再读内部 contained procedures。
4. 对涉及真实物理演化的改动，回到 `doc/physics_algorithm_crosswalk.md` 查验收指纹。

不要把内部 helper 当 public ABI。稳定边界是 Python public API、`build_extensions.py` 登记的 f2py module、以及本页列出的主 entry。

## f2py 模块和主入口

| Build module | CWD | Source closure | 主 entry | 物理/算法角色 |
| --- | --- | --- | --- | --- |
| `Dynamics_forward` | `src/Dynamics` | `Constants + dynamics_common + Dynamics_forward` | `dynamics_forward` | 正向激波动力学、ISM/wind、密度跳变和能量注入。 |
| `Dynamics_reverse` | `src/Dynamics` | `Constants + dynamics_common + reverse_jump_conditions + reverse_rhs + Dynamics_reverse` | `dynamics_reverse` | 反向激波 crossing、region-3 thermal state、磁化 jump 和次级 RS 分支。 |
| `electron_forward_fullhide_1d` | `src/Electron` | `ELECTRON_COMMON_SOURCES + electron_forward_fullhide_1d` | `fs_electron_fullhide_1d; fs_electron_fullhide_1d_coupled` | 默认 1D 电子输运和 joint feedback coupled pass。 |
| `electron_forward_fullhide_1d_hybrid` | `src/Electron` | `ELECTRON_HISTORY_SOURCES_HZ + electron_forward_fullhide_1d_hybrid` | `fs_electron_fullhide_1d_hz` | 热/非热混合谱路径。 |
| `electron_forward_charint_1d` | `src/Electron` | `ELECTRON_COMMON_SOURCES + electron_forward_charint_1d` | `fs_electron_charint_1d` | 1D 特征线输运对照。 |
| `electron_forward_dg_1d` | `src/Electron` | `ELECTRON_DG_1D_SOURCES + electron_forward_dg_1d` | `fs_electron_dg_1d` | P12 LGL-DG 正向电子输运。 |
| `electron_forward_charint_2d` | `src/Electron` | `ELECTRON_2D_SOURCES + electron_forward_transport_2d` | `fs_electron_transport_2d_core` | finite-q shell 2D 电子输运；别把 chi_grid 当强子局域坐标。 |
| `electron_forward_t2g1_1d` | `src/Electron` | `ELECTRON_COMMON_SOURCES + electron_forward_t2g1_1d` | `fs_electron_t2g1_1d` | 方法比较电子输运。 |
| `electron_forward_weno5_1d` | `src/Electron` | `ELECTRON_COMMON_SOURCES + electron_forward_weno5_1d` | `fs_electron_weno5_1d` | WENO5 方法比较电子输运。 |
| `electron_reverse_kernel` | `src/Electron` | `ELECTRON_REVERSE_SOURCES + electron_reverse_kernel` | `electron_reverse_evolve; electron_secondary_reverse_evolve` | RS primary/secondary electron transport 和 reacceleration。 |
| `electron_radiation` | `src/Electron` | `ELECTRON_RADIATION_SOURCES` | `get_syn_*; get_nu_a_*` | 电子同步/SSA/seed 低层核；通常经其他 entry 调用。 |
| `radiation_ssc_spectrum` | `src/Radiation` | `radiation_common + radiation_ssc_spectrum` | `ssc_spec; ssc_spec_nonuniform` | SSC spectrum 和 KN/Jones 积分。 |
| `radiation_gamma_gamma_absorption` | `src/Radiation` | `radiation_common + radiation_gamma_gamma_absorption` | `annihilation` | 观测侧 gamma-gamma absorption。 |
| `SED_interpolation` | `src/Interpolation` | `radiation_common + interpolation_common + SED_interpolation` | `sed_interpolation; sed_interpolation_adaptive_theta; sed_interpolation_chi; sed_interpolation_chi_electron_cached; sed_interpolation_chi_structured_axisym_ring_precomputed` | 观测者投影；Python 公开绑定集中在此模块。 |
| `SED_interpolation_structured` | `src/Interpolation` | `interpolation_common + SED_interpolation_structured` | `sed_interpolation_structured; sed_interpolation_structured_phi` | `structured_jet_1d` 内部 shell-level structured projection；不再从 `src.Interpolation` 暴露旧 Python 绑定。 |
| `hadronic_forward_1d` | `src/Hadronic` | `HADRONIC_1D_SOURCES + hadronic_forward_1d` | `fs_hadronic_formal_transport_1d; shell operators` | formal 1D proton/secondary/photon-loss shell sequence。 |
| `hadronic_reverse_1d` | `src/Hadronic` | `HADRONIC_1D_SOURCES + hadronic_reverse_1d` | `fs_hadronic_reverse_1d` | RS light proton transport + proton synchrotron。 |
| `structured_jet_1d` | `src/Structured` | `STRUCTURED_JET_1D_SOURCES + structured_jet_1d` | `structured_jet_flux_1d` | 结构化喷流 Fortran 聚合入口。 |

Fortran 改动后的最低门槛见 `doc/validation_and_benchmarks.md`。文档-only 改动不需要重建扩展；若修改本页对应源码，必须跑受影响 module 的 `build_extensions.py --force`、独立 `-Wline-truncation` 源闭包检查和最小 smoke。

## 目录级算法地图

| 目录 | 责任 | 不能混淆的边界 |
| --- | --- | --- |
| `src/Dynamics/` | 生成 `R`、`Gamma`、swept mass、RS region-3 thermal state。 | 动力学不平滑时先查这里，不在电子或投影层修光变。 |
| `src/Electron/` | 电子注入、冷却、1D/2D 输运、同步辐射 seed 和 RS 电子。 | `chi_grid` 是 finite-q 诊断坐标，不是 hadronic local transport 坐标。 |
| `src/Radiation/` | SSC、gamma-gamma absorption、polarization kernel 和共享 photon transfer。 | Observer luminosity 不能自动当 local photon density；必须有 shell volume 和 escape time。 |
| `src/Hadronic/` | formal 1D proton/secondary/photon-loss shell sequence。 | 2D/chi-resolved hadronic transport 当前未实现。 |
| `src/Interpolation/` | EATS、Doppler、redshift、SSA survival 和 structured projection。 | Observer projection 只读 transport state，不回写粒子或 photon field。 |
| `src/Structured/` | patch/axis/phi 聚合并复用现有 kernel。 | patch-local hadronic feedback 未定义为独立物理合同。 |

## 源文件与程序单元

### `src/Constants.f90`

全局物理常数和单位常量。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `Constants` | 模块命名空间；集中声明本文件共享 procedure。 |

### `src/Dynamics/Dynamics_forward.f90`

正向激波动力学入口和右端项。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 1 | `dynamics_forward` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 42 | `forward_dynamics_rhs` | ODE/PDE 右端项；定义物理源汇和动力学变量导数。 |

### `src/Dynamics/Dynamics_reverse.f90`

反向激波主入口、穿越事件、磁化 jump、次级反向激波分支。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 5 | `dynamics_reverse` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 171 | `advance_reverse_state_to_target` | 反向激波 waiting/pre-crossing/post-crossing 时间推进。 |
| `F` | 211 | `reverse_shock_pressure_ready_state` | 反向激波 pressure-balance 和 fast-mode shock 条件。 |
| `S` | 238 | `waiting_trial` | 磁化 RS 尚未成 shock 时的 waiting 分支 trial 推进。 |
| `S` | 259 | `locate_waiting_event` | 定位 waiting 到 shock-ready 的事件时间。 |
| `S` | 285 | `update_secondary_reverse_events` | 密度增强窗口内扫描次级 RS start/end。 |
| `S` | 337 | `secondary_reverse_event_root_between` | 次级 RS source 过零半径定位。 |
| `S` | 369 | `secondary_reverse_event_source` | 次级 RS mechanical source，比较新 shock 与上游绝热态。 |
| `S` | 411 | `store_secondary_branch_state` | 输出次级 RS 分支热、磁、注入和诊断量。 |
| `F` | 515 | `secondary_parent_upstream_available` | 判断上一级次级分支是否可作为当前上游。 |
| `S` | 525 | `secondary_reverse_density_branch_state` | density-jump 分支的局部密度和权重。 |

### `src/Dynamics/reverse_rhs.f90`

反向激波 RHS，显式处理 shock 注入、磁化惯性、区域 3 热演化和次级分支导数。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `reverse_dynamics_rhs` | ODE/PDE 右端项；定义反向激波主状态和次级分支导数。 |
| `F` | 201 | `secondary_parent_upstream_available` | 判断上一级次级分支是否可作为当前上游。 |
| `S` | 211 | `compute_secondary_branch_derivatives` | 次级 RS 分支质量、热能、体积和注入诊断导数。 |
| `S` | 291 | `secondary_reverse_density_branch_rhs` | density-jump 分支的局部密度和 source 权重。 |

### `src/Dynamics/dynamics_common.f90`

动力学共享介质、MHD jump、RK4 和事件分裂 primitive。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `dynamics_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 27 | `dynamics_forward_rhs_iface` | ODE/PDE 右端项；定义物理源汇和动力学变量导数。 |
| `S` | 42 | `dynamics_reverse_rhs_iface` | ODE/PDE 右端项；定义物理源汇和动力学变量导数。 |
| `S` | 63 | `dynamics_external_density_profile` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `S` | 119 | `dynamics_set_density_jump_profile` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `S` | 171 | `dynamics_density_tabulated_profile` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `F` | 205 | `rs_vegas_ud` | 有限强度 MHD jump 解析根；不要用 ultra-relativistic 近似替代。 |
| `F` | 239 | `rs_vegas_comp` | 有限强度 MHD jump 压缩比；sigma->0 极限回到 hydrodynamic baseline。 |
| `F` | 261 | `rs_mag_specific_internal` | MHD jump 下游热比内能；保持 crossing 前后和 sigma->0 极限连续。 |
| `S` | 282 | `dynamics_rk4_error_n` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 305 | `dynamics_rk4_reverse_error` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 328 | `dynamics_rk4_forward` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 413 | `dynamics_rk4_reverse_plain_step` | 单个 log-time RK4 子步；只负责给定 phase 的 RHS 采样。 |
| `S` | 459 | `dynamics_rk4_reverse_pre_step` | crossing 前以 M3 fraction 为自变量的 RK4 子步。 |
| `S` | 524 | `dynamics_rk4_reverse_pre_m3` | pressure-ready 到 crossing/目标时刻的 pre-M3 推进。 |
| `S` | 645 | `dynamics_rk4_reverse` | `select case` 展示 waiting、pre-crossing event split 和 post-crossing 推进。 |

### `src/Electron/adaptive_resampling_mod.f90`

诊断/快速谱用自适应重采样工具。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `adaptive_resampling_mod` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 9 | `adaptive_resampling_log` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 130 | `moving_average` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 152 | `unique_sorted` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 174 | `sort_integers` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 191 | `supplement_indices` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |

### `src/Electron/electron_common.f90`

电子边界向量解包、注入断点、初始谱、诊断公共工具。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 19 | `electron_unpack_boundary` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 46 | `electron_initialize_spectrum` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 76 | `electron_gamma_m_exact` | 特征 Lorentz 因子/断点诊断；用于注入、冷却或活动网格边界。 |
| `S` | 103 | `electron_gamma_c_from_loss_mean` | 特征 Lorentz 因子/断点诊断；用于注入、冷却或活动网格边界。 |
| `S` | 149 | `electron_injection_prefactor` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 158 | `electron_source_bounds` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 190 | `electron_prepare_radiation_spectrum` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 268 | `electron_max_relative_error` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 284 | `electron_initial_density` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |

### `src/Electron/electron_cooling_ssa_kernel.f90`

同步自吸收回热/冷却核，含 seed 频率缓存、SSA 几何映射和 χ 批量版本。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_cooling_ssa_kernel` | SSA 物理核模块；不承担 IC/Y 或冷却项组装分派。 |
| `S` | 24 | `ensure_ssa_geometry_workspace` | workspace/cache 管理；只服务性能和内存复用，不改变物理语义。 |
| `S` | 40 | `ensure_ssa_seed_cache` | workspace/cache 管理；只服务性能和内存复用，不改变物理语义。 |
| `S` | 74 | `advance_ssa_seed_cursor` | SSA seed 频率窗口推进 primitive。 |
| `S` | 95 | `build_ssa_geometry` | 构建 SSA 低频/高频截面区间和 prefactor。 |
| `S` | 138 | `electron_cooling_ssa_loss_batch` | χ-resolved seed photon fields 的 SSA 算子；单列调用传 `Num_chi=1`。 |
| `S` | 195 | `accumulate_ssa_batch_gamma` | `electron_cooling_ssa_loss_batch` 的 gamma/χ 局部累加。 |
| `F` | 265 | `clipped_ssa_batch_segment` | χ 批量路径中被上下限裁剪的 SSA cell 积分。 |

### `src/Electron/electron_cooling_ic_kernel.f90`

逆康普顿/KN 冷却核，含 IC 网格缓存、旧积分冷却率和 emissivity-budget 冷却率。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_cooling_ic_kernel` | IC 物理核模块；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 17 | `ensure_ic_grid_cache` | workspace/cache 管理；缓存 seed 频率中点、间距和电子能量中点。 |
| `F` | 26 | `ic_grid_cache_current` | IC cache 总命中判断。 |
| `F` | 38 | `ic_seed_grid_current` | IC seed 频率网格命中判断。 |
| `F` | 50 | `ic_gamma_grid_current` | IC 电子能量中点网格命中判断。 |
| `S` | 62 | `rebuild_ic_grid_cache` | 重建 IC 派生网格量。 |
| `S` | 83 | `electron_cooling_ic_loss` | 数值 IC/KN 冷却率双重积分。 |
| `S` | 107 | `accumulate_ic_gamma_loss` | `electron_cooling_ic_loss` 的单 gamma 积分累加。 |
| `S` | 148 | `electron_cooling_ic_loss_emissivity_budget` | 与 radiation SSC emissivity 核一致的 IC cooling budget。 |
| `S` | 183 | `accumulate_budget_gamma` | emissivity-budget 路径的单 gamma 累加；低/高 seed 分支在同一循环内显式计算。 |
| `F` | 211 | `low_seed_kernel` | Jones/KN 低 seed 侧核函数。 |

### `src/Electron/electron_cooling_y_kernel.f90`

Compton-Y 辅助量核，含 Nakar 数值积分和 Fan 解析分段模型。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_cooling_y_kernel` | Compton-Y 物理核模块；不承担 SSA/IC 回热积分或最终冷却组装。 |
| `S` | 20 | `ensure_y_nakar_workspace` | workspace/cache 管理；缓存 hat-nu、频率段 Gauss 节点和查找区间。 |
| `S` | 61 | `electron_cooling_y_nakar` | Nakar+2009 Compton-Y 数值积分。 |
| `S` | 102 | `electron_cooling_y_fan` | Fan+2008 Compton-Y 解析分段模型。 |

### `src/Electron/electron_cooling_kernel.f90`

电子冷却门面模块；只保留 batch auxiliary、assemble 和 get_forward_cooling 组装入口，具体 SSA/IC/Y 物理核从对应实现模块直接调用。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_cooling_kernel` | 冷却组装门面模块；不承载 SSA/IC/Y 具体物理核。 |
| `S` | 14 | `prepare_forward_cooling_aux_batch` | χ-resolved 冷却辅助量准备；单列调用传 `Num_chi=1`。 |
| `S` | 38 | `assemble_forward_cooling_split_batch` | χ-resolved 正向激波冷却率组装；单列调用传 `Num_chi=1`。 |
| `S` | 64 | `assemble_forward_cooling_from_terms` | 由 synch、SSA、IC/Y 辅助项组装最终 dγ/dR。 |
| `S` | 100 | `get_forward_cooling` | 正向激波 1D 冷却主入口；内部走 batch 单列路径。 |

### `src/Electron/electron_energy_coordinate_common.f90`

log-gamma 与 log-four-velocity 坐标映射。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_energy_coordinate_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `F` | 17 | `electron_coord_from_xgamma` | 按坐标类型把 `x=log10(gamma)` 映射到输运坐标；four-velocity 分支直接写出 `log10(1+u^2/u_s^2)`。 |
| `F` | 32 | `electron_xgamma_from_coord` | 按坐标类型从输运坐标恢复 `log10(gamma)`。 |
| `F` | 44 | `electron_gamma_from_coord` | 按坐标类型从输运坐标恢复电子 Lorentz factor。 |
| `F` | 56 | `electron_dxgamma_dcoord` | 计算坐标变换 Jacobian `d log10(gamma) / d coord`。 |
| `S` | 70 | `electron_build_four_velocity_grid` | 在线性 four-velocity 坐标上构造 gamma 中心和 `x=log10(gamma)` 边界。 |

### `src/Electron/electron_forward_charint_1d.f90`

1D 特征线电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `fs_electron_charint_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 235 | `prepare_characteristic_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 275 | `write_final_characteristic_diagnostics` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_forward_dg_1d.f90`

1D LGL-DG 电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `fs_electron_dg_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 95 | `initialize_forward_four_velocity_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 110 | `prepare_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 129 | `write_radiation_and_breaks` | 特征 Lorentz 因子/断点诊断；用于注入、冷却或活动网格边界。 |
| `S` | 141 | `remesh_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 171 | `advance_substep` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 204 | `prepare_forward_dg_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 223 | `limit_density_jump_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 243 | `limit_one_density_jump` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `F` | 266 | `density_jump_log_slope` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `S` | 284 | `add_density_jump_derivative` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `F` | 296 | `forward_dg_source_upper_xmax` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 304 | `forward_dg_active_xmax` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 318 | `scale_dg_state_to_grid_content` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 331 | `write_positive_output` | DG 状态投影回用户 gamma 网格，并在输出边界显式施加非负电子谱。 |

### `src/Electron/electron_forward_fullhide_1d.f90`

默认 1D 隐式有限体积电子输运入口，含 joint coupled 入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 5 | `fs_electron_fullhide_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 123 | `initialize_forward_four_velocity_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 138 | `prepare_fullhide_shell` | 单个 FS shell 的密度、磁场、注入能标、同步谱、SSA break 和冷却率准备。 |
| `S` | 372 | `fs_electron_fullhide_1d_coupled` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 461 | `prepare_coupled_shell` | joint electron-photon shell 的辐射场、SSA break 和耦合冷却率准备。 |

### `src/Electron/electron_forward_fullhide_1d_hybrid.f90`

混合热/非热谱形入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 7 | `fs_electron_fullhide_1d_hz` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 266 | `build_hybrid_or_powerlaw_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 281 | `build_hybrid_or_powerlaw_source_from_count` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |

### `src/Electron/electron_forward_slc1_1d.f90`

半拉格朗日一阶电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `fs_electron_slc1_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |

### `src/Electron/electron_forward_t2g1_1d.f90`

T2G1 方法比较电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 8 | `fs_electron_t2g1_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 93 | `prepare_t2g1_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 118 | `write_t2g1_radiation_and_cooling` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 142 | `advance_t2g1_substep` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |

### `src/Electron/electron_forward_transport_2d.f90`

2D finite-q shell 电子输运和投影预处理入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `fs_electron_transport_2d_core` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 623 | `reduce_syn_shell_from_q` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 639 | `project_q_projection_shell` | finite-q shell 几何或 chi-equivalent 投影字段。 |

### `src/Electron/electron_forward_weno5_1d.f90`

WENO5 方法比较电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `fs_electron_weno5_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 77 | `prepare_weno_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 99 | `write_weno_radiation_and_cooling` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 121 | `advance_weno_substep` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 150 | `load_weno_extended_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 163 | `compute_weno_fluxes` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 179 | `advance_weno_rk_stage` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 204 | `weno5_update_ghost_cells` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 219 | `weno5_positive_flux` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 259 | `weno5_negative_flux` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_injection_profiles.f90`

非热/热电子注入谱、log-cell 积分和源项归一化。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_injection_profiles` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 13 | `electron_exp_cutoff_factor` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 22 | `electron_initial_powerlaw_params` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 47 | `electron_dnx_powerlaw_cutoff_value` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 64 | `electron_dnx_gauss3_integral` | 积分权重或求积 primitive；影响谱积分精度。 |
| `F` | 88 | `electron_dny_gauss3_integral` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 114 | `electron_add_dnx_segment` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 135 | `electron_add_dny_segment` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 156 | `electron_profile_log_cell_edges` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `S` | 171 | `electron_initial_powerlaw_exp_cutoff` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 190 | `electron_initial_powerlaw_exp_cutoff_edges` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 228 | `electron_initial_powerlaw_exp_cutoff_coord_edges` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 267 | `electron_build_source_term_exp_cutoff_edges` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 292 | `electron_build_source_term_exp_cutoff_coord_edges` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 317 | `electron_build_kinetic_source_term_exp_cutoff_edges` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 355 | `electron_build_kinetic_source_term_exp_cutoff_coord_edges` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 391 | `electron_thermal_theta` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 399 | `electron_build_thermal_shape_dnx` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 424 | `electron_add_thermal_source_term` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 437 | `electron_add_thermal_population` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_radiation_kernel.f90`

电子同步辐射、SSA、nu_a、偏振和 seed photon 核。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `electron_radiation_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 19 | `first_greater_monotonic_from` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 40 | `first_greater_monotonic` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 59 | `first_greater_monotonic_window` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 80 | `besselk` | 同步辐射 Bessel kernel 插值/渐近 primitive。 |
| `S` | 140 | `get_syn_selected_state` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 197 | `simpson_emission_integral` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 212 | `simpson_ssa_tau_integral` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 226 | `accumulate_simpson_syn_point` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 245 | `build_reduced_log_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 274 | `project_syn_state_logbands` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 309 | `electron_syn_fx` | 同步辐射 kernel primitive。 |
| `F` | 321 | `electron_linear_interp` | 插值 primitive。 |
| `F` | 333 | `electron_syn_integrand_x` | 同步辐射 cell 积分 integrand。 |
| `F` | 344 | `electron_powerlaw_interp` | log-log/power-law 插值 primitive。 |
| `S` | 369 | `electron_log_gauss2_interval` | 积分权重或求积 primitive；影响谱积分精度。 |
| `F` | 387 | `electron_integrate_powerlaw_segment` | power-law cell 积分 primitive。 |
| `F` | 404 | `electron_ssa_segment` | SSA cell 光深积分 primitive。 |
| `F` | 436 | `electron_tau_kernel_x` | SSA optical-depth kernel primitive。 |
| `S` | 446 | `electron_syn_gauss_cell` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 470 | `electron_tau_gauss_cell` | 光深、SSA transfer 或 photon survival 相关算子。 |
| `S` | 502 | `electron_syn_cell_adaptive` | 低层单 cell adaptive diagnostic helper；public selected path 默认仍是 fixed-grid。 |
| `S` | 534 | `get_syn_polarization_selected` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 552 | `get_syn_polarization_fraction` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 587 | `get_syn_transfer` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 609 | `get_nu_a` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 681 | `evaluate_tau` | 光深、SSA transfer 或 photon survival 相关算子。 |
| `S` | 699 | `refine_nu_a_bracket` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 746 | `get_nu_a_2d_path` | 2D/chi SSA break diagnostic。 |
| `S` | 763 | `get_nu_a_2d_cell_path` | 2D/chi cell-level SSA break diagnostic。 |
| `S` | 779 | `reduce_syn_shell_from_chi` | chi-local spectra 到 shell-level baseline 的 reduction helper。 |
| `S` | 798 | `get_nu_a_from_tau_grid` | 从已计算 optical-depth grid 求 SSA break；避免重复 root search。 |
| `S` | 841 | `interpolate_log_tau_root` | 光深、SSA transfer 或 photon survival 相关算子。 |

### `src/Electron/electron_reverse_kernel.f90`

反向激波电子输运、post-crossing map、次级 RS 和 reacceleration 分支。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `electron_reverse_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 35 | `electron_reverse_evolve` | 反向激波电子演化入口；串接主 RS 注入、冷却、DG/有限体积输运和输出投影。 |
| `S` | 216 | `compute_reverse_cooling_loss` | 反向激波电子冷却率；把同步辐射、IC/Y 修正写入能量损失数组。 |
| `S` | 245 | `advance_reverse_transport_shell` | 主 RS 电子壳层输运；显式推进注入、冷却、绝热损失、DG 冷却节点插值和活动能段边界。 |
| `S` | 353 | `advance_reverse_post_cross_analytic` | crossing 后无新注入时的解析特征线重映射。 |
| `F` | 394 | `reverse_post_cross_map_value` | crossing 后特征线边界映射插值。 |
| `S` | 425 | `initialize_reverse_dg_state` | 初始化反向激波 DG 活动谱元网格和状态。 |
| `S` | 435 | `remesh_reverse_dg_state` | 按当前 cooling/injection break 重建 DG 网格并保守投影状态。 |
| `S` | 453 | `ensure_reverse_dg_work` | DG workspace 管理；只服务状态投影和子步推进暂存。 |
| `S` | 463 | `prepare_reverse_transport_substep_state` | 取半步半径处的 RS 动力学量并更新冷却/注入能标。 |
| `F` | 479 | `reverse_shell_linear_value` | 在相邻半径壳层之间线性插值反向激波动力学数组。 |
| `F` | 489 | `reverse_dg_upper_break` | DG 高能活动边界；crossing 后跟踪冷却过的高能 front。 |
| `F` | 499 | `reverse_dg_source_upper_xmax` | 由高能注入/冷却边界确定 DG 源项上界。 |
| `F` | 507 | `reverse_dg_active_xmax` | 结合尾部能量矩阈值确定 DG 活动网格上界。 |
| `F` | 521 | `reverse_dg_low_break` | DG 低能活动边界；处理近跨相对论 kinetic-source break。 |
| `S` | 531 | `advance_reverse_dg_injection_front` | 推进或重置 crossing 后的注入能标 front。 |
| `F` | 543 | `reverse_dg_injection_break` | DG 注入 break；crossing 后使用已冷却的 gamma_m front。 |
| `S` | 550 | `advance_reverse_dg_front_value` | 对单个 DG break/front 施加同步冷却和绝热漂移。 |
| `S` | 569 | `electron_secondary_reverse_evolve` | 次级反向激波电子演化入口；处理密度跳变分支和再加速分支。 |
| `S` | 651 | `advance_secondary_transport_shell` | 次级 RS 分支的电子输运推进。 |
| `S` | 692 | `electron_secondary_reverse_synchrotron` | 聚合次级 RS 分支同步辐射谱。 |
| `S` | 717 | `electron_secondary_reverse_branch_synchrotron` | 单个次级 RS 分支同步辐射输出。 |
| `S` | 749 | `electron_secondary_reverse_branch_reaccelerated` | 次级 RS 再加速电子分支输出。 |
| `S` | 809 | `build_secondary_reaccel_gamma_grid` | 构造再加速分支 gamma 网格。 |
| `S` | 834 | `transfer_reaccelerated_parent_electrons` | 把父分支电子谱转移到再加速分支网格。 |
| `S` | 873 | `advance_reaccelerated_branch_shell` | 再加速分支的输运推进。 |
| `S` | 917 | `prepare_branch_shell` | 读取单个次级分支在当前壳层的动力学和辐射参数。 |
| `S` | 945 | `compute_branch_injection` | 次级分支注入归一化和能标。 |
| `S` | 955 | `boost_log_distribution` | 对数谱的 Lorentz boost 重映射。 |
| `S` | 980 | `dsa_reaccelerate_distribution` | DSA 再加速谱构造。 |
| `F` | 997 | `distribution_energy_from_log` | 对数电子谱能量积分。 |
| `F` | 1005 | `reverse_transport_substeps` | 根据冷却步长和求解器类型选择反向输运子步数。 |
| `F` | 1017 | `reverse_dg_kinetic_break` | 反向 DG kinetic source 的低能 break。 |
| `S` | 1025 | `reverse_dg_grid_sequence` | 根据低/注入/高能 break 构造 DG 网格序列。 |
| `F` | 1080 | `reverse_interp_log_grid` | 在对数 gamma 网格上插值正谱/冷却量。 |

### `src/Electron/electron_seed_history_kernel.f90`

2D shell 历史 photon field、Doppler 映射和路径吸收。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_seed_history_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 27 | `ensure_history_cache` | workspace/cache 管理；只服务性能和内存复用，不改变物理语义。 |
| `S` | 72 | `integrate_downstream_proper_time` | 从径向网格和 bulk Lorentz factor 积分下游共动固有时间。 |
| `S` | 90 | `accumulate_comoving_history_fields` | 直接回扫历史壳层，把可由光行时连接的历史辐射叠加到当前共动 photon field。 |
| `S` | 135 | `accumulate_history_source_cell` | 计算源 chi cell 的有限光行时权重并进入目标 cell 累积。 |
| `S` | 147 | `accumulate_history_target_cell` | 沿历史源到当前目标的路径应用 Doppler 映射和壳层吸收。 |
| `S` | 188 | `advance_comoving_history_stream` | 用一阶特征线推进历史 photon stream，避免每个目标壳层回扫全部过去壳层。 |
| `S` | 233 | `accumulate_mapped_cell` | 把旧 stream 或上一壳层源场按权重、Doppler 映射和路径吸收累积到下一 stream。 |
| `S` | 262 | `build_doppler_map` | 构造目标频率到历史源频率网格的相对 Doppler 映射。 |
| `F` | 293 | `loglog_interp_mapped` | 按预计算对数分数做正谱 log-log 插值。 |
| `F` | 305 | `relative_doppler_backward` | 计算历史源区到当前目标区的相对 Doppler 因子 `D = gamma_rel*(1+beta_rel)`。 |
| `S` | 317 | `build_shell_transfer_cache` | 缓存 chi 单元宽度倒数和吸收 log-prefix，用于历史光线路径吸收。 |
| `S` | 342 | `apply_shell_path_attenuation` | 对一段壳层内光线路径累乘均匀源函数传输权重。 |
| `S` | 386 | `locate_path_cell` | 在下游面网格中定位光线路径端点所在 chi 单元。 |
| `F` | 437 | `history_transfer_weight` | 把 optical-depth segment 转换为均匀源函数逃逸/传输权重。 |

### `src/Electron/electron_shell_transport_common.f90`

shell-level fullhide/flux-split 共享 primitive。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_shell_transport_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `F` | 19 | `electron_resolve_1d_solver_id` | 1D 电子输运 solver id 解析。 |
| `S` | 28 | `electron_shell_flux_split_coord_step` | 四速度坐标上的单步 flux-split 推进；内部计算坐标 Jacobian。 |
| `S` | 46 | `electron_shell_dcoord_to_dndgamma_exp_centers` | 四速度坐标密度映射回 `dN/dgamma`。 |

### `src/Electron/electron_transport_2d_kernel.f90`

finite-q 几何、q 方向对流/扩散和 2D 能量推进。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_transport_2d_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 25 | `compute_q_geometry` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 41 | `q_geometry_point` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 67 | `compute_q_cell_geometry` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 90 | `compute_downstream_comoving_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 115 | `get_shock_transport_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 131 | `compute_q_divergence` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `F` | 153 | `q_face_transport_coeff` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 162 | `q_face_transport_coeffs` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `F` | 174 | `compute_q_step_limit` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `F` | 189 | `eta_linear_hit_time` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 203 | `eta_trace_back_faces_piecewise` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 268 | `q_split_advection_faces` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `F` | 287 | `q_depth_inverse_metric` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 309 | `q_diffusion_face_coeffs` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 334 | `build_q_advection_base_matrix` | q 方向对流矩阵 primitive。 |
| `S` | 355 | `add_q_diffusion_to_matrix` | q 方向扩散矩阵 primitive。 |
| `S` | 381 | `build_q_transport_base_matrix` | q 方向对流/扩散隐式矩阵组装。 |
| `S` | 399 | `solve_tridiagonal` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 423 | `advance_q_advection_charint` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 454 | `advance_q_diffusion_implicit` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 483 | `advance_q_implicit` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 519 | `advance_q_pwncr_implicit` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 554 | `advance_energy_loggamma_chi` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 580 | `advance_energy_loggamma_chi_pwncr` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 605 | `advance_energy_stochastic_loggamma_chi` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 639 | `advance_energy_loggamma_chi_charint` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |

### `src/Electron/electron_transport_common.f90`

1D 保守重映射、PPM、特征线和 fullhide 基础 primitive。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `electron_transport_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 21 | `electron_prepare_implicit_coeffs_common` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 33 | `electron_backward_sweep_common` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 47 | `electron_quadratic_interp3` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 60 | `electron_ppm_interfaces_nonuniform` | 非均匀电子能量网格的 PPM 界面重构。 |
| `S` | 119 | `electron_ppm_positive_cell` | PPM 单元正性限制；保持重构后的粒子数密度非负。 |
| `S` | 149 | `electron_ppm_prefix_nonuniform` | 非均匀网格 PPM 前缀积分。 |
| `F` | 168 | `electron_ppm_prefix_eval_nonuniform` | 任意位置的 PPM 保守前缀积分求值。 |
| `S` | 203 | `electron_prepare_conservative_remap_nonuniform` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 220 | `electron_prepare_exponential_source_remap` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 254 | `electron_exp_source_cell_int` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 274 | `electron_exp_source_prefix_eval` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 310 | `electron_dnx_to_dndgamma_exp_centers` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 324 | `electron_u_edges_from_x` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `F` | 338 | `electron_x_from_u` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 347 | `electron_trace_affine_u_edges_batch` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 378 | `electron_build_piecewise_affine_u` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 416 | `electron_find_u_cell_desc_hint` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 441 | `electron_trace_piecewise_affine_u_edge_from_cell` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 562 | `electron_trace_piecewise_affine_u_edges_batch` | 多滞后分段仿射 u 特征线回溯；单滞后调用传 `Num_lag=1`。 |
| `S` | 583 | `electron_characteristic_core` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 620 | `electron_characteristic_remap_edges` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 644 | `electron_characteristic_update` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 678 | `electron_ppm_cell_int` | PPM 单元积分 primitive。 |
| `S` | 695 | `electron_semi_lagrangian_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 722 | `electron_fullhide_flux_split_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 765 | `electron_fullhide_flux_split_step_nonuniform` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 814 | `electron_fullhide_flux_split_sequence_nonuniform` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 863 | `electron_fullhide_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 886 | `electron_fullhide_spacetime_sequence` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 930 | `electron_logparabola_peak_frequency` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 958 | `electron_active_gamma_hi` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 991 | `electron_active_chi_hi` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `F` | 1010 | `electron_max_xi_coeff_chi` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `F` | 1037 | `electron_max_xi_coeff_uniform` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_transport_dg_1d_kernel.f90`

DG 谱元网格、投影、正性核和特征线投影。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_transport_dg_1d_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 49 | `electron_dg1d_build_four_velocity_mesh` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 61 | `electron_dg1d_build_coord_mesh` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 94 | `electron_dg1d_initial_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 114 | `electron_dg1d_project_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 135 | `electron_dg1d_project_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 177 | `electron_dg1d_project_kinetic_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 216 | `electron_dg1d_scale_to_content` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 229 | `electron_dg1d_limit_positive_cell_preserving` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 256 | `electron_dg1d_apply_positive_kernel_filter` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 310 | `electron_dg1d_advance_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 354 | `electron_dg1d_advance_step_dense` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 376 | `electron_dg1d_assemble_transport_matrix` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 430 | `electron_dg1d_advance_characteristic_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 447 | `electron_dg1d_zero_negative_cell_means` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 462 | `electron_dg1d_project_characteristic` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 521 | `electron_dg1d_closed_low_boundary_content` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 544 | `electron_dg1d_integrate_domain_interval` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 561 | `electron_dg1d_characteristic_back_x` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 589 | `electron_dg1d_characteristic_forward_x` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 608 | `electron_dg1d_project_element` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 648 | `electron_dg1d_solve_lgl_block` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 696 | `electron_dg1d_solve_dense` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 738 | `electron_dg1d_project_to_coord_cells` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 774 | `electron_dg1d_integral` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 789 | `electron_dg1d_tail_moment_fraction` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 815 | `electron_dg1d_positive_kernel_mode` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 842 | `electron_dg1d_element_is_troubled` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 862 | `electron_dg1d_kernel_factor` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 872 | `electron_dg1d_jackson_factor` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 886 | `add_active_break` | 特征 Lorentz 因子/断点诊断；用于注入、冷却或活动网格边界。 |
| `S` | 899 | `electron_dg1d_sort_breaks` | 特征 Lorentz 因子/断点诊断；用于注入、冷却或活动网格边界。 |
| `S` | 915 | `allocate_spectral_mesh` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 929 | `ensure_reference_spectral` | workspace/cache 管理；只服务性能和内存复用，不改变物理语义。 |
| `S` | 946 | `ensure_projection_quadrature` | workspace/cache 管理；只服务性能和内存复用，不改变物理语义。 |
| `S` | 954 | `set_domain_bounds` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 982 | `fill_physical_nodes` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1000 | `lgl_nodes_weights` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1026 | `legendre_value_derivative` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1053 | `legendre_basis_values` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1067 | `gauss_legendre_nodes_weights` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1086 | `barycentric_weights` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1101 | `differentiation_matrix` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1116 | `locate_domain` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1134 | `interpolate_domain` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

### `src/Electron/hybrid_spectrum_kernel_fast.f90`

热-非热混合谱归一化和特殊函数加速。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 4 | `hybrid_spectrum_kernel_fast` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 48 | `integral_thermal1` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 145 | `integral_thermal12` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 305 | `integral_cpl` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 322 | `solve_theta` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 352 | `get_initial_theta` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 361 | `newton_method` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 396 | `normalized_hybrid_spec` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 424 | `hybrid_spec_point` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 446 | `normalized_hybrid_spec_lg` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Hadronic/hadronic_acceleration_kernel.f90`

强子加速时间、外部冷却限制和 gamma_max/injection operator。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_acceleration_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 22 | `hadronic_species_properties` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 55 | `hadronic_acceleration_timescale_s` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 77 | `hadronic_synchrotron_cooling_timescale_s` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 100 | `hadronic_external_photon_cooling_timescale_s` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 118 | `hadronic_species_injection_operator` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 165 | `hadronic_estimate_max_gamma` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 201 | `initialize_acceleration_limits` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 214 | `build_external_cooling_timescales` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 219 | `find_external_cooling_crossing` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 233 | `apply_external_cooling_limit` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 255 | `hadronic_acceleration_operator` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 284 | `hadronic_trapezoid` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Hadronic/hadronic_bethe_heitler_kernel.f90`

Bethe-Heitler 质子损失、pair source 和 photon loss kernel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_bethe_heitler_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 18 | `hadronic_bethe_heitler_operator` | Bethe-Heitler pair/source/loss 算子。 |
| `F` | 77 | `bh_proton_loss_point` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 83 | `accumulate_bh_pair_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 108 | `hadronic_bh_kernel_electron_generation` | Bethe-Heitler pair/source/loss 算子。 |
| `F` | 125 | `hadronic_bh_outer` | Bethe-Heitler pair/source/loss 算子。 |
| `F` | 145 | `hadronic_bh_inner` | Bethe-Heitler pair/source/loss 算子。 |
| `F` | 161 | `hadronic_bh_sigma_w` | Bethe-Heitler pair/source/loss 算子。 |
| `F` | 200 | `hadronic_bh_kernel_proton_loss` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 206 | `hadronic_bh_eloss_kernel_phi` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 227 | `hadronic_rk4_3` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `F` | 229 | `func` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 252 | `hadronic_rk4_4` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `F` | 254 | `func` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Hadronic/hadronic_common.f90`

强子网格、时间尺度、log-grid 校验和量子同步冷却公共工具。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 17 | `hadronic_build_gamma_p_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 42 | `hadronic_source_bounds` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 66 | `hadronic_build_gamma_edges` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `F` | 89 | `hadronic_shell_dt` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 103 | `hadronic_dynamical_time` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 113 | `hadronic_gamma_p_max` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 127 | `hadronic_validate_log_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `F` | 154 | `hadronic_quantum_syn_cooling_factor` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |

### `src/Hadronic/hadronic_decay_kernel.f90`

pi0、pi±、mu± decay、neutrino 和 electron channel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_decay_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 20 | `hadronic_pi0_to_gamma_operator` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 48 | `hadronic_pion_decay_operator` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 96 | `build_pion_decay_log_rates` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 106 | `accumulate_pion_muon_channel` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 130 | `accumulate_prompt_pion_neutrino_channel` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 151 | `hadronic_muon_decay_operator` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 199 | `build_muon_decay_log_rates` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 211 | `accumulate_muon_neutrino_channel` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 235 | `accumulate_muon_electron_channel` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 255 | `hadronic_hummer2010_decay_operator` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 313 | `hadronic_log_spacing` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 326 | `hadronic_fnu1_decay` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 336 | `hadronic_fnu2_decay` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 346 | `hadronic_log_interpolate` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

### `src/Hadronic/hadronic_forward_1d.f90`

正向激波 formal 1D 强子 f2py 入口和 shell-sequence driver。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 3 | `fs_hadronic_syn_polarization_shell` | f2py/Python 调用边界；同步偏振诊断使用的 proton synchrotron polarization wrapper。 |
| `S` | 16 | `fs_hadronic_pgamma_operator_shell` | f2py/Python 调用边界；Hummer 2010 p-gamma operator。 |
| `S` | 37 | `fs_hadronic_pair_production_shell` | f2py/Python 调用边界；gamma-gamma pair source/loss operator。 |
| `S` | 54 | `fs_hadronic_pp_delta_shell` | f2py/Python 调用边界；pp delta source/loss operator。 |
| `S` | 73 | `fs_hadronic_bethe_heitler_shell` | f2py/Python 调用边界；Bethe-Heitler pair source 和 photon/proton loss。 |
| `S` | 89 | `fs_hadronic_hadronic_ic_shell` | f2py/Python 调用边界；pion/muon/proton hadronic IC diagnostic operator。 |
| `S` | 114 | `fs_hadronic_decay_operator_shell` | f2py/Python 调用边界；pion/muon decay 到 gamma、e± 和 neutrino。 |
| `S` | 146 | `fs_hadronic_pair_cascade_sequence` | f2py/Python 调用边界；shell-sequence gamma-gamma pair/synch cascade。 |
| `S` | 172 | `fs_hadronic_1d` | light 1D proton transport + proton synchrotron shell-sequence driver。 |
| `S` | 247 | `initialize_proton_gamma_grid` | `fs_hadronic_1d` 内部网格初始化；Jacobians 不能省略。 |
| `S` | 275 | `initialize_output_grids` | `fs_hadronic_1d` 输出频率和数组初始化。 |
| `S` | 296 | `inject_protons_for_shell` | `fs_hadronic_1d` 壳层质子源项；必须同时满足粒子数和能量预算。 |
| `S` | 308 | `advance_proton_transport_for_shell` | `fs_hadronic_1d` 壳层质子输运推进。 |
| `S` | 319 | `advance_hummer_secondary_chain` | `fs_hadronic_1d` Hummer p-gamma secondary chain 推进。 |
| `S` | 337 | `emit_proton_synchrotron_for_shell` | `fs_hadronic_1d` proton synchrotron emissivity 和 seed 计算。 |
| `S` | 352 | `fs_hadronic_formal_transport_1d` | formal 1D 强子主入口；推进 proton transport、p-gamma/BH/pp、secondary、photon survival 和 secondary e± source。 |
| `S` | 394 | `fs_hadronic_positive_loglog_interp` | f2py/Python 调用边界；当前 secondary-feedback grid projection 仍直接调用。 |

### `src/Hadronic/hadronic_forward_formal_1d.f90`

formal 1D 强子壳层序列实现层；Python/f2py 只通过 `fs_hadronic_formal_transport_1d` 进入这里。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `hadronic_forward_formal_1d` | formal 1D 强子 shell-sequence 模块命名空间。 |
| `S` | 9 | `hadronic_forward_formal_transport_1d_impl` | 按半径推进 proton injection/transport、pγ/BH/pp、secondary species、secondary radiation、photon survival 和 secondary e± source。 |

### `src/Hadronic/hadronic_forward_shell_1d.f90`

formal 1D 强子底层 shell primitive 与单位/投影 helper；f2py wrapper 已收窄到运行时真正需要的入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `hadronic_forward_shell_1d` | shell-level 强子 primitive 模块命名空间。 |
| `S` | 24 | `hadronic_forward_pp_delta_shell` | pp delta source/loss operator；输出 gamma、neutrino、e± 源和 proton loss。 |
| `S` | 51 | `hadronic_forward_hadronic_ic_shell` | proton/pion/muon hadronic IC operator；输出 IC emissivity 与投影系数。 |
| `S` | 83 | `hadronic_forward_hic_projected` | hadronic IC shell emissivity 投影到 photon grid。 |
| `S` | 109 | `hadronic_forward_species_transport_step` | n、pi、mu secondary species 同壳层保守推进。 |
| `S` | 172 | `hadronic_forward_injection_content` | 壳层注入能量预算到 species source content 的归一化。 |
| `S` | 194 | `hadronic_forward_global_gamma_p_max` | 沿半径序列估计全局 proton 最大 Lorentz factor。 |
| `S` | 218 | `hadronic_forward_secondary_radiation_shell` | pion/muon synchrotron 与 IC shell emissivity。 |
| `S` | 247 | `hadronic_forward_secondary_radiation_projected` | secondary radiation 从 hadron grid 投影到 photon grid。 |
| `S` | 307 | `hadronic_forward_continuous_loss_rates` | adiabatic、synchrotron 和 quantum-synch 连续损失率。 |
| `S` | 329 | `hadronic_forward_secondary_electron_sequence` | secondary e± source 随壳层序列组装。 |
| `S` | 362 | `hadronic_forward_photon_loss_closure` | photon loss rate 到 optical-depth/survival closure。 |
| `S` | 387 | `hadronic_forward_interaction_effective_time` | interaction loss rate 的有效时间积分。 |
| `S` | 413 | `hadronic_forward_pgamma_proton_update` | pγ loss/re-injection 对 proton spectrum 的壳层更新。 |
| `S` | 433 | `hadronic_forward_proton_transport_step` | proton injection、continuous loss 和 pγ update 的单壳层推进。 |
| `S` | 454 | `hadronic_forward_exponential_sink` | 指数 sink primitive；用于已定义 interaction/loss closure。 |
| `S` | 470 | `hadronic_forward_energy_luminosity_from_rate` | rate-per-energy 到 luminosity-per-frequency/energy 的壳层换算。 |
| `S` | 483 | `hadronic_forward_project_luminosity_from_rate` | source energy grid 到目标 photon grid 的 luminosity 投影。 |
| `S` | 501 | `hadronic_forward_project_hic_luminosity` | hadronic IC 多 species emissivity 投影。 |
| `S` | 521 | `hadronic_forward_pair_source_content` | pp/BH pair source 组合成电子方程使用的 content source。 |
| `S` | 535 | `hadronic_forward_shell_density_per_gev` | shell content 到 density-per-GeV 的单位变换。 |
| `S` | 550 | `hadronic_forward_gamma_edges` | Lorentz-factor grid edges；守恒积分需要的 bin geometry。 |
| `S` | 572 | `hadronic_forward_process_power` | secondary process power diagnostic。 |
| `S` | 611 | `hadronic_forward_positive_loglog_interp` | 正值 log-log grid projection；当前 secondary-feedback Python glue 仍使用。 |
| `S` | 642 | `hadronic_forward_source_per_gamma` | source-per-energy 到 source-per-gamma 的 grid projection。 |
| `S` | 658 | `hadronic_forward_distribution_per_gev` | distribution-per-gamma 到 distribution-per-GeV 的 grid projection。 |
| `S` | 674 | `hadronic_forward_aligned_photon_grid` | hadron/photon grid 对齐 helper。 |
| `S` | 691 | `hadronic_sequence_shell_geometry` | shell dr、dt geometry helper。 |
| `S` | 712 | `hadronic_forward_shell_volumes` | 壳层体积序列。 |
| `S` | 734 | `hadronic_forward_dynamical_time` | shell dynamical time。 |
| `S` | 746 | `hadronic_forward_quantum_syn_cooling_factor` | quantum synchrotron cooling factor。 |

### `src/Hadronic/hadronic_hadronic_ic_kernel.f90`

强子 secondary inverse-Compton kernel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_hadronic_ic_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 21 | `hadronic_hadronic_ic_initialize_kernel` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 47 | `hadronic_hadronic_ic_operator_from_kernel` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 73 | `hadronic_hadronic_ic_apply_kernel` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 100 | `apply_hadronic_ic_species` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 110 | `sum_charged_pion_density` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `S` | 116 | `sum_muon_helicity_density` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `S` | 124 | `hadronic_hadronic_ic_build_species_kernel` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 149 | `hadronic_hadronic_ic_compute_channel` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `F` | 176 | `hadronic_hadronic_ic_coeff` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |

### `src/Hadronic/hadronic_interaction_kernel.f90`

Hummer 2010 p-gamma response、光子汇和 secondary family deposition。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_interaction_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 91 | `hadronic_pg_hummer2010_operator` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 165 | `accumulate_pg_photon_loss` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 195 | `hadronic_pg_deposit_family` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 214 | `hadronic_pg_deposit_baryons` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 231 | `hadronic_pg_deposit_shifted` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 247 | `hadronic_pg_family_rates_res` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 266 | `hadronic_pg_family_rates_dir` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 283 | `hadronic_pg_family_rates_mul` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 300 | `hadronic_pg_loss_coeff_res` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 310 | `hadronic_pg_loss_coeff_dir` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 320 | `hadronic_pg_loss_coeff_mul` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 332 | `hadronic_pg_family_photon_loss_res` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 340 | `hadronic_pg_family_photon_loss_dir` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 347 | `hadronic_pg_family_photon_loss_mul` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 355 | `hadronic_pg_kernel_res` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 367 | `hadronic_pg_kernel_mul` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 381 | `hadronic_pg_kernel_dir` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 397 | `hadronic_pg_idir` | p-gamma Hummer response 或 secondary family deposition。 |

### `src/Hadronic/hadronic_pair_cascade_kernel.f90`

gamma-gamma pair/synch shell-sequence cascade。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_pair_cascade_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 15 | `hadronic_cascade_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 48 | `validate_pair_cascade_inputs` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 56 | `run_pair_production_stage` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 66 | `build_pair_electron_gamma_grid` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 71 | `evolve_pair_cooling_stage` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 87 | `emit_pair_synchrotron_stage` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 107 | `hadronic_cascade_sequence` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 174 | `validate_sequence_inputs` | 系统边界校验；用于拒绝外部输入或正式 kernel contract 违反。 |
| `S` | 185 | `build_sequence_grids` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 197 | `shell_geometry` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 212 | `dynamical_time` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 217 | `pair_synchrotron_state` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 234 | `distribute_cooled_power` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 252 | `electron_loss_rates` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 262 | `advance_energy_loggamma` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `F` | 280 | `gamma_bin_index` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 294 | `advance_photon_density` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |

### `src/Hadronic/hadronic_pair_production_kernel.f90`

gamma-gamma pair injection 和 photon loss kernel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_pair_production_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 18 | `hadronic_pair_production_operator` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 65 | `hadronic_calc_pair_injection_full` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 84 | `accumulate_pair_injection_bin` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 100 | `pair_injection_energy_window` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 109 | `renormalize_pair_energy_closure` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 124 | `hadronic_build_photon_loss_kernel` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 143 | `hadronic_phibar` | 结构化喷流 patch/角向投影调度。 |
| `F` | 160 | `hadronic_phibar1` | 结构化喷流 patch/角向投影调度。 |
| `F` | 169 | `hadronic_phibar2` | 结构化喷流 patch/角向投影调度。 |
| `F` | 181 | `hadronic_phibar3` | 结构化喷流 patch/角向投影调度。 |
| `F` | 189 | `hadronic_phisum` | 结构化喷流 patch/角向投影调度。 |
| `F` | 200 | `hadronic_inner_pp` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 251 | `hadronic_outer_pp` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 258 | `hadronic_x_l` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 264 | `hadronic_gm_b` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 270 | `hadronic_fpp_m` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 276 | `hadronic_fpp_p` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 282 | `hadronic_rgg_d1` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 346 | `hadronic_xcm_u` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 358 | `hadronic_xcm_l` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 390 | `hadronic_td2` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 411 | `hadronic_a_a0_h` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 421 | `hadronic_a0` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 431 | `hadronic_a0f` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 443 | `hadronic_a_a0_hf` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 453 | `hadronic_beta` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 463 | `hadronic_sign_int` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 482 | `hadronic_grid_index_offset` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

### `src/Hadronic/hadronic_pgamma_hummer_1d.f90`

旧/诊断型 Hummer p-gamma shell 聚合 helper。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_pgamma_hummer_1d` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 12 | `hadronic_pgamma_hummer_shell` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 105 | `map_pgamma_secondary_sources` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 136 | `apply_neutron_pgamma_loss` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 152 | `interp_source_per_gamma` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 163 | `interp_positive_loglog` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

### `src/Hadronic/hadronic_pp_kernel.f90`

Delta 近似 pp 算子和 secondary source。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_pp_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 19 | `hadronic_pp_delta_operator` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 48 | `validate_pp_delta_inputs` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 59 | `set_pp_delta_options` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 84 | `build_pp_parent_collision_rate` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 90 | `emit_pp_delta_secondaries` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 101 | `hadronic_pp_sigma_inelastic_kelner2006` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 127 | `hadronic_pp_threshold_kinetic_energy_gev` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 133 | `hadronic_pp_delta_secondary_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 149 | `hadronic_pp_loglog_interp_positive` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 183 | `hadronic_pp_upper_bracket` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 211 | `hadronic_pp_loglog_linear_eval` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |

### `src/Hadronic/hadronic_pp_models_kernel.f90`

Kelner/Kafexhiu 类 pp spectral model helper。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 3 | `hadronic_pp_models_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 21 | `hadronic_pp_threshold_kinetic_energy_gev` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 28 | `hadronic_pp_spectral_shape` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 63 | `high_energy_or_geant4_shape` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 76 | `F_geant4_local` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 111 | `Egam_max` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 119 | `hadronic_pp_pi0_source_spectrum` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 150 | `sigma_pi0_model` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 172 | `Amax_model` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 196 | `sigma_1_pi_local` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 206 | `sigma_2_pi_local` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 215 | `sigma_inel_local` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 222 | `multip_pi0_SIBYLL` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 231 | `multip_pi0_QGSJET` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 240 | `multip_pi0_Geant4` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 247 | `multip_pi0_Pythia8` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |

### `src/Hadronic/hadronic_radiation_kernel.f90`

质子同步辐射和偏振。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_radiation_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 16 | `hadronic_syn_x_arg_mass` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 31 | `hadronic_syn_kernel_ultrarel_mass` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 66 | `hadronic_syn_kernel_ultrarel` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 77 | `hadronic_get_proton_syn_state` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 99 | `build_proton_syn_quadrature` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 109 | `accumulate_proton_syn_power` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 125 | `normalize_proton_syn_seed_density` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `S` | 134 | `hadronic_syn_polarization_fraction` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 162 | `accumulate_polarized_synch_integral` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |

### `src/Hadronic/hadronic_reverse_1d.f90`

反向激波轻量质子输运和质子同步辐射入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `fs_hadronic_reverse_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 35 | `initialize_reverse_proton_grid` | 反向激波局域状态；必须保持 crossing 前后和 sigma->0 极限连续。 |
| `S` | 49 | `advance_reverse_hadronic_shell` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 66 | `emit_reverse_proton_synchrotron` | 反向激波局域状态；必须保持 crossing 前后和 sigma->0 极限连续。 |

### `src/Hadronic/hadronic_secondary_radiation_kernel.f90`

pion/muon synchrotron 和 IC secondary radiation。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_secondary_radiation_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 28 | `hadronic_secondary_radiation_operator` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 67 | `hadronic_secondary_apply_radiation_kernels` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 102 | `hadronic_secondary_initialize_synchrotron_kernel` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 132 | `hadronic_secondary_initialize_inverse_compton_kernel` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 160 | `hadronic_secondary_pion_synchrotron_rate` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 175 | `hadronic_secondary_muon_synchrotron_rate` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 196 | `hadronic_secondary_pion_inverse_compton_rate` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 216 | `hadronic_secondary_muon_inverse_compton_rate` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 241 | `hadronic_secondary_syn_kernel_ultrarel` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 274 | `hadronic_secondary_build_ic_species_kernel` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 296 | `hadronic_secondary_compute_ic_channel` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `F` | 325 | `hadronic_secondary_ic_coeff` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 334 | `hadronic_secondary_validate_positive_log_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 346 | `hadronic_secondary_validate_density` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |

### `src/Hadronic/hadronic_species_transport_kernel.f90`

neutron、pion、muon explicit species transport。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_species_transport_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 24 | `hadronic_species_spherical_divergence_rate` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 37 | `hadronic_species_synchrotron_dgamma_dt` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 66 | `hadronic_species_adiabatic_dgamma_dt` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 83 | `hadronic_species_advance_operator` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 131 | `validate_species_transport_inputs` | 系统边界校验；用于拒绝外部输入或正式 kernel contract 违反。 |
| `S` | 156 | `build_species_cooling_rates` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 163 | `advance_charged_pion_pair` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 170 | `advance_muon_helicity_species` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 187 | `hadronic_species_advance_one` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 255 | `hadronic_species_validate_gamma_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 276 | `hadronic_species_validate_non_negative` | 系统边界校验；用于拒绝外部输入或正式 kernel contract 违反。 |
| `F` | 290 | `hadronic_species_log_spacing` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Hadronic/hadronic_transport_kernel.f90`

质子注入、连续损失和 log-gamma 有限体积推进。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_transport_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 10 | `hadronic_proton_injection_powerlaw` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 29 | `accumulate_powerlaw_energy_moment` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 42 | `write_powerlaw_injection` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 53 | `hadronic_proton_loss_rates` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 72 | `hadronic_advance_energy_loggamma` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 88 | `build_loss_flux_edges` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 100 | `apply_flux_divergence_with_injection` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |

### `src/Hadronic/hadronic_transport_remap_kernel.f90`

强子能量冷却重映射 primitive。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_transport_remap_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 9 | `hadronic_advance_energy_loggamma_remap` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 24 | `deposit_cooled_bin_content` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 37 | `restore_density_units` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `F` | 47 | `hadronic_remap_target` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

### `src/Interpolation/SED_interpolation.f90`

shell-level 和 chi-resolved EATS/Doppler 投影。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 9 | `sed_interpolation` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 108 | `sed_interpolation_adaptive_theta` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 168 | `integrate_theta_cell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 198 | `project_theta_sample` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 230 | `project_shell_segment` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 252 | `chi_ssa_cell_escape` | finite-q shell SSA escape factor。 |
| `F` | 266 | `lower_bound_real8` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 285 | `sed_interpolation_chi` | top-hat chi-resolved FS synchrotron+SSA lightcurve projection。 |
| `S` | 393 | `project_chi_segment_flux` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 407 | `compute_chi_segment_state` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 422 | `accumulate_chi_cell_source` | chi cell source + SSA escape 累加。 |
| `S` | 446 | `sed_interpolation_chi_electron_cached` | direct-electron chi projection；先批量计算 chi-local synchrotron/SSA 再投影。 |
| `S` | 481 | `sed_interpolation_chi_structured_axisym_ring_precomputed` | axisymmetric structured chi ring projection；输入为预计算 `F_ring/Tau_ring`。 |
| `S` | 567 | `project_precomputed_ring_segment` | structured ring EATS/Doppler segment projection。 |
| `S` | 583 | `accumulate_precomputed_ring_source` | structured ring chi cell source + SSA escape 累加。 |

### `src/Interpolation/SED_interpolation_structured.f90`

结构化喷流 theta/theta-phi shell-level 投影。该文件仍由 `structured_jet_1d` Fortran 聚合入口调用；旧 Python binding 不再从 `src.Interpolation` 导出。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 9 | `sed_interpolation_structured` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 103 | `project_structured_segment` | 结构化喷流 patch/角向投影调度。 |
| `S` | 125 | `sed_interpolation_structured_phi` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 202 | `project_structured_phi_segment` | 结构化喷流 patch/角向投影调度。 |

### `src/Interpolation/interpolation_common.f90`

SED 插值共享累加 primitive。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `interpolation_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 9 | `interpolation_accumulate_log_sed` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 46 | `interpolation_accumulate_shifted_linear_sed` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

### `src/Radiation/quantum_synchrotron_kernel.f90`

量子同步辐射修正因子。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `quantum_synchrotron_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 18 | `quantum_chi_parameter` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `F` | 30 | `quantum_syn_cooling_factor` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |

### `src/Radiation/radiation_common.f90`

辐射共享 Simpson 权重、transfer、pair cross-section 和 seed 转换。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `radiation_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 9 | `compute_simpson_weights` | 积分权重或求积 primitive；影响谱积分精度。 |
| `F` | 27 | `radiation_powerlaw_interp` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 52 | `radiation_transfer_factor` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 63 | `radiation_syn_kernel_value` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 78 | `radiation_prepare_annihilation_grid` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 94 | `radiation_pair_cross_section` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 112 | `radiation_pair_tau_headon_segment` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 142 | `radiation_syn_seed_chi_batch_core` | chi-resolved synchrotron/SSA batch kernel；单列调用传 `Num_chi=1`。 |
| `S` | 258 | `radiation_syn_flux_tau_chi_batch_core` | chi-resolved projection-only synchrotron flux/tau batch kernel。 |

### `src/Radiation/radiation_gamma_gamma_absorption.f90`

观测侧 gamma-gamma 吸收光深。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 3 | `annihilation` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 79 | `compute_mid_seed_field` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 92 | `prepare_pair_angle_weights` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 105 | `build_pair_sigma_kernel` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 123 | `set_pair_energy_window` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 145 | `clear_pair_energy_window` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `F` | 153 | `first_pair_index_gt` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `F` | 172 | `last_pair_index_lt` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |

### `src/Radiation/radiation_ssc_spectrum.f90`

SSC 谱、KN/Jones kernel、均匀/非均匀电子谱积分。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `ssc_spec` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `F` | 83 | `first_gamma_gt` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 102 | `first_gamma_ge` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 121 | `prepare_uniform_weights_and_tails` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 141 | `prepare_kn_tables` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 158 | `prepare_uniform_gamma_bounds` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 176 | `low_seed_sum` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 199 | `high_seed_tail_sum` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 220 | `accumulate_uniform_point` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 237 | `ssc_spec_nonuniform` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 289 | `validate_nonuniform_inputs` | 系统边界校验；用于拒绝外部输入或正式 kernel contract 违反。 |
| `S` | 304 | `build_nonuniform_reconstruction_and_tails` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 318 | `build_shell_reconstruction` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 336 | `set_limited_cell_slope` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 364 | `add_shell_tail_cell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 379 | `gauss2_abscissa` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 391 | `accumulate_nonuniform_point` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 432 | `first_cell_above_edge` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `F` | 461 | `ssc_minmod` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 473 | `linear_profile_value` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `F` | 483 | `linear_gamma_moment` | 特征 Lorentz 因子/断点诊断；用于注入、冷却或活动网格边界。 |
| `F` | 502 | `ssc_low_gamma_cell` | 特征 Lorentz 因子/断点诊断；用于注入、冷却或活动网格边界。 |
| `F` | 540 | `ssc_low_gamma_integral` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 565 | `ssc_high_gamma_tail` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |

### `src/Radiation/synchrotron_polarization_kernel.f90`

同步辐射偏振 F/G kernel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `synchrotron_polarization_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 12 | `synchrotron_fg_kernel` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 37 | `synchrotron_polarized_components` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 49 | `synchrotron_fg_integral` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 78 | `add_fg_node` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Structured/structured_jet_1d.f90`

结构化喷流 patch 调度、缓存和 Fortran 聚合入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 1 | `structured_jet_flux_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 57 | `run_axisymmetric` | 结构化喷流 patch/角向投影调度。 |
| `S` | 93 | `run_nonaxisymmetric` | 结构化喷流 patch/角向投影调度。 |
| `S` | 131 | `structured_solve_axisymmetric` | 结构化喷流 patch/角向投影调度。 |
| `S` | 187 | `register_axis_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 204 | `solve_axis_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 222 | `copy_axis_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 232 | `structured_solve_nonaxisymmetric` | 结构化喷流 patch/角向投影调度。 |
| `F` | 300 | `flatten_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 307 | `unflatten_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 316 | `register_phi_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 335 | `solve_phi_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 353 | `copy_phi_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 367 | `structured_solve_element` | 结构化喷流 patch/角向投影调度。 |
| `S` | 523 | `structured_apply_observer_absorption` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |

## 修改时的最短自检

- 只改网页文档：`tools/check_text_encoding.py`、`git diff --check`、`mkdocs build --strict --site-dir /tmp/asgard_mkdocs_site`。
- 改 dynamics：重建 `Dynamics_forward` / `Dynamics_reverse`，并跑 RS smoke 或受影响 lightcurve smoke。
- 改 electron：重建对应 electron module；显式从 `/tmp` 跑 `-Wline-truncation` 源闭包；看电子谱、`nu_a`、seed 和光变是否平滑。
- 改 hadronic：重建 `hadronic_forward_1d` 或 `hadronic_reverse_1d`；检查 formal grid contract、proton/secondary/photon survival 预算。
- 改 interpolation：重建 `SED_interpolation` 或 `SED_interpolation_structured`；比较 lightcurve/SED/chi projection 是否只改变预期投影层。
