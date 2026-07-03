# Fortran 程序单元索引

本文是当前工作树的 Fortran kernel 索引。它面向需要逐个进入子程序读算法的人：先看 f2py 入口和物理阶段，再进入文件内的 `module`、`subroutine`、`function`。更高层的物理到算法映射见 `doc/physics_algorithm_crosswalk.md`，运行主链见 `doc/call_chain.md`。

当前索引按 ASGARD 自有 Fortran 数值核抽取，排除第三方固定格式特殊函数依赖，共 723 个程序单元：45 个 module、497 个 subroutine、181 个 function。行号是生成本页时的源文件位置。

## 读源码顺序

1. 从 public Python API 进入 `asgard_core/asgard_state.py`，确认阶段顺序。
2. 在本页查 f2py module 和主 entry。
3. 进入对应 Fortran 文件，先读入口参数和数组形状，再读内部 contained procedures。
4. 对涉及真实物理演化的改动，回到 `doc/physics_algorithm_crosswalk.md` 查验收指纹。

不要把内部 helper 当 public ABI。稳定边界是 Python public API、`build_extensions.py` 登记的 f2py module、以及本页列出的主 entry。

## f2py 模块和主入口

| Build module | CWD | Source closure | 主 entry | 物理/算法角色 |
| --- | --- | --- | --- | --- |
| `Dynamics_forward` | `src/Dynamics` | `Constants + dynamics_density_profile + Dynamics_forward` | `dynamics_forward` | 正向激波动力学、ISM/wind、密度跳变和能量注入。 |
| `Dynamics_reverse` | `src/Dynamics` | `Constants + dynamics_density_profile + reverse_shock_state + reverse_shock_mhd_jump + reverse_jump_conditions + reverse_rhs + reverse_shock` | `dynamics_reverse` | first RS crossing/region-3 thermal state 和 density-jump multiple RS 分支。 |
| `electron_forward_fullhide_1d` | `src/Electron` | `ELECTRON_COMMON_SOURCES + electron_forward_fullhide_1d` | `fs_fullhide_1d; fs_fullhide_coupled` | 默认 1D 电子输运和 joint feedback coupled pass。 |
| `electron_forward_fullhide_1d_hybrid` | `src/Electron` | `ELECTRON_HISTORY_SOURCES_HZ + electron_forward_fullhide_1d_hybrid` | `fs_fullhide_hz` | 热/非热混合谱路径。 |
| `electron_forward_charint_1d` | `src/Electron` | `ELECTRON_COMMON_SOURCES + electron_forward_charint_1d` | `fs_charint_1d` | 1D 特征线输运对照。 |
| `electron_forward_dg_1d` | `src/Electron` | `ELECTRON_DG_1D_SOURCES + electron_forward_dg_1d` | `fs_dg_1d` | P12 LGL-DG 正向电子输运。 |
| `electron_forward_charint_2d` | `src/Electron` | `ELECTRON_2D_SOURCES + electron_forward_transport_2d` | `fs_transport_2d` | finite-q shell 2D 电子输运；别把 chi_grid 当强子局域坐标。 |
| `electron_forward_t2g1_1d` | `src/Electron` | `ELECTRON_COMMON_SOURCES + electron_forward_t2g1_1d` | `fs_t2g1_1d` | 方法比较电子输运。 |
| `electron_forward_weno5_1d` | `src/Electron` | `ELECTRON_COMMON_SOURCES + electron_forward_weno5_1d` | `fs_weno5_1d` | WENO5 方法比较电子输运。 |
| `electron_reverse_kernel` | `src/Electron` | `ELECTRON_REVERSE_SOURCES + electron_reverse_kernel` | `electron_reverse_evolve; multiple_evolve` | RS primary/secondary electron transport 和 reacceleration。 |
| `electron_radiation` | `src/Electron` | `ELECTRON_RADIATION_SOURCES` | `get_syn_*; nua_solve_*` | 电子同步/SSA/seed 低层核；通常经其他 entry 调用。 |
| `ssc_spectrum` | `src/Radiation` | `rad_common + ssc_spectrum` | `ssc_spec; ssc_spec_nonuniform` | SSC spectrum 和 KN/Jones 积分。 |
| `pair_absorption` | `src/Radiation` | `rad_common + pair_absorption` | `annihilation` | 观测侧 gamma-gamma absorption。 |
| `SED_interpolation` | `src/Interpolation` | `rad_common + interpolation_common + SED_interpolation` | `sed_interpolation; sed_adaptive_theta; sed_interpolation_chi; sed_chi_electron; sed_chi_ring` | 观测者投影；Python 公开绑定集中在此模块。 |
| `SED_interpolation_structured` | `src/Interpolation` | `interpolation_common + SED_interpolation_structured` | `sed_interpolation_structured; sed_structured_phi` | `structured_jet_1d` 内部 shell-level structured projection；不再从 `src.Interpolation` 暴露旧 Python 绑定。 |
| `hadronic_forward_1d` | `src/Hadronic` | `HADRONIC_1D_SOURCES + hadronic_forward_1d` | `formal_transport_1d; shell operators` | formal 1D proton/secondary/photon-loss shell sequence。 |
| `hadronic_reverse_1d` | `src/Hadronic` | `HADRONIC_1D_SOURCES + hadronic_reverse_1d` | `reverse_hadronic_1d` | RS light proton transport + proton synchrotron。 |
| `structured_jet_1d` | `src/Structured` | `STRUCTURED_JET_1D_SOURCES + structured_jet_1d` | `jet_flux_1d` | 结构化喷流 Fortran 聚合入口。 |

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
| `S` | 42 | `forward_rhs` | ODE/PDE 右端项；定义物理源汇和动力学变量导数。 |

### `src/Dynamics/reverse_shock.f90`

反向激波 `Dynamics_reverse.dynamics_reverse` f2py 入口：first RS 推进、等待成 shock、pre-crossing swept-mass 推进、crossing event split、post-crossing log-time RK、density-jump 分支扫描和 branch 输出都在同一个输出循环中完成。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 3 | `dynamics_reverse` | `Dynamics_reverse` f2py/Python 调用边界；解包 public `Boundary`，推进 first 分支并同步扫描 density-jump multiple 分支。 |

### `src/Dynamics/reverse_rhs.f90`

反向激波 RHS，显式处理 shock 注入、磁化惯性、区域 3 热演化和次级分支导数。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `reverse_dynamics_rhs` | ODE/PDE 右端项；定义反向激波主状态和次级分支导数。 |
| `F` | 201 | `parent_ready` | 判断上一级次级分支是否可作为当前上游。 |
| `S` | 211 | `branch_derivs` | 次级 RS 分支质量、热能、体积和注入诊断导数。 |
| `S` | 291 | `branch_density` | density-jump 分支的局部密度和 source 权重。 |

### `src/Dynamics/dynamics_density_profile.f90`

外介质密度、density jump/profile 状态和 tabulated profile 插值。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `dynamics_density_profile` | 介质密度和 density-jump/profile 状态；被 dynamics、electron、structured 路径按需引用。 |
| `S` | 20 | `density_profile` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `S` | 76 | `set_density_profile` | 从 `Boundary` 解包 density jump/profile 状态。 |
| `S` | 128 | `tab_density` | tabulated profile 的 log-log 插值。 |

### `src/Dynamics/reverse_shock_mhd_jump.f90`

有限强度反向激波 MHD jump 公式。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `reverse_shock_mhd_jump` | MHD jump 公式模块；不包含时间推进状态。 |
| `F` | 10 | `rs_vegas_ud` | 有限强度 MHD jump 解析根；不要用 ultra-relativistic 近似替代。 |
| `F` | 44 | `rs_vegas_comp` | 有限强度 MHD jump 压缩比；sigma->0 极限回到 hydrodynamic baseline。 |
| `F` | 66 | `rs_mag_internal` | MHD jump 下游热比内能；保持 crossing 前后和 sigma->0 极限连续。 |

正向 RK4 是 `Dynamics_forward.f90` 的本地 helper；反向 event-split RK 是 `reverse_shock.f90` 的本地 helper。

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
| `S` | 2 | `fs_charint_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 235 | `prepare_characteristic_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 275 | `write_final_characteristic_diagnostics` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_forward_dg_1d.f90`

1D LGL-DG 电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `fs_dg_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
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
| `S` | 5 | `fs_fullhide_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 123 | `initialize_forward_four_velocity_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 138 | `prepare_fullhide_shell` | 单个 FS shell 的密度、磁场、注入能标、同步谱、SSA break 和冷却率准备。 |
| `S` | 372 | `fs_fullhide_coupled` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 461 | `prepare_coupled_shell` | joint electron-photon shell 的辐射场、SSA break 和耦合冷却率准备。 |

### `src/Electron/electron_forward_fullhide_1d_hybrid.f90`

混合热/非热谱形入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 7 | `fs_fullhide_hz` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 266 | `build_hybrid_or_powerlaw_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 281 | `build_hybrid_or_powerlaw_source_from_count` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |

### `src/Electron/electron_forward_slc1_1d.f90`

半拉格朗日一阶电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `fs_slc1_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |

### `src/Electron/electron_forward_t2g1_1d.f90`

T2G1 方法比较电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 8 | `fs_t2g1_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 93 | `prepare_t2g1_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 118 | `write_t2g1_radiation_and_cooling` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 142 | `advance_t2g1_substep` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |

### `src/Electron/electron_forward_transport_2d.f90`

2D finite-q shell 电子输运和投影预处理入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `fs_transport_2d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 623 | `reduce_syn_shell_from_q` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 639 | `project_q_projection_shell` | finite-q shell 几何或 chi-equivalent 投影字段。 |

### `src/Electron/electron_forward_weno5_1d.f90`

WENO5 方法比较电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 2 | `fs_weno5_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 67 | `prepare_weno_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 89 | `write_weno_radiation_and_cooling` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 117 | `advance_weno_substep` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 145 | `load_weno_extended_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 158 | `compute_weno_fluxes` | 刷新零阶外推鬼点并按速度符号构造 WENO5 数值通量。 |
| `S` | 176 | `advance_weno_rk_stage` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `F` | 201 | `weno5_positive_flux` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 237 | `weno5_negative_flux` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_injection_profiles.f90`

非热/热电子注入谱、log-cell 积分和源项归一化。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_injection_profiles` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 13 | `exp_cutoff` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 22 | `pl_params` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 47 | `dnx_cutoff` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 64 | `electron_dnx_gauss3_integral` | 积分权重或求积 primitive；影响谱积分精度。 |
| `F` | 88 | `electron_dny_gauss3_integral` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 114 | `electron_dnx_segment` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 135 | `electron_dny_segment` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 156 | `log_edges` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `S` | 171 | `init_powerlaw` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 190 | `init_edges` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 228 | `init_coord` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 267 | `source_edges` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 292 | `source_coord` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 317 | `kinetic_edges` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 355 | `kinetic_coord` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 391 | `electron_thermal_theta` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 399 | `electron_build_thermal_shape_dnx` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 424 | `add_thermal` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 437 | `thermal_pop` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_radiation_kernel.f90`

电子同步辐射、SSA、nu_a、偏振和 seed photon 核。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `electron_radiation_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 19 | `greater_from` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 40 | `first_greater` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 59 | `greater_window` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 80 | `besselk` | 同步辐射 Bessel kernel 插值/渐近 primitive。 |
| `S` | 140 | `syn_state` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 197 | `simpson_emission_integral` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 212 | `simpson_ssa_tau_integral` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 226 | `accumulate_simpson_syn_point` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 245 | `reduce_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 274 | `project_syn` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 384 | `electron_syn_fx` | 同步辐射 kernel primitive。 |
| `F` | 396 | `electron_linear_interp` | 插值 primitive。 |
| `F` | 408 | `electron_syn_integrand_x` | 同步辐射 cell 积分 integrand。 |
| `F` | 419 | `pl_interp` | log-log/power-law 插值 primitive。 |
| `S` | 444 | `log_gauss2` | 积分权重或求积 primitive；影响谱积分精度。 |
| `F` | 462 | `pl_integral` | power-law cell 积分 primitive。 |
| `F` | 479 | `ssa_segment` | SSA cell 光深积分 primitive。 |
| `F` | 511 | `electron_tau_kernel_x` | SSA optical-depth kernel primitive。 |
| `S` | 521 | `electron_syn_gauss` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 545 | `electron_tau_gauss` | 光深、SSA transfer 或 photon survival 相关算子。 |
| `S` | 577 | `electron_syn_adapt` | 低层单 cell adaptive diagnostic helper；public selected path 默认仍是 fixed-grid。 |
| `S` | 609 | `syn_polarized` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 627 | `get_syn_polarization_fraction` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 662 | `syn_transfer` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 687 | `nua_solve` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 759 | `evaluate_tau` | 光深、SSA transfer 或 photon survival 相关算子。 |
| `S` | 777 | `refine_nu_a_bracket` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 824 | `nua_path` | 2D/chi cell-level SSA break diagnostic。 |
| `S` | 840 | `reduce_chi` | chi-local spectra 到 shell-level baseline 的 reduction helper。 |
| `S` | 859 | `nua_fromtau` | 从已计算 optical-depth grid 求 SSA break；避免重复 root search。 |
| `S` | 902 | `interpolate_log_tau_root` | 光深、SSA transfer 或 photon survival 相关算子。 |

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
| `F` | 499 | `reverse_dg_active_xmax` | 结合尾部能量矩阈值确定 DG 活动网格上界。 |
| `F` | 514 | `reverse_dg_low_break` | DG 低能活动边界；处理近跨相对论 kinetic-source break。 |
| `S` | 524 | `advance_reverse_dg_injection_front` | 推进或重置 crossing 后的注入能标 front。 |
| `F` | 536 | `reverse_dg_injection_break` | DG 注入 break；crossing 后使用已冷却的 gamma_m front。 |
| `S` | 543 | `advance_reverse_dg_front_value` | 对单个 DG break/front 施加同步冷却和绝热漂移。 |
| `S` | 562 | `multiple_evolve` | 次级反向激波电子演化入口；处理密度跳变分支和再加速分支。 |
| `S` | 644 | `advance_multiple` | 次级 RS 分支的电子输运推进。 |
| `S` | 685 | `multiple_synch` | 聚合次级 RS 分支同步辐射谱。 |
| `S` | 710 | `branch_synch` | 单个次级 RS 分支同步辐射输出。 |
| `S` | 742 | `branch_reaccel` | 次级 RS 再加速电子分支输出。 |
| `S` | 802 | `reaccel_grid` | 构造再加速分支 gamma 网格。 |
| `S` | 827 | `transfer_parent` | 把父分支电子谱转移到再加速分支网格。 |
| `S` | 866 | `advance_reaccel` | 再加速分支的输运推进。 |
| `S` | 910 | `prepare_branch_shell` | 读取单个次级分支在当前壳层的动力学和辐射参数。 |
| `S` | 938 | `branch_inject` | 次级分支注入归一化和能标。 |
| `S` | 948 | `boost_log` | 对数谱的 Lorentz boost 重映射。 |
| `S` | 973 | `dsa_reaccel` | DSA 再加速谱构造。 |
| `F` | 990 | `log_energy` | 对数电子谱能量积分。 |
| `F` | 998 | `reverse_transport_substeps` | 根据冷却步长和求解器类型选择反向输运子步数。 |
| `S` | 1010 | `dg_sequence` | 根据低/注入/高能 break 构造 DG 网格序列。 |
| `F` | 1067 | `log_interp` | 在对数 gamma 网格上插值正谱/冷却量。 |

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
| `F` | 154 | `compute_q_step_limit` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `F` | 170 | `eta_linear_hit_time` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 184 | `eta_trace_back_faces_piecewise` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 249 | `q_split_advection_faces` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `F` | 268 | `q_depth_inverse_metric` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 290 | `q_diffusion_face_coeffs` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 315 | `build_q_advection_base_matrix` | q 方向对流矩阵 primitive。 |
| `S` | 336 | `add_q_diffusion_to_matrix` | q 方向扩散矩阵 primitive。 |
| `S` | 362 | `build_q_transport_base_matrix` | q 方向对流/扩散隐式矩阵组装。 |
| `S` | 384 | `solve_tridiagonal` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 408 | `advance_q_advection_charint` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 442 | `advance_q_diffusion_implicit` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 471 | `advance_q_implicit` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 507 | `advance_q_pwncr_implicit` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 542 | `advance_energy_chi` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 568 | `advance_energy_chi_pwncr` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 593 | `advance_energy_stochastic_loggamma_chi` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 627 | `advance_energy_chi_charint` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |

### `src/Electron/electron_transport_common.f90`

1D 保守重映射、PPM、特征线和 fullhide 基础 primitive。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `electron_transport_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 21 | `prepare_implicit_coeffs` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 33 | `backward_sweep` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 47 | `quad_interp3` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 60 | `ppm_interfaces` | 非均匀电子能量网格的 PPM 界面重构。 |
| `S` | 119 | `ppm_positive_cell` | PPM 单元正性限制；保持重构后的粒子数密度非负。 |
| `S` | 149 | `ppm_prefix` | 非均匀网格 PPM 前缀积分。 |
| `F` | 168 | `ppm_eval_prefix` | 任意位置的 PPM 保守前缀积分求值。 |
| `S` | 203 | `prepare_remap` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 220 | `prepare_exp_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 254 | `exp_source_int` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 274 | `exp_source_prefix` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 310 | `dnx_dgamma` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 324 | `u_edges` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `F` | 338 | `x_from_u` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 347 | `trace_affine_u` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 378 | `build_piece_u` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 416 | `find_u_cell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 441 | `trace_piece_edge` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 562 | `trace_piece_u` | 多滞后分段仿射 u 特征线回溯；单滞后调用传 `Num_lag=1`。 |
| `S` | 583 | `char_core` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 620 | `remap_edges` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 644 | `char_update` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 678 | `ppm_cell_int` | PPM 单元积分 primitive。 |
| `S` | 695 | `semi_lagrangian_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 722 | `flux_split_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 765 | `flux_split_nonuniform` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 814 | `flux_seq_nonuniform` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 863 | `fullhide_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 886 | `fullhide_spacetime` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 930 | `logparabola_peak` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 958 | `gamma_active_hi` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 991 | `chi_active_hi` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `F` | 1010 | `max_xi_chi` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `F` | 1037 | `max_xi_uniform` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_dg_transport.f90`

DG 谱元网格、投影、正性核和特征线投影。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_dg_transport` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 49 | `dg_build_mesh` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 61 | `dg_build_coord` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 94 | `dg_initial_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 114 | `dg_project_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 135 | `dg_project_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 177 | `dg_kinetic_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 216 | `dg_scale_content` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 229 | `dg_limit_positive` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 256 | `dg_filter_positive` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 310 | `dg_advance_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 354 | `dg_dense_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 376 | `dg_transport_matrix` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 430 | `dg_char_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 447 | `dg_zero_bad` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 462 | `dg_project_char` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 521 | `dg_low_content` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 544 | `dg_interval_int` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 561 | `dg_back_x` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 589 | `dg_forward_x` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 608 | `dg_project_element` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 648 | `dg_solve_block` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 696 | `dg_solve_dense` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 738 | `dg_project_cells` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 774 | `dg_integral` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 789 | `dg_tail_fraction` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 815 | `dg_filter_mode` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 842 | `dg_is_troubled` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 862 | `dg_filter_factor` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 872 | `dg_jackson_factor` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 886 | `add_active_break` | 特征 Lorentz 因子/断点诊断；用于注入、冷却或活动网格边界。 |
| `S` | 899 | `sort_breaks` | 特征 Lorentz 因子/断点诊断；用于注入、冷却或活动网格边界。 |
| `S` | 915 | `allocate_spectral_mesh` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 929 | `ensure_reference_spectral` | workspace/cache 管理；只服务性能和内存复用，不改变物理语义。 |
| `S` | 946 | `ensure_projection_quadrature` | workspace/cache 管理；只服务性能和内存复用，不改变物理语义。 |
| `S` | 954 | `set_domain_bounds` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 982 | `fill_physical_nodes` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1000 | `lgl_nodes_weights` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1026 | `legendre_value_derivative` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1053 | `legendre_basis_values` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1067 | `gauss_nodes_weights` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1086 | `barycentric_weights` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1101 | `differentiation_matrix` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1116 | `locate_domain` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1134 | `interpolate_domain` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

### `src/Electron/hybrid_spectrum.f90`

热-非热混合谱归一化和特殊函数加速。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 4 | `hybrid_spectrum` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 48 | `thermal_int1` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 145 | `thermal_int12` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 305 | `cpl_integral` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 322 | `solve_theta` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 352 | `get_initial_theta` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 361 | `newton_method` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 396 | `hybrid_spec` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 424 | `spec_point` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Hadronic/hadronic_accel.f90`

强子加速时间、外部冷却限制和 gamma limit / injection rate。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_accel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 20 | `species_info` | 查询 hadronic species 的质量和电荷。 |
| `S` | 52 | `accel_time` | 费米加速时标。 |
| `S` | 75 | `syn_time` | 同步辐射冷却时标。 |
| `S` | 96 | `ext_time` | 外部光子场冷却时标。 |
| `S` | 113 | `inject_rate` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 158 | `gamma_limit` | 动力学、同步冷却和外部冷却限制给出的最大 Lorentz factor。 |
| `S` | 197 | `init_limit` | 未加入外部冷却前的动力学和同步辐射限制。 |
| `F` | 213 | `find_cross` | 查找 t_acc 与 t_ext 的交叉区间。 |
| `S` | 229 | `apply_xlimit` | 用对数插值给出外部冷却限制。 |
| `S` | 253 | `accel_calc` | 一次性计算加速、冷却、注入和最大能量限制。 |
| `F` | 276 | `trapz` | 梯形法则积分。 |

### `src/Hadronic/hadronic_bh.f90`

Bethe-Heitler 质子损失、pair source 和 photon loss kernel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_bh` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 17 | `bh_calc` | Bethe-Heitler pair/source/loss 算子。 |
| `F` | 77 | `loss_point` | 单个 photon bin 对 proton loss 的贡献。 |
| `F` | 86 | `bh_pair` | Bethe-Heitler pair 产生核；进入 pair source 与 photon sink。 |
| `F` | 104 | `bh_outer` | Bethe-Heitler 外层 omega 积分核。 |
| `F` | 124 | `bh_inner` | Bethe-Heitler 内层 ebar 积分核。 |
| `F` | 142 | `bh_sigma` | Blumenthal 1970 微分截面。 |
| `F` | 184 | `proton_loss` | Bethe-Heitler 质子能量损失核。 |
| `F` | 191 | `bh_phi` | Bethe-Heitler 能量损失 phi(x) 近似。 |
| `F` | 212 | `bh_rk3` | 三参数 RK 3/8 积分器。 |
| `F` | 237 | `bh_rk4` | 四参数 RK 3/8 积分器。 |

### `src/Hadronic/hadronic_base.f90`

强子网格、时间尺度、log-grid 校验和量子同步冷却公共工具。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_base` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 17 | `build_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 42 | `source_bounds` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 66 | `build_edges` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `F` | 89 | `shell_dt` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 103 | `dyn_time` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 113 | `proton_limit` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 127 | `check_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `F` | 154 | `quant_factor` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |

### `src/Hadronic/hadronic_decay.f90`

pi0、pi±、mu± decay、neutrino 和 electron channel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_decay` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 20 | `pi0_gamma` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 48 | `pion_decay` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 82 | `build_pion` | pion decay source 的 log-energy rate 预处理。 |
| `S` | 91 | `pion_madd` | charged-pion decay 到 helicity-resolved muon。 |
| `S` | 116 | `pion_nadd` | charged-pion prompt neutrino source。 |
| `S` | 151 | `muon_decay` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 181 | `build_muon` | muon decay source 的 log-energy rate 预处理。 |
| `S` | 191 | `muon_nadd` | muon decay 到 flavor-resolved neutrino。 |
| `S` | 217 | `muon_eadd` | muon decay 到 e± source。 |
| `S` | 255 | `decay_hummer` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 313 | `log_spacing` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 326 | `nu1_decay` | muon decay neutrino 谱函数 f_nu1(x,h)。 |
| `F` | 337 | `nu2_decay` | muon decay neutrino 谱函数 f_nu2(x,h)。 |
| `F` | 346 | `log_interp` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

### `src/Hadronic/hadronic_forward_1d.f90`

正向激波 formal 1D 强子 f2py 入口和 shell-sequence driver。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 3 | `had_syn_pol` | f2py/Python 调用边界；同步偏振诊断使用的 proton synchrotron polarization wrapper。 |
| `S` | 16 | `pg_operator` | f2py/Python 调用边界；Hummer 2010 p-gamma operator。 |
| `S` | 37 | `pair_production` | f2py/Python 调用边界；gamma-gamma pair source/loss operator。 |
| `S` | 54 | `pp_shell` | f2py/Python 调用边界；pp delta source/loss operator。 |
| `S` | 73 | `bethe_heitler` | f2py/Python 调用边界；Bethe-Heitler pair source 和 photon/proton loss。 |
| `S` | 89 | `hadronic_ic` | f2py/Python 调用边界；pion/muon/proton hadronic IC diagnostic operator。 |
| `S` | 114 | `decay_operator` | f2py/Python 调用边界；pion/muon decay 到 gamma、e± 和 neutrino。 |
| `S` | 146 | `cascade_sequence` | f2py/Python 调用边界；shell-sequence gamma-gamma pair/synch cascade。 |
| `S` | 172 | `hadronic_1d` | light 1D proton transport + proton synchrotron shell-sequence driver。 |
| `S` | 279 | `init_grid` | `hadronic_1d` 内部 proton grid 初始化；Jacobians 不能省略。 |
| `S` | 306 | `init_out` | `hadronic_1d` 输出频率和数组初始化。 |
| `S` | 326 | `inject_p` | `hadronic_1d` 壳层质子源项；必须同时满足粒子数和能量预算。 |
| `S` | 336 | `advance_p` | `hadronic_1d` 壳层质子输运推进。 |
| `S` | 345 | `advance_sec` | `hadronic_1d` Hummer p-gamma secondary chain 推进。 |
| `S` | 361 | `emit_syn` | `hadronic_1d` proton synchrotron emissivity 和 seed 计算。 |
| `S` | 352 | `formal_transport_1d` | formal 1D 强子主入口；推进 proton transport、p-gamma/BH/pp、secondary、photon survival 和 secondary e± source。 |

### `src/Hadronic/hadronic_formal.f90`

formal 1D 强子壳层序列实现层；Python/f2py 只通过 `formal_transport_1d` 进入这里。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `hadronic_formal` | formal 1D 强子 shell-sequence 模块命名空间。 |
| `S` | 9 | `formal_transport` | 按半径推进 proton injection/transport、pγ/BH/pp、secondary species、secondary radiation、photon survival 和 secondary e± source。 |

### `src/Hadronic/hadronic_shell.f90`

formal 1D 强子底层 shell primitive 与单位/投影 helper；f2py wrapper 已收窄到运行时真正需要的入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `hadronic_shell` | shell-level 强子 primitive 模块命名空间。 |
| `S` | 24 | `pp_shell` | pp delta source/loss operator；输出 gamma、neutrino、e± 源和 proton loss。 |
| `S` | 51 | `hic_shell` | proton/pion/muon hadronic IC operator；输出 IC emissivity 与投影系数。 |
| `S` | 83 | `hic_project` | hadronic IC shell emissivity 投影到 photon grid。 |
| `S` | 109 | `species_step` | n、pi、mu secondary species 同壳层保守推进。 |
| `S` | 172 | `inject_content` | 壳层注入能量预算到 species source content 的归一化。 |
| `S` | 194 | `scan_pmax` | 沿半径序列估计全局 proton 最大 Lorentz factor。 |
| `S` | 218 | `secondary_rad` | pion/muon synchrotron 与 IC shell emissivity。 |
| `S` | 247 | `secondary_project` | secondary radiation 从 hadron grid 投影到 photon grid。 |
| `S` | 307 | `loss_rates` | adiabatic、synchrotron 和 quantum-synch 连续损失率。 |
| `S` | 329 | `electron_seq` | secondary e± source 随壳层序列组装。 |
| `S` | 362 | `photon_loss` | photon loss rate 到 optical-depth/survival closure。 |
| `S` | 387 | `effective_time` | interaction loss rate 的有效时间积分。 |
| `S` | 413 | `pgamma_update` | pγ loss/re-injection 对 proton spectrum 的壳层更新。 |
| `S` | 433 | `proton_step` | proton injection、continuous loss 和 pγ update 的单壳层推进。 |
| `S` | 454 | `exp_sink` | 指数 sink primitive；用于已定义 interaction/loss closure。 |
| `S` | 470 | `rate_lum` | rate-per-energy 到 luminosity-per-frequency/energy 的壳层换算。 |
| `S` | 483 | `project_lum` | source energy grid 到目标 photon grid 的 luminosity 投影。 |
| `S` | 501 | `project_hic` | hadronic IC 多 species emissivity 投影。 |
| `S` | 521 | `pair_content` | pp/BH pair source 组合成电子方程使用的 content source。 |
| `S` | 535 | `shell_density` | shell content 到 density-per-GeV 的单位变换。 |
| `S` | 550 | `gamma_edges` | Lorentz-factor grid edges；守恒积分需要的 bin geometry。 |
| `S` | 572 | `proc_power` | secondary process power diagnostic。 |
| `S` | 611 | `pos_interp` | 正值 log-log grid projection；当前 secondary-feedback Python glue 仍使用。 |
| `S` | 642 | `gamma_source` | source-per-energy 到 source-per-gamma 的 grid projection。 |
| `S` | 658 | `dist_gev` | distribution-per-gamma 到 distribution-per-GeV 的 grid projection。 |
| `S` | 674 | `align_photon` | hadron/photon grid 对齐 helper。 |
| `S` | 691 | `shell_geom` | shell dr、dt geometry helper。 |
| `S` | 712 | `shell_volumes` | 壳层体积序列。 |

### `src/Hadronic/hadronic_ic.f90`

强子 secondary inverse-Compton kernel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_ic` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 21 | `ic_init` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 47 | `ic_operator` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 73 | `ic_apply` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 100 | `apply_species` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 110 | `sum_charged_pion_density` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `S` | 116 | `sum_muon_helicity_density` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `S` | 124 | `ic_build` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 149 | `ic_channel` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `F` | 176 | `ic_coeff` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |

### `src/Hadronic/hadronic_pg.f90`

Hummer 2010 p-gamma response、光子汇和 secondary family deposition。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_pg` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 91 | `pg_hummer` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 165 | `accumulate_pg_photon_loss` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 195 | `deposit_pions` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 214 | `deposit_baryons` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 231 | `deposit_shift` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 247 | `rates_res` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 266 | `rates_dir` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 283 | `rates_mul` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 300 | `loss_res` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 310 | `loss_dir` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 320 | `loss_mul` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 332 | `phloss_res` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 340 | `phloss_dir` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 347 | `phloss_mul` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 355 | `kernel_res` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 367 | `kernel_mul` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 381 | `kernel_dir` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 397 | `idir` | p-gamma Hummer response 或 secondary family deposition。 |

### `src/Hadronic/hadronic_cascade.f90`

gamma-gamma pair/synch shell-sequence cascade。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_cascade` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 15 | `cascade_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 48 | `validate_step` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 56 | `produce_pairs` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 66 | `build_egamma` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 71 | `cool_pairs` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 87 | `emit_syn` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 107 | `cascade_seq` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 175 | `validate_seq` | 系统边界校验；用于拒绝外部输入或正式 kernel contract 违反。 |
| `S` | 186 | `build_grids` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 198 | `shell_geom` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 213 | `pair_syn` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 232 | `cool_deposit` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 250 | `electron_loss` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 260 | `advance_energy` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `F` | 278 | `gamma_bin` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 292 | `advance_photon` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |

### `src/Hadronic/hadronic_pair.f90`

gamma-gamma pair injection 和 photon loss kernel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_pair` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 19 | `pair_operator` | gamma-gamma absorption、pair injection 和能量闭合主算子。 |
| `S` | 64 | `pair_injection` | pair 注入谱双重卷积和吸收功率重标定。 |
| `F` | 86 | `pair_bin` | 单个 electron bin 的 pair 注入积分。 |
| `S` | 104 | `energy_window` | 当前 electron energy 对应的 photon 积分窗。 |
| `S` | 115 | `energy_closure` | 将单电荷 pair source 归一到吸收功率的一半。 |
| `S` | 131 | `photon_loss` | phibar 卷积得到 photon absorption rate。 |
| `F` | 152 | `phibar` | 各向同性 gamma-gamma 截面的角度平均核。 |
| `F` | 170 | `phibar1` | 近阈区 phibar Taylor 分支。 |
| `F` | 180 | `phibar2` | 中间区 phibar 解析分支。 |
| `F` | 193 | `phibar3` | 高能区 phibar 渐近分支。 |
| `F` | 202 | `phisum` | phibar2 的辅助级数。 |
| `F` | 214 | `inner_pp` | pair-production 内层 photon 积分下界。 |
| `F` | 266 | `outer_pp` | pair-production 外层 photon 积分下界。 |
| `F` | 274 | `xlow` | pair-production kinematic 最小 photon 能量。 |
| `F` | 281 | `gmb` | f_minus/f_plus 积分分支边界。 |
| `F` | 288 | `fppm` | f_minus photon 积分边界函数。 |
| `F` | 295 | `fppp` | f_plus photon 积分边界函数。 |
| `F` | 302 | `rggd1` | Aharonian+1983 gamma-gamma 微分截面核。 |
| `F` | 368 | `xcmu` | 质心系能量上界。 |
| `F` | 381 | `xcml` | 质心系能量下界。 |
| `F` | 414 | `td2` | 微分截面中的 T 辅助函数。 |
| `F` | 436 | `a0ratio` | 截面核中的稳定 a0 ratio。 |
| `F` | 447 | `a0` | 截面核中的 a0 函数。 |
| `F` | 458 | `a0f` | a0 的基础分支实现。 |
| `F` | 471 | `a0ratiof` | a0ratio 的基础分支实现。 |
| `F` | 482 | `betag` | 相对论 beta(gamma)。 |
| `F` | 493 | `sign_int` | phibar2 级数符号。 |
| `F` | 504 | `grid_offset` | photon/electron 对数网格最小能量偏移。 |

### `src/Hadronic/hadronic_hummer.f90`

旧/诊断型 Hummer p-gamma shell 聚合 helper。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_hummer` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 17 | `hummer_shell` | 旧 Hummer p-gamma shell path，估计 photon survival 并推进 secondary species。 |
| `S` | 117 | `map_secondaries` | 将 p-gamma secondary source 映射到 species gamma 网格。 |
| `S` | 153 | `neutron_loss` | neutron p-gamma loss 插值和显式扣除。 |
| `S` | 172 | `pos_interp` | 正值 log-log 插值 helper。 |

### `src/Hadronic/hadronic_pp.f90`

Delta 近似 pp 算子和 secondary source。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_pp` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 20 | `pp_source` | pp delta source/loss kernel；输出 gamma、neutrino、e± 源和 proton loss。 |
| `S` | 51 | `validate_inputs` | pp delta 输入网格和 target density 检查。 |
| `S` | 62 | `set_options` | pp delta 能量份额和 charged/neutral pion 分支设置。 |
| `S` | 83 | `emit_secondaries` | 按 delta 能量映射输出三类 secondary source。 |
| `S` | 94 | `pp_sigma` | Kelner+2006 pp 非弹性截面。 |
| `F` | 121 | `pp_threshold` | pp 反应阈值动能。 |
| `S` | 127 | `secondary_source` | delta 近似 secondary source。 |
| `S` | 143 | `pos_interp` | 正值 log-log 插值 helper。 |
| `F` | 180 | `upper_bracket` | 单调数组 bracket 查找。 |
| `F` | 210 | `log_eval` | 单次 log-log 线性插值。 |

### `src/Hadronic/pp_models.f90`

Kelner/Kafexhiu 类 pp spectral model helper。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 4 | `pp_models` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 22 | `pp_threshold` | pp -> pi0 动能阈值。 |
| `F` | 28 | `pp_shape` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `F` | 63 | `high_shape` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 76 | `geant4_shape` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 111 | `egam_max` | pi0 衰变光子的最大能量。 |
| `S` | 118 | `pi0_source` | pi0 gamma source spectrum。 |
| `F` | 150 | `sigpi0` | 总 pi0 产生截面。 |
| `F` | 172 | `amax_model` | pp spectral model 归一化因子。 |
| `F` | 197 | `sig1pi` | Kafexhiu+2014 低能单 pion 截面分支。 |
| `F` | 207 | `sig2pi` | Kafexhiu+2014 低能双 pion 截面分支。 |
| `F` | 216 | `siginel` | Kafexhiu+2014 pp 非弹性截面 fit。 |
| `F` | 224 | `mult_sibyll` | SIBYLL pi0 multiplicity。 |
| `F` | 233 | `mult_qgsjet` | QGSJET pi0 multiplicity。 |
| `F` | 242 | `mult_geant4` | Geant4 pi0 multiplicity。 |
| `F` | 249 | `mult_pythia` | Pythia8 pi0 multiplicity。 |

### `src/Hadronic/hadronic_rad.f90`

质子同步辐射和偏振。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_rad` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 16 | `syn_xarg` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 31 | `syn_mass` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 66 | `syn_kernel` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 77 | `proton_syn` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 134 | `syn_polarization` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 162 | `polar_integral` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |

### `src/Hadronic/hadronic_reverse_1d.f90`

反向激波轻量质子输运和质子同步辐射入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 3 | `reverse_hadronic_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 40 | `init_grid` | 基于全部 RS shell 的最大加速能量建立 proton grid。 |
| `S` | 55 | `advance_shell` | 单个 RS shell 的 proton 注入和 log-gamma transport。 |
| `S` | 72 | `emit_syn` | 可选输出当前 shell 的 proton synchrotron。 |

### `src/Hadronic/hadronic_secondary.f90`

pion/muon synchrotron 和 IC secondary radiation。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_secondary` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 22 | `secondary_calc` | pion/muon secondary radiation 顶层算子。 |
| `S` | 45 | `apply_rad` | 应用已预计算的同步辐射和 IC 映射。 |
| `S` | 66 | `init_syn` | 建立 pion/muon 同步辐射核矩阵。 |
| `S` | 91 | `init_ic` | 建立 pion/muon IC 能量偏移和上界索引。 |
| `S` | 114 | `pion_syn` | π 介子同步辐射冷却率。 |
| `S` | 132 | `muon_syn` | μ 子同步辐射冷却率。 |
| `S` | 153 | `pion_ic` | π 介子 inverse-Compton 冷却率。 |
| `S` | 175 | `muon_ic` | μ 子 inverse-Compton 冷却率。 |
| `F` | 201 | `syn_kernel` | AM3 分段超相对论同步辐射核。 |
| `S` | 235 | `build_ic` | 为一个粒子种类建立 IC 偏移和上界索引。 |
| `S` | 260 | `ic_channel` | 单个 IC 通道卷积。 |
| `F` | 291 | `ic_coeff` | Thomson 截面和质量缩放的 IC 系数。 |
| `S` | 302 | `check_density` | secondary radiation 输入密度有限性检查。 |

### `src/Hadronic/hadronic_species.f90`

neutron、pion、muon explicit species transport。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_species` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 25 | `div_rate` | 球对称膨胀发散率，用于绝热冷却。 |
| `S` | 39 | `synch_loss` | charged species 同步冷却率。 |
| `S` | 70 | `ad_loss` | species 绝热冷却率。 |
| `S` | 89 | `species_advance` | 七分量 neutron/pi/mu 输运主推进。 |
| `S` | 154 | `advance_one` | 单一 species 的 decay-transport-decay 推进。 |
| `S` | 226 | `validate_gamma` | gamma 网格正值、递增和点数检查。 |
| `S` | 248 | `validate_positive` | transport/source 谱非负检查。 |
| `F` | 263 | `log_spacing` | gamma 网格平均对数间距检查。 |

### `src/Hadronic/hadronic_transport_kernel.f90`

质子注入、连续损失和 log-gamma 有限体积推进。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_transport` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 10 | `proton_inject` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 43 | `proton_loss` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 63 | `advance_loggamma` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |

### `src/Hadronic/hadronic_transport_remap_kernel.f90`

强子能量冷却重映射 primitive。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_remap` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 9 | `remap_loggamma` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `F` | 39 | `remap_target` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

### `src/Interpolation/SED_interpolation.f90`

shell-level 和 chi-resolved EATS/Doppler 投影。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 9 | `sed_interpolation` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 108 | `sed_adaptive_theta` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 168 | `integrate_theta_cell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 198 | `project_theta_sample` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 230 | `project_shell_segment` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 252 | `chi_escape` | finite-q shell SSA escape factor。 |
| `F` | 266 | `lower_bound_real8` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 285 | `sed_interpolation_chi` | top-hat chi-resolved FS synchrotron+SSA lightcurve projection。 |
| `S` | 393 | `project_chi` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 407 | `chi_state` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 422 | `accumulate_chi_cell_source` | chi cell source + SSA escape 累加。 |
| `S` | 446 | `sed_chi_electron` | direct-electron chi projection；先批量计算 chi-local synchrotron/SSA 再投影。 |
| `S` | 481 | `sed_chi_ring` | axisymmetric structured chi ring projection；输入为预计算 `F_ring/Tau_ring`。 |
| `S` | 567 | `project_ring` | structured ring EATS/Doppler segment projection。 |
| `S` | 583 | `ring_source` | structured ring chi cell source + SSA escape 累加。 |

### `src/Interpolation/SED_interpolation_structured.f90`

结构化喷流 theta/theta-phi shell-level 投影。该文件仍由 `structured_jet_1d` Fortran 聚合入口调用；旧 Python binding 不再从 `src.Interpolation` 导出。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 9 | `sed_interpolation_structured` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 103 | `project_structured_segment` | 结构化喷流 patch/角向投影调度。 |
| `S` | 125 | `sed_structured_phi` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 202 | `project_phi` | 结构化喷流 patch/角向投影调度。 |

### `src/Interpolation/interpolation_common.f90`

SED 插值共享累加 primitive。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `interpolation_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 13 | `accum_logsed` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 52 | `accum_shifted` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

### `src/Radiation/quantum_synch.f90`

量子同步辐射修正因子。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 3 | `quantum_synch` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 19 | `quantum_chi_parameter` | 量子同步辐射参数 \(\chi\)；按粒子质量缩放临界磁场。 |
| `F` | 29 | `quantum_cooling` | Landau/QED 同步冷却压制因子。 |

### `src/Radiation/rad_common.f90`

辐射共享 Simpson 权重、transfer、pair cross-section 和 seed 转换。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `rad_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 16 | `compute_simpson_weights` | 积分权重或求积 primitive；影响谱积分精度。 |
| `F` | 35 | `rad_interp` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 61 | `transfer_factor` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 72 | `syn_kernel` | 近似同步辐射核的低频、指数段和高频分支。 |
| `S` | 89 | `pair_grid` | gamma-gamma pair-production 能量网格和测度准备。 |
| `F` | 109 | `pair_sigma` | gamma-gamma pair-production 截面。 |
| `S` | 128 | `pair_tau` | 对头碰撞近似下的 gamma-gamma 光深。 |
| `S` | 159 | `syn_seed_chi` | chi-resolved synchrotron/SSA batch kernel；单列调用传 `Num_chi=1`。 |
| `S` | 279 | `syn_flux_chi` | chi-resolved projection-only synchrotron flux/tau batch kernel。 |

### `src/Radiation/pair_absorption.f90`

观测侧 gamma-gamma 吸收光深。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 7 | `annihilation` | f2py/Python 调用边界；对 photon field、角度权重和 pair 截面积分后返回 absorption factor。 |
| `S` | 85 | `mid_seed` | 把种子光子场插值到相邻频率 cell 中点。 |
| `S` | 96 | `angle_weights` | 构造各向异性 gamma-gamma absorption 的角度权重。 |
| `S` | 109 | `build_sigma` | 预计算 pair-production 截面窗口。 |
| `S` | 127 | `set_window` | 为给定角度和目标频率定位有效 seed 频率窗口。 |
| `S` | 149 | `clear_window` | 清空无贡献 pair-production 窗口。 |
| `F` | 157 | `first_gt` | 有效 seed 频率窗口的下界二分搜索。 |
| `F` | 176 | `last_lt` | 有效 seed 频率窗口的上界二分搜索。 |

### `src/Radiation/ssc_spectrum.f90`

SSC 谱、KN/Jones kernel、均匀/非均匀电子谱积分。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 6 | `ssc_spec` | 均匀电子网格 SSC f2py/Python 调用边界。 |
| `F` | 87 | `first_gt` | 电子 Lorentz factor 下界二分搜索。 |
| `F` | 106 | `first_ge` | 电子 Lorentz factor 闭区间下界二分搜索。 |
| `S` | 125 | `uniform_weights` | 均匀网格 Simpson 权重和电子尾积分预处理。 |
| `S` | 143 | `kn_tables` | 预计算 Jones/KN 散射核表。 |
| `S` | 160 | `uniform_bounds` | 均匀 SSC seed/electron 有效积分边界。 |
| `F` | 177 | `low_sum` | 低 seed 侧 SSC 积分。 |
| `F` | 200 | `high_tail` | 高 seed 侧 SSC 尾积分。 |
| `S` | 220 | `uniform_point` | 单个 shell/frequency 的均匀网格 SSC 累加。 |
| `S` | 242 | `ssc_spec_nonuniform` | 非均匀电子网格 SSC f2py/Python 调用边界。 |
| `S` | 297 | `check_inputs` | 系统边界校验；用于拒绝外部输入或正式 kernel contract 违反。 |
| `S` | 312 | `build_state` | 非均匀网格重构和尾矩预处理。 |
| `S` | 326 | `build_shell` | 单 shell 非均匀电子谱重构。 |
| `S` | 344 | `limit_slope` | minmod 斜率限制。 |
| `S` | 372 | `add_tail` | 非均匀电子谱 gamma 矩尾积分累加。 |
| `F` | 385 | `gauss_point` | 二点 Gauss-Legendre 节点。 |
| `S` | 397 | `nonuniform_point` | 单个 shell/frequency 的非均匀网格 SSC 累加。 |
| `F` | 438 | `first_cell` | 非均匀电子 cell 下界二分搜索。 |
| `F` | 469 | `ssc_minmod` | minmod 斜率限制器。 |
| `F` | 482 | `profile_value` | 线性重构剖面值。 |
| `F` | 493 | `gamma_moment` | 线性重构剖面的 gamma 矩解析积分。 |
| `F` | 513 | `low_cell` | 低 seed 侧单电子 cell SSC 积分。 |
| `F` | 552 | `low_integral` | 低 seed 侧完整电子网格 SSC 积分。 |
| `F` | 577 | `high_tail` | 高 seed 侧 SSC 尾积分。 |

### `src/Radiation/syn_polarization.f90`

同步辐射偏振 F/G kernel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `syn_polarization` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 16 | `synchrotron_fg_kernel` | 同步辐射 F/G kernel，含小/大 x 渐近式和中段积分。 |
| `S` | 42 | `synchrotron_polarized_components` | 单粒子两个线偏振分量。 |
| `S` | 55 | `fg_integral` | 有限 u 区间上的复合 4 点 Gauss 积分。 |
| `S` | 85 | `add_node` | 累加单个 Gauss 节点的 F/G 被积函数值。 |

### `src/Structured/structured_jet_1d.f90`

结构化喷流 patch 调度、缓存和 Fortran 聚合入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 1 | `jet_flux_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 57 | `run_axisymmetric` | 结构化喷流 patch/角向投影调度。 |
| `S` | 93 | `run_nonaxisymmetric` | 结构化喷流 patch/角向投影调度。 |
| `S` | 131 | `structured_solve_axisymmetric` | 结构化喷流 patch/角向投影调度。 |
| `S` | 195 | `register_axis_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 212 | `solve_axis_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 232 | `structured_solve_nonaxisymmetric` | 结构化喷流 patch/角向投影调度。 |
| `S` | 302 | `register_phi_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 322 | `solve_phi_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 340 | `copy_phi_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 354 | `structured_solve_element` | 结构化喷流 patch/角向投影调度。 |
| `S` | 515 | `apply_absorption` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |

## 修改时的最短自检

- 只改网页文档：`tools/check_text_encoding.py`、`git diff --check`、`mkdocs build --strict --site-dir /tmp/asgard_mkdocs_site`。
- 改 dynamics：重建 `Dynamics_forward` / `Dynamics_reverse`，并跑 RS smoke 或受影响 lightcurve smoke。
- 改 electron：重建对应 electron module；显式从 `/tmp` 跑 `-Wline-truncation` 源闭包；看电子谱、`nu_a`、seed 和光变是否平滑。
- 改 hadronic：重建 `hadronic_forward_1d` 或 `hadronic_reverse_1d`；检查 formal grid contract、proton/secondary/photon survival 预算。
- 改 interpolation：重建 `SED_interpolation` 或 `SED_interpolation_structured`；比较 lightcurve/SED/chi projection 是否只改变预期投影层。
