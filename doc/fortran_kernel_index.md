# Fortran 程序单元索引

本文是当前工作树的 Fortran kernel 索引。它面向需要逐个进入子程序读算法的人：先看 f2py 入口和物理阶段，再进入文件内的 `module`、`subroutine`、`function`。更高层的物理到算法映射见 `doc/physics_algorithm_crosswalk.md`，运行主链见 `doc/call_chain.md`。

当前索引按 ASGARD 自有 Fortran 数值核抽取，排除第三方固定格式特殊函数依赖，共 735 个程序单元：44 个 module、501 个 subroutine、190 个 function。行号是生成本页时的源文件位置。

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
| `electron_forward_fullhide_1d` | `src/Electron` | `ELECTRON_HISTORY_SOURCES + electron_forward_fullhide_1d` | `fs_fullhide_1d; fs_fullhide_coupled` | 默认 1D 电子输运和 joint feedback coupled pass。 |
| `electron_forward_fullhide_1d_hybrid` | `src/Electron` | `ELECTRON_HISTORY_SOURCES_HZ + electron_forward_fullhide_1d_hybrid` | `fs_fullhide_hz` | 热/非热混合谱路径。 |
| `electron_forward_charint_1d` | `src/Electron` | `ELECTRON_COMMON_SOURCES + electron_forward_charint_1d` | `fs_charint_1d` | 1D 特征线输运对照。 |
| `electron_forward_dg_1d` | `src/Electron` | `ELECTRON_DG_1D_SOURCES + electron_forward_dg_1d` | `fs_dg_1d` | P12 LGL-DG 正向电子输运。 |
| `electron_forward_charint_2d` | `src/Electron` | `ELECTRON_2D_SOURCES + electron_forward_transport_2d` | `fs_transport_2d` | finite-q shell 2D 电子输运；别把 chi_grid 当强子局域坐标。 |
| `electron_forward_t2g1_1d` | `src/Electron` | `ELECTRON_COMMON_SOURCES + electron_forward_t2g1_1d` | `fs_t2g1_1d` | 方法比较电子输运。 |
| `electron_forward_weno5_1d` | `src/Electron` | `ELECTRON_COMMON_SOURCES + electron_forward_weno5_1d` | `fs_weno5_1d` | WENO5 方法比较电子输运。 |
| `electron_reverse_kernel` | `src/Electron` | `ELECTRON_REVERSE_SOURCES + electron_reverse_kernel` | `electron_reverse_evolve; multiple_synch; branch_reaccel` | RS primary/secondary electron transport 和 reacceleration。 |
| `electron_radiation` | `src/Electron` | `ELECTRON_RADIATION_SOURCES` | `nua_solve; syn_state; syn_transfer` | 电子同步/SSA/seed 低层核；通常经其他 entry 调用。 |
| `ssc_spectrum` | `src/Radiation` | `rad_common + ssc_spectrum` | `ssc_spec; ssc_spec_nonuniform` | SSC spectrum 和 KN/Jones 积分。 |
| `pair_absorption` | `src/Radiation` | `rad_common + pair_absorption` | `annihilation` | 观测侧 gamma-gamma absorption。 |
| `SED_interpolation` | `src/Interpolation` | `rad_common + interpolation_common + SED_interpolation` | `sed_interpolation; sed_adaptive_theta; sed_interpolation_chi; sed_chi_electron; sed_chi_ring` | 观测者投影；Python 公开绑定集中在此模块。 |
| `SED_interpolation_structured` | `src/Interpolation` | `interpolation_common + SED_interpolation_structured` | `sed_interpolation_structured; sed_structured_phi` | structured_jet_1d 内部 shell-level structured projection；不再从 src.Interpolation 暴露旧 Python 绑定。 |
| `hadronic_forward_1d` | `src/Hadronic` | `HADRONIC_1D_SOURCES + hadronic_forward_1d` | `formal_transport_1d; shell operators` | formal 1D proton/secondary/photon-loss shell sequence。 |
| `hadronic_reverse_1d` | `src/Hadronic` | `HADRONIC_1D_SOURCES + hadronic_reverse_1d` | `reverse_hadronic_1d` | RS light proton transport + proton synchrotron。 |
| `structured_jet_1d` | `src/Structured` | `STRUCTURED_JET_1D_SOURCES + structured_jet_1d` | `jet_flux_1d` | 结构化喷流 Fortran 聚合入口。 |

Fortran 改动后的最低门槛见 `doc/validation_and_benchmarks.md`。文档-only 改动不需要重建扩展；若修改本页对应源码，必须跑受影响 module 的 `build_extensions.py --force`、独立 `-Wline-truncation` 源闭包检查和最小相关直接运行。

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
| `S` | 65 | `forward_rk4` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 172 | `forward_rhs` | ODE/PDE 右端项；定义物理源汇和动力学变量导数。 |

### `src/Dynamics/dynamics_density_profile.f90`

外介质密度、density jump/profile 状态和 tabulated profile 插值。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `dynamics_density_profile` | 介质密度和 density-jump/profile 状态；被 dynamics、electron、structured 路径按需引用。 |
| `F` | 24 | `uniform_density` | 纯函数；集中判定 ISM、legacy jump、modern jump 和 tabulated profile 是否满足均匀介质 fast-path 合同。 |
| `S` | 32 | `density_profile` | 唯一局部介质密度合同；匹配 `R0`、wind-to-ISM、Gaussian jump 和 tabulated profile。 |
| `S` | 91 | `density_moment` | 解析计算 `integral_0^R r^2 n(r)dr`，供 Forward/Reverse Dynamics 与 Electron 共享初始扫掠量。 |
| `S` | 167 | `set_density_profile` | 从 `Boundary` 解包互斥的 density jump/profile，并预存高精度 log-log 幂律指数。 |
| `S` | 242 | `tab_moment` | 逐段解析积分首末外推的 tabulated profile，显式处理 `s=-3` 极限和原点发散。 |
| `S` | 293 | `gauss_moments` | 解析 Gaussian 的零阶和二阶有限区间矩；同侧远尾使用稳定 `erfc` 差。 |
| `F` | 314 | `tab_logdensity` | 在 log 域执行 tabulated profile 的分段插值与首末幂律外推。 |
| `S` | 339 | `tab_density` | 将共享的 log-density 合同转换为局部数密度。 |

### `src/Dynamics/reverse_jump_conditions.f90`

动力学相关 Fortran 数值核。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `reverse_jump_conditions` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 8 | `reverse_contact` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 39 | `pressure_difference` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 52 | `solve_comp` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 84 | `unity_comp` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 107 | `shock_diff` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Dynamics/reverse_rhs.f90`

反向激波 RHS，显式处理 shock 注入、磁化惯性、区域 3 热演化和次级分支导数。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `reverse_rhs_module` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 9 | `reverse_dynamics_rhs` | ODE/PDE 右端项；定义反向激波主状态和次级分支导数。 |
| `F` | 212 | `parent_ready` | 判断上一级次级分支是否可作为当前上游。 |
| `S` | 222 | `branch_derivs` | 次级 RS 分支质量、热能、体积和注入诊断导数。 |
| `S` | 301 | `branch_density` | density-jump 分支的局部密度和 source 权重。 |

### `src/Dynamics/reverse_shock.f90`

反向激波 `Dynamics_reverse.dynamics_reverse` f2py 入口：first RS 推进、等待成 shock、pre-crossing swept-mass 推进、crossing event split、post-crossing log-time RK、density-jump 分支扫描和 branch 输出都在同一个输出循环中完成。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 5 | `dynamics_reverse` | `Dynamics_reverse` f2py/Python 调用边界；解包 public `Boundary`，推进 first 分支并同步扫描 density-jump multiple 分支。 |
| `S` | 109 | `init_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 142 | `init_dec` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 157 | `init_cross` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 173 | `advance_out` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 184 | `save_prim` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 198 | `save_r3` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 210 | `save_g34` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 220 | `advance_m3` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 353 | `rk_m3` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 425 | `advance_log` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 547 | `rk_log` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 587 | `advance_to` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 622 | `wait_trial` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 644 | `wait_root` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 673 | `pressure_ok` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 705 | `init_output` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 714 | `init_branch` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 729 | `init_total` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 740 | `init_diag` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 754 | `init_history` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 761 | `init_event` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 774 | `init_scan` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 790 | `step_scan` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 815 | `close_events` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 831 | `scan_events` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 884 | `event_root` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 918 | `event_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 962 | `save_multi` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 997 | `save_branch` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1032 | `sum_total` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1045 | `sum_diag` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1078 | `get_upstream` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1114 | `sum_weight` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1136 | `freqs` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1155 | `branch_density` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1183 | `parent_ready` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1198 | `unpack` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Dynamics/reverse_shock_mhd_jump.f90`

有限强度反向激波 MHD jump 公式。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `reverse_shock_mhd_jump` | MHD jump 公式模块；不包含时间推进状态。 |
| `F` | 10 | `rs_vegas_ud` | 有限强度 MHD jump 解析根；不要用 ultra-relativistic 近似替代。 |
| `F` | 57 | `cubic_max` | 用正判别式与 `atan2` 求原/倒数三次多项式的最大根。 |
| `F` | 71 | `cubic_disc` | Vegas 三次方程判别式的正系数因式分解，避免近双根消减。 |
| `S` | 92 | `rs_mhd_state` | 单次解析根返回下游四速度、压缩比、热比内能和快模判据。 |

### `src/Dynamics/reverse_shock_state.f90`

动力学相关 Fortran 数值核。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `reverse_shock_state` | 模块命名空间；集中声明本文件共享 procedure。 |

### `src/Electron/electron_common.f90`

电子边界向量解包、注入断点、初始谱、诊断公共工具。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 16 | `electron_unpack_boundary` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 48 | `electron_initialize_spectrum` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 80 | `electron_gm_exact` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 107 | `electron_gc_loss` | 冷却或能量损失计算；必须和 emissivity/source 单位保持一致。 |
| `S` | 156 | `electron_injection_prefactor` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 166 | `electron_source_bounds` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 201 | `electron_relerr_max` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 218 | `electron_initial_density` | 在 `R(1)` 复用共享局部密度和解析扫掠矩，初始化首列总电子数。 |

### `src/Electron/electron_cooling_ic_kernel.f90`

逆康普顿/KN 冷却核，含 IC 网格缓存、旧积分冷却率和 emissivity-budget 冷却率。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_ic_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 17 | `ensure_ic_grid` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 27 | `ic_grid_current` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 39 | `ic_seed_current` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `F` | 51 | `ic_gamma_current` | 特征 Lorentz 因子、断点或谱峰诊断。 |
| `S` | 63 | `rebuild_ic_grid` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 83 | `electron_ic_loss` | 冷却或能量损失计算；必须和 emissivity/source 单位保持一致。 |
| `S` | 106 | `accumulate_ic_loss` | 冷却或能量损失计算；必须和 emissivity/source 单位保持一致。 |
| `S` | 148 | `electron_ic_budget` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 184 | `accumulate_budget` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 212 | `low_seed_kernel` | Jones/KN 低 seed 侧核函数。 |

### `src/Electron/electron_cooling_kernel.f90`

电子冷却门面模块；只保留统一 `forward_cooling` 入口，具体 SSA/IC/Y 物理核从对应实现模块直接调用。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_cooling_kernel` | 冷却组装门面模块；不承载 SSA/IC/Y 具体物理核。 |
| `S` | 17 | `forward_cooling` | 统一准备 Compton auxiliary 并组装 χ-resolved 或单列正向激波冷却率。 |

### `src/Electron/electron_cooling_ssa_kernel.f90`

同步自吸收回热/冷却核，含 seed 频率缓存、SSA 几何映射和 χ 批量版本。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_ssa_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 20 | `ensure_ssa_cache` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 54 | `advance_ssa_cursor` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 76 | `build_ssa_geometry` | 构建 SSA 低频/高频截面区间和 prefactor。 |
| `S` | 122 | `electron_ssa_loss` | 冷却或能量损失计算；必须和 emissivity/source 单位保持一致。 |
| `S` | 177 | `accumulate_ssa_cell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 243 | `clipped_ssa_segment` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_cooling_y_kernel.f90`

Compton-Y 辅助量核，含 Nakar 数值积分和 Fan 解析分段模型。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_y_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 18 | `ensure_y_nakar` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 58 | `electron_y_nakar` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 92 | `electron_y_fan` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_dg_transport.f90`

DG 谱元网格、投影、正性核和特征线投影。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_dg_transport` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 49 | `dg_build_mesh` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 64 | `dg_build_coord` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 98 | `dg_initial_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 120 | `dg_project_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 143 | `dg_project_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 187 | `dg_kinetic_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 228 | `dg_scale_content` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 243 | `dg_limit_positive` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 272 | `dg_filter_positive` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 329 | `dg_advance_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 377 | `dg_dense_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 403 | `dg_transport_matrix` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 459 | `dg_char_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 478 | `dg_zero_bad` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 495 | `dg_project_char` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 557 | `dg_low_content` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 583 | `dg_interval_int` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 603 | `dg_back_x` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 633 | `dg_forward_x` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 654 | `dg_project_element` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 694 | `dg_solve_block` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 744 | `dg_solve_dense` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 789 | `dg_project_cells` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 829 | `dg_integral` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 846 | `dg_tail_fraction` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 875 | `dg_filter_mode` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 904 | `dg_is_troubled` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 926 | `dg_filter_factor` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 938 | `dg_jackson_factor` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 954 | `add_active_break` | 特征 Lorentz 因子/断点诊断；用于注入、冷却或活动网格边界。 |
| `S` | 969 | `sort_breaks` | 特征 Lorentz 因子/断点诊断；用于注入、冷却或活动网格边界。 |
| `S` | 987 | `allocate_spectral_mesh` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1003 | `ensure_reference_spectral` | workspace/cache 管理；只服务性能和内存复用，不改变物理语义。 |
| `S` | 1022 | `ensure_projection_quadrature` | workspace/cache 管理；只服务性能和内存复用，不改变物理语义。 |
| `S` | 1032 | `set_domain_bounds` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1064 | `fill_physical_nodes` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1084 | `lgl_nodes_weights` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1112 | `legendre_value_derivative` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1141 | `legendre_basis_values` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1157 | `gauss_nodes_weights` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1178 | `barycentric_weights` | 积分权重或求积 primitive；影响谱积分精度。 |
| `S` | 1195 | `differentiation_matrix` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1212 | `locate_domain` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1232 | `interpolate_domain` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

### `src/Electron/electron_energy_coordinate_common.f90`

log-gamma 与 log-four-velocity 坐标映射。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_coord_common` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 18 | `coord_from_xg` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 33 | `xg_from_coord` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 45 | `gamma_from_coord` | 特征 Lorentz 因子、断点或谱峰诊断。 |
| `F` | 57 | `dxg_dcoord` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 73 | `build_fourvel_grid` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_forward_charint_1d.f90`

1D 特征线电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 3 | `fs_charint_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 242 | `prepare_characteristic_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 284 | `write_final_diag` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_forward_dg_1d.f90`

1D LGL-DG 电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 3 | `fs_dg_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 115 | `init_fourvel_grid` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 138 | `prepare_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 159 | `write_rad_breaks` | 特征 Lorentz 因子、断点或谱峰诊断。 |
| `S` | 173 | `write_final` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 195 | `remesh_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 234 | `advance_substep` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 280 | `prepare_dg_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 299 | `prepare_thermal_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 308 | `limit_jump_step` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 329 | `limit_one_jump` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 353 | `jump_log_slope` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 371 | `add_jump_deriv` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 383 | `dg_source_xmax` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 391 | `dg_active_xmax` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 405 | `scale_dg_content` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 418 | `write_positive_output` | DG 状态投影回用户 gamma 网格，并在输出边界显式施加非负电子谱。 |

### `src/Electron/electron_forward_fullhide_1d.f90`

默认 1D 隐式有限体积电子输运入口，含 joint coupled 入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 7 | `fs_fullhide_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 101 | `init_electron_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 133 | `init_fourvel_grid` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 150 | `prepare_fullhide_shell` | 单个 FS shell 的密度、磁场、注入能标、同步谱、SSA break 和冷却率准备。 |
| `S` | 188 | `write_finaldiag` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 212 | `advance_fixed_shell` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 227 | `advance_uniform_shell` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 266 | `advance_general_shell` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 311 | `advance_adaptive_shell` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 412 | `fs_fullhide_coupled` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 508 | `prepare_coupled_shell` | joint electron-photon shell 的辐射场、SSA break 和耦合冷却率准备。 |
| `S` | 551 | `write_finaldiag` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 574 | `advance_coupled_shell` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |

### `src/Electron/electron_forward_fullhide_1d_hybrid.f90`

混合热/非热谱形入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 3 | `fs_fullhide_hz` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 253 | `init_fourvel_grid` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 269 | `write_finaldiag` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 293 | `build_hybrid_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 310 | `build_hybrid_count` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_forward_slc1_1d.f90`

半拉格朗日一阶电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 3 | `fs_slc1_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |

### `src/Electron/electron_forward_t2g1_1d.f90`

T2G1 方法比较电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 3 | `fs_t2g1_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 78 | `prepare_t2g1_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 105 | `write_t2g1_cooling` | 冷却或能量损失计算；必须和 emissivity/source 单位保持一致。 |
| `S` | 128 | `advance_t2g1_substep` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |

### `src/Electron/electron_forward_transport_2d.f90`

2D finite-q shell 电子输运和投影预处理入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 7 | `fs_transport_2d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 303 | `init_front` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 396 | `advance_shell` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 530 | `store_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 599 | `finish_output` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 682 | `transport_step_fullhide` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 775 | `update_bchi` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 794 | `shock_scales` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 805 | `assemble_cooling_chi` | 冷却或能量损失计算；必须和 emissivity/source 单位保持一致。 |
| `S` | 841 | `break_freqs` | 特征 Lorentz 因子、断点或谱峰诊断。 |

### `src/Electron/electron_forward_weno5_1d.f90`

WENO5 方法比较电子输运入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 3 | `fs_weno5_1d` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 75 | `prepare_weno_shell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 99 | `write_weno_cooling` | 冷却或能量损失计算；必须和 emissivity/source 单位保持一致。 |
| `S` | 131 | `write_finaldiag` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 153 | `advance_weno_substep` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 183 | `load_weno_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 198 | `compute_weno_fluxes` | 刷新零阶外推鬼点并按速度符号构造 WENO5 数值通量。 |
| `S` | 218 | `rk_weno_stage` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 244 | `weno5_positive_flux` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 281 | `weno5_negative_flux` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_injection_profiles.f90`

非热/热电子注入谱、log-cell 积分和源项归一化。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 5 | `electron_injection_profiles` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 26 | `exp_cutoff` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 36 | `pl_params` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 62 | `dnx_cutoff` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 80 | `dnx_gauss3` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 104 | `dny_gauss3` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 130 | `dnx_segment` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 152 | `dny_segment` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 174 | `log_edges` | 介质密度或 density-jump 分支；直接影响 swept mass、动力学和注入源项。 |
| `S` | 190 | `init_powerlaw` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 211 | `init_edges` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 251 | `init_coord` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 292 | `source_edges` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 319 | `source_coord` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 346 | `kinetic_edges` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 386 | `kinetic_coord` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 425 | `thermal_theta` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 435 | `thermal_shape` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 464 | `add_thermal` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 481 | `thermal_pop` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/electron_radiation_kernel.f90`

电子同步辐射、SSA、nu_a、偏振和 seed photon 核。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `electron_radiation_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 20 | `greater_from` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 43 | `greater_window` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 65 | `besselk` | 同步辐射 Bessel kernel 插值/渐近 primitive。 |
| `S` | 125 | `syn_state` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 148 | `syn_fixed` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `F` | 212 | `emit_simpson` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 227 | `tau_simpson` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `S` | 241 | `store_syn` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `S` | 260 | `reduce_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 289 | `project_syn` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 323 | `syn_cyclotron` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `S` | 343 | `add_cyclotron` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 383 | `syn_fx` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `F` | 395 | `linear_interp` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 407 | `syn_integrand` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `F` | 418 | `pl_interp` | log-log/power-law 插值 primitive。 |
| `S` | 443 | `log_gauss2` | 积分权重或求积 primitive；影响谱积分精度。 |
| `F` | 461 | `pl_integral` | power-law cell 积分 primitive。 |
| `F` | 478 | `tau_kernel` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `S` | 488 | `syn_gauss` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `S` | 511 | `tau_gauss` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `S` | 542 | `syn_adapt` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `S` | 573 | `syn_transfer` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 603 | `nua_solve` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 676 | `evaluate_tau` | 光深、SSA transfer 或 photon survival 相关算子。 |
| `S` | 694 | `refine_tau` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `S` | 740 | `nua_path` | 2D/chi cell-level SSA break diagnostic。 |
| `S` | 757 | `nua_fromtau` | 从已计算 optical-depth grid 求 SSA break；避免重复 root search。 |
| `S` | 800 | `tau_root` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |

### `src/Electron/electron_reverse_kernel.f90`

反向激波电子输运、post-crossing map、次级 RS 和 reacceleration 分支。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `electron_reverse_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 36 | `electron_reverse_evolve` | 反向激波电子演化入口；串接主 RS 注入、冷却、DG/有限体积输运和输出投影。 |
| `S` | 222 | `compute_cooling` | 冷却或能量损失计算；必须和 emissivity/source 单位保持一致。 |
| `S` | 255 | `advance_shell` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 361 | `advance_postcross` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `F` | 398 | `postcross_map` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 429 | `init_dg_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 439 | `remesh_dg_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 457 | `ensure_dg_work` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 467 | `prepare_substep_state` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `F` | 483 | `shell_linear_value` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 494 | `dg_upper_break` | 特征 Lorentz 因子、断点或谱峰诊断。 |
| `F` | 504 | `dg_active_xmax` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 519 | `dg_low_break` | 特征 Lorentz 因子、断点或谱峰诊断。 |
| `S` | 529 | `advance_dg_inject` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 541 | `dg_inject_break` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 548 | `advance_dg_front` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 568 | `multiple_synch` | 聚合次级 RS 分支同步辐射谱。 |
| `S` | 598 | `branch_reaccel` | 次级 RS 再加速电子分支输出。 |
| `S` | 667 | `reaccel_grid` | 构造再加速分支 gamma 网格。 |
| `S` | 692 | `transfer_parent` | 把父分支电子谱转移到再加速分支网格。 |
| `S` | 732 | `advance_reaccel` | 再加速分支的输运推进。 |
| `S` | 778 | `prepare_branch_shell` | 读取单个次级分支在当前壳层的动力学和辐射参数。 |
| `S` | 806 | `branch_inject` | 次级分支注入归一化和能标。 |
| `S` | 816 | `boost_log` | 对数谱的 Lorentz boost 重映射。 |
| `S` | 842 | `dsa_reaccel` | DSA 再加速谱构造。 |
| `F` | 860 | `log_energy` | 对数电子谱能量积分。 |
| `F` | 868 | `reverse_transport_substeps` | 根据冷却步长和求解器类型选择反向输运子步数。 |
| `S` | 880 | `dg_sequence` | 根据低/注入/高能 break 构造 DG 网格序列。 |
| `F` | 942 | `log_interp` | 在对数 gamma 网格上插值正谱/冷却量。 |

### `src/Electron/electron_seed_history_kernel.f90`

2D shell 历史 photon field、Doppler 映射和路径吸收。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_seed_history_kernel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 31 | `ensure_cache` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 77 | `integrate_proper_time` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 96 | `advance_history_stream` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 143 | `add_mapped_cell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 175 | `build_map` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 208 | `log_interp` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 221 | `relative_doppler` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 234 | `build_transfer_cache` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 260 | `apply_path_tau` | 辐射 emissivity、seed、SSA/transfer 或吸收诊断。 |
| `S` | 308 | `locate_path_cell` | 在下游面网格中定位光线路径端点所在 chi 单元。 |
| `F` | 361 | `history_transfer_weight` | 把 optical-depth segment 转换为均匀源函数逃逸/传输权重。 |

### `src/Electron/electron_shell_transport_common.f90`

shell-level fullhide/flux-split 共享 primitive。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_shell_transport` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 21 | `resolve_solver` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 32 | `shell_coord_step` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 55 | `coord_to_dgamma` | 特征 Lorentz 因子、断点或谱峰诊断。 |

### `src/Electron/electron_transport_2d_kernel.f90`

finite-q 几何、q 方向对流/扩散和 2D 能量推进。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `electron_transport_2d` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 28 | `q_geometry` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 48 | `q_point` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 75 | `q_cell_geometry` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 101 | `downstream_grid` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 130 | `shock_state` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 147 | `q_divergence` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 172 | `q_step_limit` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `F` | 191 | `eta_hit_time` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 207 | `trace_faces` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 276 | `split_q_faces` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 297 | `q_inv_metric` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 320 | `q_diff_faces` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 348 | `build_q_adv` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 372 | `add_q_diff` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 401 | `build_q_base` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 427 | `solve_tridiagonal` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 453 | `advance_q_charint` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 491 | `advance_q_diffusion` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 523 | `advance_q_implicit` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 561 | `advance_pwncr_q` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 599 | `advance_energy_chi` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 630 | `advance_pwncr_energy` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 660 | `advance_stoch_chi` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
| `S` | 695 | `advance_energy_charint` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |

### `src/Electron/electron_transport_common.f90`

1D 保守重映射、PPM、特征线和 fullhide 基础 primitive。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `electron_transport_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 19 | `prepare_implicit_coeffs` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 34 | `backward_sweep` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 50 | `quad_interp3` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 64 | `ppm_interfaces` | 非均匀电子能量网格的 PPM 界面重构。 |
| `S` | 125 | `ppm_positive_cell` | PPM 单元正性限制；保持重构后的粒子数密度非负。 |
| `S` | 156 | `ppm_prefix` | 非均匀网格 PPM 前缀积分。 |
| `F` | 177 | `ppm_eval_prefix` | 任意位置的 PPM 保守前缀积分求值。 |
| `S` | 216 | `prepare_remap` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 236 | `prepare_exp_source` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 273 | `exp_source_int` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `F` | 293 | `exp_source_prefix` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 332 | `dnx_dgamma` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 349 | `u_edges` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `F` | 364 | `x_from_u` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 374 | `trace_affine_u` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 409 | `build_piece_u` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 452 | `find_u_cell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 479 | `trace_piece_edge` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 602 | `trace_piece_u` | 多滞后分段仿射 u 特征线回溯；单滞后调用传 `Num_lag=1`。 |
| `S` | 626 | `char_core` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 669 | `remap_edges` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 697 | `char_update` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 738 | `ppm_cell_int` | PPM 单元积分 primitive。 |
| `S` | 756 | `semi_lagrangian_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 789 | `flux_split_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 836 | `flux_split_nonuniform` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 890 | `flux_seq_nonuniform` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 943 | `fullhide_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `F` | 969 | `logparabola_peak` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 998 | `gamma_active_hi` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1033 | `chi_active_hi` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `F` | 1054 | `max_xi_chi` | finite-q shell 几何或 chi-equivalent 投影字段。 |

### `src/Electron/hybrid_special.f90`

电子注入、冷却、输运或同步辐射相关 Fortran 数值核。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `hybrid_special` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 10 | `dgamic` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 15 | `dbsk0e` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 20 | `dbsk1e` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 30 | `gamma_uic` | 特征 Lorentz 因子、断点或谱峰诊断。 |
| `F` | 40 | `bessel_k0` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 50 | `bessel_k1` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Electron/hybrid_spectrum.f90`

热-非热混合谱归一化和特殊函数加速。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 4 | `hybrid_spectrum` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 37 | `thermal_int1` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 141 | `thermal_int12` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 307 | `cpl_integral` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 326 | `solve_theta` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 349 | `theta_objective` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 360 | `get_initial_theta` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 446 | `newton_method` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 488 | `hybrid_coord` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 528 | `add_segment` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 552 | `spec_gamma` | 特征 Lorentz 因子、断点或谱峰诊断。 |
| `S` | 575 | `hybrid_thermal_coord` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 603 | `add_segment` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Hadronic/hadronic_accel.f90`

强子加速时间、外部冷却限制和 gamma limit / injection rate。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_accel` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 22 | `species_info` | 查询 hadronic species 的质量和电荷。 |
| `S` | 56 | `accel_time` | 费米加速时标。 |
| `S` | 78 | `syn_time` | 同步辐射冷却时标。 |
| `S` | 99 | `ext_time` | 外部光子场冷却时标。 |
| `S` | 114 | `inject_rate` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 152 | `gamma_limit` | 动力学、同步冷却和外部冷却限制给出的最大 Lorentz factor。 |
| `S` | 187 | `init_limit` | 未加入外部冷却前的动力学和同步辐射限制。 |
| `F` | 201 | `find_cross` | 查找 t_acc 与 t_ext 的交叉区间。 |
| `S` | 217 | `apply_xlimit` | 用对数插值给出外部冷却限制。 |
| `S` | 234 | `accel_calc` | 一次性计算加速、冷却、注入和最大能量限制。 |
| `F` | 256 | `trapz` | 梯形法则积分。 |

### `src/Hadronic/hadronic_base.f90`

强子网格、时间尺度、log-grid 校验和量子同步冷却公共工具。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_base` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 25 | `build_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 50 | `source_bounds` | 粒子源项或注入谱归一化；必须同时满足粒子数和能量预算。 |
| `S` | 75 | `build_edges` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `F` | 98 | `dyn_time` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 108 | `proton_limit` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 122 | `check_grid` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `F` | 147 | `quant_factor` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |

### `src/Hadronic/hadronic_bh.f90`

Bethe-Heitler 质子损失、pair source 和 photon loss kernel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_bh` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 19 | `bh_calc` | Bethe-Heitler 单电荷 pair source、fractional proton loss 与 photon sink 算子。 |
| `F` | 77 | `loss_point` | 单个 photon bin 对 proton loss 的贡献。 |
| `F` | 86 | `bh_pair` | Bethe-Heitler pair 产生核；进入 pair source 与 photon sink。 |
| `F` | 103 | `bh_outer` | Bethe-Heitler 外层 omega 积分核。 |
| `F` | 124 | `bh_inner` | Bethe-Heitler 内层 ebar 积分核。 |
| `F` | 141 | `bh_sigma` | Blumenthal 1970 微分截面。 |
| `F` | 181 | `proton_loss` | Bethe-Heitler 质子能量损失核。 |
| `F` | 189 | `bh_phi` | Bethe-Heitler 能量损失 phi(x) 近似。 |
| `F` | 211 | `bh_rk3` | 三参数 RK 3/8 积分器。 |
| `F` | 213 | `func` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 237 | `bh_rk4` | 四参数 RK 3/8 积分器。 |
| `F` | 239 | `func` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Hadronic/hadronic_cascade.f90`

gamma-gamma pair/synch shell-sequence cascade。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 3 | `hadronic_cascade` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 17 | `cascade_step` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `S` | 54 | `validate_step` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 62 | `produce_pairs` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 72 | `build_egamma` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
| `S` | 77 | `cool_pairs` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 94 | `emit_syn` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 115 | `cascade_seq` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 189 | `validate_seq` | 系统边界校验；用于拒绝外部输入或正式 kernel contract 违反。 |
| `S` | 200 | `build_grids` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 212 | `shell_geom` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 227 | `pair_syn` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 248 | `cool_deposit` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 269 | `electron_loss` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `S` | 282 | `advance_energy` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `F` | 305 | `gamma_bin` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 322 | `advance_photon` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |

### `src/Hadronic/hadronic_decay.f90`

pi0、pi±、mu± decay、neutrino 和 electron channel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_decay` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 19 | `pi0_gamma` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 45 | `pion_decay` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 88 | `build_pion` | pion decay source 的 log-energy rate 预处理。 |
| `S` | 98 | `pion_madd` | charged-pion decay 到 helicity-resolved muon。 |
| `S` | 122 | `pion_nadd` | charged-pion prompt neutrino source。 |
| `S` | 144 | `muon_decay` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 185 | `build_muon` | muon decay source 的 log-energy rate 预处理。 |
| `S` | 195 | `muon_nadd` | muon decay 到 flavor-resolved neutrino。 |
| `S` | 218 | `muon_eadd` | muon decay 到 e± source。 |
| `S` | 238 | `decay_hummer` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 274 | `log_spacing` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 289 | `nu1_decay` | muon decay neutrino 谱函数 f_nu1(x,h)。 |
| `F` | 302 | `nu2_decay` | muon decay neutrino 谱函数 f_nu2(x,h)。 |
| `F` | 315 | `log_interp` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

### `src/Hadronic/hadronic_formal.f90`

formal 1D 强子壳层序列实现层；Python/f2py 只通过 `formal_transport_1d` 进入这里。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `hadronic_formal` | formal 1D 强子 shell-sequence 模块命名空间。 |
| `S` | 10 | `formal_transport` | 按半径推进 proton injection/transport、pγ/BH/pp、secondary species、secondary radiation、photon survival 和 secondary e± source。 |

### `src/Hadronic/hadronic_forward_1d.f90`

正向激波 formal 1D 强子 f2py 入口和 shell-sequence driver。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 3 | `had_syn_pol` | f2py/Python 调用边界；同步偏振诊断使用的 proton synchrotron polarization wrapper。 |
| `S` | 17 | `pg_operator` | f2py/Python 调用边界；Hummer 2010 p-gamma operator。 |
| `S` | 41 | `pair_production` | f2py/Python 调用边界；gamma-gamma pair source/loss operator。 |
| `S` | 61 | `pp_shell` | f2py/Python 调用边界；pp delta source/loss operator。 |
| `S` | 85 | `bethe_heitler` | f2py/Python 调用边界；Bethe-Heitler pair source 和 photon/proton loss。 |
| `S` | 103 | `hadronic_ic` | f2py/Python 调用边界；pion/muon/proton hadronic IC diagnostic operator。 |
| `S` | 129 | `decay_operator` | f2py/Python 调用边界；pion/muon decay 到 gamma、e± 和 neutrino。 |
| `S` | 164 | `cascade_sequence` | f2py/Python 调用边界；shell-sequence gamma-gamma pair/synch cascade。 |
| `S` | 193 | `hadronic_1d` | light 1D proton transport + proton synchrotron shell-sequence driver。 |
| `S` | 279 | `init_grid` | `hadronic_1d` 内部 proton grid 初始化；Jacobians 不能省略。 |
| `S` | 306 | `init_out` | `hadronic_1d` 输出频率和数组初始化。 |
| `S` | 326 | `inject_p` | `hadronic_1d` 壳层质子源项；必须同时满足粒子数和能量预算。 |
| `S` | 336 | `advance_p` | `hadronic_1d` 壳层质子输运推进。 |
| `S` | 345 | `advance_sec` | `hadronic_1d` Hummer p-gamma secondary chain 推进。 |
| `S` | 361 | `emit_syn` | `hadronic_1d` proton synchrotron emissivity 和 seed 计算。 |
| `S` | 374 | `formal_transport_1d` | formal 1D 强子主入口；推进 proton transport、p-gamma/BH/pp、secondary、photon survival 和 secondary e± source。 |

### `src/Hadronic/hadronic_hummer.f90`

旧/诊断型 Hummer p-gamma shell 聚合 helper。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_hummer` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 17 | `hummer_shell` | 旧 Hummer p-gamma shell path，估计 photon survival 并推进 secondary species。 |
| `S` | 117 | `map_secondaries` | 将 p-gamma secondary source 映射到 species gamma 网格。 |
| `S` | 153 | `neutron_loss` | neutron p-gamma loss 插值和显式扣除。 |
| `S` | 172 | `pos_interp` | 正值 log-log 插值 helper。 |

### `src/Hadronic/hadronic_ic.f90`

强子 secondary inverse-Compton kernel。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_ic` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 22 | `ic_init` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 47 | `ic_operator` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `S` | 75 | `apply_ic` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 105 | `apply_species` | hadronic secondary/decay/pp 过程；输出 gamma、e± 或 neutrino 源项。 |
| `S` | 119 | `build_ic` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 145 | `ic_channel` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |
| `F` | 175 | `ic_coeff` | inverse-Compton/KN 相关核；joint 模式要求 cooling 和 emissivity 使用同一 photon field。 |

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

### `src/Hadronic/hadronic_pg.f90`

Hummer 2010 p-gamma response、光子汇和 secondary family deposition。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_pg` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 45 | `pg_hummer` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 122 | `accum_phloss` | 冷却或能量损失计算；必须和 emissivity/source 单位保持一致。 |
| `S` | 153 | `deposit_pions` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 174 | `deposit_baryons` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 193 | `deposit_shift` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 210 | `rates_res` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 231 | `rates_dir` | p-gamma Hummer response 或 secondary family deposition。 |
| `S` | 250 | `rates_mul` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 269 | `loss_res` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 280 | `loss_dir` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 291 | `loss_mul` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 304 | `phloss_res` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 313 | `phloss_dir` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 321 | `phloss_mul` | 冷却/损失算子；自然时间率进入 R 坐标前必须乘 dtprime/dR。 |
| `F` | 330 | `kernel_res` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 343 | `kernel_mul` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 358 | `kernel_dir` | p-gamma Hummer response 或 secondary family deposition。 |
| `F` | 375 | `idir` | p-gamma Hummer response 或 secondary family deposition。 |

### `src/Hadronic/hadronic_pp.f90`

Delta 近似 pp 算子和 secondary source。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_pp` | 模块命名空间；集中声明本文件共享 procedure。 |
| `S` | 20 | `pp_source` | pp delta kernel；输出随 `pden` 线性的 secondary 源和与其归一化无关的单粒子 proton loss。 |
| `S` | 51 | `validate_inputs` | pp delta 输入网格和 target density 检查。 |
| `S` | 62 | `set_options` | pp delta 能量份额和 charged/neutral pion 分支设置。 |
| `S` | 87 | `emit_secondaries` | 按 delta 动能映射输出三类 secondary source。 |
| `S` | 96 | `pp_sigma` | Kelner+2006 总能量 pp 非弹性截面。 |
| `F` | 121 | `pp_threshold` | pp 反应阈值总能量。 |
| `S` | 127 | `secondary_source` | delta 近似 secondary source。 |
| `S` | 143 | `pos_interp` | 正值 log-log 插值 helper。 |
| `F` | 180 | `upper_bracket` | 单调数组 bracket 查找。 |
| `F` | 210 | `log_eval` | 单次 log-log 线性插值。 |

### `src/Hadronic/hadronic_rad.f90`

质子同步辐射和偏振。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_rad` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 17 | `syn_xarg` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 33 | `syn_mass` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 69 | `syn_kernel` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 81 | `proton_syn` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 122 | `syn_polarization` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 147 | `polar_integral` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |

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
| `S` | 23 | `secondary_calc` | pion/muon secondary radiation 顶层算子。 |
| `S` | 44 | `apply_rad` | 应用已预计算的同步辐射和 IC 映射。 |
| `S` | 62 | `init_syn` | 建立 pion/muon 同步辐射核矩阵。 |
| `S` | 85 | `init_ic` | 建立 pion/muon IC 能量偏移和上界索引。 |
| `S` | 106 | `pion_syn` | π 介子同步辐射冷却率。 |
| `S` | 122 | `muon_syn` | μ 子同步辐射冷却率。 |
| `S` | 140 | `pion_ic` | π 介子 inverse-Compton 冷却率。 |
| `S` | 160 | `muon_ic` | μ 子 inverse-Compton 冷却率。 |
| `F` | 182 | `syn_kernel` | AM3 分段超相对论同步辐射核。 |
| `S` | 214 | `build_ic` | 为一个粒子种类建立 IC 偏移和上界索引。 |
| `S` | 237 | `ic_channel` | 单个 IC 通道卷积。 |
| `F` | 267 | `ic_coeff` | Thomson 截面和质量缩放的 IC 系数。 |
| `S` | 277 | `check_density` | secondary radiation 输入密度有限性检查。 |

### `src/Hadronic/hadronic_shell.f90`

formal 1D 强子底层 shell primitive 与单位/投影 helper；f2py wrapper 已收窄到运行时真正需要的入口。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 1 | `hadronic_shell` | shell-level 强子 primitive 模块命名空间。 |
| `S` | 23 | `pp_delta` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 48 | `hic_shell` | proton/pion/muon hadronic IC operator；输出 IC emissivity 与投影系数。 |
| `S` | 81 | `hic_project` | hadronic IC shell emissivity 投影到 photon grid。 |
| `S` | 110 | `species_step` | n、pi、mu secondary species 同壳层保守推进。 |
| `S` | 176 | `scan_pmax` | 沿半径序列估计全局 proton 最大 Lorentz factor。 |
| `S` | 202 | `secondary_rad` | pion/muon synchrotron 与 IC shell emissivity。 |
| `S` | 231 | `secondary_project` | secondary radiation 从 hadron grid 投影到 photon grid。 |
| `S` | 294 | `loss_rates` | adiabatic、synchrotron 和 quantum-synch 连续损失率。 |
| `S` | 317 | `electron_seq` | secondary e± source 随壳层序列组装。 |
| `S` | 361 | `photon_loss` | photon loss rate 到 optical-depth/survival closure。 |
| `S` | 387 | `effective_time` | interaction loss rate 的有效时间积分。 |
| `S` | 414 | `pgamma_update` | pγ loss/re-injection 对 proton spectrum 的壳层更新。 |
| `S` | 434 | `proton_step` | proton injection、continuous loss 和 pγ update 的单壳层推进。 |
| `S` | 455 | `exp_sink` | 指数 sink primitive；用于已定义 interaction/loss closure。 |
| `S` | 472 | `rate_lum` | rate-per-energy 到 luminosity-per-frequency/energy 的壳层换算。 |
| `S` | 487 | `project_lum` | source energy grid 到目标 photon grid 的 luminosity 投影。 |
| `S` | 506 | `project_hic` | hadronic IC 多 species emissivity 投影。 |
| `S` | 526 | `pair_content` | pp source 与两个 BH 轻子电荷组合成电子方程使用的 rate source。 |
| `S` | 543 | `shell_density` | shell content 到 density-per-GeV 的单位变换。 |
| `S` | 557 | `gamma_edges` | Lorentz-factor grid edges；守恒积分需要的 bin geometry。 |
| `S` | 579 | `proc_power` | secondary process power diagnostic。 |
| `F` | 608 | `trapz` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 622 | `pos_interp` | 正值 log-log grid projection；当前 secondary-feedback Python glue 仍使用。 |
| `S` | 655 | `gamma_source` | source-per-energy 到 source-per-gamma 的 grid projection。 |
| `S` | 672 | `dist_gev` | distribution-per-gamma 到 distribution-per-GeV 的 grid projection。 |
| `S` | 689 | `align_photon` | hadron/photon grid 对齐 helper。 |
| `S` | 709 | `shell_geom` | 对相邻动力学状态的共动 proper time 作梯形积分；首点为零时长初态。 |
| `S` | 732 | `shell_volumes` | 壳层体积序列。 |

### `src/Hadronic/hadronic_species.f90`

neutron、pion、muon explicit species transport。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `hadronic_species` | 模块命名空间；集中声明本文件共享 procedure。 |
| `F` | 25 | `div_rate` | 球对称膨胀发散率，用于绝热冷却。 |
| `S` | 39 | `synch_loss` | charged species 同步冷却率。 |
| `S` | 70 | `ad_loss` | species 绝热冷却率。 |
| `S` | 89 | `species_advance` | 七分量 neutron/pi/mu 输运主推进。 |
| `S` | 126 | `validate_inputs` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 144 | `advance_muons` | 输运推进或离散更新 helper；改动需验证守恒量和谱形。 |
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
| `S` | 11 | `remap_loggamma` | 推进 primitive；把源项、通量或事件分裂推进到下一半径/时间步。 |
| `F` | 41 | `remap_target` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

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

### `src/Interpolation/SED_interpolation.f90`

shell-level 和 chi-resolved EATS/Doppler 投影。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 11 | `sed_interpolation` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 113 | `sed_adaptive_theta` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 180 | `integrate_theta_cell` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 210 | `project_theta_sample` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 245 | `projshellseg` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 268 | `sed_interpolation_chi` | top-hat chi-resolved FS synchrotron+SSA lightcurve projection。 |
| `S` | 384 | `project_chi` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 398 | `chi_state` | finite-q shell 几何或 chi-equivalent 投影字段。 |
| `S` | 425 | `sed_chi_electron` | direct-electron chi projection；先批量计算 chi-local synchrotron/SSA 再投影。 |
| `S` | 465 | `sed_chi_ring` | axisymmetric structured chi ring projection；输入为预计算 `F_ring/Tau_ring`。 |
| `S` | 485 | `sed_chiring_batchlum` | 观测者投影或插值 helper；只读局域状态，不回写 transport。 |
| `S` | 519 | `sed_chiring_batchlum_ray` | 观测者投影或插值 helper；只读局域状态，不回写 transport。 |
| `S` | 595 | `init_leaves` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 618 | `init_geom` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 641 | `init_tbase` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 649 | `init_tgrid` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 665 | `leaves_flux` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 680 | `sample_phi` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 725 | `sample_event` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 786 | `build_hits` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 849 | `append_hit` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 872 | `fill_radiation` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 911 | `interpescidx` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 928 | `local_escape` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 941 | `radial_range` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 960 | `pathfactor` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 980 | `cut_slab` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1004 | `outer_face` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1015 | `inner_face` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1028 | `interplogidx` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1044 | `ray_flux` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `F` | 1081 | `hitpatch` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1109 | `sed_chiring_core` | 观测者投影或插值 helper；只读局域状态，不回写 transport。 |
| `S` | 1155 | `projphisamp` | 局部 helper；语义由所在文件的算法阶段决定。 |
| `S` | 1205 | `project_ring` | structured ring EATS/Doppler segment projection。 |
| `S` | 1230 | `ring_source` | structured ring chi cell source + SSA escape 累加。 |
| `F` | 1255 | `chi_escape` | finite-q shell SSA escape factor。 |
| `F` | 1269 | `lowerreal8` | 局部 helper；语义由所在文件的算法阶段决定。 |

### `src/Interpolation/SED_interpolation_structured.f90`

结构化喷流 theta/theta-phi shell-level 投影。该文件仍由 `structured_jet_1d` Fortran 聚合入口调用；旧 Python binding 不再从 `src.Interpolation` 导出。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `S` | 9 | `sed_interpolation_structured` | f2py/Python 调用边界或主聚合入口；稳定性高于内部 helper，改动需同步 wrapper、构建和冒烟测试。 |
| `S` | 106 | `project_structured_segment` | 结构化喷流 patch/角向投影调度。 |
| `S` | 130 | `sed_structured_phi` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 212 | `project_phi` | 结构化喷流 patch/角向投影调度。 |

### `src/Interpolation/interpolation_common.f90`

SED 插值共享累加 primitive。

| Kind | Line | Program unit | 算法/物理责任 |
| --- | ---: | --- | --- |
| `M` | 2 | `interpolation_common` | 公共模块；提供多个入口复用的物理/数值 primitive。 |
| `S` | 13 | `accum_logsed` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |
| `S` | 52 | `accum_shifted` | 网格、坐标变换、插值或保守重映射 primitive；Jacobians 不能省略。 |

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
| `S` | 16 | `sampled_weights` | 实际节点 Simpson 权重；偶数样点用 Cartwright 末区间校正。 |
| `F` | 35 | `rad_interp` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `S` | 63 | `transfer_factor` | 辐射 emissivity、seed、SSA/transfer 或偏振计算。 |
| `F` | 74 | `syn_kernel` | 近似同步辐射核的低频、指数段和高频分支。 |
| `S` | 91 | `pair_grid` | gamma-gamma pair-production 能量网格和测度准备。 |
| `F` | 111 | `pair_sigma` | gamma-gamma pair-production 截面。 |
| `S` | 130 | `pair_tau` | 对头碰撞近似下的 gamma-gamma 光深。 |
| `S` | 163 | `syn_seed_chi` | chi-resolved synchrotron/SSA batch kernel；单列调用传 `Num_chi=1`。 |
| `S` | 284 | `syn_flux_chi` | chi-resolved projection-only synchrotron flux/tau batch kernel。 |

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
| `S` | 295 | `build_state` | 非均匀网格重构和尾矩预处理。 |
| `S` | 309 | `build_shell` | 单 shell 非均匀电子谱重构。 |
| `S` | 327 | `limit_slope` | minmod 斜率限制。 |
| `S` | 355 | `add_tail` | 非均匀电子谱 gamma 矩尾积分累加。 |
| `F` | 368 | `gauss_point` | 二点 Gauss-Legendre 节点。 |
| `S` | 380 | `nonuniform_point` | 单个 shell/frequency 的非均匀网格 SSC 累加。 |
| `F` | 421 | `first_cell` | 非均匀电子 cell 下界二分搜索。 |
| `F` | 452 | `ssc_minmod` | minmod 斜率限制器。 |
| `F` | 465 | `profile_value` | 线性重构剖面值。 |
| `F` | 476 | `gamma_moment` | 线性重构剖面的 gamma 矩解析积分。 |
| `F` | 496 | `low_cell` | 低 seed 侧单电子 cell SSC 积分。 |
| `F` | 535 | `low_integral` | 低 seed 侧完整电子网格 SSC 积分。 |
| `F` | 560 | `high_tail` | 高 seed 侧 SSC 尾积分。 |

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
| `S` | 62 | `run_axisymmetric` | 结构化喷流 patch/角向投影调度。 |
| `S` | 99 | `run_nonaxisymmetric` | 结构化喷流 patch/角向投影调度。 |
| `S` | 138 | `structured_solve_axisymmetric` | 结构化喷流 patch/角向投影调度。 |
| `S` | 204 | `register_axis_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 221 | `solve_axis_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 241 | `structured_solve_nonaxisymmetric` | 结构化喷流 patch/角向投影调度。 |
| `S` | 314 | `register_phi_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 334 | `solve_phi_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 352 | `copy_phi_patch` | 结构化喷流 patch/角向投影调度。 |
| `S` | 366 | `structured_solve_element` | 结构化喷流 patch/角向投影调度。 |
| `S` | 534 | `apply_absorption` | 光深、gamma-gamma absorption、pair injection 或 photon survival 相关算子。 |
