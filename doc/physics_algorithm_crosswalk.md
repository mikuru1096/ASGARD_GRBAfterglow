# 物理算法交叉指南

本文把“物理问题 -> 离散变量 -> 程序单元/调用路径 -> 验收指纹”放在同一页。它不是 public API 说明；public API 见 `doc/public_api.md`，完整程序单元表见 `doc/fortran_kernel_index.md`。

ASGARD 的第一性原理主线是：先在局域壳层坐标中闭合动力学、粒子输运、光子场和强子源汇，再做等到达时间面投影。观测光变不是状态变量，不能用后处理平滑、经验 time shift 或 fallback 修补局域物理错误。

## 1. 共同坐标和守恒量

| 对象 | 离散变量 | 核心公式 | 程序单元/调用路径 | 验收指纹 |
| --- | --- | --- | --- | --- |
| 半径主坐标 | `R(i)` | `dtprime/dR = 1/(beta Gamma c)` | `Dynamics_forward`, `Dynamics_reverse`, `shell_geom` | 所有自然时间率进入电子/强子方程前先换算到 `R` 坐标。 |
| log 能量谱 | `dN/dln(gamma)` | `N_x = gamma dN/dgamma` | electron `log_edges`/transport，hadronic `electron_seq` | 粒子数积分用 `sum N_x dx`，不能丢 Jacobian。 |
| 2D 厚壳 | `q_grid/q_face/dq` | `chi_BM=(1-q)^[-(4-k)/(3-k)]` | `compute_q_geometry`, `q_geometry_point`, `sed_interpolation_chi` | `chi_grid` 只作 BM-equivalent 诊断，投影读 `chi_radius_cm`、`chi_gamma_bulk`、`chi_dvolume_weight`。 |
| photon field | `Seed_syn(nu,R)` | `n_nu = u_nu/(h nu)` | `syn_seed_chi`, Python direct `E=h\nu` conversion | Observer luminosity 不能直接当 target density，必须经 shell volume/escape time/单位转换。 |

## 2. 物理过程到程序单元

| 物理过程 | 算法方程或决策 | 主要 Fortran 入口 | 关键内部程序单元 | 必须检查 |
| --- | --- | --- | --- | --- |
| 正向激波动力学 | RK4 推进 `Gamma, M_sw, R`，支持 ISM/wind、密度跳变和注入。 | `dynamics_forward` | `forward_rhs`, `density_profile`, `forward_rk4` | 无跳变/无注入时 `Gamma(R)` 应接近 BM 标度并保持平滑。 |
| 反向激波动力学 | crossing 前后分支，`gamma34` 注入能标，显式 `U3/V3`。 | `dynamics_reverse` | `reverse_dynamics_rhs`, `advance_pre_m3_to_target`, `compute_region3_field`, `compute_region3_thermal_state` | `M3/Mej` crossing 端点连续；`B3, gamma34, U3/V3` 无孤立尖峰。 |
| 磁化 RS jump | 有限强度 MHD jump，`E_iso/[(1+sigma) Gamma0 c^2]` baryonic mass。 | `Dynamics_reverse` + `reverse_rhs` + `reverse_shock_mhd_jump` | `rs_mhd_state` 一次返回四输出，主流程内处理有序场和磁压焓惯性 | `sigma -> 0` 回到 unmagnetized baseline，不能用 `4 gamma43 + 3` 替代有限强度极限。 |
| 多密度增强 multiple RS | 扫描 jump window，记录可触发的新 shock branch。 | `dynamics_reverse` | `update_multiple_reverse_events`, `multiple_reverse_event_source`, `store_multiple_branch_state` | 只在密度增强和 pressure 条件允许时出现，分支权重随半径平滑。 |
| 电子注入 | shock electron count 与能量预算归一化截断幂律/热分布。 | electron entry files | `electron_gamma_m_exact`, `electron_injection_prefactor`, `source_edges`, `electron_build_thermal_shape_dnx` | `sum Q dx` 和能量矩同时闭合；`p≈2` 不用渐近公式替代精确归一化。 |
| 1D fullhide 电子输运 | 隐式有限体积解 `partial_R N_x + partial_x(v_x N_x)=Q_x`。 | `fs_fullhide_1d` | `prepare_fullhide_shell`, `flux_seq_nonuniform`, `forward_cooling` | 冷却谱不出现网格锯齿；`nu_a` 来自同一次 `Tau_syn` 网格。 |
| joint coupled electron pass | 电子冷却 seed 和 secondary `Q_e,R` 进入同一 fullhide solve。 | `fs_fullhide_coupled` | `prepare_coupled_shell`, `electron_cooling_ic_loss_emissivity_budget` | 只在 `fullhide_1d + index_y=1 + am3_1d` 合同内运行；失败应暴露 grid contract。 |
| DG 电子输运 | LGL 元、正性核和闭合低能边界。 | `fs_dg_1d` | `dg_advance_step`, `dg_filter_positive`, `dg_low_content` | DG 不是 smoothing；谱形、守恒量和断点位置共同验收。 |
| WENO/T2G1/SLC/charint | 方法对照和数值诊断。 | `fs_weno5_1d`, `fs_t2g1_1d`, `fs_slc1_1d`, `fs_charint_1d` | `compute_weno_fluxes`, `advance_t2g1_substep`, `char_update` | 用于判别输运误差，不作为物理 fallback。 |
| 2D q-shell 电子输运 | `x` 能量方向 + `q` 厚壳方向的对流/扩散/源项。 | `fs_transport_2d` | `compute_q_cell_geometry`, `advance_q_implicit`, `advance_energy_chi`, `project_q_projection_shell` | 有限半径、非零体积权重和时间连续 lightcurve；不宣称 chi-local hadronic。 |
| 同步辐射与 SSA | 从电子谱计算 `P_syn`, `Seed_syn`, `Tau_syn`, `nu_a`。 | `syn_state` | `electron_syn_gauss`, `electron_tau_gauss`, `nua_fromtau` | `P_syn/Seed/Tau` 同源；避免重复 root search 造成 runtime 和语义分裂。 |
| 同步偏振 | 频率相关 F/G kernel 到 Stokes emissivity。 | `synchrotron_polarized_components` | `synchrotron_fg_kernel`, `synchrotron_fg_integral` | 只覆盖同步分支；峰时偏差优先查 dynamics/projection。 |
| SSC 与 IC cooling | 同一 photon seed 决定 cooling 与 emissivity。 | `ssc_spec`, `ssc_spec_nonuniform`, `electron_cooling_ic_loss` | `prepare_kn_tables`, `accumulate_uniform_point`, `accumulate_nonuniform_point` | joint 预算中不能只改 cooling 不改 photon source。 |
| gamma-gamma absorption | 目标 photon field 上积分 pair cross-section。 | `annihilation` | `pair_grid`, `sigma_kernel`, `set_window` | 输出是 observer attenuation；不要把它误作局域 cascade 全闭合。 |
| pair synch cascade | shell-sequence gamma-gamma pair/synch branch。 | `cascade_sequence` | `cascade_seq`, `produce_pairs`, `emit_syn` | `pair_cascade_iterations>1` 仍不是 IC-mediated electromagnetic cascade。 |
| proton injection/transport | log-gamma 质子源项、连续损失和保守推进。 | `formal_transport_1d` | `proton_inject`, `advance_loggamma`, `proton_step` | 质子数/能量源项、loss 和 photon survival 同壳层对齐。 |
| p-gamma Hummer response | proton loss、re-injection、photon loss、secondary families。 | `formal_pgamma_operator` | `pg_hummer`, `deposit_pions`, `pgamma_update` | `pgamma_scheme="hummer_2010_response"` 才是正式反馈路径。 |
| Bethe-Heitler | fractional proton loss 是单电荷 pair source 的双电荷能量矩，并共享 photon sink。 | `bethe_heitler` | `bh_calc`, `bh_pair`, `bh_quad3`, `bh_quad4` | \(\omega\le600\) 与 GL12 由 pair 核唯一拥有；当前 electron grid 内逐 proton 能格闭合，完整支撑仍由 grid owner 保证。 |
| pp 过程 | baryon target density 上的单粒子连续 loss 与分布 secondary source。 | `pp_shell` | `pp_source`, `secondary_source` | 截面用总能量，loss 与 delta 能量份额用动能；loss 与 `pden` 无关，secondary source 与其线性。 |
| secondary species/decay | n、pi、mu 输运与 decay/neutrino/electron channel。 | `formal_transport_1d`, `decay_operator` | `species_advance`, `pion_decay`, `muon_decay` | neutrino 逃逸不反馈；e± source 只有 formal 输出才能进入电子方程。 |
| hadronic secondary radiation | pion/muon synchrotron 与 IC。 | `formal_transport_1d` | `secondary_rad`, `secondary_calc`, `ic_channel` | 分量可输出/诊断，不能自动回灌 photon field。 |
| RS hadronic light path | RS proton transport + proton synchrotron。 | `reverse_hadronic_1d` | `advance_shell`, `emit_syn` | full-chain RS hadronic 复用 formal 1D kernel；不是新 2D RS transport。 |
| EATS 投影 | Doppler、redshift、SSA survival 到 `Fnu(t_obs)`。 | `sed_interpolation`, `sed_adaptive_theta`, `sed_interpolation_chi` | `project_shell_segment`, `integrate_theta_cell`, `accumulate_chi_cell_source` | projection 只读本地 emissivity，不修正上游谱/动力学。 |
| 结构化喷流 | axisymmetric 或 theta-phi patch 调度。 | `jet_flux_1d` | `run_axisymmetric`, `run_nonaxisymmetric`, `structured_solve_element`, `structured_apply_observer_absorption` | patch 权重/角向投影可变，微物理 kernel 语义不变。 |

## 3. 不可交换的算法顺序

默认 separated afterglow 主线：

```text
dynamics
-> electron transport
-> photon field build
-> hadronic formal transport
-> reverse-shock emission
-> gamma-gamma absorption / pair branch
-> observer projection
```

joint feedback 主线：

```text
dynamics
-> primary electron
-> photon field
-> fixed shell-level feedback iterations:
       formal hadronic solve
       photon survival/source update
       secondary e± source assembly
       coupled electron solve
       photon field rebuild
-> observer projection
```

这些顺序不能随意交换。例子：先把 BH e± 后处理并入 `dN_e`，再拿合并后的谱解释同一壳层 IC cooling，会破坏局域能量预算；这也是 joint 路径单独存在的原因。

## 4. 从异常结果反查源码

| 现象 | 首先检查 | 不应采取的做法 |
| --- | --- | --- |
| `Gamma(R)` 或 `R_Tobs` 不平滑 | `Dynamics_forward` / `Dynamics_reverse` 右端项、density jump、event splitting。 | 在 SED 投影层乘经验时间因子。 |
| RS 光变 crossing 附近尖峰 | `advance_pre_m3_to_target`、`reverse_dynamics_rhs`、`electron_reverse_evolve`。 | 用 smoothing 或 clipping 处理光变。 |
| `nu_a` 跳格 | `syn_state` 与 `nua_fromtau` 是否共用 `Tau_syn`。 | 对 `nu_a` 做后处理插值平滑。 |
| 2D chi lightcurve 断崖 | `chi_radius_cm`、`chi_gamma_bulk`、`chi_dvolume_weight`、SSA path attenuation。 | 裁剪负半径或用 1D 曲线替代 2D 投影。 |
| joint/RS full-chain hadronic 报 log-grid contract | formal hadronic electron-energy grid、`check_grid` 调用链。 | 放宽断言或临时重采样到“看起来”均匀。 |
| 偏振峰时偏早 | dynamics/jet evolution benchmark、EATS 几何。 | 在 polarization projection 加经验 time shift。 |
| 强子谱能量预算不闭合 | proton source、p-gamma/BH/pp loss、secondary source、photon survival 是否同壳层同单位。 | 只调 observer component 归一化。 |

## 5. 文档和源码同步规则

- 新增 Fortran entry 时，同步 `build_extensions.py`、`doc/fortran_kernel_index.md`、`doc/call_chain.md` 或专题算法页。
- 修改物理方程时，同步 `doc/project_physics_design.md` 或 `doc/physical_processes.md`，并在本页更新物理过程表。
- 修改数组形状、单位或坐标时，同步 `doc/algorithm_workflow.md`、`doc/electron_solver_algorithms.md` 或 hadronic 专题页。
- 修改后端支持边界时，同步 `doc/public_backend_limits.md`，不要把未实现路径写成可用能力。
