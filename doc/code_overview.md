# ASGARD Code Overview

本文档按当前工作树整理，目标是让开发者快速定位运行主链、接口边界和仍未闭合的物理/数值范围。更细的电子算法细节见 `doc/electron_solver_algorithms.md`，调用图见 `doc/call_chain.md`。

## 1. Public API

外部入口集中在 `ASGARD/`。

- `ASGARD/__init__.py`
  - 导出稳定用户对象：`Model`、`Medium`、`ISM`、`Wind`、各类 `Jet`、`Observer`、`Radiation`、`Setups`、`ObsSet`、`Fitter`、`Param`、`observe`、`run_fit`、`units`。
- `ASGARD/api_model.py`
  - 定义面向用户的模型对象和结果对象。
  - 主要输入对象：
    - `Medium` / `ISM` / `Wind`
    - `TophatJet` / `GaussianJet` / `PowerLawJet` / `TwoComponentJet` / `StepPowerLawJet` / `Ejecta`
    - `Observer`
    - `Radiation`
    - `Setups`
    - `Model`
  - `Model` 的常用方法：
    - `flux_density_grid(times_s, nu_hz)`：返回频率-时间二维网格上的 `FluxResult`。
    - `flux_density(times_s, nu_hz)`：支持一一对应的时间/频率采样。
    - `flux_density_exposures(...)`：按曝光窗口采样并返回曝光平均通量密度。
    - `spectrum(time_s, nu_hz)`：返回单时刻 SED。
    - `flux(time_s, nu_min_hz, nu_max_hz, num_points=64)`：对频段积分。
    - `sky_image(t_obs, nu_obs, fov, npixel=128)`：观测者平面成像。
    - `details()`：返回前向/反向激波轨迹、特征频率、粒子谱、hadronic 诊断等。
- `ASGARD/api_observe.py`
  - `observe(model, config=None, spectrum_output=None)`：高层观测入口。内部构造默认观测时间和多波段频率，调用 `Model.flux_density_grid`，再输出兼容旧拟合流程的 `FitResult`。
  - `run_fit(config=None)`：从 `FitConfig` 构造一个 `Model`，再调用 `observe(...)`。
- `ASGARD/api_fit.py`
  - `Fitter`、`Param`、`FitResult` 等拟合侧 API。

公共配置有两套：

- `ASGARD/api_model.py::Setups`：面向 `Model` 的运行设置。
- `asgard_core/asgard_config.py::FitConfig`：旧拟合/脚本路径仍在使用；`HadronicConfig`、`ReverseShockConfig`、`SpectrumOutputConfig` 也在这里。

当前 electron solver 名称：

- 正式名称：`fullhide_1d`、`slc1_1d`、`charint_1d`、`charint_2d`、`t2g1_1d`、`weno5_1d`、`fullhide_2d`
- 兼容别名：`fullhide`、`slc1`、`charint`、`t2g1`、`weno5`
- `*_1d` 强制 `num_chi = 1`；`*_2d` 在未显式指定时使用 `num_chi = 64`。

## 2. 运行主链

当前有效主链是：

```text
Model / observe / run_fit
  -> FitConfig / SimulationSetup
  -> solve_state_from_setup
  -> solve_dynamics
  -> solve_electron
  -> photon_field_stage
  -> solve_hadronic
  -> solve_reverse_shock_emission
  -> observer assembly
  -> Radiation.annihilation
  -> Interpolation.sed_interpolation
  -> API result
```

核心状态对象在 `asgard_core/asgard_types.py`：

- `DynamicsSolution`
  - `r_tobs`、`r_gamma`、`radius`、`swept_mass_g`，可附带 `ReverseShockDynamics`。
- `ElectronSolution`
  - `gam_e`、`d_n_gam_e`、`l_syn_spec`、`seed_syn`、`nu_m`、`nu_c`、`nu_a`。
  - 2D solver 额外携带 `d_n_gam_e_chi` 和 `chi_grid`。
  - BH 次级电子并回前向电子时，`d_n_gam_e_bh` 会被写入。
- `PhotonFieldState`
  - 严格区分前向同步 seed、供 hadronic 使用的目标场、供吸收使用的 seed 场。
- `HadronicSolution`
  - 1D hadronic proton/secondary/radiation 结果。
- `ObserverState`
  - 存放吸收因子、`tau_pair` 和投影前后的分量。
- `FluxComponents`
  - 统一装配 `total`、FS synch/SSC、hadronic、RS synch/SSC、cross-zone IC。

`asgard_core/asgard_state.py` 是运行链控制中心：

- `solve_state_from_setup(...)`
  - 依次调用 dynamics、electron、photon field、hadronic、reverse shock、observer assembly。
- `_build_photon_field_stage(...)`
  - 复制 electron `seed_syn`。
  - hadronic 启用且有 `epsilon_p` 时，可把前向 SSC seed 加入 hadronic target field。
  - 查询频率不改变内部 seed 频率网格，避免输出请求反向影响动力学/冷却演化。
- `_solve_hadronic_stage(...)`
  - 调用 `solve_hadronic(...)`。
  - 若有 BH 次级电子，则并回前向电子分布并重算 `l_syn_spec/seed_syn`。
  - 若有 `pγ` photon survival，则直接作用到 photon field。
- `_assemble_observer_stage(...)`
  - 组装 FS synch、FS SSC、RS synch/SSC、cross-zone IC、hadronic 分量。
  - hadronic 光子使用 electron Fortran kernel 给出的同步自吸收 transfer。
  - pair production 当前在 observer 侧生成 `tau_pair` 与 pair synchrotron 支路。
  - `Radiation.annihilation` 接收 pair-production `tau_extra`；当前 `pγ` photon loss 已先作为 local survival factor 进入 photon field 和 hadronic luminosity。
- `project_flux_grid(...)` / `project_spec(...)`
  - 把壳层谱投影到观测时间和观测频率。

拟合时的最短闭环：

```text
Fitter.loglike(...)
  -> asgard_core.asgard_fit.compile_problem(...)
  -> asgard_core.asgard_fit.eval_loglike(...)
  -> solve_state_from_setup(...)
  -> project_flux_grid(...)
  -> combine_multiband_flux(...)
  -> compute_light_curve_redchi(...)
```

## 3. Python 编排层

`asgard_core/` 中 Python 的职责是装配状态、选择后端、做观测投影和拟合，不应成为最终 AM3-derived 微观物理的归宿。

主要模块：

- `asgard_core/asgard_setup.py`
  - 从 `FitConfig` 构造 `SimulationSetup`，包括边界数组、seed frequency grid、observer time grid、 luminosity distance。
- `asgard_core/asgard_runtime.py`
  - 运行时后端选择和 Fortran 扩展调用。
  - 解析 electron solver、hadronic solver、`pgamma_scheme`。
  - 把 Fortran 返回数组封装成 `DynamicsSolution`、`ElectronSolution`、`HadronicSolution`。
- `asgard_core/asgard_state.py`
  - 主状态机和跨阶段耦合。
- `asgard_core/asgard_ssc.py`
  - 前向 SSC auxiliary grid 与 seed 计算。
- `asgard_core/asgard_coupling.py`
  - FS/RS cross-zone IC 的几何和 seed 场耦合。
- `asgard_core/asgard_postprocess.py`
  - 观测投影、band 聚合、拟合后处理。
- `asgard_core/asgard_fit.py`
  - 拟合问题编译和 likelihood 路径。
- `asgard_core/asgard_types.py`
  - 运行时 dataclass contract，开发时应优先看这里确认跨阶段字段。

hadronic Python 模块当前分成两类：

- Fortran wrapper / orchestration：
  - `hadronic_hummer.py`
  - `hadronic_bethe_heitler.py`
  - `hadronic_hadronic_ic.py`
  - `hadronic_pp.py`
  - `hadronic_pair_production.py`
  - `hadronic_species_transport.py`
  - `hadronic_secondary_radiation.py`
  - `hadronic_acceleration.py`
- 参考和 benchmark：
  - `hadronic_pgamma.py`
  - `hadronic_am3_solver.py`
  - `hadronic_am3_benchmark.py`
  - `hadronic_cascade.py`

这些模块可以负责输入检查、数组单位转换、benchmark glue 和过程编排。最终核心微观核应在 `src/Hadronic/*.f90`。

## 4. Fortran 内核层

Fortran 扩展由 `build_extensions.py` 统一登记并用 f2py 构建。Python 通过 `src/*/__init__.py` 暴露这些扩展。

### Dynamics

- `src/Dynamics/Dynamics_forward.f90`
  - 前向激波动力学。
  - 支持 ISM/wind、密度跳变、能量注入、不同 `index_dyn` 分支。
- `src/Dynamics/Dynamics_reverse.f90`
  - 反向激波动力学、RS 初始电子谱、RS 磁场。
- `src/Dynamics/dynamics_common.f90`
  - 共享动力学辅助量。

### Electron

- 主入口：
  - `src/Electron/FS_electron_fullhide_1d.f90`
  - `src/Electron/FS_electron_fullhide_2d.f90`
  - `src/Electron/FS_electron_charint_1d.f90`
  - `src/Electron/FS_electron_charint_2d.f90`
  - `src/Electron/FS_electron_slc1_1d.f90`
  - `src/Electron/FS_electron_t2g1_1d.f90`
  - `src/Electron/FS_electron_weno5_1d.f90`
- 共享内核：
  - `electron_common.f90`：log-gamma 网格、保守 remap、特征线、隐式步进、注入、特征频率辅助。
  - `electron_cooling_kernel.f90`：同步/IC/SSA 回热相关冷却组装。
  - `electron_radiation_kernel.f90`：同步谱和 seed。
  - `electron_seed_history_kernel.f90`：历史 photon field。
  - `electron_transport_2d_kernel.f90`：log-chi 几何、eta/log-chi 输运、2D energy advance。
  - `electron_injection_profiles.f90`：注入谱 profile。
  - `calling_modules.f90`：聚合导出层。

2D solver 当前在 `log10(gamma_e) x log10(chi)` 上演化下游结构，公开输出仍以 shell-reduced `P_syn/Seed_syn` 和 shell 级 `nu_m/nu_c/nu_a` 为主。

### Radiation / Interpolation

- `src/Radiation/SSC_spec.f90`
  - SSC 谱和 seed 计算；有 uniform 与 nonuniform 入口。
- `src/Radiation/Seed_reverse.f90`
  - 反向激波同步辐射和 seed。
- `src/Radiation/Annihilation.f90`
  - gamma-gamma 吸收；当前 pair-production optical depth 通过 `tau_extra` 输入。
- `src/Interpolation/SED_interpolation.f90`
  - 观测者系 EATS/Doppler 插值主路径。
- `src/Interpolation/SED_interpolation_structured.f90`
  - 结构化喷流插值入口。

### Hadronic

`src/Hadronic/FS_hadronic_1d.f90` 是 f2py 暴露入口，内部转调下列 Fortran kernel：

- `hadronic_transport_kernel.f90`
  - proton injection、adiabatic/synchrotron loss、log-gamma energy advance。
- `hadronic_radiation_kernel.f90`
  - proton synchrotron。
- `hadronic_interaction_kernel.f90`
  - Hummer 2010 response 结构的 photopion operator。
- `hadronic_decay_kernel.f90`
  - pi0 gamma、pion/muon decay、neutrino emissivity。
- `hadronic_pair_production_kernel.f90`
  - gamma-gamma pair production kernel。
- `hadronic_pp_kernel.f90`
  - pp delta channel。
- `hadronic_bethe_heitler_kernel.f90`
  - Bethe-Heitler pair source 和 proton loss。
- `hadronic_hadronic_ic_kernel.f90`
  - hadronic inverse Compton。
- `hadronic_species_transport_kernel.f90`
  - neutron、charged pion、charged muon 显式输运。
- `hadronic_acceleration_kernel.f90`
  - species-aware acceleration timescale、injection operator、`gamma_max` estimate。
- `hadronic_secondary_radiation_kernel.f90`
  - pion/muon synchrotron 和 IC。
- `hadronic_common.f90`
  - hadronic 共享常数和辅助。

`FS_hadronic_1d.f90` 暴露的 shell 级入口包括：

- `fs_hadronic_1d`
- `fs_hadronic_proton_syn_shell`
- `fs_hadronic_pgamma_operator_shell`
- `fs_hadronic_pair_production_shell`
- `fs_hadronic_pp_delta_shell`
- `fs_hadronic_bethe_heitler_shell`
- `fs_hadronic_hadronic_ic_shell`
- `fs_hadronic_species_transport_shell`
- `fs_hadronic_acceleration_shell`
- `fs_hadronic_secondary_radiation_shell`
- `fs_hadronic_decay_operator_shell`

## 5. Hadronic 当前状态

配置入口：

- `Radiation`
  - `epsilon_p`、`proton_synch`、`pg`、`bethe_heitler`、`hadronic_inverse_compton`、`pp`、`neutrino`、`eta_acc`、`pgamma_scheme`。
- `Setups`
  - `hadronic_enabled`、`hadronic_solver`、`num_gam_p`、`num_nu_nu`、`pgamma_scheme`。
- `HadronicConfig`
  - `enabled`、`solver`、`epsilon_p`、`p_p`、`eta_acc`、`include_*` flags。

当前 solver 名称：

- `legacy_1d`
  - proton transport + proton synchrotron。
  - 不支持 p-gamma/neutrino feedback；请求这些通道会报错。
- `am3_1d`
  - 当前 formal hadronic 主路径。
  - 只支持 forward shock + 1D electron solver。

`pgamma_scheme` canonical 名称：

- `disabled`
- `hummer_2010_response`
- `ka2008_reference`

兼容别名仍可解析：

- `hummer_2010`、`hummer2010`、`am3_reference`、`am3_numeric`、`am3_numerical`、`am3` -> `hummer_2010_response`
- `ka2008`、`aharonian_2008`、`kelner_aharonian_2008` -> `ka2008_reference`

当前 formal 1D hadronic 已接入：

- proton injection、adiabatic cooling、synchrotron cooling。
- proton synchrotron。
- `hummer_2010_response` photopion source/loss：
  - `alpha_p(E)` 和 `Q_p^{reinj}(E)` 反馈到 proton transport。
  - `alpha_gamma^{pγ}` 在 hadronic stage 转为 local shell survival factor：
    \[
    \tau_{p\gamma}(\nu,r)=\alpha_\gamma^{p\gamma}(\nu,r)\frac{R}{12\Gamma c}
    \]
    \[
    f_{\rm surv}=\frac{1-\exp(-\tau)}{\tau}
    \]
  - survival factor 先作用在 photon field 和 hadronic luminosity 上，再进入 observer projection。
- neutrino output。
- Bethe-Heitler：
  - proton continuous cooling。
  - secondary electron/positron transport。
  - BH 次级并回 forward electron distribution 后重算 forward synchrotron seed。
- proton-proton：
  - gamma、neutrino、pair source、proton loss。
- hadronic inverse Compton：
  - formal runtime 中 proton channel 可用。
  - pion/muon IC 依赖显式 secondary species 分布。
- explicit secondary species transport：
  - neutron。
  - pion plus/minus。
  - muon minus/plus left/right。
- secondary radiation：
  - pion synchrotron。
  - muon synchrotron。
  - pion IC。
  - muon IC。
- pair production：
  - 当前是 observer-side attenuation + pair synchrotron branch。
  - 不是完整 time-dependent photon kinetic cascade PDE。

`ka2008_reference` 当前是 reference emission backend，不提供 proton transport feedback terms。

## 6. 构建和测试入口

默认环境是 WSL Ubuntu + uv，命令必须从 WSL 内以 `rtk` 开头执行。

构建相关文件：

- `build_extensions.py`
  - 直接 f2py 编译入口，也是 native module 编译拓扑的权威来源。
- `pyproject.toml`
  - 注册 `asgard-demo = lc_spec_demo:main` 和 `asgard-build = build_extensions:main`。
- `setup.py`
  - `build_native` 会触发 native Fortran 编译。

常用构建入口：

```bash
rtk uv run python build_extensions.py --module FS_electron_fullhide_2d --force
rtk uv run python build_extensions.py --module FS_electron_charint_2d --force
rtk uv run python build_extensions.py --module FS_hadronic_1d --force
rtk /usr/bin/gfortran --version
```

基础 smoke：

```bash
rtk uv run python tests/readme_smoke_bench.py
rtk uv run python tests/fullhide_2d_smoke_bench.py
rtk uv run python tests/fullhide_2d_medium_diag.py
```

hadronic runtime 回归入口：

```bash
rtk uv run python tests/hadronic_1d_smoke.py
rtk uv run python tests/hadronic_species_transport_smoke.py
rtk uv run python tests/hadronic_secondary_radiation_smoke.py
rtk uv run python tests/hadronic_acceleration_smoke.py
rtk uv run python tests/hadronic_bethe_heitler_smoke.py
rtk uv run python tests/hadronic_hadronic_ic_smoke.py
rtk uv run python tests/hadronic_pair_production_smoke.py
rtk uv run python tests/hadronic_pp_smoke.py
rtk uv run python tests/hadronic_pg_neutrino_1d_diag.py
rtk uv run python tests/hadronic_proton_synch_1d_diag.py
rtk uv run python tests/hadronic_pgamma_benchmark_report.py
```

文档/物理对比入口：

```bash
rtk uv run python tests/vegas_afterglow_comparison.py
rtk uv run python tests/sed_electron_compare.py
```

Fortran 重要改动后必须覆盖：

- 受影响的 `build_extensions.py --module ... --force`。
- `-Wline-truncation` 编译检查。
- 最小相关 smoke 或 diagnostic test。

本文件只整理代码内容，不执行构建或物理 benchmark。

## 7. 已知边界

当前不要把以下内容描述成已完成：

- reverse-shock hadronic processes。
- 2D / chi-resolved hadronic transport。
- 完整 time-dependent pair cascade PDE。
- 完整 photon kinetic cascade 闭环。
- full-resolution Vegas benchmark figures 的最终平滑性闭环。
- proton-synchrotron peak placement 的 AM3-side 物理交叉验证。

开发时需要保持的边界：

- ASGARD 主架构是 shell-evolving blast-wave + observer projection，不是 AM3 one-zone PDE driver。
- AM3 只能作为微观核和过程拆分参考，不应替换 ASGARD 的 dynamics/electron/observer 主链。
- Python hadronic 代码可以做 orchestration、wrapping、benchmark、API glue；最终 AM3-derived microphysics 应落在 `src/Hadronic/*.f90`。
- 不要用绘图后处理隐藏非连续或非光滑的物理时间演化量；出现不光滑结构时，先按物理/数值 bug 处理。
- 不要提交 `.buildcache/`、临时 debug 脚本或失败占位图。
