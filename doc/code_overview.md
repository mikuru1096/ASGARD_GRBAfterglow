# ASGARD 代码总览

本文是源码结构、计算主链和程序单元的唯一概览。算法公式见
[`project_algorithm_design.md`](project_algorithm_design.md)，公开接口见
[`public_api.md`](public_api.md)，开发与验证方法见
[`developer_guide.md`](developer_guide.md)。

## 1. 分层原则

ASGARD 把一次计算分成三层：

1. `asgard_core` 接收公开对象，建立配置并编排阶段；
2. `src` 中的 Fortran 核推进动力学、粒子输运和辐射积分；
3. observer 投影把壳层本征光度转换为观测流量。

Python 不重复数值物理，Fortran 不解释公开字符串键。离散选择在公开边界
映射一次，再以整数和数组进入 f2py ABI。运行状态使用 dataclass 明示，不用
manager/context 类隐藏物理阶段。

## 2. 公开入口

`asgard_core/api_model.py` 定义：

- `Model` 与 `Observer`；
- `UniformMedium`、`WindMedium`、`TabulatedMedium`；
- `top_hat_jet`、`gaussian_jet`、`power_law_jet`；
- `Radiation`、`Numerics`、`SolverOptions`、`ReverseShock`、`Hadronic`；
- flux、spectrum、sky image 和 polarization 查询。

`asgard_core/api_fit.py` 定义 `Fitter`、`Param` 和 `FitResult`。`prompt/` 是内部
激波 snapshot 研究工具，不属于 afterglow `Model` 主链。

电子求解器公开名称为 `fullhide_1d`、`fullhide_1d_hz`、`slc1_1d`、
`charint_1d`、`dg_1d`、`charint_2d`、`t2g1_1d`、`weno5_1d` 和
`fullhide_2d`。后端支持矩阵以
[`public_backend_limits.md`](public_backend_limits.md) 为准。

## 3. 计算主链

```text
Model query
  -> RuntimeConfig -> SimulationSetup
  -> solve_dynamics
  -> solve_electron
  -> separated or joint photon/hadron stage
  -> reverse-shock emission
  -> observer assembly
  -> project_flux
  -> API result
```

`asgard_core/asgard_state.py::solve_setup` 是主状态机。它不改变公开配置，只把
同一半径网格上的阶段结果装配为后续输入。`project_flux` 根据查询类型选择
lightcurve、SED、adaptive-theta 或 χ-resolved 投影。

拟合路径复用同一主链：

```text
Fitter.loglike -> compile_problem -> eval_loglike
  -> solve_setup -> project_flux -> combine_flux
```

## 4. 运行状态

`asgard_core/asgard_types.py` 中的主要状态为：

- `DynamicsSolution`：半径、观测时间、bulk Lorentz factor 和扫掠质量；
- `ReverseShockDynamics`：区域 3 热状态、磁场、`gamma34` 和 crossing records；
- `ElectronSolution`：电子谱、同步光度、seed，以及可选 finite-q shell 数组；
- `PhotonFieldState`：同步 seed、强子 target 和吸收 seed；
- `HadronicSolution`：质子、次级粒子、强子辐射和反馈源；
- `ObserverState`：吸收、光度分量和最终投影输入；
- `FluxComponents`：FS/RS synchrotron、SSC、cross-IC 和 hadronic flux。

二维电子状态中的 `chi_grid` 是 q cell 的 BM 等效诊断坐标。真正的 observer
几何使用 `chi_radius_cm`、`chi_gamma_bulk` 和 `chi_dvolume_weight`；它不是
无限 `log10(chi)` 网格。

## 5. Python 文件责任

- `api_model.py`：公开对象、查询缓存和 `Model -> RuntimeConfig` 适配；
- `api_fit.py`：拟合问题编译与 likelihood；
- `asgard_config.py`：配置和 setup dataclass；
- `asgard_setup.py`：配置到网格及 kernel 参数；
- `asgard_runtime.py`：扩展解析、ABI 调用和结果包装；
- `asgard_state.py`：阶段编排及跨区耦合；
- `asgard_postprocess.py`：observer 投影、band 聚合和数据比较；
- `structured_jet_kernel.py`：结构化喷流采样、ring/patch 调度；
- `hadronic_processes.py`：pγ 包装和 log-log 插值；
- `hadronic_am3_solver.py`：formal hadronic shell-sequence 编排；
- `hadronic_cascade.py`：pair-cascade 序列编排。

## 6. Fortran 文件责任

### 动力学

- `src/Dynamics/Dynamics_forward.f90`：FS、介质、注能和密度跃变；
- `src/Dynamics/reverse_shock.f90`：RS 推进与 crossing event split；
- `dynamics_density_profile.f90`：共享密度模型；
- `reverse_shock_state.f90`：RS 状态下标和 phase；
- `reverse_shock_mhd_jump.f90`：有限强度横向场 MHD jump。

### 电子

公开扩展入口位于 `electron_forward_*.f90`。共享实现按职责拆分：

- `electron_common.f90`：常量、网格和基础 helper；
- `electron_cooling_*_kernel.f90`：同步、SSA、IC 和 Y cooling；
- `electron_radiation_kernel.f90`：同步辐射与 transfer；
- `electron_injection_profiles.f90`：注入；
- `electron_shell_transport_common.f90`：壳层推进公共量；
- `electron_transport_common.f90`：1D transport helper；
- `electron_transport_2d_kernel.f90`：有限 q-shell 输运；
- `electron_dg_transport.f90`：DG 空间离散；
- `electron_reverse_kernel.f90`：RS 电子与辐射。

### 辐射与 observer

- `src/Radiation/ssc_spectrum.f90`：SSC；
- `pair_absorption.f90`：γγ opacity；
- `rad_common.f90`：积分权重、截面、同步 seed 和 transfer；
- `syn_polarization.f90`：局域同步 Stokes emissivity；
- `quantum_synch.f90`：量子同步 helper；
- `src/Interpolation/SED_interpolation.f90`：shell、adaptive-theta 和 χ EATS；
- `SED_interpolation_structured.f90`：structured shell projection；
- `src/Structured/structured_jet_1d.f90`：结构化喷流聚合入口。

### 强子

`src/Hadronic/hadronic_forward_1d.f90::formal_transport_1d` 是正式正向壳层序列。
`hadronic_reverse_1d.f90` 提供轻量 RS 入口；完整 RS 链复用正式 1D 核。

- `hadronic_transport_kernel.f90`：质子注入、损失和推进；
- `hadronic_transport_remap_kernel.f90`：守恒 remap；
- `hadronic_rad.f90`：proton synchrotron；
- `hadronic_pg.f90`：Hummer 2010 pγ response；
- `hadronic_decay.f90`：π/μ 衰变和 neutrino；
- `hadronic_bh.f90`：Bethe--Heitler；
- `hadronic_pp.f90`：pp loss、次级粒子及可选 AM3-derived π0 gamma 模型；
- `hadronic_pair.f90`、`hadronic_cascade.f90`：γγ pair 与 cascade；
- `hadronic_ic.f90`：强子次级电子 IC；
- `hadronic_species.f90`：neutron、π± 和 μ± 输运；
- `hadronic_accel.f90`：加速时标和上限；
- `hadronic_secondary.f90`：π/μ synchrotron 与 IC；
- `hadronic_base.f90`：共享常量和网格。

`delta` 是默认 pp gamma 模型；Geant4、SIBYLL、QGSJET 和 Pythia8 是显式
opt-in。

## 7. 关键分支

- separated：电子完成后再计算强子；BH 次级电子可触发一次电子辐射重算；
- joint：同一 R 网格上联合 electron、photon、hadron 和 secondary source；
- `chi_eats_2d`：只替换 FS synchrotron+SSA 的 observer projection；
- reverse hadronic：使用 RS 磁场、seed、能量和 baryon target；
- structured：轴对称使用 theta rings，非轴对称使用 theta-phi patches。

二维强子输运、IC-mediated 完整电磁级联及其它未支持组合，不由 silent fallback
代替。公开请求在边界直接失败。

## 8. 构建图与 ABI

`build_extensions.py` 是扩展源闭包、编译顺序和 f2py 导出面的唯一依据。支持
module 先编译，公开 wrapper 文件最后进入签名生成。删除 Fortran 程序单元前必须
同时检查 Python caller、源闭包和 f2py export tuple。

典型构建：

```bash
TMPDIR=/tmp uv run python build_extensions.py --module Dynamics_forward --force
TMPDIR=/tmp uv run python build_extensions.py --module electron_forward_fullhide_1d --force
TMPDIR=/tmp uv run python build_extensions.py --module hadronic_forward_1d --force
```

Fortran 修改还需从干净 `/tmp` 按同一有序闭包执行 `-Wline-truncation`。不要用
邻近文件的临时集合替代真实闭包，也不要从仓库根目录拾取旧 `.mod`。

## 9. 阅读顺序

新开发者建议按以下顺序定位：

1. `api_model.py` 的公开查询；
2. `asgard_setup.py` 的配置展开；
3. `asgard_state.py::solve_setup` 的阶段顺序；
4. `asgard_runtime.py` 的具体 ABI 调用；
5. `build_extensions.py` 的源闭包；
6. 对应 `src/*` 数值核。

需要确认公式时从算法总纲反查源文件；需要确认公开限制时只查 backend limits；
需要确认未完成问题时只查根目录 `BUG.md` 和 `TODO.md`。不再维护逐程序单元索引，
因为它会与源码和构建图重复并快速失真。
