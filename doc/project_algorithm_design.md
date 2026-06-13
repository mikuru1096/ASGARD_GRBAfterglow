# ASGARD 项目算法设计总纲

本文档给出 ASGARD 当前工作树的项目级算法设计。它描述 public API 如何进入 runtime，Python 和 Fortran 如何分工，各类求解器如何组织，以及验证和 benchmark 如何证明结果可信。

物理总纲见 `doc/project_physics_design.md`；算法专题见 `doc/algorithm_workflow.md`；模块级实现细节见 `doc/code_overview.md`、`doc/numerical_methods.md`、`doc/call_chain.md` 和 `doc/electron_solver_algorithms.md`。

## 1. 总体架构

ASGARD 的架构分为四层：

```text
Public API
-> Python runtime/state orchestration
-> Fortran numerical kernels
-> observer projection / postprocess / fitting
```

职责边界：

- Public API 负责用户语义：`Model`、`Radiation`、`Setups`、`Fitter`。
- Python runtime 负责配置归一化、数组准备、阶段调度、跨阶段耦合、benchmark 和 diagnostics。
- Fortran 负责高代价数值物理核：动力学、电子输运、辐射积分、强子微物理、observer interpolation。
- Postprocess 负责投影、band aggregation、fit likelihood 和输出对象组装。

Python 不替代最终微物理核；Fortran 不直接承担 public API 兼容逻辑。

## 2. Public API 到配置对象

典型入口：

```python
model = Model(jet, medium, observer, radiation, setups=setups)
flux = model.flux_density_grid(times_s, nu_hz)
```

内部转换为：

```text
Model
-> FitConfig
-> SimulationSetup
-> SolveState
-> FluxResult / details / fit output
```

关键 dataclass：

- `FitConfig`：物理参数、数值网格、solver names、flags。
- `SimulationSetup`：boundary vector、seed frequency grid、observer time grid、distance。
- `SolveState`：dynamics、electron、photon field、hadronic、observer 和 diagnostics 的完整状态。
- `FluxComponents`：各辐射分量的 observer-side arrays。

## 3. 主状态机

核心入口是：

```text
asgard_core/asgard_state.py::solve_state_from_setup
```

默认 separated 主链：

```text
solve_dynamics
-> solve_electron
-> build_photon_field
-> solve_hadronic
-> solve_reverse_shock_emission
-> assemble_observer_stage
-> SolveState
```

joint feedback 主链：

```text
solve_dynamics
-> solve_joint_forward_stage
   -> electron / photon / hadronic / secondary feedback iterations
-> solve_reverse_shock_emission
-> assemble_observer_stage
-> SolveState
```

状态机规则：

- 每个阶段只接受已经明确归一化的上游状态。
- 不支持组合在进入高代价 kernel 前报错。
- 结果通过 dataclass 字段传递，不用隐式全局状态。
- observer projection 不回写 transport state。

## 4. 网格设计

主要网格：

- radius / observer-time shell grid：动力学和 shell state。
- `gamma_e`：电子 log-gamma grid。
- `gamma_p`：质子 log-gamma grid。
- `nu`：radiation / seed photon frequency grid。
- `nu_nu`：neutrino frequency grid。
- `chi`：2D electron / chi-resolved observer projection grid。
- `theta, phi`：observer angular / structured jet integration grids。

网格设计原则：

- 输运方程在自己的保守变量上推进。
- 频率/能量变量转换必须显式保留 Jacobian。
- shell-to-observer projection 可以插值，但不能改变本地输运结果。
- chi-resolved projection 是 FS synchrotron+SSA 的 observer geometry contract，不自动扩展到 hadronic。

## 5. 动力学算法

Forward dynamics：

- Fortran source：`src/Dynamics/Dynamics_forward.f90`。
- 输入 boundary vector 和 medium/jet 参数。
- 输出 `R`, `Gamma`, `t_obs`, swept mass。
- 支持 ISM、wind `k=2`、density jump、energy injection。

Reverse dynamics：

- Fortran source：`src/Dynamics/Dynamics_reverse.f90`。
- crossing 前使用 `m3_frac` 分支；crossing 后固定 shocked ejecta mass。
- 输出 region-3 thermal/magnetic state 和 shock-front diagnostics。

算法边界：

- jet spreading backend dynamics 未实现。
- 任意用户自定义介质的 Fortran dispatch 未实现。
- RS magnetization 必须保持 `sigma -> 0` 回归。

## 6. 电子输运算法

所有 1D 电子 solver 都围绕 log-gamma transport：

```text
source injection
-> cooling face speeds
-> conservative / implicit / characteristic update
-> synchrotron seed recomputation
```

当前求解器：

- `fullhide_1d`：默认稳定 baseline。
- `slc1_1d`：semi-Lagrangian family。
- `charint_1d`：characteristic integration。
- `t2g1_1d`：legacy implicit transport。
- `weno5_1d`：高阶谱解析路径。
- `fullhide_2d`：energy + chi resolved transport。
- `charint_2d`：2D characteristic path。

`fullhide_1d` 的 coupled entry：

```text
fs_electron_fullhide_1d_coupled
```

用于 joint feedback，额外输入 photon cooling seed 和 secondary electron source。

## 7. 冷却与辐射算法

Cooling assembly 将不同物理项组装为 electron transport face speed：

```text
adiabatic cooling
synchrotron cooling
IC / SSC cooling
SSA cooling
```

重要源文件：

- `src/Electron/electron_cooling_kernel.f90`
- `src/Electron/electron_radiation_kernel.f90`
- `src/Radiation/radiation_common.f90`
- `src/Radiation/radiation_ssc_spectrum.f90`
- `src/Radiation/radiation_gamma_gamma_absorption.f90`
- `src/Radiation/radiation_reverse_seed.f90`

算法原则：

- `index_syn_integr=1/2` 是固定网格快速路径。
- adaptive synchrotron integration 只作为显式诊断路径。
- IC cooling 与 SSC photon source 的预算测试必须使用同一 seed。
- SSA transfer 要保持频率 cell 间连续。

## 8. 光子场构建算法

`_build_photon_field_stage` 从 electron solution 构建：

```text
forward_syn_seed
hadronic_forward_ssc_seed
hadronic_target_seed
absorption_syn_seed
absorption_ssc_seed
```

默认 path 中 photon field 是 hadronic target 和 absorption 的输入。joint path 中 photon field 会被 pγ/BH/gamma-gamma survival 和 pair-synch seed 更新。

算法边界：

- `PhotonFieldState` 字段语义不因 joint 改变。
- observer luminosity 不能直接作为 local photon density。
- 新 photon feedback source 必须明确 luminosity-to-density 归一化。

## 9. 强子算法

入口：

```text
asgard_core/asgard_runtime.py::solve_hadronic
src/Hadronic/hadronic_forward_1d.f90
```

formal `am3_1d` path 组合：

- proton injection / transport。
- proton synchrotron。
- pγ Hummer response。
- BH pair production。
- pp source/loss。
- secondary species transport。
- secondary radiation。
- hadronic IC。
- pair production / cascade。
- neutrino luminosity。

Python wrapper 负责：

- 把 electron seed photon field 转换成 hadronic target units。
- 调用 Fortran microphysics。
- 把 shell-local rates 转换成 transport source/loss。
- 组装 `HadronicSolution`。

当前 hadronic transport 的正式 contract 是 1D shell-level。

## 10. 含时二级反馈算法

joint feedback 的算法专题见 `doc/joint_secondary_feedback_algorithm.md`。核心约束是：

```text
所有输运都在 R 坐标上推进
所有 s^-1 rate 通过 1/(beta Gamma c) 转为 per-R 系数
二级 e± 源项在电子求解前进入电子方程
photon survival/source 在观测投影前更新 joint photon field
```

其中 rate 换算因子为

\[
\frac{\mathrm{d}t'}{\mathrm{d}R}
=\frac{1}{\beta\Gamma c}.
\]

当前 fixed small iteration count 是内部实现细节，不是 public convergence parameter。

## 11. Pair cascade 算法

\(\gamma\gamma\) pair branch 有两种路径：

- single-shell / single-pass pair production diagnostics。
- `pair_cascade_iterations > 1` 时使用 shell-sequence time-dependent pair/synch cascade。

算法输出：

- pair density / pair injection。
- pair synchrotron luminosity。
- pair synchrotron seed。
- tau / photon survival diagnostics。

IC-mediated electromagnetic cascade 不在当前算法 contract 内。

## 12. Reverse-shock 和 cross-zone 算法

Reverse-shock emission 由专门 runtime stage 求解：

```text
solve_reverse_shock_emission
```

包括：

- RS electron synchrotron。
- RS SSC。
- optional cross-zone IC。
- optional RS hadronic light / formal dispatch。

Cross-zone IC 需要同时构建 FS/RS seed fields，并在 observer assembly 中作为独立分量输出。

## 13. Observer projection 算法

projection 入口：

```text
project_flux_grid
observe_components_from_setup
```

主要 Fortran interpolation：

- `src/Interpolation/SED_interpolation.f90`
- `src/Interpolation/SED_interpolation_structured.f90`

Projection 步骤：

```text
local shell luminosity
-> Doppler transform
-> redshift
-> luminosity distance
-> EATS interpolation
-> requested observer time/frequency grid
```

`projection_kind="lightcurve"` 是光变和 fitting 热路。`projection_kind="sed"` 是固定时刻谱和频段积分路径。

`chi_eats_2d` 只替换 FS synchrotron+SSA 的 lightcurve projection；非 chi 分量保持 shell-level projection。

## 14. 结构化喷流算法

结构化喷流通过 angular patches / structured Fortran backend 聚合：

```text
jet profile E_iso(theta, phi), Gamma0(theta, phi)
-> patch-local solve
-> observer projection
-> flux summation
```

算法约束：

- patch 间共享 public model semantics，但每个 patch 有自己的局域能量和初始 Lorentz factor。
- 非轴对称投影必须使用合理 `phi` 采样。
- hadronic feedback 当前不提升为 patch-local chi-resolved transport。

## 15. Polarization 算法

Polarization path 对 synchrotron emissivity 做 Stokes projection：

```text
局域同步辐射 emissivity
-> magnetic geometry model
-> Stokes I/Q/U
-> observer projection
-> polarization fraction and angle
```

非同步辐射不进入 Stokes \(I,Q,U\)。峰值幅度和峰时分开验收；projection 层不做经验 timing correction。

## 16. Fitting 算法

`Fitter` 将数据和参数路径编译成 likelihood problem：

```text
Param path -> model field update
-> solve_state_from_setup
-> project_flux_grid
-> combine_multiband_flux
-> redchi / loglike
```

算法原则：

- fitting 不改变物理求解器。
- cached query time 和 cold solve time 必须分开报告。
- 参数边界是用户/外部输入边界，可以校验；内部 stage 不添加防御性 fallback。

## 17. 构建算法

`build_extensions.py` 负责 f2py extension 构建：

- 选择 module source closure。
- 处理 f2py signature compatibility。
- 调用 Meson/f2py 构建。
- 保持真实 gfortran build 使用原始 Fortran source。

Fortran 重要改动后必须：

```text
affected build_extensions.py --force
gfortran -Wall -Wline-truncation -Werror=line-truncation -fsyntax-only
minimal smoke test
```

line-truncation check 必须使用干净 `/tmp` module directory，避免 stale `.mod` 污染。

## 18. Benchmark 与验证算法

验证分层：

- syntax / encoding / whitespace。
- build extension。
- line-truncation source closure。
- narrow smoke tests。
- physics benchmark。
- comparison benchmark。

关键 benchmark：

- VegasAfterglow comparison。
- magnetized reverse-shock sigma scan。
- chi-resolved EATS benchmark。
- time-dependent BH / joint secondary feedback benchmark。
- runtime breakdown benchmark。

Benchmark 应回答明确假设，例如：

- weak-feedback 是否回到 separated baseline。
- strong wind / strong BH 是否平滑。
- `sigma -> 0` 是否回到 unmagnetized RS。
- chi grid refinement 是否收敛。
- IC/BH/gamma-gamma 能量预算是否闭合。

不为填表穷举低信息增益实验。

## 19. 输出与诊断设计

主要输出层：

- `FluxResult` / `FluxComponents`：public flux。
- `TrackBundle` / `details()`：内部轨迹和 diagnostics。
- benchmark CSV / PNG / PDF / metadata。

诊断字段的设计原则：

- 能追踪到上游物理量。
- 不改变主求解器状态。
- 不替代正式源项。
- 不用显示裁剪或 smoothing 掩盖异常。

## 20. 开发准则

算法改动必须满足：

- 最小 public interface 变化。
- 保持 existing dataclass semantics，除非明确迁移。
- 物理核优先放 Fortran。
- Python 负责 orchestration，不写经验物理补丁。
- 不支持组合直接报错。
- 不加无法发生场景的防御性代码。
- 任何非光滑物理结果先查 bug。

当前唯一 TODO 入口是 `TODO.md`。不要在新文档中分散待办列表；文档只描述当前能力、边界和准入条件。
