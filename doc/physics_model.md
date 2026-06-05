# 物理模型

本文档描述当前 ASGARD 已实现的物理模型和边界。这里记录的是代码契约，不是完整教材式推导。

## 总体图像

ASGARD 的主线是：

```text
relativistic blast wave
  -> electron / proton distribution evolution
  -> local radiation and photon fields
  -> absorption / cascade / cross-zone coupling
  -> equal-arrival-time observer projection
```

Python 层组织状态机、配置和观测投影；Fortran 层求解电子、辐射、动力学和强相互作用微物理核。

## 正激波

正激波分支包含：

- blast-wave dynamics
- electron injection and transport
- synchrotron radiation
- synchrotron self-absorption
- SSC and Compton cooling
- gamma-gamma absorption
- optional hadronic emission and feedback
- observer-frame projection

重要微物理参数：

- `eps_e`：电子能量分数。
- `eps_B`：磁场能量分数。
- `p`：注入电子谱指数。
- `xi_N`：非热电子数分数。
- `ssc`：是否输出 SSC。
- `ssc_cooling` / `index_y`：冷却中是否包含 Compton 项以及使用哪种近似。

电子谱演化是核心物理状态，不能用后处理 smoothing 修复不连续结果。

## 电子输运

当前登记求解器：

- `fullhide_1d`：默认稳定基线，适合一般 public runtime 和拟合。
- `slc1_1d`：semi-Lagrangian / characteristic family path。
- `charint_1d`：characteristic integration path。
- `t2g1_1d`：legacy implicit transport path。
- `weno5_1d`：高阶电子谱解析路径。
- `fullhide_2d`：energy + chi resolved electron transport。
- `charint_2d`：2D characteristic path。

电子输运输出：

- `gam_e`
- `d_n_gam_e`
- `l_syn_spec`
- `seed_syn`
- `nu_m`
- `nu_c`
- `nu_a`

2D path 额外输出：

- `d_n_gam_e_chi`
- `chi_grid`

## 同步辐射、SSA 与 SSC

同步辐射核：

- `src/Electron/electron_radiation_kernel.f90`
- `src/Radiation/radiation_common.f90`

SSA cooling 和 transfer：

- `src/Electron/electron_cooling_kernel.f90`
- `src/Radiation/radiation_common.f90`

SSC：

- `src/Radiation/radiation_ssc_spectrum.f90`
- `asgard_core/asgard_ssc.py`

当前同步辐射积分选择中，`index_syn_integr=1/2` 是固定网格快速路径；adaptive path 只作为显式诊断路径使用，不作为 public 默认。

## 反激波

反激波当前基线：

- electron synchrotron
- RS SSC
- FS/RS cross-zone IC
- optional RS hadronic light path
- optional full-chain RS hadronic dispatch through formal 1D hadronic kernels

关键物理约束：

- 注入能标使用 shock-front `gamma34`。
- 区域 3 turbulent field 和 crossing 后热演化使用显式 `U3/V3` thermal state。
- `reverse_sigma` 引入 upstream magnetization；磁化 jump 使用 VegasAfterglow 的 jump-condition 形式作为来源和 comparison backend。
- `B3` 是 turbulent + ordered total field。
- `sigma -> 0` 必须回到当前非磁化 baseline。

VegasAfterglow 在当前项目中是 comparison backend，不是光变目标或 RS 物理基准。

## 强子过程

Forward-shock hadronic solver：

- `legacy_1d`：proton transport + proton synchrotron baseline。
- `am3_1d`：formal research path，覆盖 p-gamma、BH、pp、hadronic IC、secondary species transport、secondary radiation、pair production branch、neutrino。

Hadronic process switches：

- `epsilon_p`
- `proton_synch`
- `pg`
- `bethe_heitler`
- `hadronic_inverse_compton`
- `pp`
- `neutrino`
- `pair_production`
- `pgamma_scheme`
- `eta_acc`

RS hadronic：

- `hadronic_reverse_1d`：light proton injection/transport + proton synchrotron。
- Full-chain RS path 使用 `hadronic_forward_1d` formal kernels，并使用 RS seed photons、RS `B3`、shell energy 和 baryon target density。

当前 hadronic 边界：

- Hadronic transport 保持 1D shell-level。
- 2D / chi-resolved hadronic transport 有意不实现。
- 未来若扩展，必须先定义 chi-local photon field、hadron density、secondary feedback 和 observer projection。

## 对产生与级联

当前已实现路径：

- observer-side gamma-gamma attenuation
- pair-production branch
- `pair_cascade_iterations > 1` 时使用 shell-sequence time-dependent gamma-gamma pair/synch cascade

未实现：

- IC-mediated electromagnetic cascade。

原因是该扩展需要 photon/e± source-sink 方程、IC kernel 契约和 energy-budget benchmark。见 `doc/pair_cascade_extension_boundary.md`。

## 偏振

偏振路径：

- synchrotron Stokes projection
- FS/RS electron synch
- FS/RS hadronic synch

非同步辐射分支不混入 polarization Stokes。当前 Lan 2023 对比记录显示峰值幅度已匹配，峰时仍偏早；现有证据指向 dynamics/jet-evolution benchmark，而不是 projection-layer 修正。

## 观测者投影

Observer stage 组合：

- EATS/Doppler interpolation
- redshift 和 luminosity distance
- synchrotron/SSC/hadronic/RS/cross-zone components
- absorption factors
- structured jet 和 sky image 的 patch integration

主要 Fortran interpolation：

- `src/Interpolation/SED_interpolation.f90`
- `src/Interpolation/SED_interpolation_structured.f90`

Projection-layer 修正不能用来掩盖 dynamics 或 transport bug。

## 物理验收规则

- 物理 rate 的时间/空间演化应连续且平滑。
- 非光滑最终物理参数轨道在证明前都应优先视为 bug。
- 不用经验 smoothing 或 projection-layer time shift 修正物理不连续。
- 不在真实系统边界之外添加数值保护。
- Python 不替代最终数值微物理核；Python 只做 orchestration、wrapping 和 benchmark。
