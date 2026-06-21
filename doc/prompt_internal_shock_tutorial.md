# Prompt 内部激波快照教程

本页说明当前 `prompt/` 目录中的内部激波快照工作流。它是 prompt emission 的研究诊断入口，不是 `asgard_core` 顶层 public API、不是拟合接口，也不是完整 prompt emission 模块。当前目的有三点：

- 用两壳碰撞验证磁化 jump condition、FS/RS 分支和 EATS 投影是否能复用 ASGARD 已有数值核。
- 生成 formal prompt light curve、spectrum、component 和 Band-reference 诊断图。
- 暴露内部激波 prompt 模型进入正式主线前还缺少的物理契约。

## 1. 入口和文件

当前入口集中在：

| 文件 | 角色 |
| --- | --- |
| `prompt/internal_shock.py` | 两壳碰撞、接触面 Lorentz 因子、FS/RS jump、branch history。 |
| `prompt/eats.py` | prompt branch 的 EATS 投影包装，复用 `src/Interpolation/SED_interpolation.f90`。 |
| `prompt/radiation.py` | 复用现有 reverse electron kernel、SSC、\(\gamma\gamma\) absorption，输出 FS/RS sync/SSC。 |
| `prompt/run_snapshot.py` | 快速 snapshot 图件入口。 |
| `prompt/run_formal_results.py` | formal prompt 图件、Band reference、CSV/NPZ 诊断入口。 |
| `tests/internal_shock_prompt_smoke.py` | 当前最小物理和接口冒烟测试。 |

导入方式：

```python
from prompt.internal_shock import InternalShockShell, InternalShockNumerics, simulate_internal_shock
from prompt.radiation import InternalShockMicrophysics, RadiationNumerics, compute_prompt_observed_flux
from prompt.eats import EATSNumerics
```

不要把这些对象写成：

```python
from asgard_core import InternalShockShell
```

`asgard_core` 当前只公开余辉 `Model`、`Fitter`、介质、喷流、辐射和求解器对象。

## 2. 两壳碰撞几何

慢壳和快壳由 `InternalShockShell` 给出：

```python
slow = InternalShockShell(gamma=200.0, total_energy_iso_erg=4.0e52, duration_s=0.2, sigma=0.01)
fast = InternalShockShell(gamma=600.0, total_energy_iso_erg=4.0e52, duration_s=0.2, sigma=0.03)
```

壳层宽度为实验室系宽度

\[
\Delta_i=c\,\Delta t_i .
\]

若中心引擎发射间隔为 \(t_{\rm gap}\)，慢壳与快壳速度为 \(\beta_s c\) 和 \(\beta_f c\)，碰撞半径由

\[
\beta_s c\,t_{\rm col}
=
\beta_f c\,(t_{\rm col}-t_{\rm gap})
\]

推出：

\[
R_{\rm col}
=
c\,t_{\rm gap}
\frac{\beta_s\beta_f}{\beta_f-\beta_s},
\qquad
t_{\rm col}
=
\frac{R_{\rm col}}{\beta_s c}.
\]

代码对应 `simulate_internal_shock` 中的 `radius_collision` 和 `lab_collision_time`。

## 3. 磁化能量拆分和壳层密度

`total_energy_iso_erg` 是总 isotropic-equivalent 壳层能量。上游磁化为

\[
\sigma
=
\frac{B_4^2}{4\pi\Gamma_4^2\rho_4c^2}.
\]

与反向激波余辉主线相同，baryonic matter fraction 是 \(1/(1+\sigma)\)。因此壳层 baryonic mass 为

\[
M_{b,i}
=
\frac{E_{{\rm iso},i}}{(1+\sigma_i)\Gamma_i c^2}.
\]

代码路径是 `baryonic_mass_g(shell)`，实际调用 `dynamics_common.rs_shell_matter_fraction(shell.sigma)` 固定同一物理闭合。

碰撞半径处的壳层本征数密度为

\[
n'_i(R_{\rm col})
=
\frac{M_{b,i}}
{4\pi R_{\rm col}^2\,\Gamma_i\,\Delta_i\,m_p}.
\]

这是 `shell_proper_number_density(shell, radius_collision)` 的输入。增加 \(\sigma\) 不只是增加 ordered field，也会降低同一总能量下的 baryonic mass 和上游密度。

## 4. 接触面压力平衡

接触面 Lorentz 因子 \(\Gamma_c\) 在两壳 Lorentz 因子之间求根。对任一上游壳层 \(i\)，相对 Lorentz 因子是

\[
\gamma_{{\rm rel},i}
=
\Gamma_c\Gamma_i(1-\beta_c\beta_i).
\]

MHD jump 给出压缩比 \(C_i(\gamma_{\rm rel},\sigma)\) 和下游热比内能 \(\epsilon_{{\rm th},i}(\gamma_{\rm rel},\sigma)\)。代码通过

```text
dynamics_common.rs_mag_comp
dynamics_common.rs_mag_specific_internal
```

复用余辉反向激波基线。下游总压强写成热压强加 ordered magnetic pressure：

\[
p_i
=
\frac{1}{3}
C_i\epsilon_{{\rm th},i}
n'_i m_pc^2
+
\frac{1}{2}
C_i^2\sigma_i n'_i m_pc^2 .
\]

`_solve_contact_gamma` 求解

\[
p_{\rm FS}(\Gamma_c)-p_{\rm RS}(\Gamma_c)=0.
\]

若根不在慢壳和快壳 Lorentz 因子之间，代码直接暴露错误，不构造兜底接触面。

## 5. FS/RS jump 与 branch history

`_build_jump` 为 forward branch 和 reverse branch 分别建立 `BranchJump`。相对上游速度满足

\[
\beta_{\rm sh,lab}
=
\frac{\beta_c+\beta_{\rm sh,cd}}
{1+\beta_c\beta_{\rm sh,cd}}
\quad({\rm FS}),
\]

\[
\beta_{\rm sh,lab}
=
\frac{\beta_c-\beta_{\rm sh,cd}}
{1-\beta_c\beta_{\rm sh,cd}}
\quad({\rm RS}),
\]

其中 contact-frame shock speed 来自 MHD jump 返回的 downstream four-velocity：

\[
\beta_{\rm sh,cd}
=
\frac{u_d}{\sqrt{1+u_d^2}}.
\]

crossing 时间由壳层宽度除以 shock 与上游的实验室系相对扫过速度：

\[
t_{\rm cross,FS}
=
\frac{\Delta_{\rm slow}}
{c(\beta_{\rm sh,FS}-\beta_{\rm slow})},
\]

\[
t_{\rm cross,RS}
=
\frac{\Delta_{\rm fast}}
{c(\beta_{\rm fast}-\beta_{\rm sh,RS})}.
\]

`fast_shock_allowed(gamma_rel, sigma)` 仍调用余辉反向激波共享函数。若 fast magnetosonic 条件不满足，对应 branch 的 `valid_shock=False`，辐射源项为零；这不是后处理裁剪，而是 shock 不存在的物理状态。

每个有效 branch 在 `num_branch_steps` 个点上记录：

\[
R(t)
=
R_{\rm col}+\beta_{\rm sh,lab}ct,
\]

\[
t_{\rm obs,axis}
=
(1+z)
\left[
t_{\rm col}+t-\frac{R(t)}{c}
\right]
-t_{0,{\rm axis}},
\]

\[
\rho'_u(R)=m_p n'_u(R),
\qquad
u_{\rm th}
=
C\epsilon_{\rm th}\rho'_u c^2.
\]

ordered field 与 turbulent field 为

\[
B_{\rm ord}
=
C B_4,
\qquad
B_{\rm turb}
=
\sqrt{8\pi\epsilon_B u_{\rm th}},
\]

\[
B_{\rm tot}
=
\sqrt{B_{\rm ord}^2+B_{\rm turb}^2}.
\]

扫过质量率为

\[
\frac{\mathrm{d}M}{\mathrm{d}t_{\rm lab}}
=
4\pi R^2\Gamma_u n'_u m_p c
|\beta_{\rm sh}-\beta_u|.
\]

电子注入 luminosity 是

\[
L'_{e}
=
\epsilon_e
\epsilon_{\rm th}c^2
\Gamma_c
\frac{\mathrm{d}M}{\mathrm{d}t_{\rm lab}}.
\]

这组量对应 `BranchHistory` 的 `thermal_luminosity_comoving_erg_s`、`electron_luminosity_comoving_erg_s`、`ordered_b_g`、`turbulent_b_g`、`total_b_g`、`swept_mass_g`、`internal_energy_erg` 和 `comoving_volume_cm3`。

## 6. 辐射链路

`compute_prompt_observed_flux` 对 FS 和 RS branch 分别执行：

```text
BranchHistory
-> electron_reverse_kernel.electron_reverse_evolve
-> electron_secondary_reverse_synchrotron
-> Radiation.ssc_spec
-> Radiation.annihilation
-> prompt.eats.project_branch_flux
```

当前复用 reverse electron kernel，是为了先验证 shock history、电子冷却、同步辐射、SSC 和 EATS 投影的组合。它不表示 prompt branch 已经有独立、最终版的 prompt electron kernel。

电子特征 Lorentz 因子按 branch luminosity 估计。若电子谱指数 \(p>2\)，平均能量给出

\[
\gamma_m
=
1+
\frac{p-2}{p-1}
\frac{L'_e}{\dot{N}'_e m_ec^2},
\]

\[
\dot{N}'_e
=
\xi_e
\frac{\Gamma_c\,\mathrm{d}M/\mathrm{d}t_{\rm lab}}{m_p}.
\]

冷却 Lorentz 因子和最大 Lorentz 因子诊断为

\[
\gamma_c
=
\frac{6\pi m_ec}
{\sigma_TB_{\rm tot}^2t'_{\rm age}},
\]

\[
\gamma_{\max}
=
\left(
\frac{6\pi e}
{\sigma_TB_{\rm tot}\eta_{\rm acc}}
\right)^{1/2}.
\]

这些量写入 `BranchRadiation.gamma_m`、`gamma_c`、`gamma_max`，用于检查谱峰和冷却 regime。

同步与 SSC 源项都乘以 \(\gamma\gamma\) survival：

\[
F_{\nu,{\rm source}}
=
\frac{1+z}{4\pi d_L^2}
L'_{\nu'}
\exp(-\tau_{\gamma\gamma}).
\]

这里的 \(L'_{\nu'}\) 是对应分量的 comoving luminosity array。若要检查能量预算，应分别看 `fs_sync`、`fs_ssc`、`rs_sync`、`rs_ssc`，不要只看 `total`。

## 7. EATS 投影

`project_branch_flux` 将 branch luminosity 投影到观测者时间和频率。边界向量中：

```text
boundary[7] = redshift
boundary[8] = opening_angle_rad
boundary[9] = viewing_angle_rad
```

频率先排序传入 Fortran 插值核，返回后恢复用户输入顺序。若

\[
\theta_{\rm obs}\ne 0
\quad{\rm and}\quad
N_\phi=1,
\]

代码直接拒绝：

```text
off-axis EATS projection requires num_phi >= 2
```

这是轴对称 collapse 的几何边界。离轴投影必须显式解析 \(\phi\)。

## 8. 最小运行示例

```python
import numpy as np

from prompt.eats import EATSNumerics
from prompt.internal_shock import InternalShockNumerics, InternalShockShell, simulate_internal_shock
from prompt.radiation import InternalShockMicrophysics, RadiationNumerics, compute_prompt_observed_flux

slow = InternalShockShell(gamma=200.0, total_energy_iso_erg=4.0e52, duration_s=0.2, sigma=0.01)
fast = InternalShockShell(gamma=600.0, total_energy_iso_erg=4.0e52, duration_s=0.2, sigma=0.03)
microphysics = InternalShockMicrophysics(epsilon_e=0.1, epsilon_b=0.01, electron_index_p=2.3)

solution = simulate_internal_shock(
    slow,
    fast,
    engine_gap_s=0.2,
    redshift=0.5,
    luminosity_distance_cm=1.0e28,
    opening_angle_rad=0.1,
    epsilon_e=microphysics.epsilon_e,
    epsilon_b=microphysics.epsilon_b,
    numerics=InternalShockNumerics(num_branch_steps=64),
)

flux = compute_prompt_observed_flux(
    solution,
    observer_frequency_hz=np.logspace(16.0, 24.0, 64),
    observer_time_s=np.linspace(1.0e-4, 2.0, 128),
    microphysics=microphysics,
    radiation_numerics=RadiationNumerics(num_electron_gamma=121, num_photon_frequency=161, num_threads=4),
    eats_numerics=EATSNumerics(num_theta=32, num_phi=1, num_threads=4, adaptive_rtol=3.0e-3, adaptive_max_depth=6),
)

print(flux.total.shape)
print(solution.gamma_contact, solution.fs.valid_shock, solution.rs.valid_shock)
```

`flux.total` 的形状是 `(num_frequency, num_time)`。

## 9. Formal 图件入口

生成当前 formal prompt snapshot 图件：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python prompt/run_formal_results.py'
```

输出目录是 `prompt/results/`。该目录只提交 `.gitignore`，formal 图件默认作为本地诊断，不自动进入文档资产。

当前脚本会生成：

- `formal_flux.npz`：时间、频率和四个分量矩阵。
- `formal_diagnostics.json`：碰撞半径、\(\Gamma_c\)、FS/RS jump、磁场、\(\gamma_m/\gamma_c/\gamma_{\max}\)、compactness、\(\gamma\gamma\) survival。
- `formal_summary.csv`：代表性能段的峰时和峰值。
- `formal_band_reference.csv`：与 Band function 形状的参考比较。
- light curve、linear light curve、spectrum、component figures。

## 10. 验证口径

最小测试：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/internal_shock_prompt_smoke.py'
```

该测试检查：

- 离轴 `num_phi=1` 被拒绝。
- \(\sigma\rightarrow0\) 的弱磁化 case 接近非磁化 hydrodynamic baseline。
- 磁化 case 的 baryonic mass 下降，ordered field 和 magnetic fraction 上升。
- fast shock 不允许时对应分量为零。
- FS/RS sync/SSC 和 total flux 有限、非负，并且活动光变段分段平滑。

当前不要把它解释为 prompt light curve 已经与观测 GRB 拟合一致。它只是代码和物理链路冒烟。

## 11. 当前边界

- 不从 `asgard_core` 顶层导出。
- 不接 `Fitter`，不作为 afterglow posterior 的参数分支。
- 不含 shell spreading、multi-collision train、curvature-tail engine history、photosphere 或 subphotospheric dissipation。
- 不含 prompt hadronic cascade、pair photosphere 或 IC-mediated electromagnetic cascade。
- 当前电子辐射复用 reverse electron kernel；独立 prompt electron kernel 进入前要先写清方程、网格、能量预算和 smoke benchmark。
- Off-axis 投影必须 `num_phi>=2`；on-axis 可以用 `num_phi=1` 的轴对称 collapse。
