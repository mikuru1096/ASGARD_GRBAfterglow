# ASGARD Agent Notes

## 1. 项目定位

ASGARD 是一个用于 GRB 余辉建模与拟合的数值代码库。当前仓库的主目标不是重写物理核，而是在保留 Fortran 数值求解器的前提下，统一 Python 接口、构建链路、后处理和测试入口。所有程序的第一优先是精确，第二是尽量快的运行速度，太慢的测试尽量不做，要做也要先进行可行范围内的速度优化。

当前第一阶段重构的边界非常明确：

- 保留 Fortran 数值离散、辐射公式和主计算流程。
- 将 Python 入口统一为 `FitConfig -> run_fit -> FitResult`。
- 补齐能谱输出与绘图接口。
- 修复已知坏点，但不在非数值模拟代码中添加“数值保护式”逻辑。

- 默认只做与当前改动直接相关的最少编译和最小回归测试；不要无关地全量重编译或全量回归，除非当前问题本身要求那样做。
- 主进程负责构建计划、review、处理困难任务；边界清晰的其余任务优先交给子进程完成。
- 使用子进程时必须严格控制子进程上下文长度(max<10000 tokens)，只传递完成当前子任务所必需的最小上下文；每轮对话都要检查上下文占用。
- 当任一工作线程或子进程的上下文窗口达到50000 tokens时，必须自动进行一次上下文压缩，再继续后续工作。
- 给子进程分配任务时，优先使用能够完成该子任务的最快、最便宜模型；由主线程自行判断能力与成本的平衡，不默认使用高成本模型。

## 1.1 子 Agent 使用约束

- subagent_model = "gpt-5.4-mini"
- 当前项目允许使用子 agents 协助分析与实现。
- 同一轮任务中最多使用 `5` 个子 agents。
- 主线阻塞修改仍由主 agent 直接完成，不把关键路径完全外包。

## 2. 当前目录结构

### 顶层 Python 门面

- `mergered.py`
  - 当前门面模块。
  - 导出 `FitConfig`、`FitResult`、`run_fit`、`plot_light_curve`、`plot_spectrum`。
  - 保留旧 `fit(**kwargs)` 兼容包装，但不再是主接口。

- `asgard_models.py`
  - 定义配置与结果数据结构。
  - 包含：
    - `SpectrumOutputConfig`
    - `FitConfig`
    - `FitResult`
    - 波段和频率常量

- `asgard_solver.py`
  - 新的单入口编排层。
  - 负责：
    - 构建观测时间网格
    - 构建频率网格
    - 调用动力学、电子谱、辐射、插值模块
    - 组装多波段光变
    - 计算 `redchi`
    - 在启用时生成稠密能谱矩阵

- `asgard_runtime.py`
  - 运行时绑定层。
  - 负责按 `FitConfig` 选择真实 Fortran 后端。
  - 当前已经接通：
    - `weno5=False` -> `FS_electron_fullhide`
    - `weno5=True` -> `FS_electron_weno5`
  - `reverse=True` 当前显式报错，不再作为静默无效开关保留。

- `asgard_plot.py`
  - 绘图层。
  - 提供：
    - `plot_light_curve(result, ...)`
    - `plot_spectrum(result, times_s, quantity="nufnu"|"fnu", ...)`

### 运行与拟合脚本

- `hand_my.py`
  - 最小运行示例。
  - 当前已迁移到新 API，并包含能谱输出调用。

- `MCMC.py`
  - MCMC 采样入口，已迁移到新 API。

- `pymultinest_demo.py`
  - MultiNest 示例入口，已迁移到新 API。

### 拟合与观测数据

- `cal_chi2_light_curve.py`
  - 光变卡方模块。
  - 读取 `data_light_curve/` 下的观测文件。
  - 当前已修正时间单位：数据按“天”读入，统一转换为“秒”后与模型比较。

- `cal_chi2_spectrum.py`
  - 能谱卡方模块。
  - 当前已重写为可正常解析 2/3/4/6 列数据格式。

- `data_light_curve/`
  - 光变观测数据。

- `data_spectrum/`
  - 能谱观测数据。

### Fortran 数值核

- `src/Dynamics/`
  - 激波动力学演化。

- `src/Electron/`
  - 电子分布演化与同步辐射相关物理。

- `src/Radiation/`
  - SSC、对湮灭、EBL 等辐射过程。

- `src/Interpolation/`
  - 从共动系谱到观测者谱的 EATS 与多普勒插值。

- `src/Constants.f90`
  - 全局常数与单位换算。

## 3. 当前公开 API

### `FitConfig`

统一管理以下输入：

- 动力学参数：
  - `index_dyn`
  - `num_r`
  - `eta_0`
  - `e_iso`
  - `d_ne`
  - `a_star`
  - `r_tr`
  - `f_jump`
  - `f_wide`
  - `r0`
  - 注能参数 `e_inj_t1/e_inj_t2/l_inj_0/q_inj`

- 微物理参数：
  - `epsilon_e`
  - `epsilon_b`
  - `p`
  - `f_e`
  - `index_y`
  - `index_syn_integr`

- 观测与几何参数：
  - `z`
  - `opening_angle_jet`
  - `theta_v`
  - `num_theta`
  - `num_phi`

- 网格与输出参数：
  - `num_gam_e`
  - `num_nu`
  - `num_tobs`
  - `t_obs_min_log10`
  - `t_obs_max_log10`
  - `plot_lc`
  - `show_plots`

- 能谱输出子配置：
  - `spectrum_output.enabled`
  - `spectrum_output.num_nu_obs`
  - `spectrum_output.nu_min_hz`
  - `spectrum_output.nu_max_hz`

### `FitResult`

统一返回：

- `t_obs_s`
- `bands`
- `bands_flux`
- `redchi`
- `nu_m`
- `nu_c`
- `nu_a`
- `spectrum_freq_hz`
- `spectrum_fnu`

其中：

- `bands_flux` 的行顺序与 `result.bands` 一致。
- 默认仅在 `spectrum_output.enabled=True` 时填充 `spectrum_freq_hz` 与 `spectrum_fnu`。

### 主入口

```python
from mergered import FitConfig, run_fit

config = FitConfig()
result = run_fit(config)
```

### 新增能谱绘图接口

```python
from mergered import plot_spectrum

fig = plot_spectrum(
    result,
    times_s=[1e3, 1e4, 1e5],
    quantity="nufnu",
    outfile="Radiation_Spectra.pdf",
    show=False,
)
```

接口行为：

- `quantity` 只允许 `"nufnu"` 或 `"fnu"`。
- `times_s` 按最近邻观测时刻取样。
- 若 `result` 不含能谱矩阵，明确抛出异常，不做静默回退。

## 4. 数值主链与物理含义

当前主链在 `asgard_solver.py::run_fit` 中按下面顺序执行：

1. 动力学：`Dynamics.dynamics_forward`
2. 电子谱：`Electron.fs_electron_fullhide`
3. SSC：`Radiation.ssc_spec`
4. `gamma-gamma` 湮灭：`Radiation.annihilation`
5. 观测者系插值：`Interpolation.sed_interpolation`
6. 光变卡方：`cal_chi2_light_curve`

### 4.1 动力学模块

文件：

- `src/Dynamics/Dynamics_forward.f90`
- `src/Dynamics/Dynamics_reverse.f90`

当前默认路径使用 `dynamics_forward`。

`Dynamics_forward.f90` 中：

- 状态量是 `Y = (Gamma, m, U, R)`。
- 时间推进采用 `GRKT4`，本质是带误差控制的 RK4 子步细分。
- 支持三种动力学闭合：
  - `index_dyn = 1`：Huang 模型
  - `index_dyn = 2`：Pe'er 模型
  - `index_dyn = 3`：Zhang 模型

环境与额外物理：

- `A_star > 0` 时走风环境近似，密度按 `n \propto R^{-2}`。
- `A_star <= 0` 时走 ISM，并允许在 `R_tr` 附近加入对数高斯密度跃迁。
- `E_inj_t1/E_inj_t2/l_inj_0/q_inj` 控制注能项。
- `R0` 用于近端半径截断设定。

### 4.2 电子分布模块

文件：

- `src/Electron/FS_electron_fullhide.f90`
- `src/Electron/FS_electron_weno5.f90`
- `src/Electron/FS_electron_t2g1.f90`
- `src/Electron/FS_electron_t2g2.f90`

当前默认运行路径是 `fs_electron_fullhide`。

该模块做了两件核心事情：

1. 根据 `gamma_m`、`gamma_c`、`gamma_max` 初始化电子分布。
2. 在对数 `gamma_e` 网格上推进电子连续性方程。

物理上包含：

- 注入谱
- 同步冷却
- IC 冷却
- SSA 对电子冷却项的反馈
- 输出同步辐射特征频率 `nu_m`、`nu_c`、`nu_a`

当前 `index_y` 控制 Compton-Y 的处理：

- `1`：数值 IC
- `2`：Nakar 近似
- `3`：Fan 近似

当前 `index_syn_integr` 控制同步辐射积分方式：

- `1`：普通实现
- `2`：Simpson 积分实现

### 4.3 SSC 模块

文件：

- `src/Radiation/SSC_spec.f90`

特点：

- 在频率和电子能量两层上做复合 Simpson 积分。
- 散射核中显式包含 Klein-Nishina 修正变量 `q` 和 `q_gamma`。
- 输出：
  - `P_SSC_spec`
  - `seed_SSC`

这一步给出 SSC 辐射功率和后续吸收所需的 SSC seed field。

### 4.4 `gamma-gamma` 湮灭模块

文件：

- `src/Radiation/Annihilation.f90`

物理作用：

- 对高能辐射加入 `gamma-gamma -> e+e-` 吸收。
- 使用角向积分与目标光子场积分得到光深 `tau`。
- 最终返回吸收因子

```text
absorption = (1 - exp(-tau)) / tau
```

当前实现同时使用同步 seed 和 SSC seed。

### 4.5 观测者系插值

文件：

- `src/Interpolation/SED_interpolation.f90`
- `src/Interpolation/SED_interpolation_structured.f90`

物理作用：

- 将共动系谱通过 EATS 和 Doppler boosting 映射到观测者系。
- 在半径和频率两个方向上都采用对数插值。
- 几何项由 `theta`、`phi` 网格离散完成。

当前 Python 主链默认调用：

- `Interpolation.sed_interpolation`

### 4.6 光变与能谱后处理

多波段光变在 `asgard_solver.py` 中由一组固定观测频率构造：

- XRT：`0.5-10 keV`，先做频率积分得到能流
- 光学：`g/r/i/z`
- 射电：`9/5.5/3 GHz`
- 高能：`1 GeV` 与 `1 TeV`

当 `spectrum_output.enabled=True` 时：

- 额外构造一组稠密观测频率网格
- 直接计算 `spectrum_fnu(freq, time)`
- 供 `plot_spectrum` 后处理调用

## 5. 目前重构后的软件结构

这次重构的原则是“接口改造优先，不碰数值核公式”。

实际分层如下：

1. 配置层：`asgard_models.py`
2. 运行时绑定层：`asgard_runtime.py`
3. 求解编排层：`asgard_solver.py`
4. 绘图层：`asgard_plot.py`
5. 兼容门面层：`mergered.py`
6. 拟合脚本层：`hand_my.py`、`MCMC.py`、`pymultinest_demo.py`

这样做的直接结果是：

- 拟合、示例、测试都走同一个 `run_fit(config)`；
- 新增能谱输出不需要再在旧 `fit(**kwargs)` 上堆参数；
- 主链不再把调试绘图和求解流程混在一起。

## 6. 构建与运行

### 唯一真实构建入口

- `build_extensions.py`

它负责用 `numpy.f2py` 编译以下扩展：

- `Constants`
- `Dynamics_forward`
- `Dynamics_reverse`
- `FS_electron_weno5`
- `FS_electron_fullhide`
- `FS_electron_t2g1`
- `SED_interpolation`
- `SED_interpolation_structured`
- `Annihilation`
- `Seed_reverse`
- `SSC_spec`

### Linux 包装入口

- `install.sh`

当前仅作为 Linux 上调用 `build_extensions.py` 的薄包装，不再维护第二套构建逻辑。

### Windows 说明

- `src/__init__.py` 保留了 `ASGARD_MINGW_BIN` 机制。
- 在 Windows 上如果需要显式指定 MinGW DLL 路径，可设置：

```powershell
$env:ASGARD_MINGW_BIN="C:\msys64\mingw64\bin"
```

### 常用命令

编译：

```powershell
python build_extensions.py
```

运行示例：

```powershell
python hand_my.py
```

回归测试：

```powershell
python tests\regression_check.py
```

## 7. 已完成修复

### `Cal_ebl.py`

文件：

- `src/Radiation/Cal_ebl.py`

已修复：

- `Path` 导入
- 变量名错误
- 旧插值实现替换
- 保持 `cal_ebl(z, v_obs, model)` 接口

### `cal_chi2_spectrum.py`

已修复：

- 编码/语法损坏
- 多列数据解析
- 异常处理路径

### `cal_chi2_light_curve.py`

已修复：

- 观测时间单位不一致问题
- 当前明确将光变数据时间列按“天”转成“秒”

这一步是必要修正，不是保护性绕过。因为模型输出 `t_obs_s` 的单位就是秒。

### `src/Electron/calling_modules.f90`

已修复：

- `get_nu_a` 在第一档频率就满足 `tau <= 1` 时，`Tau_high` 未定义的问题
- 现在若吸收转折频率落在搜索下界以下，返回搜索下界 `1e4 Hz`

这不是保护性分支，而是把原来缺失的搜索边界条件补完整。

## 8. 编译与回归现状

截至本文件更新时，已完成以下检查：

- `python tests\ssc_aux_route_check.py` 通过
- `python tests\ssc_aux_grid_scheme_check.py` 通过
- `python tests\physical_closure_check.py` 通过
- `src/Radiation/SSC_spec.f90` 已完成 line truncation 复查

当前默认策略：

- 只做与当前改动直接相关的最少编译和最小回归测试。
- 不把无关全量重编译或全量回归作为默认检查动作。
- Windows 下若已有 `.pyd` 被当前 Python 进程占用，`python build_extensions.py` 可能因目标文件无法覆盖而失败；这类锁文件问题不应被误判为当前物理或接口改动本身出错。

## 9. Fortran 代码检查要求

当前仓库明确要求：

- Fortran 代码必须做编译检查。
- 必须检查 line truncation。

当前状态：

- `src/Interpolation/SED_interpolation.f90`
- `src/Interpolation/SED_interpolation_structured.f90`
- `src/Electron/FS_electron_t2g2.f90`
- `src/Radiation/SSC_spec.f90`

都已做超长行拆分处理，不改公式。

其中：

- `FS_electron_t2g2.f90` 已完成独立编译检查；
- 该文件目前不是默认运行路径，但代码本身已恢复到可编译状态。

## 10. 约束与后续开发原则

1. 第一阶段不改 Fortran 数值核公式与离散方案。
2. 新功能优先放在 Python 接口层、编排层和后处理层。
3. 非数值模拟代码禁止添加数值保护。
4. 若物理结果异常，应先定位单位、网格、插值或接口错误，再考虑文献对照。
5. 默认不生成稠密能谱；只有显式启用 `spectrum_output.enabled=True` 才计算。
6. 当前默认电子演化路径是 `FS_electron_fullhide`，不是 `FS_electron_t2g2`。
7. 每完成一次重要迭代，都要同步更新 `AGENTS.md`，删改已过时条目，并补上新的项目事实、约束和当前基准。
8. 每完成一次已被接受的修改后，都要及时执行一次 git commit；提交时只包含当前修改直接相关的文件，不混入无关脏树改动。

## 11. 最小调用示例

```python
from mergered import FitConfig, SpectrumOutputConfig, run_fit, plot_light_curve, plot_spectrum

config = FitConfig(
    z=0.4,
    eta_0=1.0e2,
    epsilon_e=1.0e-1,
    epsilon_b=1.0e-3,
    p=2.5,
    opening_angle_jet=1.0e-1,
    theta_v=0.0,
    f_e=1.0e-1,
    e_iso=1.0e53,
    d_ne=1.0e-1,
    a_star=-1.0,
    spectrum_output=SpectrumOutputConfig(enabled=True, num_nu_obs=200),
)

result = run_fit(config)
plot_light_curve(result)
plot_spectrum(result, times_s=[1e3, 1e4, 1e5], quantity="nufnu")
```

## 12. 当前最重要的事实

- 主入口已经统一到 `run_fit(config)`。
- `weno5` 已经在新 API 中真正接通。
- 能谱数据管线已经补齐。
- 能谱绘图接口已经可用。
- 现有重构没有改动 Fortran 数值物理公式。
- 反向激波参数接口已经接入 `ReverseShockConfig`，FS/RS 特征频率和公开示例链路已经打通。
- `slc1_mmg2` 路径下的高层 SSC 主链已经从直接 `ssc_spec_nonuniform(...)` 切换到 auxiliary uniform `log(gamma)` 网格重建后再调用旧 `ssc_spec(...)`。
- `ssc_spec_nonuniform(...)` 现在保留为诊断/对照工具，不再作为 `slc1_mmg2` 高层主链默认路径。
- `tests/ssc_aux_grid_scheme_check.py` 现在除了比较共动系 SSC 谱外，还输出观测量侧前向 SSC 在 `X-ray/TeV` 的光变平滑性，以及 `slc1_mmg2` 相对 `fullhide/slc1` 的偏差，用于直接检查 aux-grid 主链在观测量侧是否连续平滑。
- `slc1_mmg2` 的前沿保峰问题当前已进入数值层修正：优先改高能前沿定位与高能锚带，而不是在非数值层补保护；本轮 `electron_common.f90` 已做 line truncation 复查，并用手工 `gfortran -c` 验证 `electron_common/calling_modules/FS_electron_slc1/FS_electron_fullhide` 编译链通过。当前 `python build_extensions.py --module FS_electron_slc1` 在本环境仍会被 `f2py/meson` 的 `electron_common.mod` 依赖顺序问题卡住，这一项不能误判成新数值改动本身未通过编译。
- 当前后续工作重点应放在 SSC/RS/cross-zone IC 的物理一致性核查与基准收束，而不是继续堆新接口。

## 13. 2026-04 Runtime Update

- 新增 `ReverseShockConfig`，作为 `FitConfig.reverse_shock` 子配置。
- 当前反向激波链路已经接入新 API：
  - `Dynamics.dynamics_reverse`
  - `Radiation.seed_reverse`
- 当前实现方式：
  - 若启用反向激波，主动力学使用 `dynamics_reverse` 返回的 `R_Tobs/R_Gamma/R`
  - `slc1_mmg2` 的 FS-SSC 高层主链当前走 auxiliary-grid 重建后复用 `ssc_spec(...)`
  - 总观测谱当前已接入 `FS synchrotron + FS-SSC + RS synchrotron + RS-SSC + cross-zone IC`
  - `gamma-gamma` 吸收当前使用 `seed_syn,FS + seed_syn,RS + seed_ssc,FS + seed_ssc,RS + seed_cross-zone,IC`
- 新增 `asgard_observables.py`，统一管理波段、频率表和多波段后处理。
- `FitResult` 现在显式区分两套时间轴：
  - `t_obs_s`：观测光变与能谱输出时间轴。
  - `characteristic_time_s`：`nu_m/nu_c/nu_a` 与 `rs_nu_m/rs_nu_c/rs_nu_a` 的特征频率时间轴。
- 新增 `plot_characteristic_frequencies(result, include_reverse=True, ...)`，按 `characteristic_time_s` 绘制 FS/RS 特征频率，避免将内部动力学网格误当作观测输出网格。
- 新增 `hand_reverse.py`，作为反向激波 API 与 FS/RS 特征频率绘图的公开示例脚本。
- 新增 `PhysicalSolution`，作为主链物理求解与观测后处理之间的中间结果对象。
- 新增 `asgard_postprocess.py`，统一承载：
  - 观测频率投影
  - 多波段光变聚合
  - `redchi` 计算
  - 稠密能谱后处理
- `asgard_solver.py` 不再在主链中传递长位置元组，而是返回 `PhysicalSolution` 后交给后处理层。
- `tests/api_contract_check.py` 已加入回归集合，检查：
  - `characteristic_time_s` 与 `nu_m/nu_c/nu_a`、`rs_nu_*` 的 shape 一致
  - `spectrum_fnu.shape == (N_nu, N_tobs)`
  - `plot_characteristic_frequencies` 与 `plot_spectrum` 可正常输出文件
- `hand_my.py`、`hand_reverse.py`、`MCMC.py`、`pymultinest_demo.py` 已改为 `main()` 入口，导入时不再执行计算。
- 新增 `asgard_inference.py`，统一 `MCMC.py` 与 `pymultinest_demo.py` 的 `FitConfig` 装配默认值，避免推断脚本默认网格和物理开关继续漂移。
- 新增 `asgard_legacy.py`，承接旧 `fit(**kwargs)` 的参数映射；`mergered.py` 现在只保留门面导出与弃用提示。
- 新增 `tests/legacy_api_check.py`，检查旧接口仍可返回有限 `redchi`，且会发出 `DeprecationWarning`。
- 新增 `SimulationSetup` 与 `asgard_setup.py`，统一管理：
  - `luminosity_distance_cm`
  - `boundary`
  - `seed_frequency_hz`
  - `observer_time_s`
- `asgard_postprocess.py` 现在只接收 `SimulationSetup + PhysicalSolution`，不再在接口上传递多组裸数组。
- 新增 `asgard_presets.py`，统一基线、能谱示例和反向激波示例配置。
- 新增 `tests/import_smoke_check.py`，检查 `hand_my.py` 与 `hand_reverse.py` 导入时无绘图副作用。
- 当前总辐射已接入 `RS-SSC`：
  - `L_tot = L_syn,FS + L_syn,RS + L_ssc,FS + L_ssc,RS`
  - `gamma-gamma` 吸收种子场当前使用 `seed_syn,FS + seed_syn,RS + seed_ssc,FS + seed_ssc,RS`
- 新增 `asgard_coupling.py`，当前区间耦合 IC 实现采用：
  - `FS` 与 `RS` 共用接触不连续面速度 `Gamma_2`
  - 用 `m2/m3` 反推两区共动厚度 `Δ'_2` 与 `Δ'_3`
  - 用 `0.5(Δ'_2+Δ'_3)/c` 构造中心到中心的共动传播延迟
  - 在共动时间轴上对外区同步种子场做延迟插值
  - 再复用 `ssc_spec` 计算
    - `FS` 电子散射 `RS` 光子
    - `RS` 电子散射 `FS` 光子
- 当前区间耦合 IC 的角分布仍采用半球因子 `1/2`，因为现有 `seed_syn` 是球壳平均量，不包含更细的角分辨信息。
- 为支持双区几何，`Dynamics_forward.f90` 现在导出 `m2`，`Dynamics_reverse.f90` 现在导出 `m2/m3`。
- 新增 `tests/coupled_ic_check.py`，检查：
  - 双区几何量有限且非负
  - 共动时间单调
  - 两个区间耦合 IC 分量有限且非零

## 14. 后续重构路线

1. 推断脚本统一
   - 抽出 `MCMC.py` 与 `pymultinest_demo.py` 的公共参数装配逻辑
   - 避免不同推断脚本继续分叉出不同默认物理设置
2. 结果对象扩充
   - 评估是否将 FS/RS 的更多中间物理量纳入只读结果对象
   - 目标是减少脚本层再次回探内部模块
3. SSC / 反向激波物理链路继续核查
   - `RS-SSC` 与区间耦合 `IC` 已接入，但当前重点不是继续盲目扩链，而是做物理口径、平滑性与参考解一致性核查
   - `slc1_mmg2` 的 aux-grid SSC 主链已经替代 direct nonuniform 主链，后续应继续围绕这条路径做基准收束
4. 测试继续前移
   - 将更多 API 契约和绘图输出检查固定到 `tests/`
   - 保持每轮重构后先做接口回归，再做物理曲线回归

## 15. 2026-04 MMG2 Front Update

- `mmg2` 激波前沿抹峰问题已经进入数值层修正，不再接受非数值绕过。
- 本轮已修改 `src/Electron/electron_common.f90`：
  - `electron_find_high_energy_front(...)` 改为优先跟踪峰值后的最大下降梯度位置。
  - `electron_anchor_high_energy_edges(...)` 改为给高能前沿留更密的局部锚带，并增加尾部单元数。
- 新增 `tests/mmg2_front_sharpness_check.py`，输出 `fullhide/slc1/slc1_mmg2` 以及 `mmg2 public/work-grid` 的前沿指标。
- 本轮 Fortran 改动已完成 line truncation 复查，并通过手工 `gfortran -c` 链式编译检查。
- 当前 Windows + Python 3.12 + `numpy.f2py/meson` 对 `FS_electron_slc1` 仍存在 Fortran module 依赖排序问题：`calling_modules.f90` 会早于 `electron_common.mod` 被编译。
- 为保持最少必要编译与回归，本轮采用“手工最小 `f2py` wrapper + 手工链接”生成新的 `src/Electron/FS_electron_slc1.cp312-win_amd64.pyd` 运行时，只暴露当前主线所需的 `fs_electron_slc1/fs_electron_slc1_mmg2`。
- 当前主线验收文件：
  - `output/asgard_doc/mmg2_front_sharpness.json`
  - `output/asgard_doc/mmg2_front_old.npz`
  - `output/asgard_doc/mmg2_front_new.npz`
  - `output/asgard_doc/mmg2_front_observed_compare.png`

## 16. 2026-04 Order And Build Baseline

- `tests/order_convergence_check.py` 已从“只检查 `slc1_mmg2` 是否正阶”升级为正式主线诊断：
  - 电子谱链单独检查 `fullhide/slc1/slc1_mmg2`
  - 辐射谱链单独检查 `synchrotron/SSC/observer-side total` 以及 band flux、`nu_a`
  - 当前判据明确写死为：
    - 电子谱链有效阶数 `> 2`
    - 辐射谱链有效阶数 `> 2`
- `dynamics-forward` 仍单独保留为独立检查项，只要求有效阶数为正，不与电子谱链和辐射谱链混为一谈。
- `build_extensions.py` 当前对 `FS_electron_slc1` 已加入受控 fallback：
  - 先尝试原来的 `numpy.f2py -c` 直接构建
  - 若再次触发 Windows + Python 3.12 + meson 的 Fortran module 排序问题，则自动改走“按依赖顺序手工编译对象文件 + 生成最小 `pyf` + 再由 `f2py` 链接”的路径
  - fallback 仍只服务当前运行时所需的 `fs_electron_slc1/fs_electron_slc1_mmg2`
- 这轮主线目标已经明确收束为：
  - 先把阶数验证和 `FS_electron_slc1` 可复现构建入口固定下来
  - 再决定是否需要继续改 `slc1/slc1_mmg2` 的数值层以把电子谱链和辐射谱链都推到二阶以上
- 当前 `python tests/order_convergence_check.py` 已能稳定产出 `output/asgard_doc/order_convergence.json`，且第一轮真实结果表明：
  - 电子谱链当前远未整体达到二阶，最低项出现在 `slc1_mmg2-electron-support-low`
  - 辐射谱链中 `slc1_mmg2` 的 `sync/ssc/observer total/bands` 已有多项超过二阶，但 `nu_a` 仍未过线
  - `fullhide/slc1` 的 observer-side 与部分 SSC/SSA 项仍低于二阶，说明降阶不只在 moving-mesh
- 因而“算法阶数任务”当前状态应被视为：
  - 构建入口已固定
  - 诊断链已固定
  - 物理/数值改进尚未完成，不能宣称任务结束
- 当前 `plan.md` 作为主线任务清单；后续默认先以 `plan.md` 的任务顺序推进，再同步回写 `AGENTS.md`。

## 18. 2026-04 SLC1 Solver Note

- ?????????`CC` ? `fullhide` ??????????????? `CC` ??????????
- `slc1` ???????????? `conservative SL + IMEX`?

## 19. 2026-04 Charint Current Baseline

- `charint` 已接入运行时：
  - `FitConfig.electron_solver="charint"`
  - 主文件：
    - `src/Electron/FS_electron_charint.f90`
    - `src/Electron/electron_common.f90`
- 当前有效算法口径：
  - 固定网格 `x = log10(gamma)`
  - 特征变量 `u = 1/gamma`
  - `index_y=0`：
    - affine characteristic
    - 数值源项版本
    - 不再使用 exact-source 链
  - `index_y=1/2/3`：
    - piecewise-affine characteristic
  - 源项时间积分：
    - `4` 点 `Gauss-Legendre`
- 当前有效子步口径：
  - `charint_cfl_relax = 32`
  - `charint_substep_rtol_relax = 8`
  - `L1 = max(4,min(64,...))`
- 当前保留的主核优化：
  - `electron_prepare_conservative_remap_nonuniform(...)`
  - `electron_conservative_remap_nonuniform_prepared(...)`
  - `electron_ppm_prefix_sweep_nonuniform(...)`
  - affine / piecewise tracing 的 batch 路径
  - `index_y=0` 的 prepared-source 路径
- 当前不再保留的实验链：
  - `index_y=0` exact-source
  - 幂律段时间闭式

## 20. 2026-04 Charint Rules

- `charint` 理论阶数口径固定为整体二阶。
- 当前不再把以下量作为 `charint` 主线验收项：
  - `peak-gamma`
  - 低能前沿结构
  - `g_lo` 跳变
- 粗阶数诊断只保留：
  - `electron-spectrum`
  - `radiation-bands-flux`
  - `dynamics-forward-gamma/radius`
- benchmark 默认只比较：
  - `fullhide`
  - `slc1`
  - `charint`
- 非必要不跑重 benchmark。
- 涉及 `FS_electron_charint` 的编译与回归必须顺序执行，不能并行占用 `.pyd`。
- `charint` 后续 done-when 固定为：
  - 相对 `fullhide` 达到 `10x`
  - 相对 `slc1` 达到 `5x`
  - 精度不低于当前基线

## 21. 2026-04 Build Environment

- 固定编译环境：
  - `python = C:\Users\jia\AppData\Local\Programs\Python\Python312\python.exe`
  - `gfortran = C:\msys64\mingw64\bin\gfortran.exe`
  - `ASGARD_MINGW_BIN = C:\msys64\mingw64\bin`
  - `CC = C:\msys64\mingw64\bin\gcc.exe`
  - `CXX = C:\msys64\mingw64\bin\g++.exe`
  - `FC/F77/F90 = C:\msys64\mingw64\bin\gfortran.exe`
  - `AR = C:\msys64\mingw64\bin\gcc-ar.exe`
- 当前最小编译命令：
  - `python build_extensions.py --module FS_electron_charint --force`
- `FS_electron_charint` 与 `FS_electron_slc1` 当前都依赖 ordered-object fallback。
- 最近一次 `FS_electron_charint` 编译已通过，fallback 日志未见 `line truncation`。

## 22. 2026-04 Current Test Grids

- `tests/order_convergence_check.py`
  - 三组粗阶数网格固定为 `61 / 101 / 161`
- `tests/electron_exp_tail_spectrum_compare.py`
  - 当前谱图电子网格固定为统一 `121 / 121 / 121`
- `tests/electron_solver_comparison.py`
  - 当前 benchmark 电子网格固定为统一 `121 / 121 / 121`

## 23. 2026-04 Current Charint Facts

- 旧 `MC-limited` 线性 remap 核不是恢复较好结果的主因。
- 激进 `CFL/rtol/L1` 单独也不是恢复较好结果的主因。
- 把 `index_y=0` 从 exact-source 撤回到数值源项版本后，当前粗阶数结果回到：
  - `charint-electron-spectrum ≈ 0.8606`
  - `charint-radiation-bands-flux ≈ 0.7480`
- 当前最新谱图输出：
  - `output/benchmark_exp_tail/spectrum_compare.png`
  - `output/benchmark_exp_tail/spectrum_compare.pdf`

## 24. 2026-04 Git Backup Rule

- 每次发生重大代码更新后，必须执行一次 git 备份。
- “重大代码更新”至少包括：
  - Fortran 数值核主逻辑改动
  - 构建链路改动
  - 运行时绑定改动
  - 测试口径或 benchmark 口径改动
- git 备份要求：
  - 先更新 `AGENTS.md`
  - 再提交当前有效版本
  - 提交信息必须明确说明这次备份对应的主线变更



## 25. 2026-04 Repository Cleanup Baseline

- ???? `asgard_paths.py`??????
  - `ROOT`
  - `OUTPUT_ROOT`
  - `ASGARD_DOC_DIR`
  - `BENCHMARK_EXP_TAIL_DIR`
- ????????????????? `asgard_paths.py` ??????????????? `output/asgard_doc` ? `output/benchmark_exp_tail`?
- ???????????????
  - `doc.md`
  - `charint.md`
  - `version history.md`
  - `TODO lists.txt`
- `plan.md` ?????????????????????????????
- ????????????????????
  - ??? PDF ??
  - `f2py` ????? wrapper ???
  - ???? PDF ????
- ?????????????
  - ??
  - ??
  - ????
  - ????
  - ????
