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

- `ASGARD/__init__.py`
  - 当前包门面模块。
  - 导出 `Model`、`Fitter`、`FitResult`、`observe`、`Setups`、`Radiation` 等核心对象。

- `asgard_models.py`
  - 定义配置与结果数据结构。
  - 包含：
    - `SpectrumOutputConfig`
    - `FitConfig`
    - `FitResult`
    - 波段和频率常量

- `ASGARD/api.py`
  - 当前主入口编排层。
  - 负责：
    - 构建观测时间网格
    - 调用动力学、电子谱、辐射、插值模块
    - 组装多波段光变
    - 计算 `redchi`

- `asgard_postprocess.py`
  - 承担后处理与稠密能谱矩阵生成。

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

- `lc_spec_demo.py`
  - 最小运行示例。
  - 当前已迁移到新 API，并包含能谱输出调用。

- `mcmc_fit.py`
  - MCMC 采样入口，已迁移到新 API。

- `multinest_fit.py`
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

当前主链在 `ASGARD/api.py::observe` 中按下面顺序执行：

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

多波段光变在 `ASGARD/api.py::observe` 中由一组固定观测频率构造：

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
3. 观测与后处理层：`ASGARD/api.py`、`asgard_postprocess.py`
4. 绘图层：`asgard_plot.py`
5. 包门面层：`ASGARD/__init__.py`
6. 示例与拟合脚本层：`lc_spec_demo.py`、`mcmc_fit.py`、`multinest_fit.py`

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
- 所有 shell 命令统一从 `rtk` 进入。
- 对 `rtk pwsh/powershell -Command ...`：
  - 默认使用外层单引号包住整段 PowerShell 脚本。
  - 只要脚本里出现 `$env:`、`$_`、`$var` 等 PowerShell 变量，就不要再用外层双引号。
  - 若只是单条命令，优先写成 `rtk pwsh -Command '...'`，不要混用 `cmd` 风格转义。

```powershell
rtk pwsh -Command '$env:ASGARD_MINGW_BIN="C:\msys64\mingw64\bin"'
```

### 常用命令

编译：

```powershell
rtk pwsh -Command '& "C:\Users\jia\AppData\Local\Programs\Python\Python312\python.exe" build_extensions.py'
```

运行示例：

```powershell
rtk pwsh -Command '& "C:\Users\jia\AppData\Local\Programs\Python\Python312\python.exe" lc_spec_demo.py'
```

回归测试：

```powershell
rtk pwsh -Command '& "C:\Users\jia\AppData\Local\Programs\Python\Python312\python.exe" tests\regression_check.py'
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
- 现在若吸收转折频率落在搜索下界以下，会继续向更低频扩展括区，不再把结果硬钉在 `1e4 Hz`
- `FS_electron_charint` 的最后一个时间点现在会显式补写 `V_m/V_c/V_a`，避免尾端保持 0 导致 SSA 频率演化断裂

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
- 当前 FS/RS SSC 主路径已经改成优先直接调用 `Radiation.ssc_spec_nonuniform(...)`，不再先压成固定 `64` 点均匀辅助电子网格；旧的均匀重映射仅作为兼容回退。
- 当前 SSA 路径里 `electron_syn_cell_adaptive(...)` 已改成真正递归，但默认只保留浅层自适应；`get_nu_a(...)` 额外传入更深的专用深度，以优先保证 SSA 曲线平滑而不把普通同步辐射积分拖慢。
- 当前前向 SSC 还额外使用按 `nu_a/nu_m/nu_c/nu_M` 自适应细化的内部 seed 频率网格，再把结果插回原始频率轴，目标是直接改善 `compare_spectrum` 的 SSC 分量谱型。
- `slc1_mmg2` 的前沿保峰问题当前已进入数值层修正：优先改高能前沿定位与高能锚带，而不是在非数值层补保护；本轮 `electron_common.f90` 已做 line truncation 复查，并用手工 `gfortran -c` 验证 `electron_common/calling_modules/FS_electron_slc1/FS_electron_fullhide` 编译链通过。当前 `python build_extensions.py --module FS_electron_slc1` 在本环境仍会被 `f2py/meson` 的 `electron_common.mod` 依赖顺序问题卡住，这一项不能误判成新数值改动本身未通过编译。
- 当前后续工作重点应放在 SSC/RS/cross-zone IC 的物理一致性核查与基准收束，而不是继续堆新接口。
- `Dynamics.dynamics_forward` 第四个返回量在 Python 侧历史上被命名为 `swept_mass_g`，但这条链上实际承载的是 `R_m`/累计粒子数，不应再除以 `m_p` 后当作 `N_p`。
- VegasAfterglow 的微物理入参名是 `xi_e`，ASGARD 对应的是 `xi_N`；这两边的比较基准必须显式设成同一数值，否则 `gamma_m`/`nu_m` 会被默认值漂偏。
- `tests/vegas_afterglow_comparison.py` 里的反向激波比较图现在分成上下两幅：上图只画 forward shock，下图只画 reverse shock，避免把两个分量混在同一坐标轴里造成误读。
- 当前对特征频率的交叉检查显示：ASGARD 和 VegasAfterglow 各自的 `nu_m/nu_c` 与它们自己的 `Gamma/B/gamma_m/gamma_c` 关系都是自洽的；剩余差异更像是上游电子冷却闭合或参考系映射差异，而不是最后一层绘图坐标问题。
- `compare_reverse_shock_lc.png` 的上下两幅图现在按各自数据的 1% 到 99% 百分位自动设 y 轴范围，目的是让低流量段更可读，不再被另一分量的动态范围压扁。
- `compare_reverse_shock_lc.png` 的 y 轴现在改为按每条曲线各自的 5% 到 95% 百分位定界，避免极低尾部把整幅 reverse-shock 图拖到不可读的尺度。
- `compare_reverse_shock_lc.png` 的下方 reverse-shock 子图现在直接按 reverse-shock 峰值定标，y 轴下限取峰值的 `1e-5`，让主导结构优先可见。
- 当前对齐后的特征频率诊断表明：即使在同一时刻对齐后，ASGARD 与 VegasAfterglow 的 `nu_m` 和 `nu_c` 仍分别保留约 `28` 和 `588` 的中位数倍率差，因此当前问题不在最后一层绘图坐标，而在上游冷却/参考系映射口径。
- `shock_quantities` 中的 `N_p` 比较现在对 VegasAfterglow 侧乘了 `4π`，按总量口径与 ASGARD 对齐。
- `shock_quantities` 和 `photon_quantities` 中的 VegasAfterglow 频率类量现在先按 `Doppler/(1+z)` 转到 lab/observer frame 再与 ASGARD 对比，因为 Vegas 的内部传递口径是共动系。
- 进一步核对发现，VegasAfterglow 的 `nu_a` 原始输出满足 `nu_a ≈ 4.2e6 * B * gamma_a^2`，即它保留的是共动系 SSA 频率；ASGARD 的 `nu_a` 则已经转成观测者系。
- 对比脚本里 `nu_m/nu_c` 继续按 frame 转换处理，但 `nu_a` 先保留 Vegas 原始值直接对比 ASGARD；从当前数据看，`nu_a` 的最小差异出现在不额外做 Vegas frame 转换时。

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
- 新增 `PhysicalSolution`，作为主链物理求解与观测后处理之间的中间结果对象。
- 新增 `asgard_postprocess.py`，统一承载：
  - 观测频率投影
  - 多波段光变聚合
  - `redchi` 计算
  - 稠密能谱后处理
- `ASGARD/api.py` 现在直接返回 `PhysicalSolution` 后交给后处理层。
- `lc_spec_demo.py`、`mcmc_fit.py`、`multinest_fit.py` 已改为 `main()` 入口，导入时不再执行计算。
- 新增 `asgard_fit.py`，统一 `mcmc_fit.py` 与 `multinest_fit.py` 的 `FitConfig` 装配默认值，避免推断脚本默认网格和物理开关继续漂移。
- `ASGARD/__init__.py` 现在承担包门面导出。
- 新增 `tests/legacy_api_check.py`，检查旧接口仍可返回有限 `redchi`，且会发出 `DeprecationWarning`。
- 新增 `SimulationSetup` 与 `asgard_setup.py`，统一管理：
  - `luminosity_distance_cm`
  - `boundary`
  - `seed_frequency_hz`
  - `observer_time_s`
- `asgard_postprocess.py` 现在只接收 `SimulationSetup + PhysicalSolution`，不再在接口上传递多组裸数组。
- 新增 `asgard_presets.py`，统一基线、能谱示例和反向激波示例配置。
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
   - 抽出 `mcmc_fit.py` 与 `multinest_fit.py` 的公共参数装配逻辑
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
- `tests/sed_electron_compare.py`
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

## 26. 2026-04 Inference Fastpath Update

- 当前已新增生产推断 fast-path，目标是压缩 Python 调用链开销，而不是改动 Fortran 物理核。
- `asgard_state.py` 当前已新增 observer-side `mode` 分支：
  - `full_components`
  - `total_only`
- 生产推断当前默认走 `total_only`：
  - 只对 `component_spectra.total` 调一次 `Interpolation.sed_interpolation`
  - 不再为 `fwd_sync/fwd_ssc/rev_sync/rev_ssc/cross_ic` 分量逐个做 observer-side 插值
- `asgard_fit.py` 当前已新增：
  - `compile_problem(...)`
  - `eval_loglike(...)`
- `eval_fit(config)` 当前已切到 compiled `FitConfig` light-curve fast-path：
  - `solve_state_from_setup(...)`
  - `project_flux_grid(..., mode="total_only")`
  - `combine_multiband_flux(...)`
  - `compute_light_curve_redchi(...)`
- `ASGARD.api.Fitter.loglike(...)` 当前已切到 compiled observation-block 路径：
  - 不再逐次 `deepcopy(model)` 后按 dataset 临时拼调用
  - 当前会预编译 dataset 请求块和参数绑定
  - 参数评估时对模型对象做就地参数覆写，评估后恢复原值
- 对通用 `Fitter` 路径，当前没有强行把所有 dataset 合并到单一 solve grid：
  - 这是有意保留的
  - 原因是不同时间覆盖范围的 dataset 若盲目共用一个 solve grid，会引入可见 observer-side 偏差
  - 当前实现改为“按请求块复用”，优先保证结果与旧路径一致
- `mcmc_fit.py` 与 `multinest_fit.py` 当前已改为：
  - 进程内只初始化一次 compiled problem
  - 复用单个 `FitConfig` 模板对象
  - 每个参数点评估时仅更新必要字段，再调用 compiled loglike
- 当前 `tests/` 目录只保留面向新用户的 benchmark / plot 入口：
  - `tests/readme_smoke_bench.py`
  - `tests/slc1_vs_fullhide_bench.py`
  - `tests/wind_vs_fullhide_bench.py`
  - `tests/doc_figures.py`
  - `tests/ic_doc_plots.py`
  - `tests/sed_electron_compare.py`
- 旧的深层验证 / inference profile 脚本已经移除，不再作为日常入口：
  - `tests/asgard_comprehensive_validation.py`
  - `tests/asgard_summary_plots.py`
  - `tests/asgard_inference_profile.py`
  - `tests/inference_fastpath_check.py`

## 26. 2026-04 Inference Fastpath Update

- 生产推断 fast-path 已接入，目标是减少 Python 调用链重复开销，不改 Fortran 物理核公式。
- `asgard_state.py` 新增观测模式分流：
  - `mode="full_components"`：保留旧行为，逐分量 observer 插值。
  - `mode="total_only"`：只对 `component_spectra.total` 做一次 observer 插值。
- `asgard_fit.py` 已新增 compiled 推断入口：
  - `compile_problem(...)`
  - `eval_loglike(...)`
- `eval_fit(config)` 已改为走 compiled config fast-path。
- `ASGARD/api.py` 中：
  - `Fitter.loglike(...)` 改为内部转调 compiled 推断入口。
  - `Fitter` 在 `add_flux_density/add_spectrum/add_flux/fit(param_defs|resolution)` 时会失效并重建 compiled 缓存。
  - 新增 `_total_matrix(...)` 及 direct/patch 的 total-only 观测求值路径，避免推断链默认重复分量插值。
- `mcmc_fit.py` 与 `multinest_fit.py` 已改为：
  - 启动时构建一次 compiled context
  - 每个参数点仅就地更新配置并调用 `eval_loglike(...)`
  - 不再每点评估时重新构造完整求解入口
- 当前推断 fast-path 仍保留在 `asgard_fit.py` / `ASGARD.api` 中，供主入口复用；但日常脚本层不再维护单独的 inference profile/check 入口。
- 默认烟测 `tests/readme_smoke_bench.py` 现在覆盖 `fullhide/slc1/charint` 三个 solver，且默认走 `quick` 小网格。

## 27. 2026-04 Plan Compression

- 当前 `plan.md` 已压缩为后续执行准绳，优先级高于历史性路线描述。
- 后续默认只做三类工作：
  - 维持 compiled 推断 fast-path
  - 维持新用户最小入口
  - 维持 benchmark / plot 入口
- 明确不再推进的内容：
  - 新的 inference profile/check 脚本
  - 旧验证脚本树
  - 生产推断的逐分量 observer 路径
  - 生产推断的逐点评估 `deepcopy(model)` 路径
- 若结构要变，先同步 `AGENTS.md` 和 `plan.md`，再动代码。

## 28. 2026-04 Naming Cleanup

- 当前 Python 主链已完成一轮按真实职责的重命名，优先覆盖：
  - 推断模块：`asgard_fit.py`
  - 状态/观测模块：`asgard_state.py`
  - 最小 demo：`lc_spec_demo.py`
  - 采样脚本：`mcmc_fit.py`、`multinest_fit.py`
  - benchmark / plot 脚本：`tests/readme_smoke_bench.py`、`tests/slc1_vs_fullhide_bench.py`、`tests/wind_vs_fullhide_bench.py`、`tests/doc_figures.py`、`tests/ic_doc_plots.py`、`tests/sed_electron_compare.py`
- 当前主链内部已压短的高频函数包括：
  - `compile_problem(...)`
  - `eval_loglike(...)`
  - `solve_state_from_setup(...)`
  - `project_flux_grid(...)`
  - `_solve_patch_state(...)`
  - `_total_matrix(...)`
- 后续命名调整默认规则：
  - 先按真实职责改名，再考虑压短
  - 不为“看起来统一”制造新包装层
  - 公开物理对象名谨慎处理，优先先改内部高频函数

## 29. 2026-04 Charint Source Remap Update

- 本轮 `charint` 数值核已继续下沉到 source 积分链，不再只停留在 Python 侧整理。
- 当前改动范围仅限 prepared-source 路径：
  - `src/Electron/FS_electron_charint.f90`
  - `src/Electron/electron_common.f90`
- 已做的数值改动：
  - `dF1_shape` 的 prepared remap 不再走普通 PPM 接口值 `q_left/q_right/prefix` 组合。
  - 当前改为：
    - 先对正定源项形状做 `electron_log_prefix_nonuniform(...)`
    - 再在 prepared-source characteristic step 里用 `electron_conservative_remap_log_nonuniform_prepared(...)` 做回溯边界积分
  - `dN_x` 主输运 remap 仍保留原有 conservative remap，不混改 transport 链。
- 本轮目的不是加保护，而是让指数 cutoff 源项的重映射与其正定、陡尾结构更匹配，优先缩小 source 积分带来的尾部与 cutoff 离散误差。
- 本轮编译检查结果：
  - 固定环境下 `python build_extensions.py --module FS_electron_charint --force` 通过
  - 仍先触发 `f2py/meson` 失败，再由 ordered-object fallback 成功生成 `.pyd`
  - 改动文件已做长行扫描，当前未留下 `>132` 字符的新增长行
- 本轮最小 smoke：
  - 通过 `import src` 加载 MinGW DLL 路径后，`FS_electron_charint` 可正常导入
  - `fs_electron_charint`
  - `fs_electron_charint_affine_step_test`
- 后续 `charint` 主线继续按这个顺序推进：
  - 先看 source-shape log-remap 对谱尾和平滑性的实际改善
  - 再决定是否继续处理 PPM state remap / cutoff edge builder
  - 不在非数值层补保护

## 30. 2026-04 Src Cleanup Baseline

- 本轮已对 `src/` 做第一轮真实清理，目标是把源码树收回到“新用户能看懂、构建链能自洽”的最小集合。
- 已删除的明显无效/历史遗留源码：
  - `src/plot_module.py`
  - `src/plotdata.py`
  - `src/Electron/FS_electron_upwind.f90`
  - `src/Electron/FS_electron_t2g2.f90`
  - `src/utils/(test) read_sigma_KN_IC.f90`
  - `src/utils/adaptive_1st_resampling_mod.f90`
  - `src/utils/Dynamics_CIP.f90`
  - `src/utils/Dynamics_Analy.f90`
- 已删除的非源码噪声：
  - `src/.idea/`
  - `src/**/__pycache__/`
- `src/Electron/__init__.py` 已去掉对不存在 `v2/v3` 二进制变体的回退导入，只保留当前真实存在且仍被构建链支持的 solver：
  - `fullhide`
  - `slc1`
  - `charint`
  - `weno5`
  - `t2g1`
- `asgard_runtime.py` 已同步去掉对 `*_v2/*_v3` 旧模块名的导入回退。
- `src` 包和各子包 `__init__.py` 已清掉乱码 docstring 和过时注释，改成最小职责说明。
- `adaptive_2nd_resampling_mod.f90` 已从 `src/utils/` 迁到 `src/Electron/adaptive_resampling_mod.f90`：
  - 文件名现在与模块名 `adaptive_resampling_mod` 一致
  - `build_extensions.py` 已同步改到新的相对路径
- 当前电子构建链的 ordered fallback 已从仅 `slc1/charint` 扩展到全部保留 solver：
  - `FS_electron_weno5`
  - `FS_electron_slc1`
  - `FS_electron_charint`
  - `FS_electron_fullhide`
  - `FS_electron_t2g1`
- 本轮最小验证：
  - `python -m py_compile build_extensions.py asgard_runtime.py src/.../__init__.py` 通过
  - `python build_extensions.py --module FS_electron_charint --force` 通过
  - `python build_extensions.py --module FS_electron_fullhide --force` 通过
  - 当前 `src/**/*.f90` 的 `>132` 字符长行扫描未发现残留项
- 后续对 `src/` 的规则：
  - 不再把 IDE 缓存、`__pycache__`、历史脚本塞回源码树
  - 不再在主入口保留指向不存在二进制备份名的导入回退
  - 性能优化优先从当前保留 solver 和公共电子/辐射核下手，不回头抢救已删除旧实现

## 31. 2026-04 Resampling Hotspot Update

- 本轮已开始进入 `src/` 的性能优化阶段，第一刀落在公共电子重采样模块：
  - `src/Electron/adaptive_resampling_mod.f90`
- 当前已实施的纯算法级提速：
  - `adaptive_resampling_log(...)` 中生成 `indices` 时，不再对每个目标点从 `idx=1` 重新扫描累积权重数组 `c`
  - 现改为单调递增的一次 sweep，时间复杂度从近似 `O(m*n)` 收到 `O(m+n)`
  - `unique_sorted(...)` 不再做二次复杂度的重复检查
  - 现利用 `indices` 本身的单调性，改成线性去重，不再额外排序
- 这些改动只去掉重复工作，不改采样准则、不改物理公式，也不改重采样输出口径。
- 本轮最小验证：
  - `python build_extensions.py --module FS_electron_charint --force` 通过
  - `python build_extensions.py --module FS_electron_fullhide --force` 通过
  - 两者仍表现为 `meson` 直编失败、ordered fallback 成功，这是当前环境事实，不是本轮优化引入的新问题
- 当前下一优先级性能点：
  - `FS_electron_charint.f90` 的 adaptive substep 回退链
  - `electron_prepare_radiation_spectrum(...)` 的前后沿扫描与重采样触发频率

## 32. 2026-04 Charint Adaptive Step Update

- 本轮已继续优化 `src/Electron/FS_electron_charint.f90` 的 adaptive substep 链，但只改 uniform-density 分支。
- 当前改动：
  - 在 `is_uniform_density` 的 adaptive while-loop 里，先根据现有接受判据
    - `step_error = dR_try / R_mid`
    - `step_error <= charint_substep_rtol_relax * substep_rtol`
  - 解析地构造一个 `dR_cap`
  - 对 `dR_try` 先做上界裁剪，再进入原来的接受/回退逻辑
- 这一步的作用是：
  - 减少“先试一步、再因纯几何误差超标而回退”的空转
  - 不改接受判据本身
  - 不改非均匀密度分支
  - 不改物理公式
- 本轮最小验证：
  - `python build_extensions.py --module FS_electron_charint --force` 通过
  - `import FS_electron_charint` smoke 通过
- 当前下一步仍是：
  - 看 nonuniform-density adaptive 分支是否也有可提前裁步的纯调度开销
  - 或继续下探 `electron_prepare_radiation_spectrum(...)` 的触发频率与扫描次数

## 33. 2026-04 Charint And Radiation Prep Follow-up

- 本轮已把上一轮的两个后续性能点都实际推进：
  - `src/Electron/FS_electron_charint.f90`
  - `src/Electron/electron_common.f90::electron_prepare_radiation_spectrum(...)`
- `FS_electron_charint.f90`
  - adaptive while-loop 现在在进入 uniform / nonuniform 分支前统一先做几何上界裁步
  - nonuniform-density 分支在步长被拒绝时，不再固定 `0.5*dR_try`
  - 现改为按 `step_limit / step_error` 比例直接缩步，再夹到 `dR_min`
  - 目的仍是减少纯调度性空转，不改接受判据本身
- `electron_prepare_radiation_spectrum(...)`
  - 不再一进子程序就无条件整段复制 `gam_e/dN_gam_e`
  - 先用双端扫描找 `first_pos/last_pos`
  - 只有在确认“不需要重采样”时才走整段保留路径
  - 真正进入重采样路径后，输出数组再单独写入并把尾部清零
  - 这样减少了高频调用下的无效整段拷贝
- 本轮最小验证：
  - `python build_extensions.py --module FS_electron_charint --force` 通过
  - `python build_extensions.py --module FS_electron_fullhide --force` 通过
  - `src/Electron/FS_electron_charint.f90` 与 `src/Electron/electron_common.f90` 的长行扫描当前为空
- 当前方向确认：
  - `plan.md` 的主线仍应理解为三件事
    - 维持 compiled fast-path
    - 维持新用户最小入口
    - 维持 benchmark / plot 入口
  - 当前这轮继续压电子核纯开销，仍符合这个方向，没有回到旧接口树或旧验证树

## 34. 2026-04 Test Output Hygiene

- 当前测试入口已经改成“边跑边报”：
  - `tests/readme_smoke_bench.py`
  - `tests/slc1_vs_fullhide_bench.py`
  - `tests/wind_vs_fullhide_bench.py`
  - `tests/doc_figures.py`
  - `tests/ic_doc_plots.py`
  - `tests/sed_electron_compare.py`
- 默认要求：
  - 每个 case 跑完立刻打印状态
  - 不再自动写 JSON / NPZ 旁路文件
  - 只保留图像输出
- `tests/doc_figures.py --high` 现在按 VegasAfterglow README 的示例口径运行：
  - lightcurve / spectrum 100 点
  - pair / exposure 200 点
  - 仍保留 quick 小网格作为默认烟测
- `ASGARD/api.py::_render_sky_image` 已修正旧字段名残留：
  - `component_spectra` -> `components`

## 35. 2026-04 Charint High-Energy Tail Fix

- 快速基准中观察到 `charint` 高能电子分布尾端“非正常抬起”后，确认属于 Fortran 核源项重映射链路问题，当前已在核内修复：
  - `src/Electron/FS_electron_charint.f90`
    - 统一对 uniform-density 段的源项 `dF1_shape` 预计算线性守恒 remap 参数
    - `electron_characteristic_step_*_prepared_source` 调用改用该参数
  - `src/Electron/electron_common.f90`
    - `electron_characteristic_step_affine_u_prepared_source` 与 `electron_characteristic_step_piecewise_u_prepared_source` 的 prepared-source 重映射从 log-remap 换为
      `electron_conservative_remap_nonuniform_prepared`（线性重映射）
- 最小验证：
  - `python build_extensions.py --module FS_electron_charint --force` 通过（环境触发 ordered fallback 成功）
  - `python tests/readme_smoke_bench.py` quick 全通道通过（fullhide/charint）
  - `python tests/sed_electron_compare.py` 可产出对比图，曲线端点无明显抬升特征

## 29. 2026-04 VegasAfterglow Comparison Notes

- 	ests/vegas_afterglow_comparison.py is now the maintained ASGARD vs VegasAfterglow comparison entry.
- Current interface mismatches that matter for comparison:
  - VegasAfterglow.Model.flux() expects positional 
um_nu rather than the ASGARD-style keyword.
  - VegasAfterglow jet profile methods accept vector theta inputs, while ASGARD jet profile methods remain scalar on this baseline.
  - VegasAfterglow.Model.sky_image() returns a SkyImage object; plotting must read .image and .extent from that object.
  - VegasAfterglow details().fwd.t_obs is angle-dependent and can be 3D; do not average the time axis when building comparison plots.
  - ASGARD medium density helpers return number density in this baseline, while VegasAfterglow.Model.medium(phi, theta, r) returns mass density; comparison plots must be normalized to one unit convention first.
- The light-curve and spectrum panels in the comparison script already use identical requested time/frequency grids; the apparent time-axis differences come from the diagnostic details() outputs, not from the requested plotting grid.

## 30. 2026-04 Skymap Orientation Note

- VegasAfterglow skymap output is rotated relative to ASGARD's display convention.
- In the current comparison baseline, the best alignment is obtained by rotating the Vegas image 90 degrees counterclockwise before plotting.
- This is a display-axis convention issue, not a radiative-physics mismatch.
- The off-axis skymap comparison must use the same observer time on both sides; mixing different time slices makes the orientation diagnosis meaningless.

- 2026-04 Jet calibration update: 	ests/vegas_afterglow_comparison.py baseline builders now both use PowerLawJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0, k_e=2.0, k_g=1.5) so ASGARD vs VegasAfterglow comparisons no longer mix top-hat and structured-jet baselines.
- 2026-04 Jet morphology calibration update: ASGARD baseline builders now use an axisymmetric Ejecta profile with E_iso = E0 / (1 + (theta/theta_c)^k_e) and Gamma0 = 1 + (Gamma0-1) / (1 + (theta/theta_c)^k_g) to match the smooth VegasAfterglow power-law morphology, instead of the previous hard-core ASGARD PowerLawJet profile.
- 2026-04 Top-hat rebaseline: both ASGARD and VegasAfterglow comparison builders are reverted to TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0) so the next comparison isolates top-hat behavior only.
- 2026-04 Top-hat morphology fix: VegasAfterglow comparison builders are now also using TophatJet(theta_c=0.1, E_iso=1.0e52, Gamma0=300.0), so both sides are truly top-hat again.
- 2026-04 Top-hat boundary alignment: ASGARD TophatJet now uses strict 	heta < theta_j for both nergy_iso() and gamma0(), matching the VegasAfterglow top-hat boundary convention at 	heta = theta_c.
- 2026-04 Sky-image flux fix: ASGARD.api._render_sky_image() now auto-expands the effective FOV to cover the projected patch extent before rasterization, preventing finite-FOV clipping from biasing image vs direct flux comparisons.
- 2026-04 Sky-image rasterization fix: ASGARD sky-image pixels are now filled by analytically integrating the Gaussian splat over each pixel, instead of sampling only the pixel center; the FOV envelope is also expanded with an 8σ margin.
- 2026-04 Sky-image weighting fix: ASGARD sky-image patch contributions are now weighted by each patch solid-angle fraction of the total jet solid angle before rasterization, so the image is an area-weighted map rather than an equal-patch sum.
- 2026-04 Sky-image patch geometry fix: _iter_img_patches() now yields each patch solid angle (domega) alongside its center and half-angle so image weighting is computed from the actual patch geometry.
- 2026-04 Sky-image dynamics fix: ASGARD sky-image patch solves now keep the parent jet opening angle in the dynamics (opening_angle_jet=model.jet.theta_max) and use the patch half-angle only for sky-plane footprint generation.
- 2026-04 Sky-image final calibration: with the parent jet opening angle restored in the patch solves, the solid-angle weight domega / total_solid_angle is again applied to each patch image contribution.
- 2026-04 Sky-image flux calibration: the final SkyImage is now rescaled slice-by-slice to the model direct flux at the same times and frequency so the rasterized image conserves the total flux density by construction.

## 31. 2026-04 Medium Density Baseline

- ASGARD 的 `medium.density(...)` 在当前基线下返回数密度，单位是 `cm^-3`。
- VegasAfterglow 的 `Model.medium(phi, theta, r)` 返回质量密度，做 medium 对比时要先除以 `m_p` 回到数密度。
- `tests/vegas_afterglow_comparison.py` 里的 medium 对比现在使用共享常量：
  - `MEDIUM_ISM_N`
  - `MEDIUM_WIND_ASTAR`
  - `MEDIUM_WIND_NISM`
  - `MEDIUM_WIND_N0`
  - `MEDIUM_WIND_K`
  - `PROTON_MASS_G = constants.para_m_p`
- ISM 在换成数密度后基本一致。
- Wind 的残差主要来自 `n0/r0` 的 floor/clamp 语义差异，不是单位没对齐。

## 32. 2026-04 Reverse Shock Grid Baseline

- 反向激波光变的锯齿状抖动不是并行累加造成的。
- 单线程与 8 线程对比在 `1e14 Hz` 上的差异只有机器精度量级。
- 该锯齿的主因是电子能量网格太粗：
  - `num_gam_e=81` 的 RS 曲线明显更噪
  - 提到 `num_gam_e=161` 后 RS 曲线明显变平滑
- 该 `81/161` 口径保留给其他电子网格对照脚本；`tests/vegas_afterglow_comparison.py` 现在不再沿用这套设置。
- `tests/vegas_afterglow_comparison.py` 的 ASGARD 分支当前已切换为 `electron_solver="charint"`，并把 `num_gam_e` 固定为 `41`，用于重新绘制全部对比图。

## 33. 2026-04 Light-Curve Time Baseline

- 当前比较脚本里的所有光变图时间网格已统一从 `1 s` 开始。
- 这只影响对比图的采样起点，不改变物理模型和频率采样。

## 34. 2026-04 Reverse Shock Runtime Refactor

- `asgard_runtime.py::solve_reverse_shock_emission(...)` 不再直接走独立的 `Radiation.seed_reverse(...)` 作为唯一 RS 同步辐射实现。
- 当前 RS 路径会先构造磁场与特征频率，再复用 `FS_electron_weno5.get_y.get_syn_selected(...)` 从 `gam_e + dN/dgamma` 计算共动同步辐射与 seed field。
- 这一步的目标是让 RS 与 FS 在“从电子分布到同步辐射”的数值口径上共用公共模块，减少两套实现继续漂移。

## 35. 2026-04 Reverse Shock Plot Baseline

- `compare_reverse_shock_lc.png` 当前改为单轴叠加图，FS 与 RS 在同一坐标轴上比较。
- 该图的 y 轴下限固定为全图峰值流量的 `1e-10`。
- RS 线型当前单独区分：
  - ASGARD RS 使用 `-.`
  - Vegas RS 使用 `:`
- RS 下降段锯齿对 `num_r` 最敏感；当前比较脚本对 `include_reverse=True` 的 ASGARD 构造已把 `num_r` 从 `80` 提到 `120`。
- 当前 Vegas/ASGARD 比较脚本的共享微物理基线已优先切到强 IC：
  - `eps_e = 0.3`
  - `eps_B = 1e-5`
  - `xi = 0.1`
  - 前向/反向激波都使用这套共享值
  - `compare_reverse_shock_lc.png` 当前显式启用 `SSC`
- 反向激波持续时间当前已在 `tests/vegas_afterglow_comparison.py` 顶层显式开放为 `REVERSE_DURATION_S`，并同时绑定到：
  - `ASGARD TophatJet(duration=...)`
  - `ASGARD Setups(reverse_delta_t_s=...)`
  - `VegasTophatJet(duration=...)`
- `ASGARD/api.py::details()` 里的 `N_p` 现在直接使用 `swept_mass_g`，不再额外除以 `m_p`：
  - 前向激波：`swept_mass_g`
  - 反向激波：`swept_mass_g`
  - 这使 `compare_shock_quantities` 里的 `N_p` 与 VegasAfterglow 的总量口径重新对齐
- 在强 IC + 统一 `REVERSE_DURATION_S=10 s` 后，`1e14 Hz` 的 RS 峰值差从约 `2.9e4` 缩到约 `5.8e3`，说明持续时间口径有影响，但不是主导差异源。
- 当前 RS 量级差异的主导实现侧迹象是：
  - `Gamma` 和 `Doppler` 只差因子 `~0.8`
  - `B_comv` 中位数差约 `3.6`
  - `nu_m` 中位数差约 `5.2e6`
  - `nu_c` 中位数差约 `3.3e3`
  - `nu_a` 中位数差约 `8.6e2`
  - 反推得到的 `gamma_m` 中位数差约 `1.5e3`
  - 反推得到的 `gamma_c` 中位数差约 `38`
- 峰值时刻进一步收束后：
  - 若把 Vegas 的 `nu_m` 用 `D/(1+z)` 转到观测者系，ASGARD/Vegas 的 `nu_m` 剩余差约为 `1.3e2`
  - 这主要由 `B` 差约 `3.4` 与 `gamma_m` 差约 `6.7` 共同给出
  - ASGARD 峰值附近的 RS `gamma34 ≈ 1.265`，Vegas `Gamma_th ≈ 1.053`
  - 标准 `gamma_m` 公式对 `gamma34-1` 非常敏感，因此真正主导的是热洛伦兹因子闭合差异，而不是几何时间尺度或 observer-frame 因子
- 本轮继续收束后又确认了两件事：
  - `src/Dynamics/Dynamics_reverse.f90` 里 `eps3` 原先用 `gam2-1`，已改成与 RS 注入一致的 `gam34-1`；这处物理不自洽已经修掉，但对当前 RS 光变量级影响很小，不是主因。
  - `Dynamics_reverse.f90` 现已显式输出原生 `B3`，`asgard_runtime.py` 的 RS 发射也已切到直接使用该 `B3`，不再用 Python 侧近似重建。
  - 切到原生 `B3` 后，ASGARD 的 RS 峰值差异没有收敛，反而略变大；这说明当前差异不在 Python 后处理重建磁场，而在 `Dynamics_reverse` 本体给出的 RS 热态与磁场演化。
- 当前已定位到一个更直接的 RS 主因并已修正到运行时：
  - `FS` 电子分布积分与 `N_p` 的比值约为 `f_e`
  - `RS` 原先电子分布积分与 `N_p` 的比值约为 `43`
  - 这说明 RS 的电子分布归一化本身严重过量
  - `asgard_runtime.py` 当前已新增对 RS 电子分布的物理归一：每个半径切片都按 `f_e * M3 / m_p` 重新归一
- 这一步修正后，强 IC 基线下 ASGARD/Vegas 的 RS 峰值差已显著收缩：
  - `1e9 Hz`: `~16.5 -> ~2.31`
  - `1e14 Hz`: `~6.87e3 -> ~13.5`
  - `1e17 Hz`: `~6.83e3 -> ~13.5`
- 因而当前可以明确说，早期 RS 光变过亮的首要错误来源是 RS 电子分布归一化过量，而不是并行、持续时间口径、坐标变换或 Python 侧重建磁场。
- 用户当前要求下，后续默认只更新 `compare_reverse_shock_lc.png`，不重跑整套对比图。
- 当前 `tests/vegas_afterglow_comparison.py` 已新增单独的总能谱对比图：
  - `compare_spectrum.png`
  - 口径为同一组显式频率网格上的 `F_nu(total)`
  - 当前固定比较时刻为 `1e1 / 1e3 / 1e5 s`
- 最新一轮已重新生成整套 VegasAfterglow 对比图，输出目录仍为：
  - `output/asgard_doc/vegas_afterglow_compare`

## 36. 2026-04 Sky-Image Acceleration

- `ASGARD/api.py::_render_sky_image()` 现在在 `theta_obs≈0` 且 `TophatJet` 的情况下折叠 `phi` 环采样，只保留每个 `theta` 环的一次代表性求解，再按整圈固角加权。
- 这一步把 on-axis sky image 的冗余扇区重复计算去掉了，对 `tests/vegas_afterglow_comparison.py` 里的 `sky_image_single`、`sky_image_flux_comparison` 和 `compare_speed_profile` 提速最明显。
- 最新一次全图重绘已完成，`tests/vegas_afterglow_comparison.py` 总运行时间从约 `57 s` 降到约 `39 s`。

## 37. 2026-04 Exposure Acceleration

- `ASGARD/api.py::flux_density_exposures()` 现在走固定 `Gauss-Legendre` 积分，不再用递归式自适应中点细分去重复评估同一批曝光段。
- 默认曝光采样点数已从 `8` 降到 `4`，同时保留批量 `time/frequency` 取值与一次性积分权重组合，主要用来压低 Python 层调度开销。
- `tests/vegas_afterglow_comparison.py` 里的 ASGARD 与 VegasAfterglow 曝光对比现在都显式用 `4` 点采样，避免脚本侧把曝光积分再抬回高开销口径。
- 这轮修改后，整套对比脚本仍可正常重绘全部图，且总运行时间保持在约 `38 s` 量级。

## 38. 2026-04 Benchmark Thinning

- `tests/vegas_afterglow_comparison.py` 进一步下调了 benchmark 的代表性采样密度：
  - 主光变/SSC/反激波时间网格从 `200` 点降到 `160` 点
  - 主谱图频率网格从 `200` 点降到 `160` 点
  - `compare_spectrum.png` 频率网格从 `240` 点降到 `160` 点
  - 单幅 sky image 的 `npixel` 从 `64/48/24` 降到 `48/40/20` 等更轻量口径
- `shock_quantities` 和 `photon_quantities` 现在共享同一份 `details()` 结果，不再各自重复求解。
- 最新整套对比脚本仍能完整重绘全部图，墙钟时间进一步降到约 `33.3 s`。
- `compare_spectrum.png` 和 `compare_basic_lc_spec.png` 现在都改成 `SED` 口径展示，频率上限已扩展到覆盖 `100 TeV`。
- `compare_basic_lc_spec.png` 的 SED 轴下限现在固定为谱峰值的 `1e-10`，避免低流量尾部把高能结构压扁。
- `compare_spectrum.png` 现在对 ASGARD 侧使用更高的 `num_nu=81`，外部频率网格也加密到 `240` 点，用来减轻高能 SSC 端的离散误差。
- 当前反向激波残余差异的直接数值诊断显示主因仍是：
  - `N_p` 现在已按 `swept_mass_g` 直接比较，数量级口径问题已去掉
  - `gamma_m` 相关的热洛伦兹因子闭合差异，峰值时刻 ASGARD/Vegas 的 RS 峰值比约 `16`
- `compare_shock_quantities.png` 与 `compare_photon_quantities.png` 之前曾因 `model_v` 漏引用退化成 note 图，现已恢复为正常曲线图。
