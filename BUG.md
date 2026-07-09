# ASGARD defect ledger

- **P1 — Dynamics/Electron 缺少共享的真实介质初始扫掠矩。**
  症状：Forward Dynamics 以固定 `Y(2)=Boundary(2)=1e12` 启动，与介质和 `R_ini` 无关；Reverse Dynamics 只用旧 ISM/wind 公式初始化外介质扫掠质量，忽略 `R0`、density jump 和 tabulated profile；所有 Forward Electron 首列写在 `R(1)`，累计电子数却只积到 `R_ini`，且旧解析式同样看不到 modern jump/profile。
  原因：代码没有唯一的介质径向矩 `I(R)=integral_0^R r^2 n(r)dr`，Dynamics、Reverse Dynamics 和 Electron 分别使用固定哨兵或不完整的重复公式，局部 `density_profile` 与累计量没有共享同一合同。
  影响：初始动力学状态、公开 `swept_mass_g/N_p`、Electron 首列总粒子数及磁场和早期同步/SSC 归一化彼此不闭合；默认 `R_ini=1e14 cm, n=0.1 cm^-3` 的真实初始扫掠质量约为 `7.01e17 g`，已比固定 `1e12 g` 大约 `7.0e5` 倍。
  严格验收：由 `density_profile` 模块提供唯一解析 `I(R)`；Forward Dynamics 在 `R_ini` 初始化，`index_dyn=1/2` 使用 `Y(2)=I(R_ini)` 并继续以 `4*pi*m_p*Y(2)` 输出 g，`index_dyn=3` 使用 `Y(2)=4*pi*m_p*I(R_ini)` 并直接输出 g；Reverse Dynamics 使用 `Y(3)=4*pi*m_p*I(R_ini)`；Electron 首列的局部密度与累计数分别严格等于 `density_profile(R(1))` 和 `4*pi*I(R(1))`，并满足 `swept_mass_g(1)/m_p=N_e(1)` 的总粒子数闭合。uniform ISM、带 `R0` core 和实际 wind-to-ISM 切换的 wind、legacy/modern 单 jump、modern 多 jump、首末幂律外推的 tabulated profile 均使用同一合同；有限结果与独立分段高精度积分相对一致到 `1e-10`，legacy/modern 同参数单 jump 完全等价，`R0<=0` 且首段斜率 `s<=-3` 的发散输入必须明确失败。累计量随半径单调并在 `R0`、wind 切换和 profile 分段连接处连续，不允许通用数值积分、平滑或后处理补救。

- **P1 — Python 与 Fortran 的 tabulated profile 端点外推合同不一致。**
  症状：Fortran `tab_density` 在 profile 首点以前和末点以后继续使用首段/末段的 log-log 幂律斜率；`asgard_physics_utils.ambient_density` 与 `TabulatedMedium._rho` 使用 `np.interp`，端点以外保持常数。
  原因：Python 路径把 profile 半径限制到了采样表端点的隐式 `np.interp` 行为，没有实现 Fortran 已采用的首末段幂律外推；`ambient_density` 因而不能作为 Fortran profile 的独立验收基准。
  影响：profile 网格外的 coupled-shock geometry、诊断磁场、pp target density 和用户可见 `Medium.rho` 与 Dynamics/Electron 实际介质不同，早晚期结果及跨语言诊断会产生系统偏差。
  严格验收：Python 标量和数组路径必须在表内保持相同 log-log 插值，在首末点外使用与 Fortran 完全相同的首段/末段斜率；`ambient_density` 还必须匹配 Fortran 的 `R0` core，而 profile 介质不得叠加 density jump。对首点前、节点上、每个分段内、末点后及 `R0` 两侧的查询，Python/Fortran 与独立解析值相对一致到 `1e-12`；相关 coupled geometry、磁场和 pp target 使用同一密度，不允许通过夹住半径、常数端点、fallback 或后处理掩盖差异。
