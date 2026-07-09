# ASGARD defect ledger

- **P1 — Python 与 Fortran 的 tabulated profile 端点外推合同不一致。**
  症状：Fortran `tab_density` 在 profile 首点以前和末点以后继续使用首段/末段的 log-log 幂律斜率；`asgard_physics_utils.ambient_density` 与 `TabulatedMedium._rho` 使用 `np.interp`，端点以外保持常数。
  原因：Python 路径把 profile 半径限制到了采样表端点的隐式 `np.interp` 行为，没有实现 Fortran 已采用的首末段幂律外推；`ambient_density` 因而不能作为 Fortran profile 的独立验收基准。
  影响：profile 网格外的 coupled-shock geometry、诊断磁场、pp target density 和用户可见 `Medium.rho` 与 Dynamics/Electron 实际介质不同，早晚期结果及跨语言诊断会产生系统偏差。
  严格验收：Python 标量和数组路径必须在表内保持相同 log-log 插值，在首末点外使用与 Fortran 完全相同的首段/末段斜率；`ambient_density` 还必须匹配 Fortran 的 `R0` core，而 profile 介质不得叠加 density jump。对首点前、节点上、每个分段内、末点后及 `R0` 两侧的查询，Python/Fortran 与独立解析值相对一致到 `1e-12`；相关 coupled geometry、磁场和 pp target 使用同一密度，不允许通过夹住半径、常数端点、fallback 或后处理掩盖差异。
