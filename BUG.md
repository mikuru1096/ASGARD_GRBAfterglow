# ASGARD defect ledger

- **P1 — tabulated medium 在 2D chi 输运中被静默当作 `k=0`。**
  症状：公开 `fullhide_2d/charint_2d` 可接收 tabulated profile，但 Python `_chigrid` 与 Fortran `fs_transport_2d` 都只按 `a_star` 二分 `k_medium=0/2`；profile 的 `a_star=-1`，因此任意径向斜率都进入 ISM 几何。
  原因：有限-q downstream geometry、体积权重和绝热散度采用常数幂律 `n∝R^{-k}` 的 BM 参数化，Boundary 中虽有 profile 表，2D 核却没有非幂律几何合同；现有 `profile_tag` 形参未参与计算。
  影响：shock-front 注入使用真实 profile 局部密度，而 chi 半径、体积、绝热冷却、输运与投影权重使用 `k=0`，可能得到平滑但物理不自洽的 2D 电子谱和光变。
  严格验收：在推导并实现非幂律介质的 chi 几何前，Python 系统边界必须明确拒绝 tabulated profile 与 2D electron solver 的组合；若实现一般化，则必须从同一 profile/shock history 得到几何与散度，并严格退化到 `k=0/2`、保持粒子数闭合和 nchi 收敛。不得静默选 `k=0`、用局部拟合 k 或后处理修正。

- **P2 — ReverseShockCausalityDiagnostics 将 tabulated medium 标记并近似为 ISM。**
  症状：`_rsdiagnostics` 只按 `a_star>0` 区分 wind/ISM；tabulated profile 因 `a_star=-1` 被标为 `medium="ism"`，`reference_crossing_radius_cm` 使用末端 `d_ne` 的均匀介质 Sedov 标度。
  原因：解析参考诊断没有 profile 的累计扫掠矩合同，尽管 pressure ratio 已通过 `ambient_density` 使用真实局部 profile。
  影响：实际 Reverse Dynamics 可正确运行，但用户可见的介质标签、参考 crossing 半径和 criteria 对照不代表所输入的 profile，可能误导物理判断。
  严格验收：诊断必须标识 `density_profile`，并从与 Dynamics 相同的累计 profile 矩求解相应减速/参考 crossing 尺度；ISM/wind 极限须回到现有解析式。不得把 profile 近似为末端常密度、返回占位值或隐藏该字段。
