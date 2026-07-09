# ASGARD defect ledger

- **P2 — ReverseShockCausalityDiagnostics 将 tabulated medium 标记并近似为 ISM。**
  症状：`_rsdiagnostics` 只按 `a_star>0` 区分 wind/ISM；tabulated profile 因 `a_star=-1` 被标为 `medium="ism"`，`reference_crossing_radius_cm` 使用末端 `d_ne` 的均匀介质 Sedov 标度。
  原因：解析参考诊断没有 profile 的累计扫掠矩合同，尽管 pressure ratio 已通过 `ambient_density` 使用真实局部 profile。
  影响：实际 Reverse Dynamics 可正确运行，但用户可见的介质标签、参考 crossing 半径和 criteria 对照不代表所输入的 profile，可能误导物理判断。
  严格验收：诊断必须标识 `density_profile`，并从与 Dynamics 相同的累计 profile 矩求解相应减速/参考 crossing 尺度；ISM/wind 极限须回到现有解析式。不得把 profile 近似为末端常密度、返回占位值或隐藏该字段。
