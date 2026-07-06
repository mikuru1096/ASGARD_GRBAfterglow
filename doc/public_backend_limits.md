# 公开 API 与后端边界

本文固定当前公开 API 已暴露但后端不支持或只部分支持的边界。原则：不制造伪支持，不添加兜底路径；未支持边界只有在明确物理问题需要时才进入实现。

## 明确不支持

- 喷流横向扩张：
  - 公开喷流构造器接受 `spreading` 参数。
  - `_build_fit_config_for_patch(...)` 在 `model.jet.spreading=True` 时直接抛出 `NotImplementedError`。
  - 当前结构化喷流/top-hat 角向面元动力学是固定开角、独立球面楔；不包含横向扩张。
- 用户自定义 `Medium` 的数值核调度：
  - `Medium.to_kernel_params()` 对自定义介质直接抛出 `NotImplementedError`。
  - 当前数值核契约只接受 `ISM` 和 `Wind(k=2)` 映射出的参数。
- Wind `k != 2`：
  - `Wind._rho(...)` 和 `Wind.to_kernel_params()` 都拒绝 `k != 2`。
  - 当前 Fortran/Python 密度契约是 `A_star r^-2` 星风；当 `A_star 3e35 / R^2 <= n_ISM / 4` 时切到 ISM 下限，不是任意 `r^-k`。
- `fullhide_1d`/`fullhide_1d_hz`/`dg_1d` 之外的热电子：
  - `solve_electron(...)` 在 `thermal_electrons=True` 且求解器不是 `fullhide_1d`, `fullhide_1d_hz`, `dg_1d` 时直接抛出 `NotImplementedError`。
  - 现有热电子分支未定义到 `weno5_1d`, `slc1_1d`, `charint_*`, `fullhide_2d`。
- 非轴对称喷流上的环向磁场偏振：
  - `Model.polarization(..., magnetic_geometry="toroidal")` 当前要求轴对称喷流。

## 部分支持边界

- 2D 电子求解器：
  - `fullhide_2d` 和 `charint_2d` 可用作电子输运后端。
  - 强子主路径仍要求 1D 电子求解器；2D/\(\chi\) 强子输运已在 `doc/hadronic_chi_transport_decision.md` 中明确暂不实现。
- 反向激波强子：
  - 轻量路径是反向激波质子输运 + 质子同步辐射。
  - 完整链路复用正式 1D 强子核；不是 2D/\(\chi\) 反向激波强子输运。
- 对级联：
  - `pair_cascade_iterations > 1` 是壳层序列 \(\gamma\gamma\) 对产生/同步辐射级联。
  - 逆康普顿介导的电磁级联不属于当前契约，见 `doc/pair_cascade_extension_boundary.md`。
- 偏振：
  - 同步辐射 Stokes 路径覆盖 FS/RS 电子同步辐射和 FS/RS 强子同步辐射。
  - 非同步分支不混入偏振 Stokes；Lan 2023 峰时偏移当前归入动力学/喷流演化基准测试。

## 实现准入条件

要把以上边界推进为实现，必须先写清：

1. 物理动机和目标观测问题。
2. 新后端契约：输入、输出、单位、参考系、数组形状。
3. 与现有 Fortran/Python 主链的最小耦合点。
4. 连续平滑的物理诊断和能量/粒子数预算。
5. 受影响的构建、行截断检查、冒烟测试和基准测试列表。

没有这些条件，不创建新代码路径。
