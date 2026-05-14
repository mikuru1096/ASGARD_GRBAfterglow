# Public API / Backend Limits

本文固定当前 public API 暴露但 backend 不支持或只部分支持的边界。原则：不制造伪支持，不添加兜底路径；未支持边界只有在明确物理问题需要时才进入实现。

## Explicitly Unsupported

- Jet spreading:
  - Public jet constructors接受 `spreading` 参数。
  - `_build_fit_config_for_patch(...)` 在 `model.jet.spreading=True` 时直接 `NotImplementedError`。
  - 当前 structured/top-hat patch dynamics 是固定开角、独立 spherical wedge；不包含 lateral expansion。
- User-defined `Medium` kernel dispatch:
  - `Medium.to_kernel_params()` 对 custom medium 直接 `NotImplementedError`。
  - 当前 kernel contract 只接受 `ISM` 和 `Wind(k=2)` 映射出的参数。
- Wind `k != 2`:
  - `Wind._rho(...)` 和 `Wind.to_kernel_params()` 都拒绝 `k != 2`。
  - 当前 Fortran/Python density contract 是 `A_star r^-2` wind；当 `A_star 3e35 / R^2 <= n_ISM / 4` 时切到 ISM floor，不是任意 `r^-k`。
- Thermal electrons outside `fullhide_1d`:
  - `solve_electron(...)` 在 `thermal_electrons=True` 且 solver 不是 `fullhide_1d` 时直接 `NotImplementedError`。
  - 现有 thermal branch 未定义到 `weno5_1d`, `slc1_1d`, `charint_*`, `fullhide_2d`。
- Toroidal polarization on non-axisymmetric jets:
  - `Model.polarization(..., magnetic_geometry="toroidal")` 当前要求 axisymmetric jet。

## Partially Supported Boundaries

- 2D electron solvers:
  - `fullhide_2d` 和 `charint_2d` 可用作 electron transport backend。
  - Hadronic main path 仍要求 1D electron solver；2D/χ hadronic transport 已在 `doc/hadronic_chi_transport_decision.md` 中明确暂不实现。
- Reverse-shock hadronic:
  - light path 是 RS proton transport + proton synch。
  - full-chain path 复用 formal 1D hadronic kernels；不是 2D/χ RS hadronic transport。
- Pair cascade:
  - `pair_cascade_iterations > 1` 是 shell-sequence γγ pair/synch cascade。
  - IC-mediated electromagnetic cascade 不属于当前合同，见 `doc/pair_cascade_extension_boundary.md`。
- Polarization:
  - 同步辐射 Stokes 路径覆盖 FS/RS electron synch 和 FS/RS hadronic synch。
  - 非同步分支不混入 polarization Stokes；Lan 2023 峰时偏移当前归入 dynamics/jet-evolution benchmark，见 `doc/polarization_timing_diagnostic.md`。

## Implementation Entry Criteria

要把以上边界推进为实现，必须先写清：

1. 物理动机和目标观测问题。
2. 新 backend contract：输入、输出、单位、frame、数组形状。
3. 与现有 Fortran/Python 主链的最小耦合点。
4. 连续平滑的物理诊断和能量/粒子数预算。
5. 受影响的 build、line-truncation 和 smoke/benchmark 列表。

没有这些条件，不创建新代码路径。
