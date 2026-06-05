# 二维 chi-resolved 强子输运决策

本文记录当前是否应把 hadronic transport 扩展到 `fullhide_2d/charint_2d` 使用的 `chi_grid`。结论是：**暂不实现**。原因不是数值工作量，而是当前代码的 photon-field、hadron-density、secondary feedback 和 observer projection 契约还没有定义 χ 维度；直接把现有 1D hadronic kernel 外包一层 χ 循环不会得到物理自洽模型。

## 当前契约

- `ElectronSolution` 已经有 2D 电子输运输出：`d_n_gam_e_chi` 与 `chi_grid`。这组量来自 2D electron backend，用于电子谱在局部历史坐标上的输运。
- `PhotonFieldState` 仍是 shell-level 契约，包含 `seed_frequency_hz`、`forward_syn_seed`、`hadronic_forward_ssc_seed`、`hadronic_target_seed`、`absorption_*_seed`。这些 seed arrays 是 `frequency x radius`，没有 χ 轴。
- `solve_hadronic(...)` 明确拒绝非 1D electron solver，并把 target photons 作为 `seed_target_arr[:, i_r]` 传入 shell-local pγ/BH/pp/secondary/cascade operators。
- `HadronicSolution` 的 proton、secondary、radiation arrays 都是 `gamma/frequency x radius`。observer assembly 也没有 χ-resolved hadronic emissivity 的面积权重、Doppler 权重或 arrival-time projection。

## 物理问题

χ-resolved hadronic transport 不能只复用电子 2D 的坐标标签。pγ、BH、γγ pair branch、secondary synch/IC 的反应率取决于局部 photon density、baryon density、magnetic field、shell volume and residence time。若这些量没有在同一 χ 体元中闭合，则过程功率和冷却率会混用 shell 平均 photons 与局部 hadrons，能量预算没有可解释的守恒对象。

当前正式 hadronic path 的自洽对象是 shell：每个 radius shell 有一个 photon target field、一个 magnetic field、一个 proton/secondary spectrum 和一个 observer-side projection。2D electron backend 的 `chi_grid` 解决的是电子冷却/历史 photon field 的空间结构；它还没有定义 hadron injection、target baryon column、secondary pair feedback 和 photon escape 在 χ 方向上的合同。

## 当前决策

当前不进入 2D / χ-resolved hadronic transport 实现。继续使用 formal 1D hadronic path：

- FS hadronic: `am3_1d` 为正式路径，`legacy_1d` 仅覆盖 proton transport + proton synch。
- RS hadronic: light path 由 `hadronic_reverse_1d` 覆盖 proton synch；full-chain RS path 复用 formal 1D hadronic kernels，使用 RS seed photons、RS `B3`、shell energy 和 baryon target density。

进入 2D 实现前必须先完成这些契约设计：

1. `PhotonFieldState` 增加 `seed_target_chi_hz` 或等价结构，形状为 `frequency x chi x radius`。
2. 定义 χ 体积/权重，使局部 emissivity 积分回 shell luminosity 时保持能量守恒。
3. 定义 proton injection 与 hadron density per χ，明确它是否跟随 shock-front、electron χ history，或另一组 hydrodynamic residence-time 坐标。
4. 定义 pγ/BH/pp/secondary/cascade 的 χ-local source-sink 方程，以及 secondary e± 如何反馈到 `ElectronSolution` 和 observer projection。
5. 给出至少一个 energy-budget benchmark：shell-integrated χ 结果在 optically thin limit 下回到当前 1D 结果，且随 radius 的 `B`、target photon energy density、hadronic luminosity 无非物理尖跳。

在这些条件满足前，新增 χ-resolved hadronic kernel 只会制造一个接口上更复杂、物理上更弱的模型。
