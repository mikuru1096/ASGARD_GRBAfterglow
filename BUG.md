# ASGARD defect ledger

- **P1 — Forward Electron 非热初态漏乘 `f_e`。**
  症状：各 forward Electron driver 的非热初态把 `electron_initial_density` 返回的总扫掠电子数直接交给 `init_coord`、`init_edges`、`init_powerlaw` 或等价谱初始化，而后续壳层注入却按 `f_e` 缩放；thermal fullhide 还在完整非热初态上叠加 `(1-f_e)N_e`，使初始总数达到 `(2-f_e)N_e`。
  原因：`electron_initial_density` 返回的是总电子数，初态调用链没有统一区分 `f_e N_e` 的加速非热分量与 `(1-f_e)N_e` 的热分量。
  影响：首个径向点及早期时刻的电子数、同步辐射、SSC、冷却能标和 `f_e` 拟合归一化错误，并与随后正确按 `f_e` 注入的输运状态产生系统性偏置。
  严格验收：逐个 forward Electron solver 在 `f_e=0.1,0.3,1` 下证明初始非热数为 `f_e N_e`；thermal fullhide 分别闭合到 `f_e N_e` 和 `(1-f_e)N_e`、总数严格为 `N_e`；`f_e=1` 与现基线一致，首壳层粒子数闭合，电子谱和同步/SSC 输出有限、非负并随 `f_e` 连续变化。

- **P1 — `electron_initial_density` 未累计现代真实介质分布。**
  症状：tabulated density profile 或位于初始半径以内的 modern density jump 存在时，初始局部密度和累计电子数仍按无 jump 的解析 ISM/wind 公式生成。
  原因：`electron_initial_density` 只接收 `A_star/dNe_ISM/R_ini/R_start/R0`，没有使用 `density_profile` 已解包的 `jump_count/profile_count` 状态，也没有计算 `4*pi*integral_0^R_ini r^2 n(r)dr`。
  影响：Electron 初态与 Dynamics 所用介质不一致，初始粒子数、磁场、`gamma_m/gamma_c/gamma_max` 及早期辐射均产生永久归一化偏移；jump/profile 位于 `R_ini` 附近时最明显。
  严格验收：初始局部密度必须逐点等于 `density_profile(R_start)`，累计数必须等于同一介质合同的体积积分；uniform ISM/wind 解析结果保持到舍入误差，legacy/modern 单 jump 完全等价，modern 多 jump、tabulated power-law、`R0` 截断以及 `R_ini` 位于特征前/内/后的结果分别与独立高精度积分相对一致到 `1e-10`；累计数随 `R_ini` 单调，并在平滑 jump 和 profile 分段连接处连续，不允许用平滑或后处理掩盖误差。
