# ASGARD defect ledger

- 症状：多个 electron transport driver 把子步右端或中点半径传给精确壳体积因子 `((R+h)^3-R^3)/h`，使注入体积向外平移一个或半个子步。
  原因：调用方先推进 `R_loc` 后构造源项，但 `electron_injection_prefactor` 的 `R` 合同是区间左端；uniform merged 路径还精确累加了同一组右移区间。
  影响：fullhide 1D general/adaptive/coupled、DG、charint 及部分 2D transport 的电子数与能量注入归一化，density jump 附近误差更明显。
  验证：逐 driver 统一到真实子区间左端，证明总注入数等于 `4*pi*f_e*integral(r^2*n(r)dr)`；比较均匀介质和 density-jump 直接运行的粒子数闭合、谱收敛与时间平滑性。

- 症状：现代 `jump_r_cm/jump_factor/jump_width` 或 tabulated density profile 存在时，`fullhide_1d`/`charint_1d` 仍可能因 `f_jump==1` 误判为 uniform fast path，从而整壳冻结密度并漏掉跳变。
  原因：uniform 判据没有包含 `dynamics_density_profile` 的 `jump_count/profile_count`。
  影响：现代多密度跳和 tabulated-medium 的电子注入/冷却主链。
  验证：判据与 DG 已有合同一致；现代 jump 直接算例必须实际进入非均匀路径，跨 jump 的密度、电子数、谱与光变连续响应和高分辨率参考一致。
