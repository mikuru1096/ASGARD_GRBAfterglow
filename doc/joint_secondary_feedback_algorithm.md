# Joint 二级反馈算法契约

`SolverOptions(electron_photon_coupling="joint")` 在同一半径网格闭合电子、光子、
强子和二级 (e^\pm)。默认 `separated` 不改变。物理预算见
[`joint_secondary_feedback_physics.md`](joint_secondary_feedback_physics.md)。

> **已知缺陷**：pair state、source rate 与 radiation ownership 尚未唯一化，存在
> 重复注入或漏乘步长风险。正式约束记录在根目录 `BUG.md` 第 3 项；本页描述当前
> 实现，不宣称该缺陷已解决。

## 1. 支持边界

joint 当前要求：

- `fullhide_1d`、fixed substeps、`ssc_cooling_mode="numeric_ic_kn"`；
- forward 1D top-hat/Fortran chain，无 reverse shock 和 χ transport；
- hadronic enabled、`am3_1d`、Hummer pγ 与 BH；
- 正 proton energy fraction。

不满足时在公开配置边界报错，不自动退回 separated。

## 2. 两条主链

Separated：

```text
dynamics -> electron -> photon field -> hadronic
         -> optional BH merge -> observer
```

Joint：

```text
dynamics -> primary electron -> photon field
repeat fixed internal passes:
    formal hadronic transport
    photon survival
    BH/pp/gamma-gamma secondary source
    coupled electron cooling and radiation
    rebuild photon field
observer
```

固定 pass 数是内部算法常量，不是公开 convergence knob。`_mergebh` 只属于
separated，不能复用为本壳层反馈。

## 3. 状态单位

`PhotonFieldState` 分开保存 forward synch seed、hadronic target、absorption seed
和 SSC contribution。它们都是局域 photon field，不是 observer luminosity。

`HadronicSolution` 向 joint 提供：

- `secondary_electron_source_r(gamma,R)`：电子方程的半径源率；
- `tau_bh(nu,R)`：当前 shell path 的 BH optical depth；
- `bh_photon_loss_rate(nu,R)`：共动微物理损失率诊断。

源率到电子 state increment 的积分只能由输运层执行一次。这正是 BUG.md 所登记
所有权问题的修复目标。

## 4. Proper-time 步长

首输出点是零时长初态。其后使用

\[
\Delta t'_i=\frac{R_i-R_{i-1}}{2c}
\left[\frac1{\sqrt{\Gamma_{i-1}^2-1}}+
\frac1{\sqrt{\Gamma_i^2-1}}\right].
\]

`R_Tobs` 仅为既有 ABI 保留，不参与局域强子推进。每个径向区间只能推进一次。

## 5. Electron pass

`solve_coolingseed -> fs_fullhide_coupled` 接收 target seed 和 secondary source。
IC cooling 与 IC photon emissivity使用同一 Jones/KN 响应。secondary source 与
shock injection 在每个 substep 进入同一守恒方程。

## 6. Pair branch

启用 γγ 时，吸收 photon 产生 pair source；pair 继续输运并产生同步 seed。多次
iteration 是 shell-sequence pair/synch cascade，不包含 IC-mediated 完整电磁级联。

Photon survival、pair injection 和 pair synch seed 各应用一次。不得用额外 optical
depth、壳厚或 Doppler 因子修补预算。

## 7. Observer

joint 只改变 observer 前的局域状态。FS/RS/hadronic component 名称、吸收数组和
最终 projection API 不变。

## 8. 验收门槛

修复 BUG.md 第 3 项前必须证明：

1. 每个 source 的单位和 owner 唯一；
2. 冻结背景下 separated/joint 单步可逐数组比较；
3. pair energy 不超过 BH/pγ/γγ 注入预算；
4. 禁用 secondary 后回到 primary electron；
5. shell state、辐射和 light curve 连续；
6. electron 与 hadronic extension 强制构建并通过 line-truncation。

禁止删除断言、绕过 formal transport、平滑曲线或添加经验 photon sink/source。
