# 强子 `pγ` 说明

## 标准命名

- `hummer_2010_response`：正式 pγ 方案，含 transport feedback。
- `disabled`

public API 只接受上述两个标准命名。

## 当前正式耦合

主链: `dynamics → electron solver → photon target field → proton transport → hadronic source/loss operators → observer projection`

已接入:
- proton injection、adiabatic cooling、synchrotron cooling、proton synchrotron emission。
- `hummer_2010_response` photopion：\(\alpha_p(E)\) 与 \(Q_p^{\rm reinj}(E)\) 反馈到 proton transport。
- \(\alpha_\gamma^{p\gamma}\) 转为 local shell photon survival factor：

\[
\tau_{p\gamma}(\nu,r)
=\frac{\alpha_\gamma^{p\gamma}(\nu,r)R}{12\Gamma c},
\qquad
f_{\rm surv}
=\frac{1-\exp(-\tau_{p\gamma})}{\tau_{p\gamma}},
\]

并在 observer projection 前作用。
- neutrino output。
- Bethe-Heitler：proton continuous cooling + secondary e±，随后并入 forward electron 并重算 `seed_syn`。
- hadronic IC：proton channel；pion/muon IC 通过 explicit secondary species 处理。
- pp：gamma、neutrino、pair source、proton loss。
- pair production：observer-side `tau_pair` attenuation + pair synchrotron branch；`pair_cascade_iterations > 1` 选择 shell-sequence time-dependent \(\gamma\gamma\) pair/synch cascade。
- explicit secondary species transport：n、π±、μ± 的左右向输运。
- secondary radiation：pion/muon synchrotron + IC。
- hadronic acceleration/injection operators。
- reverse-shock full-chain dispatch：开启 RS hadronic full-chain flags 时，runtime 复用 formal 1D kernels，并使用 RS seed photons、RS magnetic field、RS shell energy 和 RS baryon target density。

## 当前验证状态

- `tests/hadronic_1d_smoke.py` 是当前通过的 FS 1D hadronic smoke。
- `tests/hadronic_reverse_shock_smoke.py` 的 base 和 RS light proton-synch 分支可执行；full-chain RS hadronic 分支当前在 formal kernel 的 electron-energy grid contract 处失败，错误为 `electron_energy_gev must be logarithmically uniform`。
- `tests/electron_photon_joint_secondary_feedback_smoke.py` 当前触发同一个 formal hadronic electron-energy grid contract 失败。

## 当前边界

Hadronic 未完成项集中维护在根目录 `TODO.md`。本文件只记录当前 pγ 耦合、命名和 AM3 参考边界。

Reverse-shock baseline：ASGARD 使用局部 `gamma34` 注入和显式 region-3 `U3/V3` thermal-state evolution；VegasAfterglow 只是 comparison backend。

## AM3 参考

本地: `/mnt/c/Users/jia/Documents/New project/_external/am3_reference` (HEAD: `7aba970b`)

用途: 微观核、过程拆分、benchmark 参考。不替换 ASGARD 的 dynamics/electron/observer 主链。
