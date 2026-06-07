# 强子 `pγ` 说明

## 标准命名

- `hummer_2010_response`：正式 pγ 方案，含 transport feedback。
- `ka2008_reference`：参考后端，仅 emission benchmark，不参与 proton transport feedback。
- `disabled`

兼容别名: `hummer_2010`/`hummer2010`/`am3_reference`/`am3` → `hummer_2010_response`; `ka2008`/`aharonian_2008`/`kelner_aharonian_2008` → `ka2008_reference`

## 当前正式耦合

主链: `dynamics → electron solver → photon target field → proton transport → hadronic source/loss operators → observer projection`

已接入:
- proton injection、adiabatic cooling、synchrotron cooling、proton synchrotron emission。
- `hummer_2010_response` photopion：`α_p(E)` + `Q_p^reinj(E)` 反馈到 proton transport。
- `α_γ^{pγ}` 转为 local shell photon survival factor：`τ_{pγ}(ν,r) = α_γ^{pγ}(ν,r) · R/(12 Γ c)`, `f_surv = (1-e^{-τ})/τ`，并在 observer projection 前作用。
- neutrino output。
- Bethe-Heitler：proton continuous cooling + secondary e±，随后并入 forward electron 并重算 `seed_syn`。
- hadronic IC：proton channel；pion/muon IC 通过 explicit secondary species 处理。
- pp：gamma、neutrino、pair source、proton loss。
- pair production：observer-side `tau_pair` attenuation + pair synchrotron branch；`pair_cascade_iterations > 1` 选择 shell-sequence time-dependent γγ pair/synch cascade。
- explicit secondary species transport：n、π±、μ± 的左右向输运。
- secondary radiation：pion/muon synchrotron + IC。
- hadronic acceleration/injection operators。
- reverse-shock full-chain dispatch：开启 RS hadronic full-chain flags 时，runtime 复用 formal 1D kernels，并使用 RS seed photons、RS magnetic field、RS shell energy 和 RS baryon target density。

## 当前边界

Hadronic 未完成项集中维护在根目录 `TODO.md`。本文件只记录当前 pγ 耦合、命名和 AM3 参考边界。

Reverse-shock baseline：ASGARD 使用局部 `gamma34` 注入和显式 region-3 `U3/V3` thermal-state evolution；VegasAfterglow 只是 comparison backend。

## AM3 参考

本地: `/mnt/c/Users/jia/Documents/New project/_external/am3_reference` (HEAD: `7aba970b`)

用途: 微观核、过程拆分、benchmark 参考。不替换 ASGARD 的 dynamics/electron/observer 主链。Python 侧参考实现仅作过渡层。
