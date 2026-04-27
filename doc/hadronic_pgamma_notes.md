# Hadronic `pγ` Notes

## Canonical Naming

- `hummer_2010_response` — 正式 pγ 方案，含 transport feedback
- `ka2008_reference` — 参考后端，仅 emission benchmark，不参与 proton transport feedback
- `disabled`

兼容别名: `hummer_2010`/`hummer2010`/`am3_reference`/`am3` → `hummer_2010_response`; `ka2008`/`aharonian_2008`/`kelner_aharonian_2008` → `ka2008_reference`

## Current Formal Coupling

主链: `dynamics → electron solver → photon target field → proton transport → hadronic source/loss operators → observer projection`

已接入:
- proton injection, adiabatic cooling, synchrotron cooling, proton synchrotron emission
- `hummer_2010_response` photopion: `α_p(E)` + `Q_p^reinj(E)` 反馈到 proton transport
- `α_γ^{pγ}` 转为 local shell photon survival factor: `τ_{pγ}(ν,r) = α_γ^{pγ}(ν,r) · R/(12 Γ c)`, `f_surv = (1-e^{-τ})/τ`, 在 observer projection 前作用
- neutrino output
- Bethe-Heitler: proton continuous cooling + secondary e± → merge into forward electron → recompute seed_syn
- hadronic IC: proton channel; pion/muon IC via explicit secondary species
- pp: γ, ν, pair source, proton loss
- pair production: observer-side `tau_pair` attenuation + pair synchrotron branch
- explicit secondary species transport: n, π±, μ± (left/right)
- secondary radiation: pion/muon synchrotron + IC
- hadronic acceleration/injection operators

## Current Boundaries

未完成:
- reverse-shock hadronic
- 2D / χ-resolved hadronic transport
- full time-dependent pair cascade PDE

## AM3 Reference

本地: `/mnt/c/Users/jia/Documents/New project/_external/am3_reference` (HEAD: `7aba970b`)

用途: 微观核、过程拆分、benchmark 参考。不替换 ASGARD 的 dynamics/electron/observer 主链。Python 侧参考实现仅作过渡层。
