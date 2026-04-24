# Hadronic `pγ` Notes

本文件只记录当前仓库已经落地的 `pγ` / hadronic 状态，不保留历史计划。
其中需要特别区分两层：

- 当前正式运行路径
- 当前过渡性的 Python 参考实现

按当前项目约束，所有直接从 AM3 / Hümmer / KA2008 转写的最终微观物理实现，后续都必须迁到 `src/Hadronic/` 的 Fortran 模块里。

## Canonical Naming

- `hummer_2010_response`
  - 当前正式 `pγ` 方案名
  - 当前仓库里仍有 Python 过渡实现
  - 最终目标是迁入 `src/Hadronic/` 的 Fortran 后端
- `ka2008_reference`
  - 当前参考后端
  - 只做 KA2008 stable-secondary benchmark，不参与正式质子输运反馈
- `disabled`

兼容 alias 仍允许，但文档和测试标签应优先使用 canonical 名：

- `am3_reference` / `am3_numeric` / `am3_numerical` / `am3` -> `hummer_2010_response`
- `hummer_2010` / `hummer2010` -> `hummer_2010_response`
- `ka2008` / `aharonian_2008` / `kelner_aharonian_2008` -> `ka2008_reference`

## Current Formal Coupling

当前正式 1D hadronic 主链是：

`dynamics -> electron solver -> photon target field -> proton transport -> hadronic source/loss operators -> observer projection`

已经接入正式流程的项：

- proton injection
- proton adiabatic cooling
- proton synchrotron cooling
- proton synchrotron emission
- `hummer_2010_response` photopion source
- neutrino output
- proton transport feedback:
  - `alpha_p(E)`
  - `Q_p^{reinj}(E)`
- `alpha_gamma^{pγ}` observer-side attenuation:
  \[
  \tau_{p\gamma}(\nu,r)=\alpha_\gamma^{p\gamma}(\nu,r)\,\frac{R}{12\Gamma c}
  \]
  并通过 `tau_extra` 接入 `Radiation.annihilation`
- Bethe-Heitler proton cooling
- Bethe-Heitler secondary `e^\pm` 注入
- BH 次级并回前向电子分布后重算 `seed_syn`
- hadronic inverse Compton
  - 当前正式 runtime 中 proton 通道活跃
  - pion/muon IC 需要显式演化这些分布后才会非零
- `pp` gamma / neutrino / proton-loss
- pair production observer-side attenuation 与 pair synchrotron 支路

## Strict Reference Layer

`asgard_core/hadronic_pgamma.py` 目前仍保留严格参考积分：

\[
t_{p\gamma}^{-1}(\gamma_p)=
\frac{c}{2\gamma_p^2}
\int_{\bar\epsilon_{\rm thr}}^\infty d\bar\epsilon\;
\sigma_{p\gamma}(\bar\epsilon)\kappa_{p\gamma}(\bar\epsilon)\bar\epsilon
\int_{\bar\epsilon/(2\gamma_p)}^\infty d\epsilon\;
\frac{n(\epsilon)}{\epsilon^2}.
\]

这里：

- \(\epsilon = h\nu\)
- \(\bar\epsilon_{\rm thr} \simeq 145\,{\rm MeV}\)
- `photon_density_hz_to_gev(...)` 固定执行
  \[
  n_\epsilon(\epsilon)=\frac{n_\nu(\nu)}{h_{\rm GeV}}
  \]

`ka2008_reference` 当前的作用仍然是：

- reference benchmark
- 检查阈值、分段、稳定二级粒子参数化
- 不作为正式 transport-coupled backend
- 也不应被视为最终 AM3-derived 实现

## Current Boundaries

当前还没有完成的项：

- reverse-shock hadronic
- 2D / `chi`-resolved hadronic transport
- pion / muon 显式输运与它们在正式 runtime 中的 IC 非零耦合
- 全 hadronic cascade 的自洽 PDE 闭环

## Local AM3 Reference

本地参考库：

- `/mnt/c/Users/jia/Documents/New project/_external/am3_reference`

用途边界：

- 用 AM3 做微观核、过程拆分和 benchmark 参考
- 不用 AM3 的 one-zone PDE 结构替换 ASGARD 现有 electron solver / observer chain
- Python 侧参考移植只允许作为过渡层，不能作为最终完成态
