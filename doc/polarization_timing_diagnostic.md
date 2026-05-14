# Polarization Timing Diagnostic

本文记录 Lan, Wu & Dai 2023 polarization overlay 中“峰值幅度基本一致、峰时偏早”的定位结果。当前结论：**不改 Stokes 累加、不做后处理平滑；下一步若要推进，应先做 dynamics / jet-evolution benchmark**。

## Reproduction

运行命令：

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/polarization_literature_overlay.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/polarization_smoke.py'
```

输出图和 CSV：

- `output/asgard_doc/polarization_lan2023_overlay.png`
- `output/asgard_doc/polarization_lan2023_overlay_data.csv`

当前 CSV 的峰值：

- Lan 2023 digitized curve: `t_peak = 2.515752e4 s`, `PD_peak = 46.12245%`
- ASGARD 32-point overlay: `t_peak = 1.34596032e4 s`, `PD_peak = 46.65525%`
- Ratio: `t_peak,ASGARD / t_peak,Lan = 0.535`, `PD_peak,ASGARD / PD_peak,Lan = 1.012`

更密 64 点 ASGARD polarization-only 采样给出 `t_peak ~ 1.55e4 s`，所以峰时偏早的稳健量级是约 `0.54--0.62` 倍，而不是单点采样误差。

## Projection Checks

当前 polarization 使用 `sed_interpolation_surface_element(...)` 做每个真实球面面元的 EATS + Doppler 投影：

```text
R_Tobs_theta = R_Tobs1 + R * (1 - cos theta_view_patch) * (1 + z) / c
flux weight = dOmega / 4pi
```

已有 `tests/polarization_smoke.py::test_surface_element_projection_scales_with_solid_angle` 验证面元 flux 对 `dOmega` 严格线性。`tests/polarization_smoke.py` 当前通过，说明没有发现 patch solid-angle 权重错误。

动力学返回的 on-axis time 已包含红移：`Dynamics_forward.f90` 中 `R_Tobs = T * (1 + z)`；surface-element projection 只额外加入 off-axis angular delay。这个结构与 Lan 2023 摘要强调的 EATS 会把 off-axis PD bumps 推迟的物理方向一致。

## Dynamics Clue

Lan overlay 参数中 `theta_j = 0.1`、`theta_obs = 0.2`。用 ASGARD 当前 dynamics 插值：

- ASGARD peak `1.35e4 s`: `Gamma ~ 4.02`, `1/Gamma ~ 0.249`
- ASGARD dense-grid peak `1.55e4 s`: `Gamma ~ 3.83`, `1/Gamma ~ 0.261`
- Lan peak `2.52e4 s`: `Gamma ~ 3.24`, `1/Gamma ~ 0.309`

偏移对应的是 blast-wave deceleration / jet-evolution time scale，而不是局部 polarization fraction 的幅度问题。ASGARD 当前 top-hat patch 是固定开角、独立 spherical wedge 的 dynamics；项目边界中 jet spreading 仍未支持。文献中对 polarization 的讨论也明确指出 EATS 会影响 off-axis PD bump time，而包含 lateral expansion 的 hydrodynamic jet 与无 lateral expansion 的 analytical wedge 会给出不同的 polarization peak morphology and timing。

## Decision

当前不修改 polarization Stokes 累加、不调整 `theta_obs`、不乘经验时间因子，也不做 smoothing。峰值幅度已经匹配，直接动 polarization normalization 会破坏已通过的物理量。

后续若要解决峰时，应先做一个 dynamics benchmark：

1. 固定 Lan 2023 overlay 参数，输出 `Gamma(t)`, `R(t)`, on-axis `R_Tobs`, angular-delay term, total EATS time。
2. 与 Blandford-McKee / jet-break analytic scale 对照，明确 ASGARD 的 `E_iso`, `n`, `z`, `theta_j`, `theta_obs` 是否与文献图完全同义。
3. 若 mismatch 来自缺失 jet spreading，则新增 jet spreading dynamics 任务；不要在 polarization projection 层补偿。
4. 若 mismatch 来自文献参数映射或能量定义，则修正 benchmark setup，而不是修改内核。

## References

- Lan, Wu & Dai 2023, ApJ 952, 31: https://doi.org/10.3847/1538-4357/acd6ef
- MNRAS 523, 4583--4592: https://academic.oup.com/mnras/article/523/3/4583/7202343
