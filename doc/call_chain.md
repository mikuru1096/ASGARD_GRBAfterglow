# ASGARD 调用链

本页给出 Python 编排和 Fortran 数值核的主链。若要逐个进入 Fortran `subroutine` / `function`，使用 `doc/fortran_kernel_index.md`；若要按物理过程反查算法阶段，使用 `doc/physics_algorithm_crosswalk.md`。

## Python 编排层

```mermaid
flowchart TD
    A["asgard_core/api_model.py\nModel query methods"] --> B["Model -> RuntimeConfig\n_build_fit_config_for_patch"]
    X["asgard_core/api_observe.py\nobserve / run_fit"] --> B
    Y["asgard_core/api_fit.py\nFitter.loglike / eval_loglike"] --> S
    B --> S["asgard_setup.py\nSimulationSetup"]
    S --> D["asgard_state.py\nsolve_state_from_setup"]
    D --> E["求解动力学"]
    D --> F["求解电子谱"]
    D --> G["求解反向激波辐射\nRS 电子 + 可选 RS 质子同步"]
    D --> H["构建光子场阶段"]
    D --> I["求解强子过程"]
    I --> J["BH 次级电子并入正向激波电子谱\n并重算 seed_syn"]
    I --> K["pγ 光子生存因子\n写回光子场"]
    I --> Q["pair production / shell-sequence cascade"]
    H --> L["Radiation.annihilation\ngamma-gamma 吸收"]
    Q --> L
    K --> L
    J --> L
    G --> L
    L --> M["project_flux_grid\nprojection_kind"]
    M --> N["top-hat chi_eats_2d lightcurve:\nInterpolation.sed_interpolation_chi"]
    M --> T["axisymmetric structured chi_eats_2d:\nring solve + sed_chi_ring"]
    M --> R["sed / shell components:\nInterpolation.sed_interpolation\nEATS + Doppler + redshift"]
    N --> O["combine_multiband_flux"]
    T --> O
    R --> O
    O --> P["compute_light_curve_redchi"]
```

## Fortran 数值核层

```mermaid
flowchart TD
    A["solve_dynamics"] --> B["Dynamics_forward / Dynamics_reverse"]
    B --> C["solve_electron"]
    C --> D["electron_forward_fullhide_1d / charint_1d / weno5_1d / t2g1_1d / slc1_1d"]
    C --> E["electron_forward_transport_2d / charint_2d"]
    D --> F["electron_common + electron_cooling_kernel assembly\nSSA / IC / Y cooling kernels + radiation kernel"]
    E --> G["electron_transport_2d_kernel + electron_seed_history_kernel"]
    F --> H["rad_common -> ssc_spectrum / pair_absorption"]
    G --> H
    H --> I["SED_interpolation\nsed / adaptive / chi / structured ring-precomputed"]
    H --> U["SED_interpolation_structured\ninternal structured_jet_1d shell projection"]
    C --> J["solve_hadronic -> hadronic_forward_1d"]
    J --> K["hadronic_transport + hadronic_rad + hadronic_pg + decay + pp + BH + IC + pair + species + accel + secondary"]
```

## 有效主线

```text
Model.flux_density_grid
  -> _build_fit_config_for_patch / _solve_patch_state
  -> SimulationSetup
  -> solve_state_from_setup
  -> solve_dynamics -> solve_electron -> photon_field_stage
  -> solve_hadronic -> solve_reverse_shock_emission
  -> pair-production branch / Radiation.annihilation -> project_flux_grid
  -> projection_kind="lightcurve" or "sed"
  -> combine_multiband_flux -> FluxResult
```

拟合最短路径：

```text
Fitter.loglike -> eval_loglike -> solve_state_from_setup
  -> project_flux_grid -> combine_multiband_flux -> compute_light_curve_redchi
```
