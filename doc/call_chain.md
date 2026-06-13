# ASGARD 调用链

## Python 编排层

```mermaid
flowchart TD
    A["ASGARD/api_model.py\nModel query methods"] --> B["Model -> FitConfig\n_build_fit_config_for_patch"]
    X["ASGARD/api_observe.py\nobserve / run_fit"] --> B
    Y["ASGARD/api_fit.py\nFitter.loglike"] --> C["asgard_fit.py\ncompile_problem / eval_loglike"]
    B --> S["asgard_setup.py\nSimulationSetup"]
    C --> S
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
    M --> N["lightcurve:\nInterpolation.sed_interpolation_chi\nfor chi_eats_2d FS synch+SSA"]
    M --> R["sed / shell components:\nInterpolation.sed_interpolation\nEATS + Doppler + redshift"]
    N --> O["combine_multiband_flux"]
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
    D --> F["electron_common + electron_cooling_kernel + electron_radiation_kernel"]
    E --> G["electron_transport_2d_kernel + electron_seed_history_kernel"]
    F --> H["radiation_common -> radiation_ssc_spectrum / radiation_gamma_gamma_absorption / radiation_reverse_seed"]
    G --> H
    H --> I["SED_interpolation / SED_interpolation_structured"]
    C --> J["solve_hadronic -> hadronic_forward_1d"]
    J --> K["hadronic_transport + radiation + interaction + decay + pp + BH + IC + pair_prod + species_transport + acceleration + secondary_radiation"]
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
