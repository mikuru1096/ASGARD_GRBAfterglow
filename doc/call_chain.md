# ASGARD 调用链

## Python 编排层

```mermaid
flowchart TD
    A["ASGARD/api_model.py\nModel / observe / run_fit"] --> B["FitConfig / SimulationSetup"]
    B --> C["asgard_fit.py\ncompile_problem / eval_loglike"]
    B --> D["asgard_state.py\nsolve_state_from_setup"]
    D --> E["求解动力学"]
    D --> F["求解电子谱"]
    D --> G["求解反激波辐射\nRS 电子 + 可选 RS 质子同步"]
    D --> H["构建光子场阶段"]
    D --> I["求解强子过程"]
    I --> J["BH 次级电子并入正激波电子谱\n并重算 seed_syn"]
    I --> K["pγ 光子生存因子\n写回光子场"]
    I --> Q["pair production / 可选迭代 cascade"]
    H --> L["Radiation.annihilation\ngamma-gamma 吸收"]
    Q --> L
    K --> L
    J --> L
    G --> L
    L --> M["project_flux_grid / project_spec"]
    M --> N["Interpolation.sed_interpolation\nEATS + Doppler + redshift"]
    N --> O["combine_multiband_flux"]
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
  -> solve_state_from_setup
  -> solve_dynamics -> solve_electron -> photon_field_stage
  -> solve_hadronic -> solve_reverse_shock_emission
  -> pair-production branch / Radiation.annihilation -> Interpolation.sed_interpolation
  -> combine_multiband_flux -> FluxResult
```

拟合最短路径：

```text
Fitter.loglike -> eval_loglike -> solve_state_from_setup
  -> project_flux_grid -> combine_multiband_flux -> compute_light_curve_redchi
```
