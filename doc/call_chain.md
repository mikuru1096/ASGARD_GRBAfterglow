# ASGARD Call Chain

## Python Orchestration Layer

```mermaid
flowchart TD
    A["ASGARD/api_model.py\nModel / observe / run_fit"] --> B["FitConfig / SimulationSetup"]
    B --> C["asgard_fit.py\ncompile_problem / eval_loglike"]
    B --> D["asgard_state.py\nsolve_state_from_setup"]
    D --> E["solve_dynamics"]
    D --> F["solve_electron"]
    D --> G["solve_reverse_shock_emission\nRS electron + optional RS proton synch"]
    D --> H["photon_field_stage"]
    D --> I["solve_hadronic"]
    I --> J["BH e± merge → recompute seed_syn"]
    I --> K["pγ photon survival → photon field"]
    I --> Q["pair production / optional iterative cascade"]
    H --> L["Radiation.annihilation\ngamma-gamma absorption"]
    Q --> L
    K --> L
    J --> L
    G --> L
    L --> M["project_flux_grid / project_spec"]
    M --> N["Interpolation.sed_interpolation\nEATS + Doppler + redshift"]
    N --> O["combine_multiband_flux"]
    O --> P["compute_light_curve_redchi"]
```

## Fortran Kernel Layer

```mermaid
flowchart TD
    A["solve_dynamics"] --> B["Dynamics_forward / Dynamics_reverse"]
    B --> C["solve_electron"]
    C --> D["FS_electron_fullhide_1d / charint_1d / weno5_1d / t2g1_1d / slc1_1d"]
    C --> E["FS_electron_fullhide_2d / charint_2d"]
    D --> F["electron_common + electron_cooling_kernel + electron_radiation_kernel"]
    E --> G["electron_transport_2d_kernel + electron_seed_history_kernel"]
    F --> H["radiation_common → SSC_spec / Annihilation / Seed_reverse"]
    G --> H
    H --> I["SED_interpolation / SED_interpolation_structured"]
    C --> J["solve_hadronic → FS_hadronic_1d"]
    J --> K["hadronic_transport + radiation + interaction + decay + pp + BH + IC + pair_prod + species_transport + acceleration + secondary_radiation"]
```

## Effective Mainline

```text
Model.flux_density_grid
  → solve_state_from_setup
  → solve_dynamics → solve_electron → photon_field_stage
  → solve_hadronic → solve_reverse_shock_emission
  → pair-production branch / Radiation.annihilation → Interpolation.sed_interpolation
  → combine_multiband_flux → FluxResult
```

Fit shortest loop:

```text
Fitter.loglike → eval_loglike → solve_state_from_setup
  → project_flux_grid → combine_multiband_flux → compute_light_curve_redchi
```
