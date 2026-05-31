# Public API Reference

本文档记录当前 public API 的稳定入口和语义边界。实现文件主要是 `ASGARD/api_model.py`, `ASGARD/api_observe.py`, `ASGARD/api_fit.py`。

## 导入入口

```python
from ASGARD import (
    Model,
    ISM,
    Wind,
    TophatJet,
    GaussianJet,
    PowerLawJet,
    TwoComponentJet,
    StepPowerLawJet,
    Ejecta,
    Observer,
    Radiation,
    Setups,
    Fitter,
    Param,
    Scale,
    observe,
    run_fit,
    units,
)
```

## Medium

### `ISM`

```python
ISM(n_ism=1.0)
ISM(n0=1.0)
```

Kernel parameters:

- `d_ne = n_ism`
- `a_star = -1`

### `Wind`

```python
Wind(A_star=1.0, n_ism=0.1, n0=None, k=2.0)
```

Current backend only supports `k=2.0`. If `n0` is set, wind density is capped by the corresponding transition radius.

### `Medium`

`Medium(rho=callable)` can evaluate density in Python, but current Fortran kernel dispatch does not support arbitrary user-defined media. See `doc/public_backend_limits.md`.

## Jet

### `TophatJet`

```python
TophatJet(E_iso=1.0e52, Gamma0=300.0, theta_j=0.1)
TophatJet(E_iso=1.0e52, lf=300.0, theta_j=0.1)
```

Supported fields:

- `E_iso`
- `lf` / `Gamma0`
- `theta_j` / `theta_c`
- `duration`
- `magnetar`
- `spreading`

`spreading=True` is accepted at the object level but current backend does not implement jet spreading dynamics.

### Structured jets

```python
GaussianJet(E_iso=1.0e52, Gamma0=300.0, theta_c=0.05, theta_max=0.5)
PowerLawJet(E_iso=1.0e52, Gamma0=300.0, theta_c=0.05, k=2.0, theta_max=0.5)
TwoComponentJet(E_iso_c=1.0e52, Gamma0_c=300.0, E_iso_outer=1.0e50, Gamma0_outer=30.0, theta_c=0.05, theta_o=0.2)
StepPowerLawJet(...)
Ejecta(...)
```

The public model evaluates patches through `energy_iso(phi, theta)` and `gamma0(phi, theta)`.

## Observer

```python
Observer(z=0.1, theta_obs=0.0, phi_obs=0.0)
Observer(lumi_dist_cm=1.0e28, z=0.1, theta_obs=0.0)
```

If luminosity distance is not provided, cosmology utilities determine it from redshift in the runtime path.

## Radiation

```python
Radiation(
    eps_e=0.1,
    eps_B=1.0e-3,
    p=2.3,
    xi_N=0.1,
    ssc=True,
)
```

Main electron switches:

- `eps_e`
- `eps_B`
- `p`
- `xi_N`
- `thermal_electrons`
- `ssc`
- `kn`
- `epsilon_b_floor`
- `magnetic_decay_alpha_t`
- `magnetic_decay_t0_s`

Hadronic switches:

- `epsilon_p`
- `proton_synch`
- `pg`
- `bethe_heitler`
- `hadronic_inverse_compton`
- `pp`
- `neutrino`
- `eta_acc`
- `pgamma_scheme`
- `pair_production`
- `reverse_epsilon_p`

## Setups

`Setups` controls numerical grids, solver selection, RS flags, hadronic flags and observer time range.

Important fields:

- `num_r`, `num_theta`, `num_phi`, `num_tobs`
- `num_gam_e`, `num_nu`, `num_chi`
- `observer_time_min_s`, `observer_time_max_s`
- `electron_solver`
- `index_y`
- `index_syn_integr`
- `ssc_cooling`
- `rvs_shock`, `rvs_ssc`, `include_cross_zone_ic`
- `reverse_delta_t_s`, `reverse_sigma`
- `hadronic_enabled`, `hadronic_solver`, `num_gam_p`, `num_nu_nu`
- `pair_cascade_iterations`
- `num_threads`

Solver aliases:

- `fullhide` -> `fullhide_1d`
- `slc1` -> `slc1_1d`
- `charint` -> `charint_1d`
- `t2g1` -> `t2g1_1d`
- `weno5` -> `weno5_1d`

Registered solver names:

- `fullhide_1d`
- `slc1_1d`
- `charint_1d`
- `t2g1_1d`
- `weno5_1d`
- `fullhide_2d`
- `charint_2d`

## Model

Construction:

```python
model = Model(jet=jet, medium=medium, observer=observer, fwd_rad=fwd_rad, rvs_rad=rvs_rad, setups=setups)
```

Supported positional compatibility forms:

- `Model(jet, medium, observer, radiation)`
- `Model(medium, jet, observer, radiation)`

### `flux_density_grid(times_s, nu_hz)`

Full grid projection. `times_s` and `nu_hz` are one-dimensional arrays; output `total` has shape `(num_nu, num_time)`.

### `flux_density(times_s, nu_hz)`

Pointwise projection for matched `(time_i, nu_i)` arrays, or extraction from a grid when shapes differ.

### `flux_density_exposures(times_s, nu_hz, exposures_s, num_subsamples=4)`

Exposure-averaged flux density. Input arrays must share the same shape.

### `spectrum(time_s, nu_hz)`

Single-time spectrum convenience wrapper.

### `flux(time_s, nu_min_hz, nu_max_hz, num_points=64)`

Band-integrated flux using log frequency sampling and trapezoid integration.

### `sky_image(t_obs, nu_obs, fov, npixel=128)`

Observer-plane sky image. `fov` is in angular units matching the `units` conversion used by caller.

### `polarization(times_s, nu_hz, magnetic_geometry="shock_random", local_emissivity="analytic_then_kernel")`

Returns Stokes and polarization fraction for synchrotron components.

### `details(t_min=None, t_max=None)`

Returns `ModelDetails`: dynamics tracks, branch diagnostics, characteristic frequencies and internal states.

## Fitter

```python
from ASGARD import Fitter, Param, Scale

fitter = Fitter(model=model)
fitter.add_flux_density(times_s, nu_hz, flux, err)
fitter.params = [
    Param("logE", "jet.E_iso", 50.0, 54.0, Scale.LOG10),
]
```

`Param` forms:

- `Param(name, lower, upper, scale)`
- `Param(name, path, lower, upper, scale)`

`Scale.LOG` and `Scale.LOG10` transform by `10**value`; `Scale.FIXED` always uses `lower`.

Sampling helpers:

- `Fitter.run_emcee(initial, nwalkers, nsteps)`
- `Fitter.run_multinest(...)`

## observe / run_fit

`observe(model, config=..., spectrum_output=...)` is the lower-level execution entry used by demos and scripts. `Model` methods are preferred for ordinary interactive use.

`run_fit(config)` is retained for compatibility with config-style workflows.

## Return Types

Common result objects:

- `ModelFluxResult`: `total`, `fwd`, `rev`, `cross_ic`
- `BranchView`: `sync`, `ssc`
- `SkyImage`: image array and coordinate grids
- `PolarizationResult`: Stokes and polarization diagnostics
- `ModelDetails`: dynamics/electron/hadronic/observer diagnostic state

## Public Boundary

Public API accepts some options before the backend supports them. Unsupported options should fail explicitly rather than silently falling back. Current fixed boundaries are in `doc/public_backend_limits.md`.
