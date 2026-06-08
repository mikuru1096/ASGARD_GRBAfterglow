# fullhide_2d PWN/CR Transport Extension Notes

## Scope

This note fixes the physics contract for the optional `fullhide2d_transport_model="pwn_cr_v1"` path. The default
`"legacy"` path keeps the existing 27-element `Boundary` array and the existing fullhide_2d numerical operators.
The PIC backend is outside this implementation.

## Boundary Layout

The existing common unpacker uses the modern 27-element boundary layout with `Boundary(27)=R0`. The transport
extension therefore appends new fields after `R0`:

- `Boundary(28)`: `transport_model_selector`, `0=legacy`, `1=pwn_cr_v1`.
- `Boundary(29)`: `stochastic_accel_norm`, dimensionless momentum-diffusion strength, default `0`.
- `Boundary(30)`: `escape_mode_selector`, `0=closed`, `1=free_outer`.

This deliberately differs from the original draft numbering because occupying `Boundary(27)` would overwrite the
external-density radius scale.

## Transport Equation Basis

PWN multizone models solve a Fokker-Planck transport equation with spatial transport, radiative and adiabatic
losses, and injection. Van Rensburg et al. describe a time-dependent multizone leptonic PWN model with spatially
dependent magnetic field, convection, diffusion, radiative losses, and adiabatic losses. The newer 3D PWN/CR work
uses the same transport-equation structure in a fully spatial setting. GALPROP uses the same cosmic-ray transport
taxonomy, including spatial diffusion, energy losses, source terms, optional momentum diffusion/reacceleration, and
free-escape halo boundaries.

Implementation decision: `pwn_cr_v1` keeps the existing fullhide_2d source and χ advection structure, then adds the
PWN/CR pieces that are meaningful for the current electron-only forward-shock solver: χ-space diffusion boundary
control, χ-local adiabatic losses, χ-local magnetic cooling/diffusion, and optional log-gamma diffusion.

## BM Local Adiabatic Loss

The legacy energy operator used a uniform-expansion coefficient `d log10(gamma) / dR = 1 / (R ln 10)`. The PWN/CR
path computes the local radial divergence from the BM downstream velocity profile instead of differencing adjacent
cell volumes.

For a shock radius `R`, shock Lorentz factor `Gamma_sh`, and BM coordinate `chi`,
`Gamma_2(chi)=Gamma_sh/sqrt(2 chi)` and `beta_2(chi)=sqrt(1-1/Gamma_2^2)` in the BM-valid ultrarelativistic
downstream region. With `r/R = 1 - (chi-1)/(8 Gamma_sh^2)`, the radial velocity divergence is evaluated as
`div(v)/c = 2 beta_2/r_ratio/R + d beta_2/dr`, using `d beta_2/dr = 8/(beta_2 R)` from the BM coordinate
definition. The energy-loss coefficient advanced per shock radius is then `div(v)/(3 beta_sh c ln 10)`.

The BM expression is used only where the local downstream flow remains ultrarelativistic, implemented as
`Gamma_2 >= 2`. Farther downstream the BM ultrarelativistic velocity-gradient expression is outside its domain and
would diverge near `Gamma_2 -> 1`; those cells use the uniform isotropic expansion closure `1/(R ln 10)`. This also
keeps the required `div(v)/3 = 1/R` limit explicit.

## Eta-Space Diffusion

The physical diffusion flux is `F_x = -kappa dN/dx`. The existing fullhide_2d geometry maps physical downstream
distance to `eta=log10(chi)` through the same BM shock metric used by the χ advection operator. `pwn_cr_v1` keeps
the tridiagonal solve in eta space, with coefficients derived from the metric factor already used in
`advance_eta_logchi_implicit`, avoiding physical-grid remapping between steps.

Boundary decisions:

- `closed`: zero diffusive flux at the outer eta face.
- `free_outer`: zero outside density at the outer eta face, producing an outward sink term in the outermost cell.

## Microturbulence Closure

`pwn_cr_v1` reuses the existing `epsilon_b_floor`, `magnetic_decay_alpha_t`, and `magnetic_decay_t0_s` closure.
The same local `B_chi` is used for synchrotron/SSA diagnostics, radiative cooling, and Bohm-like diffusion
coefficient evaluation. No PIC hit probability, patchy dual-field emission, or PIC scattering operator is added.

This follows the GRB-afterglow motivation that post-shock microturbulence can decay behind the shock, as in
Lemoine's microturbulence-decay treatment, while keeping the implementation tied to existing fullhide_2d fields.

## Stochastic Acceleration

`stochastic_accel_norm=0` disables this operator. For positive values, the code applies Strang splitting:

1. half-step conservative diffusion in log gamma,
2. full-step radiative plus adiabatic cooling,
3. half-step conservative diffusion in log gamma.

The log-gamma diffusion substep has zero flux at the energy-grid boundaries. Its implicit coefficient is
`dR * stochastic_accel_norm / (R * dlog10(gamma)^2)`, so `stochastic_accel_norm` is a dimensionless research-mode
normalization of diffusion per dynamical radius, not a standalone physical `D_pp`. With `stochastic_accel_norm=0`,
this path is not called and the result is identical to the no-stochastic `pwn_cr_v1` path.

## References

- Van Rensburg, Krüger & Venter, 2018, MNRAS 477, 3853, spatially dependent PWN transport:
  https://academic.oup.com/mnras/article/477/3/3853/4956543
- 3D PWN/CR transport equation context:
  https://academic.oup.com/mnras/article-abstract/528/2/2749/7529207
- Lemoine, 2013, GRB afterglow microturbulence decay:
  https://academic.oup.com/mnras/article/428/1/845/1062346
- GALPROP theory page, CR diffusion/reacceleration/escape taxonomy:
  https://galprop.stanford.edu/code.php?option=theory
- BM/chi afterglow profile context:
  https://academic.oup.com/mnras/article/442/4/3495/1339065
