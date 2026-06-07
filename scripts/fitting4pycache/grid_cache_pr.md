# Reviewer Note: GRB 180720B fitting-grid cache in `grb180720b_multinest.py`

## Summary

This change adds a Python-side cache for the **GRB 180720B fitting data grid** used by `scripts/fitting/grb180720b_multinest.py`.

The cache is intentionally source-specific: it is built from the GRB 180720B data files, GRB 180720B redshift, fixed band definitions, XRT absorption treatment, VHE EBL treatment, and the article-style fitting setup in this script. It is **not** a general ASGARD cache API and should not be assumed to apply to other GRBs without refactoring the data-loading and grid-building layer.

## What is cached

For every likelihood call, the physical model still solves a new light-curve/SED grid from the sampled parameters. What is cached is the observation-side fitting scaffold that does not change between likelihood calls:

- loaded GRB 180720B data rows from XRT, LAT, VHE, optical, and radio files;
- the unique observer frequency grid required by these data;
- the unique observer-time grid required by these data;
- per-row sampling specifications for narrow-band and wide-band measurements;
- XRT absorption and VHE EBL transmission values evaluated on the fixed bandpass sample grids;
- precomputed fit arrays:
  - direct frequency/time indices for monochromatic flux points;
  - frequency index arrays, time indices, and integration weights for wide-band flux points;
  - observed flux and uncertainty arrays used by the likelihood.

In the code this is centered around:

```python
build_observation_grid()
cached_observation_grid()
_build_sample_specs()
_build_fit_data()
_collect_fit_arrays()
model_flux_from_total_grid()
```

## How likelihood evaluation changes

Before this cache, each likelihood evaluation had to rebuild or re-derive the observation-side mapping from the model grid to the fitting data.

Now each likelihood call does the following:

1. Reuse `cached_observation_grid()` for the fixed GRB 180720B observation grid and fit metadata.
2. Build a new `FitConfig` from the sampled parameters.
3. Run the physical solver normally through `solve_state_from_setup()`.
4. Project the newly computed model onto the cached observer grid.
5. Use `model_flux_from_total_grid()` to extract the fitted model flux vector through cached indices and weights.
6. Compute chi-square/log-likelihood against cached observed flux and uncertainty arrays.

This reduces repeated Python bookkeeping during PyMultiNest sampling while leaving the physical model calculation unchanged.

## Scope and limitations

- This is a **per-process in-memory cache** via `functools.lru_cache`, not a cross-job on-disk cache.
- It only persists for the lifetime of the Python process or MPI rank.
- It is valid because the GRB 180720B fitting data, band definitions, absorption tables, EBL table, and observer sampling grid are fixed during one run.
- It does not cache sampled physical parameters, model light curves, spectra, or Fortran solver outputs.
- It does not change priors, parameter dimensionality, physical kernels, Fortran code, or numerical protection behavior.
- It currently applies only to the GRB 180720B fitting driver. Generalizing it to other sources would require moving the data-grid cache behind a source-aware interface.

## Timing instrumentation

The script also records lightweight timing labels around the cached grid and projection stages, for example:

```text
fit.cached_observation_grid
fit.config_from_params
fit.build_simulation_setup
fit.solve_state_from_setup
fit.project_flux_grid
fit.model_flux_from_total_grid
fit.loglike_from_model_flux
```

These labels are intended to make likelihood bottlenecks easier to inspect during local, Slurm, or MPI runs.

## Suggested PR wording

Title:

```text
Cache GRB 180720B fitting grid across likelihood calls
```

Body:

```text
This PR adds a Python-side in-memory cache for the GRB 180720B observation/fitting grid in scripts/fitting/grb180720b_multinest.py.

The cached objects are source-specific fitting scaffolds: loaded GRB 180720B data rows, observer frequencies/times, bandpass sample specs, absorption/EBL transmission values, and precomputed indices/weights used to map model grids onto the fit vector.

Each likelihood call still computes a fresh physical model from the sampled parameters. The cache only avoids rebuilding invariant observation-side bookkeeping. This keeps the Fortran solver, physical model, priors, and numerical behavior unchanged.

This cache is currently specific to GRB 180720B and should not be treated as a general ASGARD source cache without a separate source-aware refactor.
```

Validation:

```bash
python -m py_compile scripts/fitting/grb180720b_multinest.py
python -u scripts/fitting/grb180720b_multinest.py fullhide_1d_pic --benchmark-fit --repeat 3
```
