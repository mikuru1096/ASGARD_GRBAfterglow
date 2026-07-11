# PION Eta Car CSM at 1869.597 AD

This directory stores the 1024-cell spherical one-dimensional reduction of the public PION Eta Car wind-history example at $t=5.90\times10^{10}$ s (1869.597 AD), 11.597 yr after the 1858 transition out of the eruptive wind. No 2048-cell convergence claim is made.

The density passed to ASGARD combines:

- `pion_eta_car_sph1d_n1024_1870_analytic_inner.csv`: the unrenormalized free wind from $R_0=10^{13}$ cm to the PION injection radius $1.32\times10^{16}$ cm;
- `pion_eta_car_sph1d_n1024_1870_external_raw.csv`: 984 unchanged external PION cells with density, velocity, temperature, and tracer;
- `pion_eta_car_sph1d_n1024_1870_combined_96knots.csv`: two analytic anchors plus 94 unchanged PION cells used by `TabulatedMedium`;
- `pion_eta_car_sph1d_n1024_1870_combined_96indices.txt` and `pion_eta_car_sph1d_n1024_1870_metrics.json`: selection and physical diagnostics.

The first external-cell density differs from the analytic free wind at the same radius by -0.01496%; no normalization or smoothing is applied. See `PION_PROVENANCE_1870.txt` for code commits, wind parameters, run commands, and the deterministic knot-selection rule. `SHA256SUMS` records the tracked bytes.
