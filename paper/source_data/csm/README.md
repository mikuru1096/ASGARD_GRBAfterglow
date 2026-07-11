# PION Eta Car CSM at 1869.597 AD

This directory stores the 1024-cell spherical one-dimensional reduction of the public PION Eta Car wind-history example at $t=5.90\times10^{10}$ s. No 1024/2048 convergence claim is made.

The tracked density source contains:

- `pion_eta_car_sph1d_n1024_1870_analytic_inner.csv`: unrenormalized post-eruption free wind from $10^{13}$ to $1.32\times10^{16}$ cm;
- `pion_eta_car_sph1d_n1024_1870_external_raw.csv`: 984 unchanged PION cells;
- `pion_eta_car_sph1d_n1024_1870_analytic_outer.csv`: pre-eruption steady-wind anchor at $3\times10^{19}$ cm;
- `pion_eta_car_sph1d_n1024_1870_extended_96knots.csv`: two inner analytic anchors, 93 unchanged raw cells, and one outer analytic anchor;
- `pion_eta_car_sph1d_n1024_1870_extended_96indices.txt` and `pion_eta_car_sph1d_n1024_1870_extended_metrics.json`: selection and boundary diagnostics.

The inner and outer raw-to-analytic density mismatches are -0.01496% and -0.165693%, respectively. No normalization or smoothing is applied. See `PION_PROVENANCE_1870.txt` and `SHA256SUMS`.
