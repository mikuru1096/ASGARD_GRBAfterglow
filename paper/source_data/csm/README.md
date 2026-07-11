# PION Eta Car circumstellar-medium source data

This directory preserves the 1D spherical conversion of the public PION Eta Car wind-history example and the exact tabulated profile passed to ASGARD Figure 2. The raw table contains 1968 external-domain cell centers from the 2048-cell simulation. The 96-knot table is a deterministic subset of unchanged raw cells selected by maximum log-density interpolation residual; it is not smoothed or fitted.

ASGARD uses `proton_equivalent_n_cm3=rho/m_p` for the upstream mass density. The `nH_X033_cm3` column is retained only as a composition diagnostic. See `PION_PROVENANCE.txt` for commits, run commands, crop, and the upstream angle-patch dependency boundary. The 1024--2048 comparison establishes convergence of the main structure positions and shell mass, but not full cellwise convergence.
