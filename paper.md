# ASGARD paper restructuring execution plan

## 1. Goal and fixed claim

The paper must become an ApJ methods paper, not a kernel manual.  The main line
is:

1. observational degeneracy;
2. radius-ordered local state;
3. multi-physics modules;
4. benchmark evidence;
5. discussion;
6. hard boundaries;
7. summary.

The fixed contribution sentence is:

ASGARD records the radius-ordered local state of GRB afterglow dynamics,
particles, photons, optional reverse shocks, optional hadronic/pair feedback,
absorption, and observer components, so that complex multi-physics light-curve
features can be traced back to their local physical causes before observer
projection mixes them.

This claim is bounded.  ASGARD exposes and validates the implemented state
sequence.  It does not claim support for physical couplings that lack a local
state, projection operator, and validation evidence.

## Current progress ledger

- Phase 0: completed.  The planning document has been expanded into a full
  execution plan.
- Phase 1: completed.  The manuscript skeleton, benchmark placeholders,
  Discussion entry, and technical appendix entry points have been added.  The
  encoding check, focused diff check, Windows TeX Live compile, and LaTeX log
  scan pass; the only reported warning is the accepted AASTeX small-caps font
  substitution.
- Phase 2: completed.  Main-text compression and appendix migration have reached
  the end of `paper/main.tex` in increasing manuscript order.
  Completed windows:
  - `paper/main.tex` Dynamics opening through Jet angular profiles.
  - `paper/main.tex` Forward-shock dynamical variables.
  - `paper/main.tex` Dynamical branches.
  - `paper/main.tex` Forward-dynamics integration, with runtime container and
    output-grid details moved to Appendix~\ref{app:projection-runtime}.
  - `paper/main.tex` Electron continuum state and finite-shell source, including
    the natural-log \(x=\ln\gamma\) grid and no base-10 Jacobian.
  - `paper/main.tex` Initial electron spectrum and coordinate projection, with
    quadrature and four-velocity coordinate-map details moved to
    Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Shared conservative shell-transport primitives, with
    diagnostic reconstruction, triangular sweep, remap interpolation, and
    positive-source reconstruction moved to Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Characteristic remap primitive, with quadrature weights,
    semi-Lagrangian split, flux-split matrix, legacy triangular path, and
    centralized support diagnostics moved to Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Implicit finite-volume `fullhide_1d` solver, with fixed
    substep grid, uniform merged solve, and non-uniform or thermal substep
    details moved to Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Coupled `fullhide` electron-photon pass, with repeated
    triangular-solve coefficients moved to Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Hybrid thermal `fullhide` branch, with special-function,
    thermal-tail, root-seeding, and pointwise density details moved to
    Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Semi-Lagrangian log-gamma solver, with source-split remap
    details kept in Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Three-level log-gamma solver, with startup and BDF2
    triangular coefficients moved to Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` WENO finite-volume solver, with WENO-Z candidate flux,
    smoothness-indicator, boundary-stencil, and SSP-RK stage details moved to
    Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Direct characteristic-integration solver, with piecewise
    affine cooling-map and conservative source-integral remap details moved to
    Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Discontinuous-Galerkin state, active support, radial
    substeps, finite-shell source, and output projection, with active-edge,
    density-jump, gradient-limiter, and source-renormalization details moved to
    Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Discontinuous-Galerkin source projection, implicit
    transport, conservative remap, and characteristic projector, with source
    moment projection, nodal block entries, remap/positivity formulas, and
    affine characteristic maps moved to Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Cooling assembler, Nakar analytic Compton closure,
    Fan Klein--Nishina Compton closure, and cooling-branch boundary, with the
    fast/slow Fan \(\eta_i^{\rm KN}\) integrals moved to
    Appendix~\ref{app:photon-kernels}.
  - `paper/main.tex` Synchrotron-self-absorption heating kernel, with target
    photon state, single-particle absorption fit, and heating integral retained
    in main text, and prefix-sum/high-frequency quadrature details moved to
    Appendix~\ref{app:photon-kernels}.
  - `paper/main.tex` Grid-level inverse-Compton cooling kernel, with IC energy
    variables, Thomson/Jones kernel selection, and \(A_i^{\rm IC}\) retained in
    main text, and seed-grid centering plus coupled-budget Simpson accumulation
    moved to Appendix~\ref{app:photon-kernels}.
  - `paper/main.tex` Synchrotron emissivity quadrature, with the physical
    emissivity integral and fixed-grid log-\(\gamma\) quadrature retained in
    main text, and the single-particle \(\Phi(x)\) fit moved to
    Appendix~\ref{app:photon-kernels}.
  - `paper/main.tex` Synchrotron self-absorption shell transfer, with the shell
    optical depth, finite-cell transfer factor, escaping synchrotron luminosity,
    and IC seed photon density retained in main text.
  - `paper/main.tex` Cyclotron-bin synchrotron extension, with the low-frequency
    shell luminosity deposit kept in main text and the electron spectrum boundary
    stated explicitly.
  - `paper/main.tex` Radiation interpolation and log-Gauss quadrature helpers,
    with main text reduced to their radiation-field role and the interpolation
    plus quadrature formulas moved to Appendix~\ref{app:photon-kernels}.
  - `paper/main.tex` SSA heating segment quadrature, with the absorbed-photon
    energy-return integral retained in main text and the segment cross-section
    branch moved to Appendix~\ref{app:photon-kernels}.
  - `paper/main.tex` Self-absorption frequency solver, with the
    \(\tau_\nu(\nu_a)=1\) definition, diagnostic output, and no-crossing boundary
    kept in main text, and adaptive quadrature plus root-update details moved to
    Appendix~\ref{app:photon-kernels}.
  - `paper/main.tex` Synchrotron polarization kernel, with the local diagnostic
    role and unsupported full observer-polarization boundary stated in main text,
    and the intrinsic kernel integrals, asymptotic branches, and thermal Bessel
    table moved to Appendix~\ref{app:photon-kernels}.
  - `paper/main.tex` Adaptive radiation-grid resampling, with the main text
    reduced to its quadrature-support role and electron-spectrum boundary, and
    the active-support scan, curvature estimator, and point-selection formulas
    moved to Appendix~\ref{app:photon-kernels}.
  - `paper/main.tex` \(\chi\)-resolved downstream electron state, with the
    downstream \(q\)-coordinate, cell average, local magnetic field, and cell
    radiation outputs retained in main text, and initialization, finite-speed
    geometry, pattern-speed age, and Bohm-diffusion details moved to
    Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` \(\chi\)-resolved photon-history cooling, with the causal
    light-cone field, Doppler frequency transform, attenuated history
    contribution, and cooling-field output retained in main text, and overlap,
    positive interpolation, path-product, and streamed-update details moved to
    Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Downstream \(q\)-space advection--diffusion, with the
    conservation law, shock-front source, representative face flux, branch
    output, and outer-boundary statement retained in main text, and
    characteristic face tracing, PPM remap, diffusion metric, and matrix-boundary
    details moved to Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Downstream energy-space update, with the four-velocity
    cooling displacement, natural-log \(\ln\gamma_e\) branch, characteristic
    cooling law, local adiabatic term, stochastic-transport boundary, and module
    output retained in main text, and backward-sweep coefficients plus stochastic
    diffusion matrix details moved to Appendix~\ref{app:electron-solvers}.
  - `paper/main.tex` Returned \(\chi\)-resolved state and electron-transport
    figure/boundary, with the returned projection arrays, shell-vs-\(\chi\)
    projection choice, unsupported \(q\)-local channel boundary, figure
    diagnostic role, and thermal/adaptive substep limits stated in main text.
  - `paper/main.tex` Photon Fields and Absorption through Comoving target photon
    state, with the photon-state role, local comoving-grid boundary, separated
    target-field equation, SSC target switch, joint-branch survival factors, and
    observer-frame boundary stated in main text.
  - `paper/main.tex` Photon ABI binding and radiation interpolation reuse, with
    the main text reduced to the radiation-kernel boundary, public photon
    products, physical-kernel ownership, and positive-endpoint interpolation
    rule, and compiled binding names moved to
    Appendix~\ref{app:photon-kernels}.
  - `paper/main.tex` Composite Simpson radiation quadrature, with the helper
    reduced to the uniform-grid SSC quadrature object, logarithmic integration
    variables, named Simpson weights, branch output, and unsupported grid-layout
    boundary.
  - `paper/main.tex` Uniform-grid SSC emissivity, with the local SSC physical
    object, Jones--Klein--Nishina variables, kinematic electron limit,
    low-seed and high-seed luminosity integrals, stored SSC photon density, and
    no-electron-update boundary stated in main text.
  - `paper/main.tex` Nonuniform-grid SSC emissivity, with the same SSC kernel
    tied to unequal \(\ln\gamma\) cells, linear reconstructed \(N_x\), low-seed
    cell integrals, high-seed electron moments, returned SSC products, and
    strictly increasing electron-edge boundary.
  - `paper/main.tex` Finite-cell photon transfer, with the cell-average
    escaping-luminosity factor, local SSA optical depth, local \(\gamma\gamma\)
    optical depth, photon-output role, and no-change boundary for photon
    production and electron transport stated in main text.
  - `paper/main.tex` \(\gamma\gamma\) shell absorption, with the target photon
    cell measure, angle-dependent pair-production invariant, physical threshold
    and finite kernel-domain boundary, shell optical-depth assembly, finite-cell
    survival factor, and head-on path-opacity helper boundary stated in main
    text.
  - `paper/main.tex` External EBL attenuation, with the cosmological absorption
    role, table frequency/redshift interpolation, zero-distance and table-range
    boundaries, observer-side attenuation factor, and no-local-state-change
    boundary stated in main text.
  - `paper/main.tex` Observer-side photon component assembly, with active
    component accounting, exported FS/RS/hadronic channel boundaries, the local
    distance-and-absorption prefactor, EBL separation, and the projection handoff
    stated in main text.
  - `paper/main.tex` Reverse-shock primary ejecta state and activation, with the
    finite-shell activation boundary, baryonic and magnetized ejecta density
    split, shock-front \(\gamma_{34}\), magnetized pressure/fast-mode criteria,
    waiting-state inertia, and no-heat/no-injection boundary stated in main text.
  - `paper/main.tex` Finite-strength MHD jump solver, with the main text reduced
    to the cold-upstream variables, finite positive downstream four-speed root,
    compression output, low-\(\sigma\) hydrodynamic branch, primary-RS-only
    boundary, and the full cubic coefficients moved to
    Appendix~\ref{app:reverse-shock-technical}.
  - `paper/main.tex` Comoving region-3 field and heat state, with the upstream
    ordered field, turbulent-plus-ordered total \(B'_3\), post-crossing
    ordered-field advection, magnetic enthalpy inertia, magnetized heat
    \(\zeta_3\), thermal-response coupling, and waiting/post-crossing boundaries
    stated in main text.
  - `paper/main.tex` Primary reverse-shock ODE, with the normalized contact-state
    vector, pre-crossing mass flux, reverse-shock front speed from the jump state,
    active pre-crossing Lorentz equation, reverse radiative fraction using total
    \(B'_3\), thermal/volume update, and post-crossing no-new-mass boundary stated
    in main text.
  - `paper/main.tex` Crossing and activation event splitting, with the waiting,
    active pre-crossing, and post-crossing branch changes stated as physical
    events, \(q=M_3/M_{\rm ej}\) and \(s=\ln(t/t_0)\) used as readable step
    coordinates, the event-state register reduced to a diagnostic record, the
    frozen crossing vector \(\mathbf C_\times\) stated, and the single-branch
    refinement boundary kept in main text.
  - `paper/main.tex` Reverse electron injection contract and reverse
    radius-output sweep, with newly shocked ejecta mass used as the radial source,
    the transported injection break tied to shock-front \(\gamma_{34}\), the
    natural-log \(x=\ln\gamma_e\) source normalization and Jacobian stated,
    fullhide/DG source sharing stated, and the no-new-injection post-crossing
    boundary stated in main text.
  - `paper/main.tex` Reverse cooling assembly and finite-volume transport, with
    the positive radial cooling coefficient \(D_i\) defined, synchrotron-only,
    local IC/KN, and Nakar-style Compton modes separated, \(B'_3\) kept as the
    total reverse-shock field, \(Y_i^{\rm N}\) tied to the reverse seed field,
    the \(y\)-coordinate face speed written with \(dx/dy\), and the fixed
    fullhide/DG substep controls stated as transport controls rather than fit
    parameters.
  - `paper/main.tex` Post-crossing characteristic remap, with the no-new-source
    boundary stated first, \(N_x=dN_e/dx\) on the natural-log \(x=\ln\gamma_e\)
    grid used as the conserved spectrum, synchrotron cooling written as an
    analytic inverse-energy map, IC/KN cooling written as a piecewise affine
    inverse-energy map, and the final conservative remap and \(dN_e/d\gamma_e\)
    diagnostic stated.
  - `paper/main.tex` Secondary-reservoir electron transport, with each
    density-jump branch stated as a diagnostic reservoir carrying its own mass,
    volume, field, injection scale, and spectrum, the radial source and volume
    coefficient written in main text, inactive-mass and nonpositive-volume
    boundaries stated, branch synchrotron/SSA output stated, and unsupported
    secondary-reservoir SSC/hadronic/pair-feedback products stated directly.
  - `paper/main.tex` Reacceleration diagnostic, with the execution order stated as
    parent-mass transfer, compression/relative-motion boost, DSA remap, and fresh
    kinetic source insertion, the dissipative-branch and positive-parent-mass
    boundaries stated, \(p_r>2\) kept as the hard DSA boundary, and the seed and
    reaccelerated electron-energy diagnostics stated as diagnostics rather than
    radiative components.
  - `paper/main.tex` Density-jump secondary-shock solver, with the diagnostic
    branch role stated, Gaussian rising-side event density, secondary inertia
    contribution, first-branch and chained-parent region-4 upstream states,
    contact-pressure root, positive excess heat source \(S_j\), compressive
    shock-frame closure, dissipated mass/heat/volume source terms, the
    \(\gamma_{m,j}\) newly heated electron scale, and hot-gas/cutoff boundaries
    stated in main text.
  - `paper/main.tex` Secondary event windows and observer diagnostics, with the
    \(S_j>0\) activation rule stated before the scanner, start/end event records
    tied to adjacent primary reverse-shock states, dimensional branch outputs
    restored from ODE variables, dissipation-weighted shock diagnostics separated
    from old inactive mass, combined secondary state and break frequencies stated
    as diagnostics, and the unsupported SSC/hadronic/cross-zone feedback boundary
    kept explicit.
  - `paper/main.tex` Hadronic and Pair-Feedback Boundary opening, with the
    one-dimensional shell-averaged activation boundary stated, supported
    legacy/formal hadronic operators separated from inactive switches, the ABI
    package binding limited to kernel selection, the hadronic base grid corrected
    to the implemented \(\log_{10}\gamma\) coordinate, its separation from the
    natural-log electron continuity grid stated, geometric cell faces and shell
    time scales retained, and the legacy proton upper-bound estimate tied to
    acceleration, dynamical, and synchrotron limits.
  - `paper/main.tex` Quantum-synchrotron helper, with the helper stated as a
    synchrotron-cooling multiplier, the QED invariant \(\chi\) defined from the
    local comoving field and particle mass, the AM3-style Landau factor kept in
    main text, and the \(B'\le0\), \(\gamma\le1\), and disabled-option boundaries
    stated as \(f_q=1\) cases.
  - `paper/main.tex` Hadronic species state, with the returned proton, neutron,
    pion, muon, neutrino, and secondary-pair source arrays stated as shell-local
    objects, the compact proton operator split into source, nucleon reinjection,
    interaction losses, and continuous adiabatic/synchrotron losses, species
    masses and charges reduced to \((m_s,Z_s)\), the incorrect swept-volume
    energy formula removed in favor of the supplied shell energy \(E_{s,i}\),
    formal and legacy source normalizations separated by relativistic versus
    kinetic energy moments, and charged-species acceleration, synchrotron, and
    external-cooling limits stated with neutral-species exclusion.
  - `paper/main.tex` Bethe-Heitler pair-production kernel, with the
    \(p\gamma\to pe^+e^-\) physical channel named first, the kinematic threshold
    tied to zero pair source/photon sink/proton loss below threshold, the
    Blumenthal--Gould differential pair kernel kept in compact integral form,
    closed phase-space and finite cross-section-domain boundaries stated, the
    shell-local secondary-pair source, target-photon sink, and positive proton
    cooling rate separated, and the Chodorowski loss fit kept as named coefficient
    sets rather than raw decimal coefficients.
  - `paper/main.tex` Proton transport operator, with the formal feedback path
    corrected to conservative cooled-content remap followed by exponential
    \(p\gamma\) sink and same-step reinjection, the legacy upwind branch separated
    from the coupled feedback path, baseline adiabatic/synchrotron/BH/pp losses
    kept as positive cooling rates, and the two-pass \(p\gamma\) photon-survival
    operator tied to shell-local one-dimensional hadronic outputs before observer
    projection.
  - `paper/main.tex` Shell-sequence pair-cascade algorithm, with the supported
    time-dependent cascade branch stated first, the pair grid tied to the photon
    log-energy spacing and \(m_ec^2\) lower bound, the shell escape time written
    as the photon residence scale, the \(\gamma\gamma\) operator outputs separated
    into photon loss, pair injection, and local energy-closure diagnostics, the
    finite-grid pair cooling written as a cooled-content remap plus current
    source, and the photon-density update closed with pair synchrotron seed
    photons before observer projection.
  - `paper/main.tex` Single-shell pair-synchrotron diagnostic, with the branch
    stated as a local kernel check rather than the formal cascade path, the
    pair-operator outputs limited to photon loss and absorbed-power diagnostics,
    the synchrotron cooling track written with \(a_B\), \(\gamma_i^{\rm c}\), and
    \(t_{\rm cool}\), the grid deposition tied to the local diagnostic pair
    source, and the diagnostic emissivity written with \(C_\nu\), \(C_P\), and
    \(\Phi_{\rm d}\) instead of cgs numerical constants.
  - `paper/main.tex` Observer Projection and Component Assembly opening, with
    observer projection stated as the final read-only mapping from local
    radius-ordered spectra to \(F_\nu(t_{\rm obs})\), the inverse Doppler factor
    \(\eta_\delta=\delta^{-1}\) separated from the stored kernel variable, the
    comoving-frequency shift and angular weight written in the same convention,
    the positive-endpoint log-frequency interpolation rule stated, and the
    structured-jet active-patch register corrected so it no longer claims
    arbitrary exact-state-class reuse.
  - `paper/main.tex` Shell prefactor and component assembly, with local component
    luminosities separated from observer-projected fluxes, the prefactor corrected
    to the combined photon--photon absorption factor times
    \((1+z)/(4\pi d_L^2)\), the active component labels expanded to the actual
    FS/RS/SSC/hadronic/pair/cross-zone component set, inactive components stated
    as omitted, and the structured diagnostic radius track identified as a local
    one-dimensional patch history rather than an angular flux average.
  - `paper/main.tex` Projection symbol binding, with the
    \texttt{src.Interpolation} initializer stated as a Python-to-compiled-kernel
    ABI binding layer, the exported projection names tied to the compiled
    \texttt{SED\_interpolation} module, non-exported names stated as absent, and
    the physical radiation process, optical depth, Doppler factor, arrival-time
    map, and interpolation weights left inside the selected Fortran kernel.
  - `paper/main.tex` Direct top-hat equal-arrival-time interpolation, with
    midpoint angular cells tied to observer-facing surface elements, the
    collapsed azimuthal half-plane handled by the named \(N_{\phi,c}\) control,
    the radial arrival-time bracket written through \(\lambda_R\), \(R_R\), and
    \(\Gamma_R\), the physical Doppler factor separated from the stored inverse
    Doppler kernel variable, the shared positive-endpoint log-frequency
    interpolation rule stated, and the adaptive polar quadrature limited to
    sampling density rather than changing the arrival map or shell spectrum.
  - `paper/main.tex` \(\chi\)-resolved finite-cell projection, with downstream
    cells treated as finite emitting zones, the arrival-time map written with both
    shock-front and cell radii, the radial search corrected to a monotonic
    arrival-time bracket plus segment scan for nonmonotonic cells, the finite-cell
    escape factor written from \(\tau_{\rm front}\), \(\tau_\chi\), and
    \(\psi(x)\), the direct-electron diagnostic tied to the same projection
    kernel, and the ray-traced SSA path described as a geometric opacity operator
    over stored luminosity and SSA grids rather than a new radiation solve.
  - `paper/main.tex` Legacy structured-jet projection, with structured angular
    quadrature written as observer weighting of thin-shell jet patches, the
    axisymmetric and non-axisymmetric entries separated by their azimuthal
    measures and mirror factor, the internal structured interpolation source
    stated as a non-exported companion used by \texttt{structured\_jet\_1d}, the
    patch arrival-time map kept identical in physics to the top-hat map but with
    patch-dependent radial histories, and the output identified as the component
    flux summed over the structured angular grid.
  - `paper/main.tex` Shared frequency-grid accumulator and public projection
    boundaries, with \texttt{accum\_logsed} described as the final remap of an
    already shifted source spectrum, \texttt{accum\_shifted} separated as the
    scalar-shift version that evaluates \(S(X'_j)\), the positive-endpoint
    log-log interpolation rule and linear nonpositive-endpoint rule stated
    without adding any logarithmic floor, the accepted geometry kernels listed,
    and the \(\chi\)-resolved projection boundary limited to forward-shock
    synchrotron plus SSA while SSC, hadronic, and pair-cascade components remain
    shell-level.
  - `paper/main.tex` Ray-traced synchrotron self-absorption boundary, with the
    ray path stated as a foreground-opacity operation over already computed
    \(\chi\)-cell luminosities and SSA depths, the comoving optical-depth
    relation written from \(\tau_r=\alpha'_{\nu'}\Delta r'\), foreground hits
    limited to cells that cover the emitting sky position and lie in front of the
    source event, the final flux written with \(\tau_{{\rm front},e}\),
    \(\psi(\tau_e\ell_e)\), \(\widetilde L_{\nu'_e,e}\), and \(W_{\chi,e}\), and
    convergence limited to global resolution changes rather than local angular
    additivity tests.
  - `paper/main.tex` Fitting, Runtime, and Reproducibility opening, with the
    fitting wrapper separated from the physical kernels, the active-variable map
    limited to \({\cal P}_{\log}\), \({\cal P}_{\rm lin}\), and
    \({\cal P}_{\rm fix}\), observation-block residuals tied to one solved
    total-flux matrix, band fluxes written as trapezoid frequency integrals,
    exposure products written as normalized Gauss--Legendre time averages, the
    repeated-pair cache corrected to a runtime query cache, and all caches stated
    as cost controls that do not change the physical state.
  - `paper/main.tex` Verification and Benchmark Evidence plus Discussion, with
    source-data figures, formal cross-code benchmark subsets, ASGARD-only complex
    state diagnostics, unsupported branches excluded from error metrics, and the
    Discussion organized around degeneracy, radius-ordered state, and explicit
    method boundaries.
  - `paper/main.tex` Model Boundaries, with runtime limits rewritten as explicit
    supported and unsupported input domains, coupled-physics limits tied to the
    need for both local state and observer projection, the formal hadronic path
    stated as shell-averaged rather than \(\chi\)-resolved, and IC-mediated
    electromagnetic cascades stated as unsupported.
  - `paper/main.tex` Summary, rewritten as bounded contribution, source-data
    backed evidence, and constrained use cases, with numerical claims limited to
    rows present in the tracked source tables.
  - `paper/main.tex` Appendix Electron Solvers opening, with the appendix role
    clarified, all electron energy grids tied to the natural-log reference
    coordinate \(x_g=\ln\gamma_e\), initial spectra converted as conservative
    cell averages, four-velocity storage described as a coordinate map, and the
    diagnostic \(dN_e/d\gamma_e\) reconstruction separated from the transported
    cell averages.
  - `paper/main.tex` Appendix Electron Solvers remap details, with the local
    remap coordinate defined as \(x=x_g=\ln\gamma_e\), point and parabolic
    reconstructions described as integrals over traced faces, source
    reconstruction separated from transported spectra, and the characteristic
    update stated as old-electron transport plus source injection along the
    radial step.
  - `paper/main.tex` Appendix Electron Solvers semi-Lagrangian, three-level
    natural-log, and WENO branch details, with the semi-Lagrangian update written
    as traced-face integration plus split source injection, the three-level
    source label changed to \(Q_{x,j}^{\rm 3}\) to avoid a new uppercase \(T\)
    symbol, and WENO stated as an explicit natural-log spectrum update.
  - `paper/main.tex` Appendix Electron Solvers DG active mesh and implicit block,
    with the active high-energy domain tied to the stored high-energy moment
    rather than a physical cutoff, density-jump substeps stated as radial work
    controls, source normalization and Legendre projection separated, DG remaps
    described as same-\(\gamma_e\) coordinate transport, and the natural-log
    `fullhide` triangular solve tied back to the shared cooling sweep.
  - `paper/main.tex` Appendix Electron Solvers diagnostic ranges and
    `fullhide_1d` fixed-substep details, with diagnostic support ranges stated as
    radiation-kernel scan limits rather than physical cutoffs, fixed substeps
    written as shell source and cooling-displacement accumulation, non-uniform
    density handled by explicit radial substeps, and the density response stated
    as the supported rule rather than an open-ended approximation.
  - `paper/main.tex` Appendix Electron Solvers hybrid thermal details, with the
    thermal/nonthermal hybrid source framed as source normalization rather than a
    new transport equation, special functions tied to cutoff power-law and
    Maxwell--Juttner moments, \(\Theta\) solved from an explicit root
    \(F(\Theta)\), and the output stated as the hybrid source profile
    \(\phi_{\rm h}\) on the electron grid.
  - `paper/main.tex` Appendix Electron Solvers \(\chi\)-resolved downstream
    state details, with the \(q\)-coordinate initialization, BM/planar geometry
    bridge, downstream speed, magnetic age, local \(B'_a\), and Bohm diffusion
    written as local state construction rather than a new dynamical model.
  - `paper/main.tex` Appendix Electron Solvers \(\chi\)-resolved photon-history
    details, with causal source/target cell support, relative Doppler mapping,
    positive-field interpolation, SSA plus \(\gamma\gamma\) path attenuation, and
    streamed light-cone recurrence stated as the source of the local seed-photon
    history field.
  - `paper/main.tex` Appendix Electron Solvers downstream \(q\)-space
    advection--diffusion details, with \(A_q=(3-k)(q_\ast-q)/R\), first-cell
    source placement, characteristic face traceback, conservative PPM remap,
    Bohm-to-\(q\) diffusion metric, implicit tridiagonal solve, and \(N_\chi\)
    outer-boundary notation stated explicitly.
  - `paper/main.tex` Appendix Electron Solvers downstream energy-space update
    details, with natural-log \(x=\ln\gamma_e\) cooling, default and PWN/CR
    adiabatic coefficients, backward-sweep coefficients, nonnegative solved-cell
    projection, and stochastic zero-flux energy diffusion stated as separate
    operations.
  - `paper/main.tex` Appendix Photon, SSC, and SSA Kernels opening details, with
    the appendix framed around local synchrotron luminosity, seed photon density,
    SSC luminosity, self-absorption heating, photon survival, and ABI entry points
    stated as runtime bindings rather than new photon states.
  - `paper/main.tex` Appendix Photon, SSC, and SSA Kernels Fan-style
    Klein--Nishina closure details, with \(\hat\gamma_i\), \(\gamma_M\), fast/slow
    cooling support, \(\eta_i^{\rm KN}\), and the analytic \(Y_i^{\rm F}\) use
    stated before the piecewise formulas.
  - `paper/main.tex` Appendix Photon, SSC, and SSA Kernels SSA heating
    quadrature details, with \(A_i^{\rm abs}\) kept as the physical object,
    low-frequency complete cells accumulated by a cached logarithmic-Gauss prefix
    sum, high-frequency complete cells accumulated from the exponential branch,
    clipped boundary cells stated separately, and interpolation rules tied to the
    local seed photon density \(n_\nu\).
  - `paper/main.tex` Appendix Photon, SSC, and SSA Kernels grid-level
    inverse-Compton cooling path, with the seed/scattered photon cell-centered
    measure, Jones/KN branch kernel, discrete \(A_i^{\rm IC}\) loss accumulation,
    Simpson-weighted joint-budget coefficient \(A_i^{\rm IC,b}\), and
    no-separate-cooling-law boundary stated explicitly.
  - `paper/main.tex` Appendix Photon, SSC, and SSA Kernels synchrotron emissivity
    analytic-fit details, with the single-particle \(\Phi(x)\) approximation,
    positive-field interpolation rule, logarithmic Gauss cell measure, and reduced
    cooling-table support described as radiation-evaluation helpers rather than
    new physical electron or photon states.
  - `paper/main.tex` Appendix Photon, SSC, and SSA Kernels SSA segment and
    self-absorption-frequency quadrature details, with the segment absorption
    cross section, absorbed-energy contribution \(A_{\rm seg}^{\rm abs}\),
    \(\tau_\nu\) quadrature on \(x=\ln\gamma_e\), and logarithmic
    \(\tau_\nu(\nu_a)=1\) root solve stated as diagnostics for \(A_i^{\rm abs}\)
    and \(\nu_a\).
  - `paper/main.tex` Appendix Photon, SSC, and SSA Kernels intrinsic polarization
    and auxiliary Bessel-table details, with local \(F/G\) synchrotron
    polarization kernels, split local luminosities, the no-observer-polarization
    boundary, analytic/quadrature special-function branches, and the \(K_2(z)\)
    table stated as a Maxwell--Juttner normalization helper.
  - `paper/main.tex` Appendix Photon, SSC, and SSA Kernels adaptive
    radiation-grid helper details, with positive-electron support selection,
    natural-log \(\gamma_e\) and \(N_e\) coordinates, finite-difference
    \(dy/dx\) and \(d^2y/dx^2\) curvature weights, cumulative-weight sampling,
    endpoint/neighbor preservation by the caller, and the output stated as
    quadrature support rather than a changed electron spectrum.
  - `paper/main.tex` Appendix Hadronic and Pair Feedback opening figure and
    transport-boundary details, with the optional one-dimensional shell-averaged
    high-energy state, local inputs and outputs, the unsupported
    \(\chi\)-resolved hadronic boundary, the \(N_{\rm cas}>1\)
    \(\gamma\gamma\) pair/synchrotron cascade meaning, and the unsupported
    IC-mediated electromagnetic cascade boundary stated before the figure
    caption.
  - `paper/main.tex` Appendix Reverse-Shock Technical Details finite-strength MHD
    jump cubic, with \(g=\gamma_{34}\), \(y=u_d^2\), the low-\(\sigma\)
    hydrodynamic-strength root, the finite-\(\sigma\) cubic coefficients,
    depressed-cubic root, upstream shock-frame four-speed, compression ratio, and
    \(\sigma_{\rm cut}\) stated as an algorithmic branch parameter rather than a
    physical magnetization boundary.
  - `paper/main.tex` Appendix Projection, Fitting, and Runtime Details opening
    and forward-dynamics runtime record, with the appendix scope made definite,
    the forward boundary vector stated as an ABI-ordered container, optional
    energy-injection and density-jump slots explained, integration-start
    parameters separated from physical inputs, the source-frame on-axis output
    grid written in clear steps, and public output arrays tied to on-axis time,
    Lorentz factor, radius, and swept mass.
  - `paper/main.tex` Appendix Software and Data Availability, with the cited
    software record, manuscript code baseline, tracked source tables,
    formal figure-generation entry point, figure products, exploratory-output
    exclusion, and software dependencies stated as the evidence path for the
    method-paper claims.
- Phase 3: completed as a bounded cross-code source-data layer.  ASGARD,
  afterglowpy, jetsimpy, VegasAfterglow, and PyBlastAfterglowMag have numerical
  benchmark rows where public assumptions are matched.  PyBlastAfterglowMag is
  installed, source-pinned, cloned under `_external/PyBlastAfterglowMag`, and
  compiled into a repo-local `src/pba.out`; its rows are quoted for the matched
  FS and Gaussian-geometry subsets.
  Completed setup:
  - pinned the comparison extra in `pyproject.toml` and `uv.lock`:
    afterglowpy 0.8.1, VegasAfterglow 2.0.5, jetsimpy commit
    `53d07610830f2247d8c41364a2b469491fa22eb2`, and
    PyBlastAfterglowMag commit
    `b0a39d170930908a398de56089fde8cfab16883d`;
  - added the PyBlastAfterglowMag import dependencies `pandas==2.3.3` and
    `uncertainties==3.2.3` to the same compare extra because the package metadata
    does not declare them;
  - verified minimal imports for afterglowpy, jetsimpy, VegasAfterglow, and
    PyBlastAfterglowMag in the project `uv` environment;
  - added `tests/generate_cross_code_benchmarks.py` as the formal Phase 3
    source-data setup entry;
  - generated `paper/source_data/figB1_cross_code_fs.csv`,
    `paper/source_data/figB2_rs_ssc_geometry.csv`,
    `paper/source_data/figB3_asgard_complex_state.csv`, and
    `paper/source_data/table_cross_code_capabilities.csv` with fixed versions,
    explicit unsupported boundaries, and empty numerical fields where no formal
    benchmark run has been made.
  - filled `paper/source_data/figB1_cross_code_fs.csv` with measured
    forward-shock synchrotron flux densities for ASGARD, afterglowpy, jetsimpy,
    and VegasAfterglow on the shared high-frequency top-hat, on-axis,
    uniform-medium scenario. The table records cgs flux density, ASGARD/code
    ratios, wall time, package versions, and unit-conversion notes for packages
    that return mJy.
  - widened `paper/source_data/figB1_cross_code_fs.csv` to a formal
    broad-window grid. The rows now span \(1\)--\(10^8\,{\rm s}\) and
    \(10^8\)--\(2.4\times10^{26}\,{\rm Hz}\), including low radio, radio,
    mm/IR, optical, X-ray, GeV, and VHE observer bands. The ASGARD run uses
    \(N_R=192\), \(N_{t_{\rm obs}}=48\), \(N_{\gamma_e}=128\), and
    \(N_\nu=96\). The jetsimpy run uses its supported highest public
    calibration level, stricter relative tolerance, and an internal evolution
    endpoint beyond the plotted observer-time window. These rows are still
    assumption-audit diagnostics, not a statement that every code implements the
    same absorption, cutoff, or high-energy physics.
  - filled `paper/source_data/figB2_rs_ssc_geometry.csv` with an ASGARD/Vegas
    matched top-hat RS+SSC component subset. The table records FS synchrotron,
    FS SSC, RS synchrotron, and RS SSC flux-density components at two
    frequencies and three observer times, while afterglowpy, jetsimpy, and
    PyBlastAfterglowMag remain explicit unsupported entries for the RS/SSC subset.
    These rows demonstrate component visibility and shared-interface execution;
    they are not yet a broad RS/SSC agreement claim because weak SSC components
    have large ratios.
  - filled `paper/source_data/figB3_asgard_complex_state.csv` with ASGARD-only
    radius-ordered dashboard rows for the full state chain, the proton-synch
    hadronic state, the activated pair-feedback state, and the \(\chi\)-resolved
    projection state. The numeric rows
    record local \(\Gamma\), comoving magnetic field, electron number and energy,
    proton number and energy, hadronic synchrotron peak values, chi-volume
    weights, chi optical-depth maxima, chi seed-photon maxima, BH-pair number and
    energy, p-gamma and Bethe--Heitler peak values, p-gamma and BH optical-depth
    maxima, and minimum p-gamma photon survival. The activated pair-feedback rows
    use the one-dimensional shell-averaged `am3_1d` path with joint coupling;
    \(\chi\)-resolved hadronic feedback remains unsupported.
  - widened `paper/source_data/figB2_rs_ssc_geometry.csv`. The RS/SSC subset now
    spans \(1\)--\(10^7\,{\rm s}\) and \(10^8\)--\(2.4\times10^{26}\,{\rm Hz}\).
    The Gaussian structured/off-axis FS subset now spans
    \(10\)--\(10^8\,{\rm s}\) and \(10^8\)--\(2.4\times10^{23}\,{\rm Hz}\) for
    ASGARD, afterglowpy, jetsimpy, and VegasAfterglow. The rows use
    \(\theta_c=0.08\) rad, \(\theta_{\rm obs}=0.12\) rad, \(\Gamma_0=120\),
    uniform medium, and no spreading. ASGARD uses a \(32\times64\) structured
    angular grid, \(N_R=96\), \(24\times24\) EATS quadrature, and
    \(N_{\gamma_e}=96\). PyBlastAfterglowMag now contributes a compiled
    repo-local `pba.out` Gaussian-geometry row set. These rows are
    geometry-interface evidence, not a broad structured-jet agreement claim,
    because the public Gaussian wing treatment is not identical across packages.
  - increased `paper/source_data/figB3_asgard_complex_state.csv` from 3093 to
    5733 rows by using denser radius, electron, photon, proton, and \(\chi\)
    grids for the ASGARD-only state dashboard. The dashboard remains a
    radius-ordered state-completeness and component-visibility figure, not a
    cross-code agreement figure.
  Boundary kept explicit:
  - PyBlastAfterglowMag is used only for matched FS synchrotron and Gaussian
    off-axis FS rows.  Its RS/SSC and ASGARD complex-state entries remain
    explicit unsupported rows in the manuscript source tables;
  - \(\chi\)-resolved hadronic feedback remains an explicit unsupported boundary
    until a local hadronic state, local photon target, local secondary feedback,
    and observer projection contract exist.
- Phase 4: completed.
  - extended `tests/generate_paper_figures.py` so the formal figure entry point
    now builds `fig6_cross_code_fs_benchmark`,
    `fig7_multiphysics_geometry_benchmark`, and
    `fig8_asgard_complex_state_dashboard` from tracked source data;
  - generated PDF, SVG, and TIFF outputs for these three figures under
    `paper/figures/`;
  - set Fig. 6 to show ASGARD radio-to-VHE light curves from \(1\) to
    \(10^8\,{\rm s}\), broad-time SEDs across the full frequency grid, X-ray
    cross-code ratios on a logarithmic ratio axis, and same-machine wall time;
  - set Fig. 7 to show an optical RS/SSC component slice, a broad-grid median
    RS/SSC ratio summary, Gaussian off-axis multi-band FS light curves, and
    structured-geometry ratios;
  - set Fig. 8 to show the denser ASGARD hadronic/pair state, \(\chi\)-resolved
    SSA state, and full state chain.
  - replaced the main-text cross-code placeholder table with a capability and
    evidence map linked to `paper/source_data/table_cross_code_capabilities.csv`;
  - replaced the Fig. 6--8 placeholders in `paper/main.tex` with formal
    `\includegraphics` figures, labels, and captions that state the widened
    source-data coverage and explicit unsupported boundaries.
- Phase 5: completed.  Verification now reads as evidence interpretation rather
  than method repetition.  Discussion contains Degeneracy, State, and Boundaries
  subsections, and Model Boundaries remains an independent hard-list section.
  Each claim is tied to figures, source tables, comparison rows, or explicit
  unsupported boundaries.
- Phase 6: completed.  The Introduction, Summary, and Abstract now use the final
  paper argument: observational degeneracy, radius-ordered local state,
  comparison-code context, source-data-backed evidence, PyBlast compiled
  benchmark rows on matched subsets, and constrained supported use cases.
- Phase 7: completed.  The final validation pass ran the formal figure
  generator, text-encoding check, global `git diff --check`, Windows TeX Live
  `latexmk`, LaTeX log scan, and rendered-page inspection for the title page and
  benchmark block.  The only LaTeX warnings are the accepted AASTeX small-caps
  font substitutions.

## 2. Non-negotiable writing rules

- Main manuscript prose is English.  Process reports are Chinese.
- Every paragraph starts with the physical object or algorithmic object it
  discusses.
- Main text keeps physical equations, source terms, loss terms, projection
  equations, one representative numerical equation per solver family, module
  outputs, and hard boundaries.
- Main text removes or compresses solver-by-solver matrices, repeated ABI
  binding details, long quadrature/remap/interpolation instructions, code-name
  prose, and raw implementation controls that do not define the physical model.
- Unsupported features are written as unsupported.  Do not use fallback, repair,
  smoothing, or heuristic patch language to hide unsupported behavior.
- No benchmark number is written as fact until it exists in tracked source data
  under `paper/source_data/`.
- Historical files under `output/` are not manuscript evidence unless they are
  regenerated through the formal paper workflow.
- Do not change ASGARD physics code or the public runtime API during manuscript
  restructuring.

## 3. Final main-text architecture

The final main text should use this section order:

1. `Introduction`
2. `Model Definition and Radius-Ordered State`
3. `Dynamics and Jet/Medium Setup`
4. `Particle and Radiation Transport`
5. `Reverse Shocks`
6. `Hadronic and Pair-Feedback Boundary`
7. `Observer Projection and Component Assembly`
8. `Fitting, Runtime, and Reproducibility`
9. `Verification and Benchmark Evidence`
10. `Discussion`
11. `Model Boundaries`
12. `Summary`

The Introduction uses three blocks:

1. Observational degeneracy: the same \(F_\nu(t_{\rm obs})\) can be changed by
   dynamics, cooling, geometry, reverse shocks, SSC, hadronic losses, secondary
   pairs, and absorption.
2. High-energy and multi-messenger pressure: VHE detections, SSC
   interpretation, hadronic scenarios, and structured/off-axis jets require
   models that expose local state before projection.
3. Code ecosystem: afterglowpy, jetsimpy, VegasAfterglow, and PyBlastAfterglow
   solve important parts of the afterglow problem.  ASGARD is complementary
   because it emphasizes multi-physics state visibility and component
   accounting.

## 4. Phase 0: planning document

Input:

- current `paper.md`;
- current `paper/main.tex`;
- repository paper workflow rules.

Actions:

- replace `paper.md` with this complete execution plan;
- record phase order, artifacts, validation commands, and boundaries;
- do not edit `paper/main.tex`.

Output:

- complete `paper.md`.

Validation:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
rtk git diff --check -- paper.md
```

Acceptance:

- `paper.md` is decision complete for every later phase;
- no LaTeX compile is required if `paper/main.tex` is unchanged.

## 5. Phase 1: main-text skeleton and placeholders

Input:

- current `paper/main.tex`;
- this `paper.md`.

Actions:

- update section titles and route text to match the final architecture;
- merge the electron and photon main-line headings as `Particle and Radiation
  Transport` while leaving detailed formulas in place for later passes;
- keep reverse shocks as a core main-text section;
- change the hadronic section into a short main-text boundary section in this
  pass, with detailed transport still present until Phase 2 moves it;
- add the cross-code capability/evidence table placeholder;
- add three benchmark figure placeholders:
  - `fig:benchmark-fs`;
  - `fig:benchmark-rs-ssc-geometry`;
  - `fig:benchmark-complex-state`;
- add `Discussion` after `Verification and Benchmark Evidence` and before
  `Model Boundaries`;
- add five technical appendix entry points:
  - Electron solvers;
  - Photon, SSC, and SSA kernels;
  - Hadronic and pair feedback;
  - Reverse-shock technical details;
  - Projection, fitting, and runtime details.

Output:

- a compiling manuscript skeleton with explicit placeholders for pending
  benchmark data.

Validation:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
rtk git diff --check -- paper/main.tex paper.md
rtk cmd /c latexmk -cd -pdf -interaction=nonstopmode -halt-on-error -outdir=build paper/main.tex
```

Acceptance:

- LaTeX compiles;
- new labels are defined;
- pending benchmark values are visibly marked as pending;
- no fake benchmark number is introduced.

## 6. Phase 2: method compression and appendix migration

Input:

- compiling Phase 1 manuscript;
- code paths and formulas already described in the current manuscript.

Actions:

- work in increasing manuscript order;
- compress one physical/module block per pass;
- keep only the physical model, key formulas, one representative numerical
  update per solver family, outputs, diagnostics, and boundaries in main text;
- move detailed finite-volume, characteristic, DG, WENO, quadrature,
  interpolation, remap, and ABI material to the correct appendix;
- preserve labels that are already referenced, or update all references in the
  same pass;
- keep FS and RS as core physics in main text;
- keep SSC/IC core cooling and emissivity equations in main text;
- reduce hadronic/pair main text to activation boundary, physical objects,
  state products, and figure evidence.

Output:

- shorter main text;
- populated technical appendices;
- no loss of implemented physics claims.

Validation after each pass:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
rtk git diff --check -- paper/main.tex
rtk cmd /c latexmk -cd -pdf -interaction=nonstopmode -halt-on-error -outdir=build paper/main.tex
```

Acceptance:

- every main-text paragraph starts with a physical or algorithmic object;
- every changed paragraph ends with an output, diagnostic, constraint, or
  boundary;
- no unsupported capability is implied by moved text.

## 7. Phase 3: cross-code benchmark environment and source data

Input:

- compiling Phase 2 manuscript;
- benchmark scripts already present in `tests/` or `scripts/benchmarks/`;
- public APIs of ASGARD and comparison codes.

Actions:

- install and pin benchmark versions for:
  - afterglowpy;
  - jetsimpy;
  - VegasAfterglow;
  - PyBlastAfterglow;
- record versions in the benchmark source data;
- use common minimum physics for cross-code error metrics:
  - forward-shock synchrotron;
  - matched jet/medium/observer/microphysics parameters;
  - common frequencies and observer times;
- compare RS, SSC, and structured/off-axis geometry only when physical
  assumptions can be matched;
- write unsupported for unsupported paths, never a numerical failure value;
- report wall time as cost information, not as a speed ranking.

Required source-data outputs:

- `paper/source_data/figB1_cross_code_fs.csv`;
- `paper/source_data/figB2_rs_ssc_geometry.csv`;
- `paper/source_data/figB3_asgard_complex_state.csv`;
- `paper/source_data/table_cross_code_capabilities.csv`.

Each CSV must include:

- scenario;
- code;
- code version;
- physical assumptions;
- frequency or observer-time coordinate where relevant;
- flux or diagnostic value where relevant;
- ASGARD/code ratio where relevant;
- wall time where measured;
- unsupported boundary text where no fair comparison exists.

Validation:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
rtk git diff --check -- paper/source_data
```

Acceptance:

- every benchmark claim has tracked source data;
- every unsupported entry is explicit;
- historical `output/` files are not cited as evidence.

## 8. Phase 4: figures and comparison table

Input:

- Phase 3 CSV source data.

Actions:

- extend `tests/generate_paper_figures.py`;
- generate three new main figures:
  - `fig6_cross_code_fs_benchmark`;
  - `fig7_multiphysics_geometry_benchmark`;
  - `fig8_asgard_complex_state_dashboard`;
- export each figure as PDF, SVG, and TIFF under `paper/figures/`;
- add or update the main-text comparison table from
  `paper/source_data/table_cross_code_capabilities.csv`;
- update captions so each panel states what it tests and what it does not test.

Figure contracts:

- Fig. 6: forward-shock light curve, SED, ASGARD/code ratios, and wall-time
  cost.
- Fig. 7: left block for RS/SSC, right block for structured/off-axis geometry.
- Fig. 8: hadronic/pair state, \(\chi\)-resolved projection or downstream
  state, and full state-chain diagnostic.

Validation:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tests/generate_paper_figures.py'
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
rtk git diff --check -- paper/main.tex paper/source_data tests/generate_paper_figures.py
rtk cmd /c latexmk -cd -pdf -interaction=nonstopmode -halt-on-error -outdir=build paper/main.tex
```

Acceptance:

- figures rebuild from the formal entry point;
- all figure files in the manuscript have matching source data;
- figure captions do not overclaim beyond the panel evidence.

## 9. Phase 5: verification, discussion, and boundaries

Input:

- Phase 4 manuscript and figures.

Actions:

- rewrite Verification as evidence interpretation, not method repetition;
- write Discussion in three subsections:
  - Degeneracy;
  - State;
  - Boundaries;
- connect every Discussion claim to a figure, table, source data file, code
  path, or explicit boundary;
- keep Model Boundaries as an independent hard-list section after Discussion;
- state that new couplings require a local state, projection operator, and
  validation evidence.

Output:

- full Verification, Discussion, and Model Boundaries sequence.

Validation:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
rtk git diff --check -- paper/main.tex
rtk cmd /c latexmk -cd -pdf -interaction=nonstopmode -halt-on-error -outdir=build paper/main.tex
```

Acceptance:

- Discussion contains no generic future-work promise;
- boundaries are not softened by vague wording;
- every decision-relevant sentence has an evidence anchor.

## 10. Phase 6: Introduction, Summary, and Abstract

Input:

- stable methods, benchmark evidence, Discussion, and Boundaries.

Actions:

- rewrite Introduction after the evidence is fixed;
- rewrite Summary with three closing moves:
  - bounded contribution;
  - reproducible software/source-data evidence;
  - constrained open use cases;
- rewrite Abstract last;
- keep claims narrower than the benchmark evidence.

Output:

- final high-level manuscript argument.

Validation:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
rtk git diff --check -- paper/main.tex
rtk cmd /c latexmk -cd -pdf -interaction=nonstopmode -halt-on-error -outdir=build paper/main.tex
```

Acceptance:

- Abstract does not claim benchmark values before they are in source data;
- Introduction cites prior tools only for specific roles;
- Summary does not introduce new claims.

## 11. Phase 7: final validation and release-readiness check

Run:

```bash
rtk bash -lc 'source ~/.wsl_env && cd "/mnt/c/Users/jia/Documents/New project/ASGARD_GRBAfterglow" && uv run python tools/check_text_encoding.py'
rtk git diff --check
rtk cmd /c latexmk -cd -pdf -interaction=nonstopmode -halt-on-error -outdir=build paper/main.tex
```

Scan `paper/build/main.log` for:

- undefined citations;
- undefined references;
- overfull boxes;
- underfull boxes;
- fatal errors;
- emergency stops;
- LaTeX warnings.

Acceptance:

- no undefined citations;
- no undefined references;
- no overfull/underfull boxes unless explicitly inspected and accepted;
- only the known AASTeX Computer Modern small-caps font substitution warning is
  accepted by default.

## 12. Artifact policy

Allowed manuscript artifacts:

- `paper.md`;
- `paper/main.tex`;
- `paper/source_data/*.csv`;
- `paper/figures/*.{pdf,svg,tiff}`;
- `tests/generate_paper_figures.py`;
- benchmark scripts under `tests/` or `scripts/benchmarks/` only if they feed
  formal source data.

Forbidden as manuscript evidence:

- untracked `output/` plots;
- temporary debug scripts;
- failed placeholder figures;
- benchmark numbers not backed by `paper/source_data/`;
- undocumented manual edits to generated figures.
