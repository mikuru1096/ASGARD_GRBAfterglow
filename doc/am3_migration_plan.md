# AM3 Migration Plan For ASGARD

## Scope

This note compares:

- ASGARD's current runtime architecture
- ASGARD's actual target use cases
- AM3's one-zone hadronic process architecture

The goal is not to port AM3 wholesale. The goal is to decide:

- what must still migrate into ASGARD
- what should be delayed
- what should be explicitly abandoned


## Current Reality

ASGARD is not a generic homogeneous lepto-hadronic box solver.

Its real driver is:

`dynamics -> electron solver -> photon target field -> hadronic add-on -> observer projection`

The practical application scenarios are:

1. Forward-shock GRB afterglow broadband modeling
2. Observer-side light curves and SEDs under shell evolution
3. Parameter scans and fitting, where runtime cost matters
4. Hadronic diagnostics on top of an existing blast-wave/electron pipeline

This means ASGARD fundamentally differs from AM3 in one way:

- AM3 is organized around a self-contained time-dependent one-zone PDE closure.
- ASGARD is organized around shell-evolving blast-wave dynamics with observer mapping.

So the migration target must be:

- AM3 microphysics kernels and process structure
- not AM3's global driver architecture


## Architecture Gap

### What AM3 has structurally

AM3 exposes a one-zone process stack with:

- transport/state evolution
- acceleration/injection
- photopion
- Bethe-Heitler
- pair production / cascade
- pion decay
- muon decay
- proton-proton
- hadronic synchrotron
- hadronic inverse Compton
- escape
- expansion
- extra cooling
- a unified photon / secondary solver tree

### What ASGARD currently has structurally

ASGARD already has active hadronic kernels for:

- photopion
- decay
- Bethe-Heitler
- pair production
- proton-proton
- hadronic inverse Compton
- explicit secondary-species transport
- secondary radiation
- acceleration / injection

But the remaining closure is uneven:

- some AM3-derived logic is still wrapped through Python orchestration
- photon-loss / cascade closure is only partial
- hadronic support is still effectively 1D forward-shock only
- benchmark closure is incomplete


## Scenario-Driven Decision

### Must migrate

These are necessary for ASGARD's real use cases.

1. `photopion photon-loss` as a local hadronic-stage operator  
Reason:
- AM3 treats it as part of the physical interaction closure
- ASGARD currently only converts it to observer-side optical depth
- this is insufficient for strict in-chain closure

2. `remaining AM3-derived core kernels` out of Python wrappers and into `src/Hadronic/*.f90`  
Reason:
- current project rule already requires final microphysics to live in Fortran
- Python must remain orchestration only

3. `AM3-consistent benchmark closure`  
Reason:
- without regenerated benchmark figures and smoothness checks, the current hadronic chain is still not validated

4. `proton-synchrotron peak validation` against AM3-side physical expectations  
Reason:
- this directly affects interpretation of the hadronic SED

### Should migrate, but after the above

1. `full gamma-gamma / pair-production cascade PDE`
Reason:
- physically important for the highest-energy closure
- but not the first blocker for the main forward-shock workflow

2. `pp pi0-into-cascade treatment`
Reason:
- present in AM3
- useful for a fuller gamma-ray closure
- lower priority than photopion and benchmark closure

3. `reverse-shock hadronic`
Reason:
- scientifically useful
- but not on the critical path for the current forward-shock production workflow

4. `2D / chi-resolved hadronic transport`
Reason:
- desirable long-term
- but expensive and not required for the current 1D formal closure

### Explicitly do not port as architecture

These may be referenced, but should not be migrated as a framework layer.

1. `SimulationManager / PhysicsHandler / RunParams / AM3Arrays` as a top-level architecture  
Reason:
- that would overwrite ASGARD's shell-evolution design with a one-zone driver

2. `Escape / Expansion / ExtraCooling` as one-to-one standalone framework modules  
Reason:
- ASGARD already handles the relevant global evolution through blast-wave dynamics and shell history
- direct structural port would duplicate responsibilities

3. `AM3's full photon-solver tree` as the master runtime abstraction  
Reason:
- ASGARD needs process kernels, not a replacement for its observer-driven chain


## Priority Plan

### P0

1. Move remaining AM3-derived hadronic core logic into Fortran
2. Add local photopion photon-loss closure inside the hadronic stage
3. Re-run hadronic benchmark figures and smoothness checks

### P1

1. Add direct AM3 comparison smokes for:
- Bethe-Heitler
- acceleration / gamma max
- secondary transport

2. Close proton-synchrotron peak validation

### P2

1. Extend pair-production toward a real cascade PDE
2. Evaluate whether `pp pi0 into cascade` is needed for the intended science cases

### P3

1. Reverse-shock hadronic
2. 2D / chi-resolved hadronic transport


## Immediate Next Steps

The next practical implementation sequence should be:

1. Keep the new AM3 reference smoke working
2. Add AM3 comparison smoke for Bethe-Heitler
3. Add AM3 comparison smoke for acceleration / `gamma_p,max`
4. Move photopion photon-loss from observer-side closure into a local hadronic-stage operator
5. Re-run benchmark plots and reject any non-smooth result as a bug


## Abandon List

For the current ASGARD roadmap, the following should be treated as intentionally not in scope:

1. Replacing ASGARD's shell driver with AM3's one-zone driver
2. Porting AM3 framework classes one-to-one
3. Building a generic plugin architecture around hadronic kernels
4. Porting low-value framework abstractions that do not improve the forward-shock production path
