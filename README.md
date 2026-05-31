# ASGARD: GRB Afterglow Analysis Tool

## A Standard GRB afterglow Radiation Diagnoser (ASGARD) is a state-of-the-art simulation code for GRB afterglow. 

The code is entirely based on numerical partial differential equation methods to solve the evolution of the afterglow electron spectrum, while precisely handling the cooling process of electrons via Compton scattering. It self-consistently computes synchrotron radiation and synchrotron self-Compton (SSC) radiation using high-order integration schemes and fully accounts for observational effects. The calculated afterglow radiation spectrum covers the entire electromagnetic range from radio to very high energies (VHE), with proper treatment of synchrotron self-absorption (SSA, based on radiative transfer) and $\gamma\gamma$ annihilation (based on scattering cross-sections).

The code's greatest strengths lie in its exceptional computational efficiency and accuracy. When employing the default first-order fully implicit scheme to solve the electron continuity equation, the code can rapidly and stably generate results even under extreme conditions, such as very strong magnetic fields ($\epsilon_B \sim 1$) and high ambient medium particle densities ($n \gg 10^6 \, \text{cm}^{-3}$). When using the WENO5 scheme, the code achieves fifth-order accuracy in resolving the electron spectrum. The calculation of synchrotron and SSC radiation is based on the composite Simpson's method (O($N^4$)).

**ASGARD** is written in `Fortran`, and its computational processes are highly parallelized using `OpenMP`. When combined with `MPI` parallelization schemes employing `emcee` or `pymultinest` samplers, the code can operate with extremely high efficiency on personal computers, workstations, and even computing clusters. For instance, in the case of on-axis-viewed top-hat jet synchrotron radiation, sampling a million times on a single-node dual-socket EPYC 9754 system requires only a few hours.

## Highlights
1. Numerical PDE solvers for forward-shock electron transport, including multiple 1D schemes and active 2D paths.
2. Self-consistent synchrotron, SSC, SSA, and $\gamma\gamma$ attenuation in the formal observer chain.
3. Reverse-shock electron synchrotron, RS SSC, cross-zone IC, and hadronic dispatch in the current runtime.
4. Active 1D hadronic paths for forward shock and reverse-shock comparison workflows.

## Documentation

The documentation entry point is `doc/index.md`.

Recommended reading:

- `doc/installation.md` — environment, native build, demo commands.
- `doc/user_guide.md` — common `Model` workflows: light curves, spectra, RS, hadronic, polarization, sky image and fitting.
- `doc/public_api.md` — current public API contracts.
- `doc/physics_model.md` — implemented physics modules and explicit boundaries.
- `doc/numerical_methods.md` — Fortran kernels, solver families and numerical validation targets.
- `doc/validation_and_benchmarks.md` — build gates, smoke tests, benchmark refresh and artifact policy.
- `doc/developer_guide.md` — development workflow and review checklist.
- `doc/code_overview.md`, `doc/source_tree.md`, `doc/call_chain.md` — implementation map.

Current development baseline:

- `AGENTS.md`
- `PLAN.md`
- `doc/code_overview.md`
- `doc/public_backend_limits.md`

Benchmark and comparison scripts live under `tests/` and `scripts/benchmarks/`. Generated benchmark figures under `output/` are artifacts; commit them only when they are reproducible documentation assets following `doc/benchmark_refresh_protocol.md`.

## License
**Copyright (c) 2025 Jia Ren**  

This source code is governed by the **BSD 3-Clause License**.

## Attribution Requirement
If you use, adapt, or reference the core algorithms from this project in other software projects (whether open-source or proprietary), you are required to provide explicit attribution to this original code project in your project's documentation, 'About' section, or any publicly published papers.

### Recommended Citation Format
```bibtex
@ARTICLE{2024ApJ...962..115R,
       author = {{Ren}, Jia and {Wang}, Yun and {Dai}, Zi-Gao},
       title = "{Jet Structure and Burst Environment of GRB 221009A}",
       journal = {\apj},
       keywords = {Gamma-ray bursts, 629, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2024,
        month = feb,
       volume = {962},
       number = {2},
          eid = {115},
        pages = {115},
          doi = {10.3847/1538-4357/ad1bcd},
       archivePrefix = {arXiv},
       eprint = {2310.15886},
       primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024ApJ...962..115R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

```

This project name is **ASGARD**, Retrieved from
<https://github.com/mikuru1096/ASGARD_GRBAfterglow>

### Quick Start
The usage of this code is very simple.
Ensure you have `GNU` compilers, `python >=3.8`, `numpy`, `astropy`, `scipy`, and `matplotlib` installed on your system.
For Ubuntu/Debian systems:
```shell
sudo apt install gcc g++ gfortran
```

Clone this repository to your local machine:
```shell
git clone https://github.com/mikuru1096/ASGARD_GRBAfterglow
cd ASGARD_GRBAfterglow
```
Install dependencies
```shell
pip install -r Requirements.txt
```
Run the automatic installation script. It detects Windows vs Ubuntu/Linux and builds the native modules accordingly:
```shell
python install.py
```
On Ubuntu/Linux you can also use:
```shell
bash install.sh
```
On Windows PowerShell:
```powershell
.\install.ps1
```
After compilation completes, run:
```shell
python lc_spec_demo.py
```

If you already have the matplotlib package installed, the program should generate the first multi-band afterglow light curve image for you.
### Documentation
The main public entry point is `lc_spec_demo.py` for the simplest end-to-end demo. Full project documentation starts at `doc/index.md`.
### Current Status
The public runtime is usable for forward-shock afterglow calculations and benchmark workflows.
Reverse shock uses local shock-front `gamma34` for new electron injection, while region-3 magnetic field and post-crossing thermal evolution are closed by the explicit `U3/V3` thermal state. VegasAfterglow is a comparison backend, not the physical target for RS closure.

The forward-shock hadronic branch is a formal 1D research path. Reverse-shock hadronic has a light proton-synchrotron path plus a full-chain path for pγ/BH/pp/secondary/cascade coupling through the formal 1D hadronic kernels. The current pair cascade path is a shell-sequence time-dependent γγ pair/synch cascade. Chi-resolved hadronic transport and inverse-Compton-mediated electromagnetic cascades remain open work.
### Web Interface
We have a website available at <https://hetools.xyz>  
that requires no installation, for comparing the results of **ASGARD** and **jetsimpy**. Feel free to give it a try!
