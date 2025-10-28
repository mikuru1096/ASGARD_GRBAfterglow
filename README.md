# ASGARD: GRB Afterglow Analysis Tool

## A Standard GRB afterglow Radiation Diagnoser (ASGARD) is a state-of-art simulation code for GRB afterglow. 

The code is entirely based on numerical partial differential equation methods to solve the evolution of the afterglow electron spectrum, while precisely handling the cooling process of electrons via Compton scattering. It self-consistently computes synchrotron radiation and synchrotron self-Compton (SSC) radiation using high-order integration schemes and fully accounts for observational effects. The calculated afterglow radiation spectrum covers the entire electromagnetic range from radio to very high energies (VHE), with proper treatment of synchrotron self-absorption (SSA, based on radiative transfer) and $\gamma\gamma$ annihilation (based on scattering cross-sections).

The code's greatest strengths lie in its exceptional computational efficiency and accuracy. When employing the default first-order fully implicit scheme to solve the electron continuity equation, the code can rapidly and stably generate results even under extreme conditions, such as very strong magnetic fields ($\epsilon_B \sim 1$) and high ambient medium particle densities ($n \gg 10^6 \, \text{cm}^{-3}$). When using the WENO5 scheme, the code achieves fifth-order accuracy in resolving the electron spectrum. The calculation of synchrotron and SSC radiation is based on the composite Simpson's method (O($N^4$)).

**ASGARD** is written in `Fortran`, and its computational processes are highly parallelized using `OpenMP`. When combined with `MPI` parallelization schemes employing `emcee` or `pymultinest` samplers, the code can operate with extremely high efficiency on personal computers, workstations, and even computing clusters. For instance, in the case of on-axis-viewed top-hat jet synchrotron radiation, sampling a million times on a single-node dual-socket EPYC 9754 system requires only a few hours.

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
Ensure you have `GNU` compilers, `python <= 3.11`, `numpy`, `astropy`, `scipy`, and `matplotlib` installed on your system.

For Ubuntu/Debian systems:
```shell
sudo apt install gcc g++ gfortran
```
Clone this repository to your local machine:
```shell
git clone https://github.com/mikuru1096/ASGARD_GRBAfterglow
cd ASGARD_GRBAfterglow
```
Run the installation script:
```shell
bash install.sh
```
After compilation completes, run:
```shell
python hand_my.py
```
The program should generate the first multiband afterglow light curve image for you.
### Documentation
In `mergered.py`, we have provided the basic invocation method of the program, along with simple comments for the keywords.
### Current Status
Due to current progress limitations, we are not yet able to provide a complete demonstration of the afterglow fitting workflow. 
However, please start exploring and try to integrate it into your own fitting framework!
### Web Interface
We have a website available at <https://hetools.xyz>  
that requires no installation, for comparing the results of **ASGARD** and **jetsimpy**. Feel free to give it a try!
