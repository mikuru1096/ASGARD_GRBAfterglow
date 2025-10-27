# ASGARD: GRB Afterglow Analysis Tool

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
This project name is ASGARD, Retrieved from
       <https://github.com/mikuru1096/ASGARD_GRBAfterglow>
### Quick Start
The usage of this code is very simple.
Ensure you have GNU compilers installed on your system.

For Ubuntu/Debian systems:
    sudo apt install gcc g++ gfortran

Clone this repository to your local machine:
```shell
    git clone https://github.com/mikuru1096/ASGARD_GRBAfterglow
    cd ASGARD_GRBAfterglow
```
Run the installation script:
```shell
    bash install.sh
```
After compilation completes, generate your first multi-band afterglow light curve:
```shell
    python hand_my.py
```
If you already have the matplotlib package installed, the program should generate the first multi-band afterglow light curve image for you.
### Documentation
In merger.py, we have provided the basic invocation method of the program, along with simple comments for the keywords.
### Current Status
Due to current progress limitations, we are not yet able to provide a complete demonstration of the afterglow fitting workflow. 
However, please start exploring and try to integrate it into your own fitting framework!
### Web Interface
We have a website available at
       <https://hetools.xyz>
that requires no installation, for comparing the results of ASGARD and jetsimpy. Feel free to give it a try!
