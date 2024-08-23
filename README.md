# BHPTNRSurrogate(s)

Copyright (c) 2020 Black Hole Perturbation Toolkit Team

The BHPTNRSurrogate(s) package provides access to a family of surrogate 
gravitational waveform models built on waveforms generated with point-particle 
black hole perturbation theory (ppBHPT) framework. These models extend from 
comparable mass-ratio to large mass-ratio regimes and are calibrated 
to numerical relativity (NR) waveforms in the comprable-mass-ratio regime.
Many harmonic modes are included, for example, the BHPTNRSur1dq1e4 model
includes up to up to $\\ell=10$. Please see the 
BHPTNRSurrogate package [landing page](https://bhptoolkit.org/BHPTNRSurrogate/)
for more details. 

## Getting the package

The latest development version will always be available from the project git
repository:
```bash
git clone https://github.com/BlackHolePerturbationToolkit/BHPTNRSurrogate.git
cd BHPTNRSurrogate
git submodule init
git submodule update
```

## Available Models

#### 1. BHPTNRSur2dq1e3

This model can generate waveforms from non-spinning black hole binary 
systems with mass ratios varying from 3 to 1000 and spins 
from −0.8≤χ1≤0.8 on the larger black hole.  The waveforms include all
spin-weighted spherical harmonic modes up to ℓ=4, except the (4,1) and m=0 modes,
and their m<0 counterparts. The uncalibrated raw ppBHPT waveforms are ~13,500M long.

Model details can be found in the following paper:
[Gravitational wave surrogate model for spinning, intermediate mass
ratio binaries based on perturbation theory and numerical relativity](https://arxiv.org/abs/2407.18319)

#### 2. BHPTNRSur1dq1e4

This model can generate waveforms from non-spinning black hole binary 
systems with mass ratios varying from 2.5 to 10000. It supports a total of 50 
modes : [(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),
(5,4),(5,5),(6,4),(6,5),(6,6),(7,5),(7,6),(7,7),(8,6),(8,7),(8,8),(9,7),(9,8),(9,9),(10,8),(10,9)]
and their m<0 counterparts, with modes up to $\ell=5$ calibrated to NR. 
The uncalibrated raw ppBHPT waveforms are ~30,500M long.

Model details can be found in the following paper:
[Surrogate model for gravitational wave signals from non-spinning, comparable- to
large-mass-ratio black hole binaries built on black hole perturbation theory waveforms
calibrated to numerical relativity](https://arxiv.org/pdf/2204.01972.pdf)

#### 3. EMRISur1dq1e4 (deprecated)

EMRISur1dq1e4 is the predecessor of the BHPTNRSur1dq1e4 model, which includes numerous
numerous upgrades (please see [table 1](https://arxiv.org/pdf/2204.01972.pdf)). EMRISur1dq1e4 is
a non-spinning model trained for mass ratios $q=3$ to $q=10000$ and the dominant $(2,2)$ 
mode was calibrated to NR in the comparable mass ratios. The EMRISur1dq1e4 model is NOT supported in 
this package but can be accessed from [EMRISurrogate](https://bhptoolkit.org/EMRISurrogate/).
**CAUTION :** This model is outdated and we advise for using BHPTNRSurrogate(s).

# Requirements

This package requires Python 3, sklearn, hashlib, and gwtools.

```bash
pip install scikit-learn hashlib gwtools
```

Parts of the accompanying Jupyter notebook will require gwsurrogate, 
which can be installed with either pip

```bash
pip install gwsurrogate
```

or conda

```bash
conda install -c conda-forge gwsurrogate
```

Note that you do not need gwsurrogate to evaluate the EMRI surrogate model or 
run most parts of the notebook.

# Installation

1. Clone the repository
2. Download the datafiles from Zenodo

```bash
wget https://zenodo.org/records/13340319/BHPTNRSur1dq1e4.h5
wget https://zenodo.org/records/13340319/BHPTNRSur2dq1e3.h5
```

3. Simply move these files into the data directory `BHPTNRSurrogate/data/`.

# Examples

Example tutorial notebooks for the **BHPTNRSur1dq1e4** and **BHPTNRSur2dq1e3** models are available here `BHPTNRSurrogate/tutorials`.

### 1. BHPTNRSur1dq1e4

For example, 

```bash
jupyter notebook BHPTNRSur1dq1e4.ipynb
```

# Known problems

Known bugs are recorded in the project bug tracker:

https://github.com/BlackHolePerturbationToolkit/BHPTNRSurrogate/issues


# License

This code is distributed under the MIT License. Details can
be found in the LICENSE file.


# Authors

Ritesh Bachhar, Scott Field, Tousif Islam, Gaurav Khanna, Nur Rifat, Katie Rink, Vijay Varma

# Citation guideline

If you make use of any module from the [Black Hole Perturbation Toolkit](https://bhptoolkit.org/) in your research, please acknowledge using:

> This work makes use of the Black Hole Perturbation Toolkit.

If you make use of the BHPTNRSur models, please cite the following papers:

```
@article{Islam:2022laz,
    author = "Islam, Tousif and Field, Scott E. and Hughes, Scott A. and Khanna, Gaurav and Varma, Vijay and Giesler, Matthew and Scheel, Mark A. and Kidder, Lawrence E. and Pfeiffer, Harald P.",
    title = "{Surrogate model for gravitational wave signals from nonspinning, comparable-to large-mass-ratio black hole binaries built on black hole perturbation theory waveforms calibrated to numerical relativity}",
    eprint = "2204.01972",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.106.104025",
    journal = "Phys. Rev. D",
    volume = "106",
    number = "10",
    pages = "104025",
    year = "2022"
}

@article{Rink:2024swg,
    author = "Rink, Katie and Bachhar, Ritesh and Islam, Tousif and Rifat, Nur E. M. and Gonzalez-Quesada, Kevin and Field, Scott E. and Khanna, Gaurav and Hughes, Scott A. and Varma, Vijay",
    title = "{Gravitational wave surrogate model for spinning, intermediate mass ratio binaries based on perturbation theory and numerical relativity}",
    eprint = "2407.18319",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "7",
    year = "2024"
}
```
