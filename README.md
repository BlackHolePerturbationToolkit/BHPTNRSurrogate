# BHPTNRSurrogate(s)

Copyright (c) 2020 Black Hole Perturbation Toolkit Team

The BHPTNRSurrogate(s) package provides access to a family of surrogate 
gravitational waveform models built on waveforms generated with point-particle 
black hole perturbation theory (ppBHPT) framework. These models extend from 
comparable mass-ratio to large mass-ratio regimes and are tuned 
to numerical relativity (NR) waveforms at the comprable-mass-ratio regime.
Many harmonic modes are included, for example the BHPTNRSur1dq1e4 model
includes up to up to $\\ell=10$. Please see the the 
BHPTNRSurrogate package [landing page](https://bhptoolkit.org/BHPTNRSurrogate/)
for more detail. 

## Getting the package

The latest development version will always be available from the project git
repository:
```bash
git clone https://github.com/BlackHolePerturbationToolkit/BHPTNRSurrogate.git
```

## Available Models

#### 1. BHPTNRSur1dq1e4

This model can generate waveforms from a merging non-spinning black hole binary 
systems with mass-ratios varying from 2.5 to 10000. It supports a total of 50 
modes : [(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),
(5,4),(5,5),(6,4),(6,5),(6,6),(7,5),(7,6),(7,7),(8,6),(8,7),(8,8),(9,7),(9,8),(9,9),(10,8),(10,9)]
and their m<0 counterparts. Uncalibrated raw ppBHPT waveforms are ~30,500M long.
Modes up to $\ell=5$ are NR-calibrated. The model has been further tested against
state-of-art NR simulations at mass ratio $q=15,16,30,32$.

Model details can be found in the following paper:
[Surrogate model for gravitational wave signals from non-spinning, comparable- to
large-mass-ratio black hole binaries built on black hole perturbation theory waveforms
calibrated to numerical relativity](https://arxiv.org/pdf/2204.01972.pdf)

#### 2. EMRISur1dq1e4 (deprecated)

EMRISur1dq1e4 is the predecessor of the BHPTNRSur1dq1e4 model, which includes numerous
numerous upgrades (please see [table 1](https://arxiv.org/pdf/2204.01972.pdf)). EMRISur1dq1e4 is
a non-spinning model trained for mass ratios $q=3$ to $q=10000$ and the dominant $(2,2)$ 
mode was calibrated to NR in the comparable mass ratios. The EMRISur1dq1e4 model is NOT supported in 
this package but can be accessed from [EMRISurrogate](https://bhptoolkit.org/EMRISurrogate/).
**CAUTION :** This model is outdated and we advise for using BHPTNRSurrogate(s).

# Requirements

This package requires Python 3 and the sklearn package. Parts of the accompanying
Jupyter notebook will require gwsurrogate, which can be installed with either
pip

```bash
pip install gwsurrogate
```

or conda

```bash
conda install -c conda-forge gwsurrogate
```

Note that you do not need gwsurrogate to evalulate the EMRI surrogate model or 
run most parts of the notebook.

# Installation

1. Clone the repository
2. Download the datafile from Google Drive (this data will eventually be hosted in Zenodo)

```bash
https://drive.google.com/file/d/1oDs6zwT8oHsDMDx0t9ip-FF6JiQANpct/view?usp=sharing
```
3. Simply move the BHPTNRSur1dq1e4.h5 file into the data directory `BHPTNRSurrogate/data/`.

# Example Usage

### 1. BHPTNRSur1dq1e4
Example tutorial notebooks for the **BHPTNRSur1dq1e4** are available here `BHPTNRSurrogate/tutorials/BHPTNRSur1dq1e4/`
Please see the accompanying Jupyter notebook for example use of the model

```bash
jupyter notebook BHPTNRSur1dq1e4_notebook.ipynb
```
Please see the notebook for a comparison between **BHPTNRSur1dq1e4** and **NRHybSur3dq8** in the 
comparable mass ratio regime. 
```bash
jupyter notebook comparison_BHPTNRSur1dq1e4_NRHybSur3dq8.ipynb
```

# Known problems

Known bugs are recorded in the project bug tracker:

https://github.com/BlackHolePerturbationToolkit/BHPTNRSurrogate/issues


# License

This code is distributed under the MIT License. Details can
be found in the LICENSE file.


# Authors

Scott Field, Tousif Islam, Gaurav Khanna, Nur Rifat, Vijay Varma

# Citation guideline

If you make use of any module from the Toolkit in your research please acknowledge using:

> This work makes use of the Black Hole Perturbation Toolkit.

If you make use of the BHPTNRSur models please cite the following paper:

```
@article{Islam:2022laz,
    author = "Islam, Tousif and Field, Scott E. and Hughes, Scott A. and Khanna, Gaurav and Varma, Vijay and Giesler, Matthew and Scheel, Mark A. and Kidder, Lawrence E. and Pfeiffer, Harald P.",
    title = "{Surrogate model for gravitational wave signals from non-spinning, comparable- to large-mass-ratio black hole binaries built on black hole perturbation theory waveforms calibrated to numerical relativity}",
    eprint = "2204.01972",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "4",
    year = "2022"
}
```
