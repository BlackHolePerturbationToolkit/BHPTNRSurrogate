# BHPTNRSurrogate(s)

Copyright (c) 2020 Black Hole Perturbation Toolkit Team

The BHPTNRSurrogate(s) package provides access to a family of surrogate 
gravitational waveform models built on waveforms generated with point-particle 
black hole perturbation theory (ppBHPT) framework. These models extend from 
comparable mass-ratio to large mass-ratio regime and are tuned 
to numerical relativity (NR) waveforms at the comprable-mass-ratio regime.
These models has many higher order modes (at the least up to $\\ell=5$) apart 
from the dominant [(2,2)] mode of radiation. The m<0 modes are deduced from 
the m>0 modes.

A previous version of the model (which was kknown as EMRISurrogate) can be accessed
from the following repository which we keep it for record
```bash
git clone https://github.com/BlackHolePerturbationToolkit/EMRISurrogate.git
```
We, however, advise for using BHPTNRSurrogate(s) from now on.

# Available Models

### 1. BHPTNRSur1dq1e4

This model can generate waveforms from a merging non-spinning black hole binary 
systems with mass-ratios varying from 2.5 to 10000. It supports a total of 50 
modes : [(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),
(5,4),(5,5),(6,4),(6,5),(6,6),(7,5),(7,6),(7,7),(8,6),(8,7),(8,8),(9,7),(9,8),(9,9),(10,8),(10,9)]
and their m<0 counterparts.

Model details can be found in the following paper:
[Surrogate model for gravitational wave signals from non-spinning, comparable- to
large-mass-ratio black hole binaries built on black hole perturbation theory waveforms
calibrated to numerical relativity] (https://arxiv.org/pdf/2204.01972.pdf)

## Getting the package

The latest development version will always be available from the project git
repository:

```bash
git clone https://github.com/BlackHolePerturbationToolkit/EMRISurrogate.git
```

# Requirements

This package requires Python 3 and the sklearn package. Parts of the accompanying
Jupyter notebook will require gwsurrogate, which can be installed with 

```bash
pip install gwsurrogate
```

Note that you do not need gwsurrogate to evalulate the EMRI surrogate model or 
run most parts of the notebook.

# Installation

1. Clone the repository
2. Download the datafile (hosted on [Zenodo](https://zenodo.org/record/3612600#.YYPdG3VKg5k))

```bash
wget https://zenodo.org/record/3612600/files/EMRISur1dq1e4.h5
```

3. Create a static link from in the main directory to the EMRISur1dq1e4.h5 file.
Or simply move the EMRISur1dq1e4.h5 file into the main directory.

# Usage

Please see the accompanying Jupyter notebook

```bash
jupyter notebook EMRISur1dq1e4.ipynb
```

# Examples

Examples are included in the Jupyter notebook.

# Known problems

Known bugs are recorded in the project bug tracker:

https://github.com/BlackHolePerturbationToolkit/EMRISurrogate/issues


# License

This code is distributed under the MIT License. Details can
be found in the LICENSE file.


# Authors

Tousif Islam 
Scott Field   
Gaurav Khanna    
Vijay Varma

# Citation guideline

If you make use of any module from the Toolkit in your research please acknowledge using:

> This work makes use of the Black Hole Perturbation Toolkit.

If you make use of the EMRI surrogate model, EMRISur1dq1e4, please cite Paper 1:

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

@article{rifat2019surrogate,
  title={A Surrogate Model for Gravitational Wave Signals from Comparable-to Large-Mass-Ratio Black Hole Binaries},
  author={Rifat, Nur EM and Field, Scott E and Khanna, Gaurav and Varma, Vijay},
  journal={arXiv preprint arXiv:1910.10473},
  year={2019}
}
```
