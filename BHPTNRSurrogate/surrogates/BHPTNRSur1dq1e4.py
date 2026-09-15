##==============================================================================
## BHPTNRSur1dq1e4 : arXiv:2204.01972
## Description : generates calibrated ppBHPT surrogate model BHPTNRSur1dq1e4
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import warnings

import numpy as np

from .._data import get_data_dir
from .model_utils import load_surrogates as load
from .model_utils import eval_surrogates as eval_sur
from .common_utils import utils, fits
from .common_utils import nr_calibration as nrcalib
from .common_utils import doc_string as docs

# h5 data directory
h5_data_dir = get_data_dir()

# lazy-loaded surrogate data cache
_surrogate_data = {}

_SURROGATE_KEYS = (
    'time', 'fit_data_dict_1', 'fit_data_dict_2',
    'B_dict_1', 'B_dict_2', 'alpha_coeffs', 'beta_coeffs',
)

def _ensure_loaded():
    """Load the surrogate data from H5 file on first access."""
    if not _surrogate_data:
        time, fit_data_dict_1, fit_data_dict_2, B_dict_1, B_dict_2, \
            alpha_coeffs, beta_coeffs = load.load_BHPTNRSur1dq1e4_surrogate(h5_data_dir)
        _surrogate_data['time'] = time
        _surrogate_data['fit_data_dict_1'] = fit_data_dict_1
        _surrogate_data['fit_data_dict_2'] = fit_data_dict_2
        _surrogate_data['B_dict_1'] = B_dict_1
        _surrogate_data['B_dict_2'] = B_dict_2
        _surrogate_data['alpha_coeffs'] = alpha_coeffs
        _surrogate_data['beta_coeffs'] = beta_coeffs

def __getattr__(name):
    """Lazy access to surrogate data attributes at the module level."""
    if name in _SURROGATE_KEYS:
        _ensure_loaded()
        return _surrogate_data[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

#----------------------------------------------------------------------------------------------------
# add docstring from utility
@docs.copy_doc(docs.generic_doc_for_models,docs.BHPTNRSur1dq1e4_doc)
def generate_surrogate(q, spin1=None, spin2=None, ecc=None, ano=None, modes=None, M_tot=None, \
                       dist_mpc=None, orb_phase=None, inclination=None, neg_modes=True, \
                       mode_sum=False, lmax=5, calibrated=True, mass_scale='M'):

    _ensure_loaded()

    # modes modelled in the surrogate
    modes_available = [(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),(5,4),(5,5),
               (6,4),(6,5),(6,6),(7,5),(7,6),(7,7),(8,6),(8,7),(8,8),(9,7),(9,8),
               (9,9),(10,8),(10,9)]

    # Warning to user if inputs include spin or eccentricity
    ignored = {k for k, v in [("spin1", spin1), ("spin2", spin2), ("ecc", ecc), ("ano", ano)] if v is not None}
    if ignored:
        warnings.warn(
            "Model only takes [q] as input. Ignoring extra params: %s" % ", ".join(sorted(ignored)),
            stacklevel=2,
        )

    # validate mass_scale parameter
    if mass_scale not in ('M', 'm1'):
        raise ValueError("mass_scale must be 'M' or 'm1', got %r" % mass_scale)

    if mass_scale == 'm1' and M_tot is not None and dist_mpc is not None:
        raise ValueError(
            "mass_scale='m1' is only valid for geometric waveforms; physical waveforms require mass_scale='M'"
        )

    if calibrated and mass_scale != 'M':
        raise ValueError("mass_scale='m1' cannot be used when calibrated=True; calibrated waveforms use total mass M")

    # Raw ppBHPT data use m1 units. Rescale the default uncalibrated output to total-mass units;
    # mass_scale='m1' preserves the behavior from versions before 0.2.0.
    if not calibrated and mass_scale == 'M':
        mass_factor = 1.0 / (1.0 + 1.0/q)
    else:
        mass_factor = 1.0

    # modes requested
    if modes==None:
        modes = modes_available

    # define the parameterization for surrogate
    X_sur = np.log10(q)

    # define parameterization for nr calibration
    X_calib = 1/q

    # normalization parameter to be multiplied with the surrogate waveform
    norm = 1/q

    # domain of validity
    X_min = [np.log10(2.5)]
    X_max = [np.log10(10000)]
    X_bounds = [X_min, X_max]

    # fit type
    fit_func = 'spline_1d'

    # data decomposition functions for 22 mode and HMs
    decomposition_funcs = [utils.amp_ph_to_comp, utils.re_im_to_comp]

    # nr calibratiin function
    alpha_beta_functional_form = nrcalib.alpha_beta_BHPTNRSur1dq1e4

    # tell whether the higher modes needed to be transformed from coorbital
    # to inertial frame
    CoorbToInert = True

    # generate surrogate waveform
    t_surrogate, h_surrogate = eval_sur.evaluate_surrogate(X_sur, X_calib, X_bounds,
                                        _surrogate_data['time'], modes,
                                        modes_available, _surrogate_data['alpha_coeffs'],
                                        _surrogate_data['beta_coeffs'], alpha_beta_functional_form,\
                                        calibrated, M_tot, dist_mpc, orb_phase, inclination,
                                        _surrogate_data['fit_data_dict_1'], \
                                        _surrogate_data['fit_data_dict_2'],
                                        _surrogate_data['B_dict_1'], _surrogate_data['B_dict_2'],
                                        fit_func, decomposition_funcs,\
                                        norm, mode_sum, neg_modes, lmax, CoorbToInert,
                                        mass_factor=mass_factor)

    return t_surrogate, h_surrogate
