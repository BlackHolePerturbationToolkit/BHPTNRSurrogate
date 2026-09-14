##==============================================================================
## BHPTNRSur2dq1e3 : arXiv:2407.18319
## Description : generates calibrated ppBHPT surrogate model BHPTNRSur2dq1e3
## Author : Katie Rink, Dec 2022 [krink@utexas.edu]
## Modified : Tousif Islam, Jul 2023
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
    'times_dict', 'fit_data_dict_1_sign', 'fit_data_dict_2_sign',
    'B_dict_1_sign', 'B_dict_2_sign', 'alpha_coeffs', 'beta_coeffs',
)

def _ensure_loaded():
    """Load the surrogate data from H5 file on first access."""
    if not _surrogate_data:
        times_dict, fit_data_dict_1_sign, fit_data_dict_2_sign, B_dict_1_sign, B_dict_2_sign, \
            alpha_coeffs, beta_coeffs = load.load_BHPTNRSur2dq1e3_surrogate(h5_data_dir)
        _surrogate_data['times_dict'] = times_dict
        _surrogate_data['fit_data_dict_1_sign'] = fit_data_dict_1_sign
        _surrogate_data['fit_data_dict_2_sign'] = fit_data_dict_2_sign
        _surrogate_data['B_dict_1_sign'] = B_dict_1_sign
        _surrogate_data['B_dict_2_sign'] = B_dict_2_sign
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
@docs.copy_doc(docs.generic_doc_for_models,docs.BHPTNRSur2dq1e3_doc)
def generate_surrogate(q, spin1=0.0, spin2=None, ecc=None, ano=None, modes=None, M_tot=None, dist_mpc=None,
                       orb_phase=None, inclination=None, neg_modes=True, mode_sum=False, lmax=4,
                       calibrated=True, mass_scale='M'):

    _ensure_loaded()

    # list the modes modelled in BHPTNRSur2dq1e3
    modes_available = [(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4)]
    if modes==None:
        modes = modes_available

    # Warning to user if inputs include secondary spin or eccentricity
    ignored = {k for k, v in [("spin2", spin2), ("ecc", ecc), ("ano", ano)] if v is not None}
    if ignored:
        warnings.warn(
            "Model only takes [q, spin1] as input. Ignoring extra params: %s" % ", ".join(sorted(ignored)),
            stacklevel=2,
        )

    # validate mass_scale parameter
    if mass_scale not in ('M', 'm1'):
        raise ValueError("mass_scale must be 'M' or 'm1', got %r" % mass_scale)

    if calibrated and mass_scale != 'M':
        warnings.warn(
            "mass_scale is ignored when calibrated=True (NR calibration already uses total mass M)",
            stacklevel=2,
        )

    # compute mass_factor for uncalibrated waveforms
    if not calibrated and mass_scale == 'M':
        mass_factor = 1.0 / (1.0 + 1.0/q)
    else:
        mass_factor = 1.0

    # this model provide fits for the positive spin and negative spin cases differently
    # choose appropriate fit params here depending on the input spin value
    if spin1 < 0.0:
        times = _surrogate_data['times_dict']['negative_spin']
        fit_data_dict_1 = _surrogate_data['fit_data_dict_1_sign']['negative_spin']
        fit_data_dict_2 = _surrogate_data['fit_data_dict_2_sign']['negative_spin']
        B_dict_1 = _surrogate_data['B_dict_1_sign']['negative_spin']
        B_dict_2 = _surrogate_data['B_dict_2_sign']['negative_spin']
    else:
        times = _surrogate_data['times_dict']['positive_spin']
        fit_data_dict_1 = _surrogate_data['fit_data_dict_1_sign']['positive_spin']
        fit_data_dict_2 = _surrogate_data['fit_data_dict_2_sign']['positive_spin']
        B_dict_1 = _surrogate_data['B_dict_1_sign']['positive_spin']
        B_dict_2 = _surrogate_data['B_dict_2_sign']['positive_spin']

    # define the parameterization for surrogate
    X_sur = [np.log10(q), spin1]

    # define parameterization for nr calibration
    X_calib = [q, spin1]

    # normalization parameter to be multiplied with the surrogate waveform
    norm = 1/q

    # domain of validity
    X_min_q = np.log10(3)
    X_max_q = np.log10(1000)
    X_min_chi = -0.8
    X_max_chi = 0.8
    X_bounds = [[X_min_q, X_min_chi],[X_max_q, X_max_chi]]

    # fit type
    fit_func = 'GPR_fits'

    # data decomposition functions for each mode
    decomposition_funcs = [utils.amp_ph_to_comp, utils.amp_ph_to_comp]

    # nr calibratiin function
    alpha_beta_functional_form = nrcalib.alpha_beta_BHPTNRSur2dq1e3

    # tell whether the higher modes needed to be transformed from coorbital
    # to inertial frame
    CoorbToInert = False

    # generate surrogate waveform
    t_surrogate, h_surrogate = eval_sur.evaluate_surrogate(X_sur, X_calib, X_bounds, times, modes,\
            modes_available, _surrogate_data['alpha_coeffs'], _surrogate_data['beta_coeffs'],
            alpha_beta_functional_form,\
            calibrated, M_tot, dist_mpc, orb_phase, inclination, fit_data_dict_1,\
            fit_data_dict_2, B_dict_1, B_dict_2, fit_func, decomposition_funcs,\
            norm, mode_sum, neg_modes, lmax, CoorbToInert,
                                        mass_factor=mass_factor)

    return t_surrogate, h_surrogate
