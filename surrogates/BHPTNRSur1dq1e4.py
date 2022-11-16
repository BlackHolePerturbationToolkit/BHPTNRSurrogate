##==============================================================================
## BHPTNRSur1dq1e4 : arXiv:2204.01972
## Description : generates calibrated ppBHPT surrogate model BHPTNRSur1dq1e4
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np
import os
from os import path

import model_utils.BHPTNRSur1dq1e4.load_surrogate_fits as load
from common_utils import utils, fits
import common_utils.nr_calibration as nrcalib
import common_utils.check_inputs as checks
import common_utils.doc_string as docs

# h5 data directory
h5_data_dir = os.path.dirname(os.path.abspath(__file__)) + '/../data'

# load all fits data
time, fit_data_dict_1, fit_data_dict_2, B_dict_1, B_dict_2, \
                            alpha_coeffs, beta_coeffs = load.load_surrogate(h5_data_dir)

#----------------------------------------------------------------------------------------------------
# add docstring from utility
@docs.copy_doc(docs.generic_doc_for_models,docs.BHPTNRSur1dq1e4_doc)
def generate_surrogate(q, spin1=None, spin2=None, ecc=None, ano=None, modes=None, M_tot=None, \
                       dist_mpc=None, orb_phase=None, inclination=None, neg_modes=True, \
                       mode_sum=False, lmax=5, calibrated=True):
    
    # check inputs
    checks.check_user_inputs(modes, M_tot, dist_mpc, orb_phase, inclination, mode_sum)
    
    # modes
    if modes==None:
        modes=[(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),(5,4),(5,5),
               (6,4),(6,5),(6,6),(7,5),(7,6),(7,7),(8,6),(8,7),(8,8),(9,7),(9,8),
               (9,9),(10,8),(10,9)]
        
    
    # define the parameterization for surrogate
    X_sur = np.log10(q)
    
    # define parameterization for nr calibration
    X_calib = 1/q
    
    # normalization parameter to be multiplied with the surrogate waveform
    norm = 1/q
    
    # domain of validity
    X_min = [2.5]
    X_max = [10000]
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
    
    # uncalibrated waveforms in geometric units
    hsur_raw_dict = fits.all_modes_surrogate(modes, X_sur, fit_data_dict_1, fit_data_dict_2, \
                           B_dict_1, B_dict_2, lmax, fit_func, decomposition_funcs, norm)
    
    # process the raw surrogate output depending on the user inputs
    t_surrogate, h_surrogate = utils.obtain_processed_output(X_calib, time, hsur_raw_dict, alpha_coeffs, 
                                    beta_coeffs, alpha_beta_functional_form, calibrated, M_tot, dist_mpc, 
                                    orb_phase, inclination, mode_sum, neg_modes, lmax, CoorbToInert)
    
    return t_surrogate, h_surrogate
