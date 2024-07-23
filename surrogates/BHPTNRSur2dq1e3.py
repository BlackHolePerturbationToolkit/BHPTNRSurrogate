##==============================================================================
## BHPTNRSur2dq1e3 : arXiv:XXXX.YYYYY
## Description : generates calibrated ppBHPT surrogate model BHPTNRSur2dq1e3
## Author : Katie Rink, Dec 2022 [krink@utexas.edu]
## Modified : Tousif Islam, Jul 2023
##==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import h5py
from gwtools import gwtools as _gwtools
import os
from os import path
import subprocess

import model_utils.load_surrogates as load
import model_utils.eval_surrogates as eval_sur
from common_utils import utils, fits
import common_utils.nr_calibration as nrcalib
import common_utils.doc_string as docs

# h5 data directory
h5_data_dir = os.path.dirname(os.path.abspath(__file__)) + '/../data'

# load all fits data
# Here each of the data file are contains two separate spins
times_dict, fit_data_dict_1_sign, fit_data_dict_2_sign, B_dict_1_sign, B_dict_2_sign, \
    alpha_coeffs, beta_coeffs = load.load_BHPTNRSur2dq1e3_surrogate(h5_data_dir)

print("**** Surrogate loaded ****")

#---------------------------------------------------------------------------------------------------- 
# add docstring from utility
@docs.copy_doc(docs.generic_doc_for_models,docs.BHPTNRSur2dq1e3_doc)
def generate_surrogate(q, spin1=0.0, spin2=None, ecc=None, ano=None, modes=None, M_tot=None, dist_mpc=None, orb_phase=None, inclination=None, neg_modes=False, mode_sum=False, lmax=4, calibrated=True):

    modes_available = [(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4)]

    if modes==None:
        modes = modes_available

    if spin1 < 0.0:
        times = times_dict['negative_spin']
        fit_data_dict_1 = fit_data_dict_1_sign['negative_spin']
        fit_data_dict_2 = fit_data_dict_2_sign['negative_spin']
        B_dict_1 = B_dict_1_sign['negative_spin']
        B_dict_2 = B_dict_2_sign['negative_spin']
    else:
        times = times_dict['positive_spin']
        fit_data_dict_1 = fit_data_dict_1_sign['positive_spin']
        fit_data_dict_2 = fit_data_dict_2_sign['positive_spin']
        B_dict_1 = B_dict_1_sign['positive_spin']
        B_dict_2 = B_dict_2_sign['positive_spin']

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
            modes_available, alpha_coeffs,  beta_coeffs, alpha_beta_functional_form,\
            calibrated, M_tot, dist_mpc, orb_phase, inclination, fit_data_dict_1,\
            fit_data_dict_2, B_dict_1, B_dict_2, fit_func, decomposition_funcs,\
            norm, mode_sum, neg_modes, lmax, CoorbToInert)

    return t_surrogate, h_surrogate
