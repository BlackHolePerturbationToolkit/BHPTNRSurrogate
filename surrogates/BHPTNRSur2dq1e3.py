##==============================================================================
## BHPTNRSur2dq1e3 : arXiv:XXXX.YYYYY
## Description : generates calibrated ppBHPT surrogate model BHPTNRSur2dq1e3
## Author : Katie Rink, Dec 2022 [krink@utexas.edu]
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
times_sign, fit_data_dict_1_sign, fit_data_dict_2_sign, B_dict_1_sign, B_dict_2_sign, \
    alpha_coeffs, beta_coeffs = load.load_BHPTNRSur2dq1e3_surrogate(h5_data_dir)

print("SURROGATE LOADED")
#---------------------------------------------------------------------------------------------------- 
# add docstring from utility
@docs.copy_doc(docs.generic_doc_for_models,docs.BHPTNRSur2dq1e3_doc)
def generate_surrogate(q, spin1=0.0, spin2=None, ecc=None, ano=None, modes=None, M_tot=None, dist_mpc=None, orb_phase=None, inclination=None, neg_modes=False, mode_sum=False, lmax=4, calibrated=False):

    """ 
    Description : Top-level function to generate surrogate waveform in either geometric or physical units
    
    Input
    =====
    q:      mass ratio
    
    spin1: spin on primary object, ranging from [-0.6, 0.6] for spinning 2D surrogate.
           Default: 0.0 (uses poisitve spin data)
    
    modes:  list of modes
            Default (None) corresponds to all available modes in the model

            [(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),
             (5,3),(5,4),(5,5),(6,4),(6,5),(6,6),(7,5),(7,6),
             (7,7),(8,6),(8,7),(8,8),(9,7),(9,8),(9,9),(10,8),(10,9)]
            
    M_tot:     total mass of the binary in solar masses
               Default: None (in which case a geometric waveform is returned)
    
    dist_mpc:  distance of the binary from the observer in Mpc
               Default: None (in which case geometric wf is returned)
               
    orb_phase: orbital phase
    
    inclination: inclination angle
    
    lmax:  5 (default)
           Modes upto l=5 are NR calibrated. Modes l>6 are uncalibrated
           Note: If one provides a list of modes, modes beyond lmax 
                 will not be returned
            
    mode_sum:  If true, all modes are summed. If false all modes are returned 
               in a dictionary. Default: false
               Note: Only works when orb_phase and inclination are not None
               
    calibrated:  Whether you want NR-calibrated waveform or not
                 When set to True, it applies a scaling to the uncalibrated
                 surrogate waveform. This scaling has been obtained by calibrating
                 the ppBHPT waveforms to NR in comparable mass ratio
                 regime (3<=q<=9). If set to False, the raw (uncalibrated)
                 ppBHPT waveforms are returned.
                 Default: True

                 
    Output
    ======
    t : time
    h : waveform modes as a dictionary


    Example Uses:
    =============
    1. to obtain uncalibrated (raw) geometric waveform
            t, h = generate_surrogate(q=8, modes=[(2,1),(2,2),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),(5,4),(5,5)], calibrated=False)
    2. to obtain NR Calibrated geometric waveform
            t, h = generate_surrogate(q=8, modes=mode_list)       
    3. to obtain NR calibrated physical waveform
            t, h = generate_surrogate(q=q, modes=mode_list, M_tot=50, dist_mpc=100)
    4. to obtain NR calibrated physical waveform on a sphere
            t, h = generate_surrogate(q=q, modes=mode_list, M_tot=50, dist_mpc=100, orb_phase=np.pi/3, inclination=np.pi/4)
    5. to obtain NR calibrated physical waveform on a sphere for modes up to l=5
            t, h = generate_surrogate(q=q, modes=mode_list, M_tot=50, dist_mpc=100, orb_phase=np.pi/3, inclination=np.pi/4, lmax=5)
    6. to obtain mode-summed NR calibrated physical waveform on a sphere
            t, h = generate_surrogate(q=8, M_tot=60, dist_mpc=100, orb_phase=np.pi/3, inclination=np.pi/4, lmax=3, mode_sum=True)
         
            
    """

    modes_available = [(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4)]

    if modes==None:
        modes = modes_available
        
    if spin1 < 0.0:
        times = times_sign['negative_spin']
        fit_data_dict_1 = fit_data_dict_1_sign['negative_spin']
        fit_data_dict_2 = fit_data_dict_2_sign['negative_spin']
        B_dict_1 = B_dict_1_sign['negative_spin']
        B_dict_2 = B_dict_2_sign['negative_spin']
    else:
        times = times_sign['positive_spin']
        fit_data_dict_1 = fit_data_dict_1_sign['positive_spin']
        fit_data_dict_2 = fit_data_dict_2_sign['positive_spin']
        B_dict_1 = B_dict_1_sign['positive_spin']
        B_dict_2 = B_dict_2_sign['positive_spin']
        
    # define the parameterization for surrogate
    X_sur = [np.log10(q), spin1]
    
    # define parameterization for nr calibration
    X_calib = 1/q
    
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
