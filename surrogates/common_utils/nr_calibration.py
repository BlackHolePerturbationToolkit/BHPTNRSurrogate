##==============================================================================
## BHPTNRSurrogate module
## Description : applies NR calibration to raw ppBHPT waveforms
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
## Modified: Katie Rink, Mar 2023 [krink@utexas.edu]
##==============================================================================

import numpy as np

#----------------------------------------------------------------------------------------------------
def alpha_beta_BHPTNRSur1dq1e4(x,a,b,c,d):
    """
    functional form of alpha and beta scaling factors used in BHPTNRSur1dq1e4 model
    """
    
    return 1 + a*x + b*x**2 + c*x**3 + d*x**4

#----------------------------------------------------------------------------------------------------
def alpha_beta_BHPTNRSur2dq1e3(x,a,b,c,d):
    """
    functional form of alpha and beta scaling factors used in BHPTNRSur2dq1e3 model
    """

    return 1 + a*x + b*x**2 + c*x**3 + d*x**4

#----------------------------------------------------------------------------------------------------
def evaluate_alpha(X, l, coefs_alpha, alpha_beta_functional_form):
    """ Implements alpha-beta-scaling to match NR 
        Computes alpha value at a given point in the paprameter space
    """
    
    # find out the max value of \ell upto which modes are nr calibrated
    lmax_nrcalib = max([x[0] for x in coefs_alpha.keys()])
    
    if l<=lmax_nrcalib:
        alpha = alpha_beta_functional_form(X, coefs_alpha[(l,l)][0], coefs_alpha[(l,l)][1], 
                                           coefs_alpha[(l,l)][2], coefs_alpha[(l,l)][3])
    else:
        alpha = 1.0

    return alpha


#----------------------------------------------------------------------------------------------------
def evaluate_beta(X, coefs_beta, alpha_beta_functional_form):
    """ Implements alpha-scaling to match NR 
        Computes beta value at a given point in the paprameter space
    """
    
    beta = alpha_beta_functional_form(X, coefs_beta[0], coefs_beta[1], coefs_beta[2], 
                                      coefs_beta[3])
    return beta

#----------------------------------------------------------------------------------------------------
def alpha_scaling_h(h_raw, alpha):
    """ Implements alpha-beta-scaling to the strain
        Rescales the strain
    """
    
    return np.array(h_raw)*alpha


#----------------------------------------------------------------------------------------------------
def beta_scaling_time(t_raw, beta):
    """ Implements alpha-scaling to match NR 
        Rescales the time axis
    """

    return np.array(t_raw)*beta


#----------------------------------------------------------------------------------------------------
def generate_calibrated_ppBHPT(X_input, raw_time, h_raw_dict, coefs_alpha, coefs_beta, 
                               alpha_beta_functional_form):
    """
    rescales all raw ppBHPT waveform modes to match NR
    
    Inputs
    ======
        
        X_input : array of nr calibration parameterization e.g. [1/q, spin]
        time : array of uncalibrated time on which surrogate has been trained on 
        h_raw_dict : dictiornary of uncalibrated modes
        coeffs_alpha : dictionary of alpha values obtained from calibration mode-by-mode
        coeffs_beta : beta value obtain from calibration - used in time rescaling
        alpha_beta_functional_form : function to use for nr calibration - must come from 
                                     common_utils.nr_calibration.py
    
    Outputs
    =======
    
        t_calib : rescaled time array
        hcal_dict : dictiornary of rescaled modes 
        
    """
    
    hcal_dict = {}
    # scale all modes
    for mode in h_raw_dict.keys():
        (l,m)=mode
        # evaluate alpha
        alpha = evaluate_alpha(X_input, l, coefs_alpha, alpha_beta_functional_form) 
        # scale the strain
        hcal_dict[mode] = alpha_scaling_h(h_raw_dict[mode], alpha) 
        
    # evaluate beta
    beta = evaluate_beta(X_input, coefs_beta, alpha_beta_functional_form)
    # scale the time
    t_calib = beta_scaling_time(raw_time, beta)
    
    return t_calib, hcal_dict
    
