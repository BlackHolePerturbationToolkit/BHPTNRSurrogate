##==============================================================================
## BHPTNRSur1dq1e4 : arXiv:2204.01972
## Description : applies NR calibration to raw ppBHPT waveforms
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np

#----------------------------------------------------------------------------------------------------
def alpha_beta_BHPTNRSur1dq1e4(x,a,b,c,d):
    """
    functional form of alpha and beta scaling factors used in BHPTNRSur1dq1e4 model
    """
    
    return 1 + a*x + b*x**2 + c*x**3 + d*x**4

#----------------------------------------------------------------------------------------------------
def evaluate_alpha(X, l, coefs_alpha, alpha_beta_functional_form):
    """ Implements alpha-beta-scaling to match NR 
        Computes alpha value at a given point in the paprameter space
    """
    
    if l<=5:
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
        The model is only calibrated up to l=5
        For modes with l>5, we provide 'uncalibrated perturbation waveform'
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
    