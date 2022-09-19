##==============================================================================
## BHPTNRSur1dq1e4 : arXiv:2204.01972
## Description : applies NR calibration to raw ppBHPT waveforms
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np

#----------------------------------------------------------------------------------------------------
def alpha_beta_functional_form(x,a,b,c,d):
    """
    functional form of alpha and beta scaling factors
    """
    return 1 + a*x + b*x*x + c*x*x*x + d*x*x*x*x

#----------------------------------------------------------------------------------------------------
def alpha_scaling_h(q, h_raw, l, coefs_alpha):
    """ Implements alpha-beta-scaling to match NR 
        The model is only calibrated up to l=5
        For modes with l>5, we provide 'uncalibrated perturbation waveform'
    """
    if l<=5:
        alpha = alpha_beta_functional_form(1/q, coefs_alpha[(l,l)][0], coefs_alpha[(l,l)][1], coefs_alpha[(l,l)][2], coefs_alpha[(l,l)][3])
    else:
        alpha = 1.0
    h_scaled = np.array(h_raw)*alpha

    return h_scaled


#----------------------------------------------------------------------------------------------------
def beta_scaling_time(q, t_raw, coefs_beta):
    """ Implements alpha-scaling to match NR """
    beta = alpha_beta_functional_form(1/q, coefs_beta[0], coefs_beta[1], coefs_beta[2], coefs_beta[3])
    t_scaled=np.array(t_raw)*beta
    
    return t_scaled

#----------------------------------------------------------------------------------------------------
def generate_calibrated_ppBHPT(q_input, raw_time, h_raw_dict, coefs_alpha, coefs_beta):
    """
    rescales all raw ppBHPT waveform modes to match NR
    """
    hcal_dict = {}
    for mode in h_raw_dict.keys():
        (l,m)=mode
        hcal_dict[mode] = alpha_scaling_h(q_input, h_raw_dict[mode], l, coefs_alpha) 
        
    t_calib = beta_scaling_time(q_input, raw_time, coefs_beta)
    
    return t_calib, hcal_dict
    
