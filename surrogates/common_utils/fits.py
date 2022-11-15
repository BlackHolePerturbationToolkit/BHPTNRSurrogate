##==============================================================================
## BHPTNRSur1dq1e4 : arXiv:2204.01972
## Description : utility file for surrogate model
## Author : Tousif Islam, Nov 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np
import scipy
from scipy.interpolate import splrep, splev
from . import utils

#----------------------------------------------------------------------------------------------------
def evaluate_GPR_at_EIM_nodes(X, fit_data):
    """ Evaluate the spline at one EIM node """
    return None

#----------------------------------------------------------------------------------------------------
def evaluate_splines_at_EIM_nodes(X, fit_data):
    """ Evaluate the spline at one EIM node """
    
    [h_eim_spline, eim_indicies] = fit_data
    return np.array([splev(X, h_eim_spline[j]) for j in range(len(eim_indicies))])


#----------------------------------------------------------------------------------------------------
def EIM_B_to__waveform_datapiece(B, eim_vals):
    """ Compute the interpolated waveform for a single mode : dominant 22 mode """
    
    approx_datapiece = np.dot(B.transpose(), eim_vals)
    return approx_datapiece


#----------------------------------------------------------------------------------------------------
def evaluate_datapiece(X, fit_data, B, fit_func):
    """ Compute the datapiece for the input parameters"""
    
    # evaluates spline fits at eim nodes
    if fit_func == 'spline_1d':
        h_eim_datapiece = evaluate_splines_at_EIM_nodes(X, fit_data)
    # evaluates GPR fits at eim nodes
    elif fit_func == 'GPR_fits':
        h_eim_datapiece = evaluate_GPR_at_EIM_nodes(X, fit_data)
    # combine h_eim and  eim basis matrix to give full datapiece
    h_approx_datapiece = EIM_B_to__waveform_datapiece(B, h_eim_datapiece) 
    
    return h_approx_datapiece


#----------------------------------------------------------------------------------------------------
def evaluate_surrogate_mode(X, fit_data_1, fit_data_2, B_datapiece_1, B_datapiece_2, 
                            fit_func, norm, mode=None, orbital_phase=None):
    """ Compute the interpolated waveform for a single mode : higher order modes"""
    
    # evaluate first datapiece e.g amplitude / real part of wf
    h_approx_datapiece_1 = evaluate_datapiece(X,  fit_data_1, B_datapiece_1, fit_func)
    # evaluate second datapiece e.g phase / imag part of wf
    h_approx_datapiece_2 = evaluate_datapiece(X,  fit_data_2, B_datapiece_2, fit_func)
    
    # if orbital phase is given, we combine amp/phase to full wf
    if orbital_phase is None:
        h_approx =  utils.amp_ph_to_comp(h_approx_datapiece_1, h_approx_datapiece_2) 
    # combine real/imag to full wf, requires orbital phase information
    else:
        (l,m) = mode
        h_approx =  utils.re_im_to_comp(h_approx_datapiece_1, h_approx_datapiece_2,  m, 
                                        orbital_phase) 
        
    # needed to match convention of other surrogate models
    h_approx = np.conj(np.array(h_approx))*norm
    
    return h_approx
