##==============================================================================
## BHPTNRSurrogate module
## Description : evaluates the (spline or GPR) surrogate model
## Author : Tousif Islam, Nov 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np
import scipy
from scipy.interpolate import splrep, splev
from . import utils
from .eval_pysur import evaluate_fit as evaluate_GPR

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
                            fit_func, decomposition_func, norm):
    """ Compute the interpolated waveform for a single mode : higher order modes"""
    
    # evaluate first datapiece e.g amplitude / real part of wf
    h_approx_datapiece_1 = evaluate_datapiece(X,  fit_data_1, B_datapiece_1, fit_func)
    # evaluate second datapiece e.g phase / imag part of wf
    h_approx_datapiece_2 = evaluate_datapiece(X,  fit_data_2, B_datapiece_2, fit_func)
    
    # combine datapieces to obtain full wf either in the inertial frame or in the
    # coorbital frame; at this stage, the waveforms are rertured in their respective
    # frames where models have been built e.g. inertial for 22 or coorbital for HMs
    # in case of BHPTNRSur1dq1e4
    h_approx =  decomposition_func(h_approx_datapiece_1, h_approx_datapiece_2)
    
    # needed to match convention of other surrogate models
    # multiply surrogate amplitude with overall normalization factor
    h_approx = np.conj(np.array(h_approx))*norm
    
    return h_approx


#----------------------------------------------------------------------------------------------------
def all_modes_surrogate(modes, X_input, fit_data_dict_1, fit_data_dict_2, \
                        B_dict_1, B_dict_2, lmax, fit_func, decomposition_funcs, norm):

    """ Takes the fit data (either from splines or GPR), matrix B and computes the 
        interpolated waveform for all modes 
        
    Inputs
    ======
    
        modes : list of modes to evaluate
        
        X_input :  array of surrogate parameterization e.g. [log(q), spin1, spin2]

        fit_data_dict_1, fit_data_dict_2 : dictionary of fit data obtained for two datapieces from 
                                           the h5 file.
                                           Keys are the modes.
                                           Structure may depend on ether the data comes from spline 
                                           fits or GPR fits. However, they should always be packed 
                                           in fit_data_dict_1 and fit_data_dict_2. 
                                           Make sure to modify your data loading script to achieve
                                           this if necessary.

        B_dict_1, B_dict_2 : dictionary of the basis matrices obtained from h5 file.
                             Modes used as keys.

        fit_func : form of fitting function. options : 'spline_1d' or 'GPR_fits'

        decomposition_funcs : form of data decomposition function to combine datapieces for 22 and
                              higher modes respectively. e.g. Amp/Phase to full or real/imag to full
                              etc. These functions are available at common_utils.utils.py

        norm : overall normalization factor to be multiplied to final waveform. This depends on the 
              way the surrogate have been constructed. Mostly norm=1/q or norm=1. 
    
    Outputs
    =======
    
        t_surrogate : time array
        h_surrogate : dictiornary of modes  
        
    """
    
    # dictionary to save waveform
    h_approx_dict={}
    # evaluate all the modes     
    for mode in modes:
        (l,m) = mode
        # load modes only upto l=lmax
        if l<=lmax:
            # get the fit data for specific mode for both the datapieces
            fit_data_1 = fit_data_dict_1[mode]
            fit_data_2 = fit_data_dict_2[mode]
            
            # read the decompositon function for the modes; special treatment for the
            # 22 mode and higher order modes
            if mode==(2,2):
                decomposition_func = decomposition_funcs[0]
            else:
                decomposition_func = decomposition_funcs[1]
    
            # evaluate surrogate modes
            h_approx_dict[(mode)] = evaluate_surrogate_mode(X_input, fit_data_1, fit_data_2, 
                                                            B_dict_1[(mode)], B_dict_2[(mode)], 
                                                            fit_func, decomposition_func, norm)
                
    return h_approx_dict
