##==============================================================================
## BHPTNRSur1dq1e4 : arXiv:2204.01972
## Description : generates raw ppBHPT surrogate
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np
import scipy
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.interpolate import splrep, splev
from common_utils import fits

#----------------------------------------------------------------------------------------------------
def all_modes_surrogate(modes, X_input, eim_indicies_amp_dict, eim_indicies_ph_dict, 
                       B_amp_dict, B_ph_dict, h_eim_amp_spline_dict, h_eim_ph_spline_dict,
                       eim_indicies_re_dict, eim_indicies_im_dict, B_re_dict, B_im_dict, 
                        h_eim_re_spline_dict, h_eim_im_spline_dict, lmax, fit_func, norm):

    """ Takes the interpolation indices, spline nodes, matrix B and computes the 
        interpolated waveform for all modes """
    
    # evaluate all the modes
    h_approx_dict={}
            
    for mode in modes:
        
        # 22 mode is generated from amp/phase surrogate
        if mode==(2,2):
            
            # have the fit data for 1d splines
            fit_data_1 = [h_eim_amp_spline_dict, eim_indicies_amp_dict]
            fit_data_2 = [h_eim_ph_spline_dict, eim_indicies_ph_dict]
            # evaluate surrogate
            h_approx_dict[(mode)] = fits.evaluate_surrogate_mode(X_input, fit_data_1, fit_data_2,
                                                       B_amp_dict, B_ph_dict, fit_func, norm)
            # compute the orbital phase as it is needed to compute the higher order modes
            orbital_phase = np.unwrap(np.angle(h_approx_dict[mode]))/2
             
        # higher order modes are generated using re/im surrogates 
        else:
            (l,m)=mode
            # load modes only upto l=lmax
            if l<=lmax:
                
                # have the fit data for 1d splines
                fit_data_1 = [h_eim_re_spline_dict[(mode)], eim_indicies_re_dict[(mode)]]
                fit_data_2 = [h_eim_im_spline_dict[(mode)], eim_indicies_im_dict[(mode)]]
                # evaluate surrogate
                h_approx_dict[(mode)] = fits.evaluate_surrogate_mode(X_input, fit_data_1, fit_data_2, 
                                                         B_re_dict[(mode)], B_im_dict[(mode)], 
                                                         fit_func, norm, mode, orbital_phase)
                
    return h_approx_dict


