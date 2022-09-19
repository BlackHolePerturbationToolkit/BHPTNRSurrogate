##==============================================================================
## BHPTNRSur1dq1e4 : arXiv:2204.01972
## Description : generates raw ppBHPT surrogate
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np
import scipy
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.interpolate import splrep, splev

#----------------------------------------------------------------------------------------------------
def amp_ph_to_comp(amp,phase):
    """ Takes the amplitude and phase of the waveform and
    computes the compose them together"""
    
    full_wf = amp*np.exp(1j*phase)
    return full_wf

#----------------------------------------------------------------------------------------------------
def re_im_to_comp(re, im, m, orbital_phase):
    """ Takes the real and imaginary part of the waveform and
    combine them to obtain the coorbital frame waveform;
    then transform the wf into the inertial frame"""
    
    full_wf = (re+1j*im)*np.exp(1j*m*np.array(orbital_phase))
    return full_wf

#----------------------------------------------------------------------------------------------------
def log_surrogate_22_mode(q, h_eim_amp_spline, h_eim_ph_spline, eim_indicies_amp, eim_indicies_ph, B_amp, B_ph, calibrated):
    """ Compute the interpolated waveform for a single mode : dominant 22 mode """
    
    h_eim_amp = np.array([splev(np.log10(q), h_eim_amp_spline[j])  for j in range(len(eim_indicies_amp))])
    h_eim_ph = np.array([splev(np.log10(q), h_eim_ph_spline[j]) for j in range(len(eim_indicies_ph))])
    h_approx_amp = np.dot(B_amp.transpose(), h_eim_amp)
    h_approx_ph = np.dot(B_ph.transpose(), h_eim_ph)
    h_approx = amp_ph_to_comp(h_approx_amp, h_approx_ph)
        
    return np.array(h_approx)*(1/q) 


def log_surrogate_hm_mode(q, h_eim_re_spline, h_eim_im_spline, eim_indicies_re, eim_indicies_im, B_re, B_im, 
                   orbital_phase, mode, calibrated):
    """ Compute the interpolated waveform for a single mode : higher order modes"""
    
    (l,m) = mode
    h_eim_re = np.array([splev(np.log10(q), h_eim_re_spline[j])  for j in range(len(eim_indicies_re))])
    h_eim_im = np.array([splev(np.log10(q), h_eim_im_spline[j]) for j in range(len(eim_indicies_im))])
    h_approx_re = np.dot(B_re.transpose(), h_eim_re)
    h_approx_im = np.dot(B_im.transpose(), h_eim_im)
    h_approx = re_im_to_comp(h_approx_re, h_approx_im, m, orbital_phase) 
        
    return np.array(h_approx)*(1/q) 

#----------------------------------------------------------------------------------------------------
def all_modes_surrogate(modes, q_input,
              eim_indicies_amp_dict, eim_indicies_ph_dict, B_amp_dict, B_ph_dict, h_eim_amp_spline_dict, h_eim_ph_spline_dict,
              eim_indicies_re_dict, eim_indicies_im_dict, B_re_dict, B_im_dict, h_eim_re_spline_dict, h_eim_im_spline_dict,
              lmax, calibrated):

    """ Takes the interpolation indices, spline nodes, matrix B and computes the interpolated waveform for all modes"""
    
    h_approx_dict={}
    
    for mode in modes:
        
        # 22 mode is generated from amp/phase surrogate
        if mode==(2,2):
            
            h_approx_dict[(mode)] = log_surrogate_22_mode(q_input, 
                                                     h_eim_amp_spline_dict, h_eim_ph_spline_dict, 
                                                     eim_indicies_amp_dict, eim_indicies_ph_dict, 
                                                     B_amp_dict, B_ph_dict, 
                                                     calibrated)
            
            # compute the orbital phase as it is needed to compute the higher order modes
            orbital_phase = np.unwrap(np.angle(h_approx_dict[mode]))/2
            # needed to match convention of other surrogate models
            h_approx_dict[(mode)] = np.array(np.conj(h_approx_dict[(mode)])) 
             
        # higher order modes are generated using re/im surrogates 
        else:
            (l,m)=mode
            # load modes only upto l=lmax
            if l<=lmax:
                h_approx_dict[(mode)] = log_surrogate_hm_mode(q_input, 
                                                         h_eim_re_spline_dict[(mode)], h_eim_im_spline_dict[(mode)], 
                                                         eim_indicies_re_dict[(mode)], eim_indicies_im_dict[(mode)], 
                                                         B_re_dict[(mode)], B_im_dict[(mode)], 
                                                         orbital_phase, mode,
                                                         calibrated)
                # needed to match convention of other surrogate models
                h_approx_dict[(mode)] = np.array(np.conj(h_approx_dict[(mode)]))
        
    return h_approx_dict


