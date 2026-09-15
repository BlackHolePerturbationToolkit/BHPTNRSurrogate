##==============================================================================
## BHPTNRSurrogate module
## Description : utility file for surrogate model
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import warnings

import numpy as np
from gwtools import gwtools as _gwtools
from gwtools.harmonics import sYlm as _sYlm
from . import nr_calibration as nrcalib
from gwtools.gwtools import geo_to_SI

#----------------------------------------------------------------------------------------------------
def chars_to_string(chars):
    """ Function converts array of  unicode characters to string.
            - Copied from gwsurrogate: surrogatIO.py
            - Needed for reading in hdf5 file data.
    """
    return "".join(chr(cc) for cc in chars)

#----------------------------------------------------------------------------------------------------
def amp_ph_to_comp(amp,phase):
    """ Takes the amplitude and phase of the waveform and
    computes the compose them together"""
    
    full_wf = amp*np.exp(1j*phase)
    return full_wf

#----------------------------------------------------------------------------------------------------
def re_im_to_comp(re, im):
    """ Takes the real and imaginary part of the waveform and
    combine them to obtain the coorbital frame waveform"""
    
    full_wf = (re+1j*im)
    return full_wf

#----------------------------------------------------------------------------------------------------
def coorbital_to_inertial(h_coorb):
    """ Transform the coorbital frame wf into the inertial frame"""
    
    h_inertial = {}
    # 22 mode is in inertial frame and HMs are in coorbital phase
    for mode in h_coorb.keys():
        (l,m)=mode
        if mode==(2,2):
            # compute orbital phase
            orbital_phase = np.unwrap(np.angle(h_coorb[mode]))/2
            h_inertial[mode] = h_coorb[mode]
        else:
            # transform HMs to inertial frame
            h_inertial[mode] = h_coorb[mode]*np.exp(1j*m*np.array(orbital_phase))
            
    return h_inertial

#---------------------------------------------------------------------------------------------------- 
def phase_rotation(h, delta_orb_phase):
    """
    performs an orbital phase rotation
    """
    h_rotated = {}
    
    for mode in h.keys():
        (l,m) = mode
        # calculate the corresponding phase rotation in respective modes
        phase_rot = m*delta_orb_phase
        # apply phase rotation
        h_rotated[mode] = h[mode] * np.exp(1j*phase_rot)
        
    return h_rotated

#---------------------------------------------------------------------------------------------------- 
def evaluate_on_sphere(theta, phi, h_dict):
    """evaluate on the sphere"""

    if theta is not None:
        if phi is None: raise ValueError('phi must have a value')
            
        hdict_sphere = {}
        for mode in h_dict.keys():
            (ell,m)=mode
            # compute spherical harmonics 
            sYlm_value =  _sYlm(-2,ll=ell,mm=m,theta=theta,phi=phi)
            # compute modes
            hdict_sphere[mode] = sYlm_value*h_dict[mode]
            
    return hdict_sphere

#---------------------------------------------------------------------------------------------------- 
def sum_modes(h_dict):
    """sum all the modes on a point in the sky"""
    
    h = np.zeros(len(h_dict[(2,2)]))
    for mode in h_dict.keys():
        h = h + h_dict[mode]
    return h


#---------------------------------------------------------------------------------------------------- 
def generate_negative_m_mode(h_dict):
    """ 
    For m>0 positive modes hp_mode,hc_mode use h(l,-m) = (-1)^l h(l,m)^* to compute the m<0 mode.
    See Eq. 78 of Kidder,Physical Review D 77, 044016 (2008), arXiv:0710.0614v1 [gr-qc].
    """

    h_dict_all_modes = {}
    
    for mode in h_dict.keys():
        (l,m)=mode
        # obtain postive m mode values
        h_dict_all_modes[(l,m)] = h_dict[mode]
        # sanity checks
        if (m==0):
            raise ValueError('m must be nonnegative. m<0 will be generated for you from the m>0 mode.')
        elif (m<0):
            raise ValueError('m must be nonnegative. m<0 will be generated for you from the m>0 mode.')
        # calculate negative m mode
        else:
            h_dict_all_modes[(l,-m)] = np.power(-1,l) * np.conjugate(h_dict[mode])

    return h_dict_all_modes


#----------------------------------------------------------------------------------------------------
def obtain_processed_output(X_calib, time, hsur_raw_dict, alpha_coeffs, beta_coeffs,
                            alpha_beta_functional_form, calibrated, M_tot, dist_mpc,
                            orb_phase, inclination, mode_sum, neg_modes, lmax, CoorbToInert=False,
                            mass_factor=1.0):
    """
    Function to process the output of raw surrogate to apply :
    (i) NR calibration;
    (ii) Conversion to SI units from geometric units;
    (iii) mode summation;
    
    Inputs
    ======
        
        X_calib :--: array of nr calibration parameterization e.g. [1/q, spin]
        time :--: array of time on which surrogate has been trained on - read from h5 file
        hsur_raw_dict :--: raw surrogate waveform dictionary
        alpha_coeffs :--: dictionary of alpha values obtained from calibration mode-by-mode
        beta_coeffs :--: beta value obtain from calibration - used in time rescaling
        alpha_beta_functional_form :--: function to use for nr calibration - must come from 
                                        common_utils.nr_calibration.py
        calibrated :--: True/False - whether requested waveforms should be nr calibrated or not
        M_tot :--: total mass of the binary
        dist_mpc :--: distance in mpc for the binary
        orb_phase :--: orbital phase at the start of the waveform
        inclination :--: inclination angle wrt the observer
        mode_sum :--: indicate whether modes should be summed up.
        neg_modes :--: indicate whether negative modes should be retured using orbital plane symmetry.
        lmax :--: maximum value of l upto which modes should be returned.
        CoorbToInert :--: indicate whether higher modes have been modelled in coorbital frame. In that
                          case, additional processing will be performed.
    
    Outputs
    =======
    
        t_surrogate : time array 
        h_surrogate : dictiornary of modes
    
     
    """
    
    # transform higher modes from coorbital to inertial frame if asked
    if CoorbToInert==True:
        hsur_raw_dict = coorbital_to_inertial(hsur_raw_dict)
        
    # when nr calibration is applied
    if calibrated==True:
        t_sur, hsur_dict = nrcalib.generate_calibrated_ppBHPT(X_calib, time, hsur_raw_dict, alpha_coeffs, 
                                                                  beta_coeffs, alpha_beta_functional_form)
        if lmax>5:
            warnings.warn('Only modes up to ell=5 are NR calibrated', stacklevel=2)
    # when no nr calibration is applied
    else:
        t_sur = np.array(time) * mass_factor
        hsur_dict = {mode: np.array(h) * mass_factor for mode, h in hsur_raw_dict.items()}
        warnings.warn('Modes are NOT NR calibrated - waveforms only have 0PA contribution', stacklevel=2)

    # get all the negative m modes from postive m modes using symmetry
    if neg_modes:
        hsur_dict = generate_negative_m_mode(hsur_dict)

    # relevant for obtaining physical waveforms
    if M_tot is not None and dist_mpc is not None:
        t_sur, hsur_dict = geo_to_SI(t_sur, hsur_dict, M_tot, dist_mpc)
        # evaluate on the sphere
        if orb_phase is not None and inclination is not None:
            hsur_dict = evaluate_on_sphere(inclination, orb_phase, hsur_dict)

    # sum up the modes if it is asked
    if mode_sum==False:
        return t_sur, hsur_dict
    else:
        if M_tot is not None and dist_mpc is not None and orb_phase is not None and inclination is not None:
            h_summed = sum_modes(hsur_dict)
            return t_sur, h_summed
        
        