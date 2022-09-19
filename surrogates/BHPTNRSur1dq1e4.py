##==============================================================================
## BHPTNRSur1dq1e4 : arXiv:2204.01972
## Description : generates calibrated ppBHPT surrogate
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.interpolate import splrep, splev
import h5py
import hashlib
from gwtools import gwtools as _gwtools
import os
from os import path
import subprocess

import surrogate_utils.BHPTNRSur1dq1e4.load_surrogate_fits as load
import surrogate_utils.BHPTNRSur1dq1e4.generate_raw_surrogate as raw_sur
import surrogate_utils.BHPTNRSur1dq1e4.nr_calibration as nrcal
import surrogate_utils.BHPTNRSur1dq1e4.utils as utils


# h5 data directory
h5_data_dir = os.path.dirname(os.path.abspath(__file__)) + '/../data'
# load all fits data
time, eim_indicies_amp, eim_indicies_ph, B_amp, B_ph, h_eim_amp_spline, h_eim_ph_spline, \
eim_indicies_re_dict, eim_indicies_im_dict, B_re_dict, B_im_dict, h_eim_re_spline_dict, h_eim_im_spline_dict, alpha_coeffs, beta_coeffs = load.load_surrogate(h5_data_dir)

#---------------------------------------------------------------------------------------------------- 
def generate_surrogate(q_input, modes=None, M_tot=None, dist_mpc=None, orb_phase=None, inclination=None, neg_modes=True, mode_sum=False, lmax=5, calibrated=True):
    """ 
    Description : Top-level function to generate surrogate waveform in either geometric or physical units
    
    Inputs
    ====================
    q_input : mass ratio
    
    modes : list of modes
            default is all available modes in the model i.e. [(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),(5,4),(5,5),(6,4),(6,5),(6,6),(7,5),(7,6),(7,7),(8,6),(8,7),(8,8),(9,7),(9,8),(9,9),(10,8),(10,9)]
            
    M_tot : total mass of the binary in solar unit
              default: None (in which case geometric wf is returned)
    
    dist_mpc : distance of the binary from the observer in Mpc
               default: None (in which case geometric wf is returned)
               
    orb_phase : orbital phase
    
    inclination : inclination angle
    
    lmax : 5 (deafult)
           Modes upto l=5 are NR calibrated
           Note : If one provides a list of modes, modes beyond lmax will not be returned
            
    mode_sum : False
               Only works when orb_phase and inclination are not None
               
    calibrated : tell whether you want NR calibrated waveform or not
                 When set to True, it applies a scaling to the raw surrogate waveform 
                 This scaling has been obtained by calibrating the ppBHPT waveforms to NR in comparable mass ratio regime (3<=q<=9)
                 If set to False, the raw (uncalibrated)  ppBHPT waveforms are returned.
                 default: True
                 
    Output
    ====================
    t : time
    h : waveform modes
                 
    Example Uses:
    ====================
    1. to obtain raw geometric waveform
            t, h = generate_surrogate(q_input=8, modes=[(2,1),(2,2),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),(5,4),(5,5)], calibrated=False)
    2. to obtain NR Calibrated geometric waveform
            t, h = generate_surrogate(q_input=8, modes=mode_list)       
    3. to obtain NR calibrated physical waveform
            t, h = generate_surrogate(q_input=q, modes=mode_list, M_tot=50, dist_mpc=100)
    4. to obtain NR calibrated physical waveform on a sphere
            t, h = generate_surrogate(q_input=q, modes=mode_list, M_tot=50, dist_mpc=100, orb_phase=np.pi/3, inclination=np.pi/4)
    5. to obtain NR calibrated physical waveform on a sphere for modes up to l=5
            t, h = generate_surrogate(q_input=q, modes=mode_list, M_tot=50, dist_mpc=100, orb_phase=np.pi/3, inclination=np.pi/4, lmax=5)
    6. to obtain mode-summed NR calibrated physical waveform on a sphere
            t, h = generate_surrogate(q_input=8, M_tot=60, dist_mpc=100, orb_phase=np.pi/3, inclination=np.pi/4, lmax=3, mode_sum=True)
         
            
    """
    
    if modes==None:
        modes=[(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),(5,4),(5,5),(6,4),(6,5),(6,6),(7,5),(7,6),(7,7),(8,6),(8,7),(8,8),(9,7),(9,8),(9,9),(10,8),(10,9)]
        
    # raw geometric waveforms
    hsur_raw_dict = raw_sur.all_modes_surrogate(modes, q_input,
              eim_indicies_amp, eim_indicies_ph, B_amp, B_ph, h_eim_amp_spline, h_eim_ph_spline,
              eim_indicies_re_dict, eim_indicies_im_dict, B_re_dict, B_im_dict, h_eim_re_spline_dict, h_eim_im_spline_dict,
              lmax, calibrated)
    
    if calibrated==True:
        t_sur, hsur_dict = nrcal.generate_calibrated_ppBHPT(q_input, time, hsur_raw_dict, alpha_coeffs, beta_coeffs)
        if lmax>5:
            print('WARNING : only modes up to \ell=5 are NR calibrated')
    else:
        t_sur=np.array(time)
        hsur_dict = hsur_raw_dict
        print('WARNING : modes are NOT NR calibrated - waveforms only have 0PA contribution')
    
    # get all the negative m modes using symmetry
    if neg_modes:
        hsur_dict = utils.generate_negative_m_mode(hsur_dict)
            
        
    # relevant for obtaining physical waveforms
    if M_tot is not None and dist_mpc is not None:
        
        t_sur, hsur_dict = utils.geo_to_SI(t_sur, hsur_dict, M_tot, dist_mpc)
        # evaluate on the sphere
        if orb_phase is not None and inclination is not None:
            hsur_dict = utils.evaluate_on_sphere(inclination, orb_phase, hsur_dict)
                
        # add checks
        elif orb_phase is not None and inclination is None:
                raise ValueError("Both orb_phase and inclination should be None! Or both should have physical values to generate physical waveform")    
        elif orb_phase is None and inclination is not None:
                raise ValueError("Both orb_phase and inclination should be None! Or both should have physical values to generate physical waveform")
                         
    # add checks    
    elif M_tot is not None and dist_mpc is None:
        raise ValueError("Both M_tot and dist_mpc should be None! Or both should have physical values to generate physical waveform")
    elif M_tot is None and dist_mpc is not None:
        raise ValueError("Both M_tot and dist_mpc should be None! Or both should have physical values to generate physical waveform")
    
    
    # whether to output mode contribution summed up
    if mode_sum==False:
        return t_sur, hsur_dict
    else:
        if M_tot is not None and dist_mpc is not None and orb_phase is not None and inclination is not None:
            h_summed = utils.sum_modes(hsur_dict)
            return t_sur, h_summed
        else:
            raise ValueError("M_tot, dist_mpc, orb_phase and inclination should NOT be None")
            