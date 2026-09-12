##==============================================================================
## BHPTNRSurrogate module
## Description : loads spline surrogate fits data
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np
import h5py

#----------------------------------------------------------------------------------------------------
def read_amplitude_fits(f, lmode, mmode):
    """
    Read fit info for the amplitude of 22 mode
    """
    # EIM nodes
    eim_indicies_amp=f['l%s_m%s/eim_indices'%(lmode,mmode)][:]
    # B matrix
    B_amp=np.transpose(f['l%s_m%s/B'%(lmode,mmode)][:])
    # spline fit info : degree of the spline
    degree=f['l%s_m%s/degree'%(lmode,mmode)][:]
    # spline fit info : knots used in spline fit
    knots_amp=f['l%s_m%s/spline_knots_amp'%(lmode,mmode)][:]
    # spline fit info : values at the knots
    h_spline_amp=f['l%s_m%s/fitparams_amp'%(lmode,mmode)][:]
    # combine spline info
    h_eim_amp_spline=[(knots_amp[flag], h_spline_amp[flag],int(degree.item())) for flag
                      in range(len(eim_indicies_amp))]
    
    return eim_indicies_amp, B_amp, h_eim_amp_spline
 
#----------------------------------------------------------------------------------------------------
def read_phase_fits(f, lmode, mmode):
    """
    Read fit info for the phase of 22 mode
    """

    # EIM nodes
    eim_indicies_ph=f['l%s_m%s/eim_indices_phase'%(lmode,mmode)][:]
    # B matrix
    B_ph=np.transpose(f['l%s_m%s/B_phase'%(lmode,mmode)][:])
    # spline fit info : degree of the spline
    degree=f['l%s_m%s/degree'%(lmode,mmode)][:]
    # spline fit info : knots used in spline fit
    knots_ph=f['l%s_m%s/spline_knots_phase'%(lmode,mmode)][:]
    # spline fit info : values at the knots
    h_spline_ph=f['l%s_m%s/fitparams_phase'%(lmode,mmode)][:]
    # combine spline info
    h_eim_ph_spline=[(knots_ph[flag], h_spline_ph[flag],int(degree.item())) for flag
                     in range(len(eim_indicies_ph))]
    
    return eim_indicies_ph, B_ph, h_eim_ph_spline

#----------------------------------------------------------------------------------------------------
def read_real_part_fits(f, lmode, mmode):
    """
    Read fit info for the real part of higher modes (in coorbital phase)
    """

    # EIM nodes
    eim_indicies_re=f['l%s_m%s/eim_indices'%(lmode,mmode)][:]
    # B matrix
    B_re=np.transpose(f['l%s_m%s/B'%(lmode,mmode)][:])
    # spline fit info : degree of the spline
    degree=f['l%s_m%s/degree'%(lmode,mmode)][:]
    # spline fit info : knots used in spline fit
    knots_re=f['l%s_m%s/spline_knots_re'%(lmode,mmode)][:]
    # spline fit info : values at the knots
    h_spline_re=f['l%s_m%s/fitparams_re'%(lmode,mmode)][:]
    # combine spline info
    h_eim_re_spline=[(knots_re[flag], h_spline_re[flag],int(degree.item())) for
                                  flag in range(len(eim_indicies_re))]
    
    return eim_indicies_re, B_re, h_eim_re_spline


#----------------------------------------------------------------------------------------------------
def read_imag_part_fits(f, lmode, mmode):
    """
    Read fit info for the real part of higher modes (in coorbital phase)
    """

    # EIM nodes
    eim_indicies_im=f['l%s_m%s/eim_indices_im'%(lmode,mmode)][:]
    # B matrix
    B_im=np.transpose(f['l%s_m%s/B_im'%(lmode,mmode)][:])
    # spline fit info : degree of the spline
    degree=f['l%s_m%s/degree'%(lmode,mmode)][:]
    # spline fit info : knots used in spline fit
    knots_im=f['l%s_m%s/spline_knots_im'%(lmode,mmode)][:]
    # spline fit info : values at the knots
    h_spline_im=f['l%s_m%s/fitparams_im'%(lmode,mmode)][:]
    # combine spline info
    h_eim_im_spline=[(knots_im[flag], h_spline_im[flag],int(degree.item())) for
                                  flag in range(len(eim_indicies_im))]
    
    return eim_indicies_im, B_im, h_eim_im_spline


#----------------------------------------------------------------------------------------------------
def read_times(f, lmode, mmode):
    """
    Read fit info / time
    """
    # time
    time=f['l%s_m%s/times'%(lmode,mmode)][:]
    
    return time

#----------------------------------------------------------------------------------------------------
def read_nrcalib_info(f, nrcalib_modes):
    """
    Read alpha/beta nr calibration information
    """
    # read the coefficients for alpha
    alpha_coeffs = {}
    for mode in nrcalib_modes:
        alpha_coeffs[mode] = f["nr_calib_params/(%d,%d)"%(mode[0],mode[1])]['alpha'][:]
    # read the coefficients for beta
    beta_coeffs = f["nr_calib_params/(2,2)"]['beta'][:]
    
    return alpha_coeffs, beta_coeffs

#----------------------------------------------------------------------------------------------------
def load_surrogate(h5_data_dir, h5File, wf_modes, nrcalib_modes):
    
    """ 
    Loads all spline interpolation data for a given set of modes
    """

    # read all fit info
    with h5py.File('%s/%s'%(h5_data_dir,h5File), 'r') as f:

        # we will save the fit info for two datapieces
        # in dictionary for all modes
        B_dict_1, B_dict_2  = {}, {}
        fit_data_dict_1, fit_data_dict_2 = {}, {}
        
        # loop over all modes
        for mode in wf_modes:

            lmode,mmode=mode
            # read time data
            time = read_times(f, lmode, mmode)
            
            # special treatment for 22 mode; we model the mode with amp/phase decompositon
            # so we will use arrays to store the fit info
            if mode==(2,2):
                # read amplitude fits
                eim_indicies_1, B_dict_1[(mode)], h_eim_spline_1 \
                                                                = read_amplitude_fits(f, lmode, mmode)
                # read phase fits
                eim_indicies_2, B_dict_2[(mode)], h_eim_spline_2 \
                                                                = read_phase_fits(f, lmode, mmode)

            # read other modes
            # these modes have been modelled with real/imag decompositon
            # we will use dictionary to save the fit info for all modes
            else:
                # read real part fits
                eim_indicies_1, B_dict_1[(mode)], h_eim_spline_1 \
                                                                = read_real_part_fits(f, lmode, mmode)
                # read imag part fits
                eim_indicies_2, B_dict_2[(mode)], h_eim_spline_2 \
                                                                = read_imag_part_fits(f, lmode, mmode)
            
            # now combine eim_spline and eim_indicies values to get the full fit data for 
            # all modes - this cleans up the code significantly
            fit_data_dict_1[mode] = [h_eim_spline_1, eim_indicies_1]
            fit_data_dict_2[mode] = [h_eim_spline_2, eim_indicies_2]
                
        # nr calibration info
        alpha_coeffs, beta_coeffs = read_nrcalib_info(f, nrcalib_modes)

    return time, fit_data_dict_1, fit_data_dict_2, B_dict_1, B_dict_2, alpha_coeffs, beta_coeffs

