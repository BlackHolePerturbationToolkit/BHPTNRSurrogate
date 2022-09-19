##==============================================================================
## BHPTNRSur1dq1e4 : arXiv:2204.01972
## Description : loads surrogate fits data
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np
import h5py
import os

#----------------------------------------------------------------------------------------------------
def load_surrogate(h5_data_dir):
    
    """ 
    Loads all interpolation data for the following modes
    modes=[(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),(5,4),(5,5),(6,4),(6,5),(6,6),(7,5),(7,6),(7,7),(8,6),(8,7),(8,8),(9,7),(9,8),(9,9),(10,8),(10,9)]
    Assumes the file BHPTNRSur1dq1e4.h5 is located in the same directory as this file.
    """
    
    with h5py.File('%s/BHPTNRSur1dq1e4.h5'%h5_data_dir, 'r') as f:

        modes = [(2,2),(2,1),
                   (3,1),(3,2),(3,3),
                   (4,2),(4,3),(4,4),
                   (5,3),(5,4),(5,5),
                   (6,4),(6,5),(6,6),
                   (7,5),(7,6),(7,7),
                   (8,6),(8,7),(8,8),
                   (9,7),(9,8),(9,9),
                   (10,8),(10,9)]


        # fit info
        h_eim_re_spline_dict = {}
        h_eim_im_spline_dict = {}
        B_re_dict = {}
        B_im_dict = {}
        eim_indicies_im_dict = {}
        eim_indicies_re_dict = {}

        for mode in modes:

            lmode,mmode=mode

            if mode==(2,2):

                eim_indicies_amp_dataset=f['l%s_m%s/eim_indicies'%(lmode,mmode)]
                eim_indicies_amp=eim_indicies_amp_dataset[:]
                eim_indicies_ph_dataset=f['l%s_m%s/eim_indicies_phase'%(lmode,mmode)]
                eim_indicies_ph=eim_indicies_ph_dataset[:]
                B_ph_dataset=f['l%s_m%s/B_phase'%(lmode,mmode)]
                B_ph=np.transpose(B_ph_dataset[:])
                B_amp_dataset=f['l%s_m%s/B'%(lmode,mmode)]
                B_amp=np.transpose(B_amp_dataset[:])
                time_dataset=f['l%s_m%s/times'%(lmode,mmode)]
                time=time_dataset[:]
                degree_dataset=f['l%s_m%s/degree'%(lmode,mmode)]
                degree=degree_dataset[:]
                knots_dataset_amp=f['l%s_m%s/spline_knots_amp'%(lmode,mmode)]
                knots_amp=knots_dataset_amp[:]
                knots_dataset_ph=f['l%s_m%s/spline_knots_ph'%(lmode,mmode)]
                knots_ph=knots_dataset_ph[:]
                h_spline_amp_dataset=f['l%s_m%s/fitparams_amp'%(lmode,mmode)]
                h_spline_amp=h_spline_amp_dataset[:]
                h_spline_ph_dataset=f['l%s_m%s/fitparams_ph'%(lmode,mmode)]
                h_spline_ph=h_spline_ph_dataset[:]

                h_eim_amp_spline=[(knots_amp[flag], h_spline_amp[flag],int(degree)) for flag in range(len(eim_indicies_amp))]
                h_eim_ph_spline=[(knots_ph[flag], h_spline_ph[flag],int(degree)) for flag in range(len(eim_indicies_ph))]

            else:

                eim_indicies_re_dataset=f['l%s_m%s/eim_indicies'%(lmode,mmode)]
                eim_indicies_re_dict[(mode)]=eim_indicies_re_dataset[:]
                eim_indicies_im_dataset=f['l%s_m%s/eim_indicies_im'%(lmode,mmode)]
                eim_indicies_im_dict[(mode)]=eim_indicies_im_dataset[:]
                B_im_dataset=f['l%s_m%s/B_im'%(lmode,mmode)]
                B_im_dict[(mode)]=np.transpose(B_im_dataset[:])
                B_re_dataset=f['l%s_m%s/B'%(lmode,mmode)]
                B_re_dict[(mode)]=np.transpose(B_re_dataset[:])
                time_dataset=f['l%s_m%s/times'%(lmode,mmode)]
                time=time_dataset[:]
                degree_dataset=f['l%s_m%s/degree'%(lmode,mmode)]
                degree=degree_dataset[:]
                knots_dataset_re=f['l%s_m%s/spline_knots_re'%(lmode,mmode)]
                knots_re=knots_dataset_re[:]
                knots_dataset_im=f['l%s_m%s/spline_knots_im'%(lmode,mmode)]
                knots_im=knots_dataset_im[:]
                h_spline_re_dataset=f['l%s_m%s/fitparams_re'%(lmode,mmode)]
                h_spline_re=h_spline_re_dataset[:]
                h_spline_im_dataset=f['l%s_m%s/fitparams_im'%(lmode,mmode)]
                h_spline_im=h_spline_im_dataset[:]

                h_eim_re_spline_dict[(mode)]=[(knots_re[flag], h_spline_re[flag],int(degree)) for flag in range(len(eim_indicies_re_dict[(mode)]))]
                h_eim_im_spline_dict[(mode)]=[(knots_im[flag], h_spline_im[flag],int(degree)) for flag in range(len(eim_indicies_im_dict[(mode)]))]
           
        # nr calibration info
        alpha_coeffs = {}
        for mode in [(2,2),(3,3),(4,4),(5,5)]:
            alpha_coeffs[mode] = f["nr_calib_params/(%d,%d)"%(mode[0],mode[1])]['alpha'][:]
        beta_coeffs = f["nr_calib_params/(2,2)"]['beta'][:]

    return time, \
        eim_indicies_amp, eim_indicies_ph, B_amp, B_ph, h_eim_amp_spline, h_eim_ph_spline, \
        eim_indicies_re_dict, eim_indicies_im_dict, B_re_dict, B_im_dict, h_eim_re_spline_dict, h_eim_im_spline_dict,\
        alpha_coeffs, beta_coeffs
    
