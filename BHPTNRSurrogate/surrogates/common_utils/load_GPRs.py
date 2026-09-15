##==============================================================================
## BHPTNRSurrogate module
## Description : loads GPR surrogate fits data
## Author: Katie Rink, Mar 2023 [krink@utexas.edu]
## Modified: Tousif Islam, July 2024
##
##
## NOTE: despite this file being in the common_utils directory, it is specific 
## to the BHPTNRSur2dq1e3 model and its associated h5 data file.
##
##==============================================================================

import numpy as np
import h5py
import copy
from . import load_splines
from . import utils

#----------------------------------------------------------------------------------------------------
def read_times(file):
    """
    Read fit info / time keeping in mind that there are two sub-surrogates
    for the negative spin and positive spin cases
    file: h5file opened in h5py with read mode
    """
    times = {}
    for spin_sign in ['negative_spin', 'positive_spin']:
        # same times for all modes
        times[spin_sign] = copy.deepcopy(file[spin_sign]["l2_m2"]["times"][()])
    return times

#----------------------------------------------------------------------------------------------------
def mk_deepcopy_dictionary(h5dict, targetdict, keys):
    """
    Copy dictionary entries from the h5 file dictionary into an empty surrogate fit dictionary
    for future evaluation
    h5dict: dictionary in the h5 file opened through h5py
    targetdict: empty target directory
    keys: keys to be copied
    """
    for key in keys:
        if key=='name':
            targetdict['name'] = utils.chars_to_string(h5dict['name'][()])
        elif key=='fitType':
            targetdict['fitType'] = utils.chars_to_string(h5dict['fitType'][()])
        else:
            targetdict[key] = copy.deepcopy(h5dict[key][()])
            
#----------------------------------------------------------------------------------------------------
def extract_h5filegprsettings_to_emptydict(h_gpr, h_file, nnodes):
    """
    Copy all GPR fit related settings for any node in the amplitude/phase of real/imag fit
    h_gpr: empty dictionary to copy all GPR setttings
    h_file: dictionary located in the h5 file containing all GPR settings
    nnodes: number of EIM nodes
    """
    # iterate through time interpolation nodes (EIM indicies)
    for eim_indx in range(nnodes):

        h_gpr['node%s'%eim_indx] = {}
        h_gpr['node%s'%eim_indx]['GPR_params'] = {}

        h_gpr['node%s'%eim_indx]['lin_reg_params'] = {}
        h_gpr['node%s'%eim_indx]['GPR_params']['kernel_'] = {}
        for key in ['k1__k1', 'k2', 'k1__k2', 'k1']:
            h_gpr['node%s'%eim_indx]['GPR_params']['kernel_'][key] = {}
        for key in ['k1', 'k2']:
            h_gpr['node%s'%eim_indx]['GPR_params']['kernel_']['k1'][key] = {}

        h_gpr['node%s'%eim_indx]['fitType'] = utils.chars_to_string(h_file['fitType'][()])

        mk_deepcopy_dictionary(h_file['node%s'%eim_indx], 
                               h_gpr['node%s'%eim_indx],
                               ['data_mean', 'data_std'])

        mk_deepcopy_dictionary(h_file['node%s'%eim_indx]['lin_reg_params'], 
                               h_gpr['node%s'%eim_indx]['lin_reg_params'], 
                               ['coef_', 'intercept_'])

        mk_deepcopy_dictionary(h_file['node%s'%eim_indx]['GPR_params'], 
                               h_gpr['node%s'%eim_indx]['GPR_params'], 
                               ['X_train_', 'alpha_', '_y_train_mean', 'L_'])

        mk_deepcopy_dictionary(h_file['node%s'%eim_indx]['GPR_params']['kernel_'], 
                               h_gpr['node%s'%eim_indx]['GPR_params']['kernel_'], 
                               ['name', 'k2__noise_level', 'k2__noise_level_bounds', 'k1__k2__length_scale', 
                                'k1__k2__length_scale_bounds', 'k1__k1__constant_value', 'k1__k1__constant_value_bounds'])

        mk_deepcopy_dictionary(h_file['node%s'%eim_indx]['GPR_params']['kernel_']['k1'], 
                               h_gpr['node%s'%eim_indx]['GPR_params']['kernel_']['k1'], 
                               ['name', 'k1__constant_value', 'k1__constant_value_bounds', 'k2__length_scale', 
                                'k2__length_scale_bounds'])

        mk_deepcopy_dictionary(h_file['node%s'%eim_indx]['GPR_params']['kernel_']['k1']['k1'], 
                               h_gpr['node%s'%eim_indx]['GPR_params']['kernel_']['k1']['k1'], 
                               ['name', 'constant_value', 'constant_value_bounds'])

        mk_deepcopy_dictionary(h_file['node%s'%eim_indx]['GPR_params']['kernel_']['k1']['k2'], 
                               h_gpr['node%s'%eim_indx]['GPR_params']['kernel_']['k1']['k2'], 
                               ['name', 'length_scale', 'length_scale_bounds'])

        mk_deepcopy_dictionary(h_file['node%s'%eim_indx]['GPR_params']['kernel_']['k1__k1'], 
                               h_gpr['node%s'%eim_indx]['GPR_params']['kernel_']['k1__k1'], 
                               ['name', 'constant_value', 'constant_value_bounds'])

        mk_deepcopy_dictionary(h_file['node%s'%eim_indx]['GPR_params']['kernel_']['k1__k2'], 
                               h_gpr['node%s'%eim_indx]['GPR_params']['kernel_']['k1__k2'], 
                               ['name', 'length_scale', 'length_scale_bounds'])

        mk_deepcopy_dictionary(h_file['node%s'%eim_indx]['GPR_params']['kernel_']['k2'], 
                               h_gpr['node%s'%eim_indx]['GPR_params']['kernel_']['k2'], 
                               ['name', 'noise_level', 'noise_level_bounds'])
                        

#----------------------------------------------------------------------------------------------------
def load_surrogate(h5_data_dir, fname, wf_modes, nrcalib_modes):
    """ Loads all GPR interpolation data
            - Included modes = [(2,1),(3,1),(2,2),(3,2),(4,2),(3,3),(4,3),(4,4)]
            - NOTE: Requires h5 file to be in the same directory as this script.
            - Returns:
                - times, eim_indicies_amp, eim_indicies_ph, b_amp, b_ph, h_amp_gpr, h_ph_gpr
                - NOTE: times is dictionary with times.keys() = ['negative_spin', 'positive_spin']
    """
    
    with h5py.File('%s/%s'%(h5_data_dir,fname), 'r') as file:
        
        # obtain training time values
        times = read_times(file)
        
        # dicts to copy .h5 file data into (needed for surrogate generation)
        # Now for 2d surrogae each of dictionary will contain data for positive and negative spins 
        B_dict_amp, B_dict_ph = {}, {}
        fit_data_dict_amp, fit_data_dict_ph = {}, {}
        
        for spin_sign in ['negative_spin', 'positive_spin']:
            
            # empty dictionaries for holding basis vectors, EIM info and GPR settings
            B_dict_amp[spin_sign] = {}
            B_dict_ph[spin_sign] = {}
            fit_data_dict_amp[spin_sign] = {}
            fit_data_dict_ph[spin_sign] = {}
    
            # Copy data groups we need to access from hdf5 file into output dicts
            for mode in wf_modes:
                
                # splice out the relevant dictionary from h5 file for each mode
                f_mode = file[spin_sign]['l%s_m%s'%(mode[0], mode[1])]
                # obtain data groups for the amplitude and phases
                h_amp_file, h_ph_file = dict(f_mode['gpr_amp']), dict(f_mode['gpr_phase'])

                # basis matrix
                B_dict_amp[spin_sign][mode] = copy.deepcopy(f_mode["B_amp"][()])
                B_dict_ph[spin_sign][mode] = copy.deepcopy(f_mode["B_phase"][()])
                
                # EIM indicies
                eim_indicies_amp = copy.deepcopy(f_mode["eim_indicies_amp"][()])
                eim_indicies_ph = copy.deepcopy(f_mode["eim_indicies_phase"][()])
                
                # GPR settings
                h_eim_gpr_amp, h_eim_gpr_ph = {}, {}
                extract_h5filegprsettings_to_emptydict(h_eim_gpr_amp, h_amp_file, len(eim_indicies_amp))
                extract_h5filegprsettings_to_emptydict(h_eim_gpr_ph, h_ph_file, len(eim_indicies_ph))
                        
                # construct fit data for each mode
                fit_data_dict_amp[spin_sign][mode] = [h_eim_gpr_amp, eim_indicies_amp]
                fit_data_dict_ph[spin_sign][mode] = [h_eim_gpr_ph, eim_indicies_ph]

        # nr calibration info
        alpha_coeffs, beta_coeffs = load_splines.read_nrcalib_info(file, nrcalib_modes)

    return times, fit_data_dict_amp, fit_data_dict_ph, B_dict_amp, B_dict_ph, alpha_coeffs, beta_coeffs
                        


