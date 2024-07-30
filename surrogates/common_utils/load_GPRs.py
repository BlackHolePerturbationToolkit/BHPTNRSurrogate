##==============================================================================
## BHPTNRSurrogate module
## Description : loads GPR surrogate fits data
## Author: Katie Rink, Mar 2023 [krink@utexas.edu]
##==============================================================================

import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import copy
from . import load_splines
from . import utils


def load_surrogate(h5_data_dir, fname, wf_modes, nrcalib_modes):
    """ Loads all interpolation data.
            - Included modes = [(2,1),(3,1),(2,2),(3,2),(4,2),(3,3),(4,3),(4,4)]
            - NOTE: Requires h5 file to be in the same directory as this script.
            - Returns:
                - times, eim_indicies_amp, eim_indicies_ph, b_amp, b_ph, h_amp_gpr, h_ph_gpr
                - NOTE: times is dictionary with times.keys() = ['negative_spin', 'positive_spin']
    """
    with h5py.File('%s/%s'%(h5_data_dir,fname), 'r') as file:
        
        # dicts to copy .h5 file data into (needed for surrogate generation)
        # Now for 2d surrogae each of dictionary will contain data for positive and negative spins 
        b_dict_amp, b_dict_ph = {}, {}
        fit_data_dict_amp, fit_data_dict_ph = {}, {}
        times = {}

        for spin_sign in ['negative_spin', 'positive_spin']:
            times[spin_sign] = copy.deepcopy(file[spin_sign]["l2_m2"]["times"][()]) # same times for all modes
            b_dict_amp[spin_sign] = {}
            b_dict_ph[spin_sign] = {}
            fit_data_dict_amp[spin_sign] = {}
            fit_data_dict_ph[spin_sign] = {}

            for mode in wf_modes:
                # Copy data groups we need to access from hdf5 file into output dicts.
                f_mode = file[spin_sign]['l%s_m%s'%(mode[0], mode[1])]
                h_amp_file, h_ph_file = dict(f_mode['gpr_amp']), dict(f_mode['gpr_phase'])

                h_eim_gpr_amp, h_eim_gpr_ph = {}, {}
                eim_indicies_amp, eim_indicies_ph = {}, {}

                b_dict_amp[spin_sign][mode] = copy.deepcopy(f_mode["B_amp"][()])
                b_dict_ph[spin_sign][mode] = copy.deepcopy(f_mode["B_phase"][()])
                eim_indicies_amp = copy.deepcopy(f_mode["eim_indicies_amp"][()])
                eim_indicies_ph = copy.deepcopy(f_mode["eim_indicies_phase"][()])

                # iterate through GPR fit params once per data piece (amp, ph)

                h_params_file = [h_amp_file, h_ph_file]
                h_params_gpr = [h_eim_gpr_amp, h_eim_gpr_ph]
                eim_indicies = [eim_indicies_amp, eim_indicies_ph]

                for i in range(2): # run through once for amp data, then again for phase data
                    h_file = h_params_file[i] # hdf5 source data
                    h_gpr = h_params_gpr[i] # empty dictionary

                    # iterate through time interpolation nodes (EIM indicies)
                    for eim_indx in range(len(eim_indicies[i])):

                        h_file_node = h_file['node%s'%eim_indx]
                        h_file_node_params = h_file_node['GPR_params']

                        h_gpr['node%s'%eim_indx] = {}
                        h_gpr_node = h_gpr['node%s'%eim_indx]

                        h_gpr_node['GPR_params'] = {}
                        h_gpr_node_params = h_gpr_node['GPR_params']

                        h_gpr_node['lin_reg_params'] = {}
                        h_gpr_node_params['kernel_'] = {}
                        h_gpr_node_params['kernel_']['k1__k1'] = {}
                        h_gpr_node_params['kernel_']['k2'] = {}
                        h_gpr_node_params['kernel_']['k1__k2'] = {}
                        h_gpr_node_params['kernel_']['k1'] = {}
                        h_gpr_node_params['kernel_']['k1']['k1'] = {}
                        h_gpr_node_params['kernel_']['k1']['k2'] = {}

                        h_gpr_node['fitType'] = utils.chars_to_string(h_file['fitType'][()])
                        h_gpr_node['data_mean'] = copy.deepcopy(h_file_node['data_mean'][()])
                        h_gpr_node['data_std'] = copy.deepcopy(h_file_node['data_std'][()])

                        h_gpr_node['lin_reg_params']['coef_'] = copy.deepcopy(h_file_node['lin_reg_params']['coef_'][()])
                        h_gpr_node['lin_reg_params']['intercept_'] = copy.deepcopy(h_file_node['lin_reg_params']['intercept_'][()])

                        h_gpr_node_params['X_train_'] = copy.deepcopy(h_file_node_params['X_train_'][()])
                        h_gpr_node_params['alpha_'] = copy.deepcopy(h_file_node_params['alpha_'][()])
                        h_gpr_node_params['_y_train_mean'] = copy.deepcopy(h_file_node_params['_y_train_mean'][()])
                        h_gpr_node_params['L_'] = copy.deepcopy(h_file_node_params['L_'][()])

                        h_gpr_node_params['kernel_']['name'] = utils.chars_to_string(h_file_node_params['kernel_']['name'])
                        h_gpr_node_params['kernel_']['k2__noise_level'] = copy.deepcopy(h_file_node_params['kernel_']['k2__noise_level'][()])
                        h_gpr_node_params['kernel_']['k2__noise_level_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k2__noise_level_bounds'][()])
                        h_gpr_node_params['kernel_']['k1__k2__length_scale'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k2__length_scale'][()])
                        h_gpr_node_params['kernel_']['k1__k2__length_scale_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k2__length_scale_bounds'][()])
                        h_gpr_node_params['kernel_']['k1__k1__constant_value'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k1__constant_value'][()])
                        h_gpr_node_params['kernel_']['k1__k1__constant_value_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k1__constant_value_bounds'][()])

                        h_gpr_node_params['kernel_']['k1']['name'] = utils.chars_to_string(h_file_node_params['kernel_']['k1']['name'])
                        h_gpr_node_params['kernel_']['k1']['k1__constant_value'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k1__constant_value'][()])
                        h_gpr_node_params['kernel_']['k1']['k1__constant_value_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k1__constant_value_bounds'][()])
                        h_gpr_node_params['kernel_']['k1']['k2__length_scale'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k2__length_scale'][()])
                        h_gpr_node_params['kernel_']['k1']['k2__length_scale_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k2__length_scale_bounds'][()])

                        h_gpr_node_params['kernel_']['k1']['k1']['name'] = utils.chars_to_string(h_file_node_params['kernel_']['k1']['k1']['name'])
                        h_gpr_node_params['kernel_']['k1']['k1']['constant_value'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k1']['constant_value'][()])
                        h_gpr_node_params['kernel_']['k1']['k1']['constant_value_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k1']['constant_value_bounds'][()])
                        h_gpr_node_params['kernel_']['k1']['k2']['name'] = utils.chars_to_string(h_file_node_params['kernel_']['k1']['k2']['name'])
                        h_gpr_node_params['kernel_']['k1']['k2']['length_scale'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k2']['length_scale'][()])
                        h_gpr_node_params['kernel_']['k1']['k2']['length_scale_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k2']['length_scale_bounds'][()])

                        h_gpr_node_params['kernel_']['k1__k1']['name'] = utils.chars_to_string(h_file_node_params['kernel_']['k1__k1']['name'])
                        h_gpr_node_params['kernel_']['k1__k1']['constant_value'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k1']['constant_value'][()])
                        h_gpr_node_params['kernel_']['k1__k1']['constant_value_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k1']['constant_value_bounds'][()])

                        h_gpr_node_params['kernel_']['k1__k2']['name'] = utils.chars_to_string(h_file_node_params['kernel_']['k1__k2']['name'])
                        h_gpr_node_params['kernel_']['k1__k2']['length_scale'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k2']['length_scale'][()])
                        h_gpr_node_params['kernel_']['k1__k2']['length_scale_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k2']['length_scale_bounds'][()])

                        h_gpr_node_params['kernel_']['k2']['name'] = utils.chars_to_string(h_file_node_params['kernel_']['k2']['name'])
                        h_gpr_node_params['kernel_']['k2']['noise_level'] = copy.deepcopy(h_file_node_params['kernel_']['k2']['noise_level'][()])
                        h_gpr_node_params['kernel_']['k2']['noise_level_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k2']['noise_level_bounds'][()])

                fit_data_dict_amp[spin_sign][mode] = [h_eim_gpr_amp, eim_indicies_amp]
                fit_data_dict_ph[spin_sign][mode] = [h_eim_gpr_ph, eim_indicies_ph]

        # nr calibration info
        alpha_coeffs, beta_coeffs = load_splines.read_nrcalib_info(file, nrcalib_modes)

    return times, fit_data_dict_amp, fit_data_dict_ph, b_dict_amp, b_dict_ph, alpha_coeffs, beta_coeffs
