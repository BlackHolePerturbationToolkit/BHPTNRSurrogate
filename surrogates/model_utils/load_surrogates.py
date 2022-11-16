##==============================================================================
## BHPTNRSurrogate module
## Description : loads BHPTNRSur1dq1e4 surrogate fits data
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
## Modified : 
##==============================================================================

import numpy as np
import h5py
import os
from os import path
import hashlib
from common_utils import load_splines as load_spl
from common_utils import load_GPRs
from common_utils import filehash

#----------------------------------------------------------------------------------------------------
def load_BHPTNRSur1dq1e4_surrogate(h5_data_dir):
    
    """ 
    Loads all interpolation data for the following modes
    modes=[(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),
          (5,4),(5,5),(6,4),(6,5),(6,6),(7,5),(7,6),(7,7),(8,6),
          (8,7),(8,8),(9,7),(9,8),(9,9),(10,8),(10,9)]
    Assumes the file BHPTNRSur1dq1e4.h5 is located in the h5_data_dir directory.
    """
    
    # h5 file name
    fname = 'BHPTNRSur1dq1e4.h5'
    # provide current zenodo hash
    zenodo_current_hash = "58a3a75e8fd18786ecc88cf98f694d4a"
    # zenodo url
    url = 'https://zenodo.org/record/7125742'
    # obtain zenodo ID
    zenodo_ID = url.rsplit("/")[-1]
    # obtain the hash for the current file; also downloads the file
    # if it doesn't exist in h5_data_dir
    file_hash = filehash.md5(fname, h5_data_dir, zenodo_ID)
    # check hash is the most recent
    filehash.check_current_hash(file_hash, zenodo_current_hash, url, fname)
    
    # modes to read fit data for
    wf_modes = [(2,2),(2,1),(3,1),(3,2),(3,3), (4,2),(4,3),(4,4),
             (5,3),(5,4),(5,5),(6,4),(6,5),(6,6),(7,5),(7,6),(7,7),
             (8,6),(8,7),(8,8),(9,7),(9,8),(9,9),(10,8),(10,9)]
    
    # modes used in nr calibration
    nrcalib_modes = [(2,2),(3,3),(4,4),(5,5)]
    
    # obtain all fit data
    time, fit_data_dict_1, fit_data_dict_2, B_dict_1, B_dict_2, alpha_coeffs, beta_coeffs \
                    = load_spl.load_surrogate(h5_data_dir, fname, wf_modes, nrcalib_modes)

    return time, fit_data_dict_1, fit_data_dict_2, B_dict_1, B_dict_2, alpha_coeffs, beta_coeffs


#----------------------------------------------------------------------------------------------------
def load_BHPTNRSur2dq1e3_surrogate(h5_data_dir):
    NotImplemented
    