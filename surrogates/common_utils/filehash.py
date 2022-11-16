##==============================================================================
## BHPTNRSurrogate module
## Description : checks hash of h5 files
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np
import os
from os import path
import hashlib

#----------------------------------------------------------------------------------------------------
def md5(fname, h5_data_dir, zenodo_ID):
    """ Compute hash from file. code taken from 
    https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file"""
    
    # download file if not already there
    if path.isfile('%s/%s'%(h5_data_dir,fname))==False:
        print('%s file is not found in the directory - PATH-TO/BHPTNRSurrogate/data/'%(fname))
        print('... downloading h5 file from zenodo')
        print('... this might take some time')
        os.system('wget https://zenodo.org/record/%s/files/BHPTNRSur1dq1e4.h5 -P %s'%(zenodo_ID,h5_data_dir))
        print('... downloaded')
    
    hash_md5 = hashlib.md5()
    with open('%s/%s'%(h5_data_dir,fname), "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

#----------------------------------------------------------------------------------------------------
def check_current_hash(file_hash, zenodo_current_hash, url, fname):
    # chech if the h5file is the most recent one or if it is corrupted
    if file_hash != zenodo_current_hash:
        raise AttributeError("%s out of date.\n \
                             Please download new version from %s"\
                             %(fname,url))
    else:
        pass
        