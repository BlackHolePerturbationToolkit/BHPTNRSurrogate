##==============================================================================
## BHPTNRSurrogate module
## Description : checks whether the user inputs are correct to generate a waveform
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np

#---------------------------------------------------------------------------------------------------- 
def check_user_inputs(modes, M_tot, dist_mpc, orb_phase, inclination, mode_sum):
    """ 
        Checks whether the user inputs are valid 
         
    """
        
    # geometric or SI units
    if (M_tot is None) ^ (dist_mpc is None):
        raise ValueError("Either specify both M_tot and dist_mpc, or neither")
           
    # check phase and inclination
    if (orb_phase is None) ^ (inclination is None):
        raise ValueError("Either specify both orb_phase and inclination, or neither")   
        
    # for mode summation, check mass, distance and angles have been specified
    if mode_sum==True:
        if M_tot is None and dist_mpc is None and orb_phase is None and inclination is None:
            raise ValueError("M_tot, dist_mpc, orb_phase and inclination should NOT be None")
            
    # TOD0 : add a check for domain of validity
            