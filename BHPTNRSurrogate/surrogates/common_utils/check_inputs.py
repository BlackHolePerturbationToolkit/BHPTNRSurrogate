##==============================================================================
## BHPTNRSurrogate module
## Description : checks whether the user inputs are correct to generate a waveform
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import warnings

import numpy as np

#---------------------------------------------------------------------------------------------------- 
def check_input_modes(modes_requested, modes_available):
    """
        Checks whether modes requested by the user are avilable in h5 file for the
        surrogate 
    """
    if set(modes_available).issuperset(modes_requested) == True:
        pass
    else:
        raise ValueError("Set of the requested modes is NOT a subset of the \
                         available modes for the surrogate")

#---------------------------------------------------------------------------------------------------- 
def check_domain_of_validity(X_in, X_bounds):
    """
        Checks whether the user input parameters are within surrogate training
        space
    """
    # if the input X is just a number, make sure to array it up
    if isinstance(X_in,(list))==False and isinstance(X_in,(float,int))==True:
        X_in= [X_in]
    # raise error for all other scenarios
    #else:
    #    raise ValueError("param types are not matching with bound types")

    for param_indx in range(len(X_in)):
        if X_in[param_indx]<X_bounds[0][param_indx] or X_in[param_indx]>X_bounds[1][param_indx]:
            warnings.warn(
                "Input parameter is outside bounds for parameter value at index %d" % param_indx,
                stacklevel=3,
            )
        
#---------------------------------------------------------------------------------------------------- 
def check_extrinsic_params(M_tot, dist_mpc, orb_phase, inclination, mode_sum):
    """ 
        Checks whether the user inputs realted to extrinsic params are valid  
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
            
            
#---------------------------------------------------------------------------------------------------- 
def check_user_inputs(X_in, X_bounds, modes_requested, modes_available, M_tot, dist_mpc, 
                      orb_phase, inclination, mode_sum):
    """ 
        Checks whether the user inputs are valid   
        
    Inputs
    ======
        X_in : array of surrogate parameterization e.g. [log(q), spin1, spin2]
        X_bounds : bounds on all parameters - 2D array
        modes_requested : list of modes to evaluate
        modes_available : available mode in the surrogate
        M_tot : total mass of the binary
        dist_mpc : distance in mpc for the binary
        orb_phase : orbital phase at the start of the waveform
        inclination : inclination angle wrt the observer
        mode_sum : indicate whether modes should be summed up
    
    """
        
    # check extrinsic param inputs make sense
    check_extrinsic_params(M_tot, dist_mpc, orb_phase, inclination, mode_sum)
    
    # check input params lie within the surrogate training space
    check_domain_of_validity(X_in, X_bounds)
    
    # check requested modes exits
    check_input_modes(modes_requested, modes_available)
    
