##==============================================================================
## BHPTNRSur1dq1e4 : arXiv:2204.01972
## Description : utility file for surrogate model
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import numpy as np
from gwtools import gwtools as _gwtools
from gwtools.harmonics import sYlm as _sYlm

#---------------------------------------------------------------------------------------------------- 
def geo_to_SI(t_geo, h_geo, M_tot, dist_mpc):
    """
    transforms the waveform from geomeric unit to physical unit
    given geoemtric time, geometric waveform, total mass M, distance dL
    """    
    # Physical units
    G=_gwtools.G
    MSUN_SI = _gwtools.MSUN_SI
    PC_SI = _gwtools.PC_SI
    C_SI = _gwtools.c
    M = M_tot * MSUN_SI
    dL = dist_mpc * PC_SI
    
    # scaling of time and strain
    t_SI = t_geo * (G*M/C_SI**3)
    
    strain_geo_to_SI = (G*M/C_SI**3)/dL
    h_SI={}
    for mode in h_geo.keys():
        h_SI[(mode)] = np.array(h_geo[mode])*strain_geo_to_SI
    
    return t_SI, h_SI


#---------------------------------------------------------------------------------------------------- 
def phase_rotation(h, delta_orb_phase):
    """
    performs an orbital phase rotation
    """
    h_rotated = {}
    
    for mode in h.keys():
        (l,m) = mode
        phase_rot = m*delta_orb_phase
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
            sYlm_value =  _sYlm(-2,ll=ell,mm=m,theta=theta,phi=phi)
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
        h_dict_all_modes[(l,m)] = h_dict[mode]
        if (m==0):
            raise ValueError('m must be nonnegative. m<0 will be generated for you from the m>0 mode.')
        elif (m<0):
            raise ValueError('m must be nonnegative. m<0 will be generated for you from the m>0 mode.')
        else:
            h_dict_all_modes[(l,-m)] = np.power(-1,l) * np.conjugate(h_dict[mode])

    return h_dict_all_modes