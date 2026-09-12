##==============================================================================
## BHPTNRSurrogate module
## Description : docs for different surrogate model wrappers
## Author : Tousif Islam, Nov 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import inspect

def copy_doc(base_func, model_func):
    """
    function to copy doc string from the base function and the model function
    and join them to create docstring fot the target function
    """
    def wrapper(target_func):
        target_func.__doc__ = base_func.__doc__ + model_func.__doc__
        return target_func
    return wrapper
    
    
def generic_doc_for_models() -> None:
    """
    ## -------------------------------------------------------------------------- ##
    ## general overview
    ## -------------------------------------------------------------------------- ##
    
    Description : wrapper to generate BHPT surrogate waveforms
    
    Input
    =====
    q: mass ratio (with q:=m_1/m_2>=1 where m_1 is the mass of the larger black hole 
                   and m_2 is the mass of the smaller black hole)
    
    chi1: dimensionless spin of the primary black hole where -1 <= chi1 <= 1
          Default: None
    
    chi2: dimensionless spin of the secondary black hole where -1 <= chi2 <= 1
          Default: None
          Not Implemented in any model so far
    
    ecc: eccentricity / Default: None
         Not Implemented in any model so far
    
    ano: mean anomaly / Default: None
         Not Implemented in any model so far
    
    modes:  list of modes
            Default (None) corresponds to all available modes in the model

            e.g. [(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),
                  (5,3),(5,4),(5,5),(6,4),(6,5),(6,6),(7,5),(7,6),
                  (7,7),(8,6),(8,7),(8,8),(9,7),(9,8),(9,9),(10,8),(10,9)]
            for BHPTNRSur1dq1e4
            
    M_tot: total mass of the binary in solar masses
           Default: None (in which case a geometric waveform is returned)
    
    dist_mpc:  distance of the binary from the observer in Mpc
               Default: None (in which case geometric wf is returned)
               
    orb_phase: orbital phase at the start of the waveform
    
    inclination: inclination angle
    
    lmax:  5 (default for BHPTNRSur1dq1e4)
           modes are only calibrated to NR up to some maximum value of l. 
           Please consult a specific model's documentation for information on how modes 
           are calibrated to NR.
           Default value changes depending on the model.
           Note: If one provides a list of modes, modes beyond lmax will not be returned
           Computes both positive and negative m modes upto lmax.
            
    mode_sum:  If true, all modes are summed. If false all modes are returned 
               in a dictionary. Default: false
               Note: Only works when orb_phase and inclination are not None
               
    calibrated:  Whether you want NR-calibrated waveform or not
                 When set to True, it applies a scaling to the uncalibrated
                 surrogate waveform. This scaling has been obtained by calibrating
                 the ppBHPT waveforms to NR in comparable mass ratio
                 regime.
                 Please consult a specific model's documentation for information on
                 how modes are calibrated to NR.
                 If set to False, the raw (uncalibrated) ppBHPT waveforms are returned.
                 Default: True

    mass_scale:  Mass convention for the waveform output. Options: 'M' or 'm1'.
                 Default: 'M'
                 When calibrated=True, this is ignored (NR calibration already uses
                 total mass M).
                 When calibrated=False and mass_scale='M', the raw ppBHPT waveform
                 (which uses m1 as the mass scale) is rescaled to total mass M by
                 applying a factor of q/(q+1) to both time and strain.
                 When calibrated=False and mass_scale='m1', the raw ppBHPT waveform
                 is returned without any rescaling.

    Output
    ======
    t : time
    h : waveform modes as a dictionary


    Example Uses:
    =============
    1. to obtain uncalibrated (raw) geometric waveform
            t, h = generate_surrogate(q=8, modes=[(2,1),(2,2),(3,1),(3,2),(3,3),
                        (4,2),(4,3),(4,4),(5,3),(5,4),(5,5)], calibrated=False)
    2. to obtain NR Calibrated geometric waveform
            t, h = generate_surrogate(q=8, modes=[(2,1),(2,2)])       
    3. to obtain NR calibrated physical waveform
            t, h = generate_surrogate(q=8, modes=[(2,1),(2,2)], M_tot=50, dist_mpc=100)
    4. to obtain NR calibrated physical waveform on a sphere
            t, h = generate_surrogate(q=8, modes=[(2,1),(2,2)], M_tot=50, dist_mpc=100, 
                                     orb_phase=np.pi/3,inclination=np.pi/4)
    5. to obtain NR calibrated physic al waveform on a sphere for modes up to l=5
            t, h = generate_surrogate(q=8, modes=[(2,1),(2,2)], M_tot=50, dist_mpc=100, 
                                      orb_phase=np.pi/3, inclination=np.pi/4, lmax=5)
    6. to obtain mode-summed NR calibrated physical waveform on a sphere
            t, h = generate_surrogate(q=8, M_tot=60, dist_mpc=100, orb_phase=np.pi/3, 
                                      inclination=np.pi/4, lmax=3, mode_sum=True)
              
    """
    return


def BHPTNRSur1dq1e4_doc() -> None:
    """
    ## -------------------------------------------------------------------------- ##
    ## BHPTNRSur1dq1e4
    ## -------------------------------------------------------------------------- ##
    
    This is a surrogate for non-spinning black hole binary systems with mass-ratios 
    varying from 2.5 to 10000. 
    
    This surrogate model is trained on waveform data generated by point-particle black 
    hole perturbation theory (ppBHPT), and tuned to NR simulations in the comparable 
    mass ratio regime (q=3 to q=10). 
    
    Available modes are: (2,1),(2,2),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),
                         (5,4),(5,5),(6,4),(6,5),(6,6),(7,5),(7,6),(7,7),(8,6),
                         (8,7),(8,8),(9,7),(9,8),(9,9),(10,8),(10,9)].                     
    The m<0 modes are deduced from the m>0 modes. 
    
    Model details can be found in Islam et al. 2022, arXiv:2204.01972.      
    """
    return


def BHPTNRSur2dq1e3_doc() -> None:
    """
    ## -------------------------------------------------------------------------- ##
    ## BHPTNRSur2dq1e3
    ## -------------------------------------------------------------------------- ##
    
    This is a surrogate for aligned-spinning black hole binary systems (where spin is
    on the primary black hole) with mass-ratios varying from 3 to 1000 and dimensionless
    spin on the primary black hole varying from -0.8 to 0.8. 
    
    This surrogate model is trained on waveform data generated by point-particle black 
    hole perturbation theory (ppBHPT), and tuned to NR simulations in the comparable 
    mass ratio regime (q=3 to q=10). 
    
    Available modes are: (2,1),(2,2),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4)].                     
    The m<0 modes are deduced from the m>0 modes. 
    
    Model details can be found in arXiv:2407.18319. 
    """
    return

