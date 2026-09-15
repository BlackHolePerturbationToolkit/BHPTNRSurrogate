##==============================================================================
## BHPTNRSurrogate module
## Description : generates ppBHPT surrogate model waveforms
## Author : Tousif Islam, Nov 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

from ..common_utils import utils, fits
from ..common_utils import check_inputs as checks

#----------------------------------------------------------------------------------------------------
def evaluate_surrogate(X_sur, X_calib, X_bounds, time, modes, modes_available, alpha_coeffs,\
                       beta_coeffs, alpha_beta_functional_form, calibrated, M_tot, dist_mpc,
                       orb_phase, inclination, fit_data_dict_1, fit_data_dict_2, B_dict_1, \
                       B_dict_2, fit_func, decomposition_funcs, norm, mode_sum, neg_modes, \
                       lmax, CoorbToInert, mass_factor=1.0):
    """
    Inputs
    ======
        X_sur : array of surrogate parameterization e.g. [log(q), spin1, spin2]
        
        X_calib : array of nr calibration parameterization e.g. [1/q, spin]

        X_bounds : bounds on all parameters - 2D array
        
        time : array of time on which surrogate has been trained on - read from h5 file

        modes : list of modes to evaluate
        
        modes_available : available mode in the surrogate

        alpha_coeffs : dictionary of alpha values obtained from calibration mode-by-mode

        beta_coeffs : beta value obtain from calibration - used in time rescaling

        alpha_beta_functional_form : function to use for nr calibration - must come from 
                                     common_utils.nr_calibration.py

        calibrated : True/False - whether requested waveforms should be nr calibrated or not

        M_tot : total mass of the binary

        dist_mpc : distance in mpc for the binary

        orb_phase : orbital phase at the start of the waveform

        inclination : inclination angle wrt the observer

        fit_data_dict_1, fit_data_dict_2 : dictionary of fit data obtained for two datapieces from 
                                           the h5 file.
                                           Keys are the modes.
                                           Structure may depend on ether the data comes from spline 
                                           fits or GPR fits. However, they should always be packed 
                                           in fit_data_dict_1 and fit_data_dict_2. 
                                           Make sure to modify your data loading script to achieve
                                           this if necessary.

        B_dict_1, B_dict_2 : dictionary of the basis matrices obtained from h5 file.
                             Modes used as keys.

        fit_func : form of fitting function. options : 'spline_1d' or 'GPR_fits'

        decomposition_funcs : form of data decomposition function to combine datapieces for 22 and
                              higher modes respectively. e.g. Amp/Phase to full or real/imag to full
                              etc. These functions are available at common_utils.utils.py

        norm : overall normalization factor to be multiplied to final waveform. This depends on the 
              way the surrogate have been constructed. Mostly norm=1/q or norm=1. 

        mode_sum : indicate whether modes should be summed up.

        neg_modes : indicate whether negative modes should be retured using orbital plane symmetry.

        lmax : maximum value of l upto which modes should be returned.

        CoorbToInert : indicate whether higher modes have been modelled in coorbital frame. In that
                       case, additional processing will be performed.
    
    Outputs
    =======
    
        t_surrogate : time array
        
        h_surrogate : dictiornary of modes if mode_sum not requested
                      full waveform if mode_sum is requested
    
     
    """
    
    # check inputs
    checks.check_user_inputs(X_sur, X_bounds, modes, modes_available, M_tot, dist_mpc, 
                      orb_phase, inclination, mode_sum)
    
    # uncalibrated waveforms in geometric units
    hsur_raw_dict = fits.all_modes_surrogate(modes, X_sur, fit_data_dict_1, fit_data_dict_2, \
                           B_dict_1, B_dict_2, lmax, fit_func, decomposition_funcs, norm)
    
    # process the raw surrogate output depending on the user inputs
    t_surrogate, h_surrogate = utils.obtain_processed_output(X_calib, time, hsur_raw_dict, alpha_coeffs,
                                    beta_coeffs, alpha_beta_functional_form, calibrated, M_tot, dist_mpc,
                                    orb_phase, inclination, mode_sum, neg_modes, lmax, CoorbToInert,
                                    mass_factor=mass_factor)
    
    return t_surrogate, h_surrogate
