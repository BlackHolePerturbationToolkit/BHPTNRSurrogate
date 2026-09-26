"""Unit tests for BHPTNRSurrogate/surrogates/common_utils/check_inputs.py."""

import warnings

import pytest
from BHPTNRSurrogate.surrogates.common_utils.check_inputs import (
    check_input_modes,
    check_domain_of_validity,
    check_extrinsic_params,
)


# ---------------------------------------------------------------------------
# check_input_modes
# ---------------------------------------------------------------------------
class TestCheckInputModes:
    def test_valid_subset_passes(self):
        available = [(2, 2), (2, 1), (3, 3)]
        requested = [(2, 2), (3, 3)]
        check_input_modes(requested, available)  # should not raise

    def test_invalid_mode_raises(self):
        available = [(2, 2), (2, 1)]
        requested = [(2, 2), (4, 4)]
        with pytest.raises(ValueError, match="NOT a subset"):
            check_input_modes(requested, available)


# ---------------------------------------------------------------------------
# check_domain_of_validity
# ---------------------------------------------------------------------------
class TestCheckDomainOfValidity:
    def test_within_bounds_no_warning(self):
        X_bounds = [[0.0], [10.0]]
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            check_domain_of_validity(5.0, X_bounds)

    def test_below_lower_bound_warns(self):
        X_bounds = [[2.0], [10.0]]
        with pytest.warns(UserWarning, match="outside bounds"):
            check_domain_of_validity(1.0, X_bounds)

    def test_above_upper_bound_warns(self):
        X_bounds = [[0.0], [5.0]]
        with pytest.warns(UserWarning, match="outside bounds"):
            check_domain_of_validity(6.0, X_bounds)

    def test_list_input_within_bounds(self):
        X_bounds = [[0.0, -1.0], [10.0, 1.0]]
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            check_domain_of_validity([5.0, 0.0], X_bounds)

    def test_list_input_one_param_out(self):
        X_bounds = [[0.0, -1.0], [10.0, 1.0]]
        with pytest.warns(UserWarning, match="outside bounds"):
            check_domain_of_validity([5.0, 2.0], X_bounds)


# ---------------------------------------------------------------------------
# check_extrinsic_params
# ---------------------------------------------------------------------------
class TestCheckExtrinsicParams:
    def test_all_none_ok(self):
        check_extrinsic_params(None, None, None, None, mode_sum=False)

    def test_all_specified_ok(self):
        check_extrinsic_params(50.0, 100.0, 0.0, 0.5, mode_sum=False)

    def test_mtot_without_dist_raises(self):
        with pytest.raises(ValueError, match="both M_tot and dist_mpc"):
            check_extrinsic_params(50.0, None, None, None, mode_sum=False)

    def test_dist_without_mtot_raises(self):
        with pytest.raises(ValueError, match="both M_tot and dist_mpc"):
            check_extrinsic_params(None, 100.0, None, None, mode_sum=False)

    def test_phase_without_inclination_raises(self):
        with pytest.raises(ValueError, match="both orb_phase and inclination"):
            check_extrinsic_params(None, None, 0.0, None, mode_sum=False)

    def test_inclination_without_phase_raises(self):
        with pytest.raises(ValueError, match="both orb_phase and inclination"):
            check_extrinsic_params(None, None, None, 0.5, mode_sum=False)

    def test_mode_sum_without_params_raises(self):
        with pytest.raises(ValueError, match="should NOT be None"):
            check_extrinsic_params(None, None, None, None, mode_sum=True)
