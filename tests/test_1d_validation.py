"""Validation tests for BHPTNRSur1dq1e4 (mocked — no h5 data needed)."""

import warnings
from unittest.mock import patch, MagicMock

import numpy as np
import pytest

MODULE = "BHPTNRSurrogate.surrogates.BHPTNRSur1dq1e4"


def _make_dummy_result():
    """Return (time_array, mode_dict) matching evaluate_surrogate output."""
    t = np.linspace(-100, 0, 500)
    h = {(2, 2): np.ones(500, dtype=complex)}
    return t, h


@pytest.fixture(autouse=True)
def _mock_loading():
    """Patch _ensure_loaded and evaluate_surrogate so no h5 file is needed."""
    with patch(f"{MODULE}._ensure_loaded"), \
         patch(f"{MODULE}._surrogate_data", {
             'time': None, 'fit_data_dict_1': None, 'fit_data_dict_2': None,
             'B_dict_1': None, 'B_dict_2': None,
             'alpha_coeffs': None, 'beta_coeffs': None,
         }), \
         patch(f"{MODULE}.eval_sur.evaluate_surrogate", return_value=_make_dummy_result()):
        yield


class TestExtraParamsWarning:
    def test_spin1_ignored(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4
        with pytest.warns(UserWarning, match="Ignoring extra params.*spin1"):
            BHPTNRSur1dq1e4.generate_surrogate(q=10, spin1=0.5)

    def test_ecc_ignored(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4
        with pytest.warns(UserWarning, match="Ignoring extra params.*ecc"):
            BHPTNRSur1dq1e4.generate_surrogate(q=10, ecc=0.1)

    def test_ano_ignored(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4
        with pytest.warns(UserWarning, match="Ignoring extra params.*ano"):
            BHPTNRSur1dq1e4.generate_surrogate(q=10, ano=1.0)


class TestMassScaleValidation:
    def test_invalid_mass_scale_raises(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4
        with pytest.raises(ValueError, match="mass_scale must be"):
            BHPTNRSur1dq1e4.generate_surrogate(q=10, mass_scale='invalid')

    def test_calibrated_with_m1_warns(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4
        with pytest.warns(UserWarning, match="mass_scale is ignored"):
            BHPTNRSur1dq1e4.generate_surrogate(q=10, calibrated=True, mass_scale='m1')

    def test_valid_mass_scale_M(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4
        t, h = BHPTNRSur1dq1e4.generate_surrogate(q=10, mass_scale='M')
        assert isinstance(h, dict)

    def test_valid_mass_scale_m1(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4
        t, h = BHPTNRSur1dq1e4.generate_surrogate(q=10, calibrated=False, mass_scale='m1')
        assert isinstance(h, dict)

    def test_physical_waveform_rejects_mass_scale_m1(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4
        with pytest.raises(ValueError, match="only valid for geometric waveforms"):
            BHPTNRSur1dq1e4.generate_surrogate(q=10, M_tot=50, dist_mpc=100, calibrated=False, mass_scale='m1')
