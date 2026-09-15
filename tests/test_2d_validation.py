"""Validation tests for BHPTNRSur2dq1e3 (mocked — no h5 data needed)."""

from unittest.mock import patch

import numpy as np
import pytest

MODULE = "BHPTNRSurrogate.surrogates.BHPTNRSur2dq1e3"


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
             'times_dict': {'positive_spin': None, 'negative_spin': None},
             'fit_data_dict_1_sign': {'positive_spin': None, 'negative_spin': None},
             'fit_data_dict_2_sign': {'positive_spin': None, 'negative_spin': None},
             'B_dict_1_sign': {'positive_spin': None, 'negative_spin': None},
             'B_dict_2_sign': {'positive_spin': None, 'negative_spin': None},
             'alpha_coeffs': None, 'beta_coeffs': None,
         }), \
         patch(f"{MODULE}.eval_sur.evaluate_surrogate", return_value=_make_dummy_result()):
        yield


class TestExtraParamsWarning:
    def test_spin2_ignored(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
        with pytest.warns(UserWarning, match="Ignoring extra params.*spin2"):
            BHPTNRSur2dq1e3.generate_surrogate(q=10, spin1=0.3, spin2=0.1)

    def test_ecc_ignored(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
        with pytest.warns(UserWarning, match="Ignoring extra params.*ecc"):
            BHPTNRSur2dq1e3.generate_surrogate(q=10, spin1=0.3, ecc=0.1)

    def test_ano_ignored(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
        with pytest.warns(UserWarning, match="Ignoring extra params.*ano"):
            BHPTNRSur2dq1e3.generate_surrogate(q=10, spin1=0.3, ano=1.0)


class TestMassScaleValidation:
    def test_invalid_mass_scale_raises(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
        with pytest.raises(ValueError, match="mass_scale must be"):
            BHPTNRSur2dq1e3.generate_surrogate(q=10, mass_scale='bad')

    def test_calibrated_with_m1_raises(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
        with pytest.raises(ValueError, match="cannot be used when calibrated=True"):
            BHPTNRSur2dq1e3.generate_surrogate(q=10, calibrated=True, mass_scale='m1')

    def test_valid_mass_scale_M(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
        t, h = BHPTNRSur2dq1e3.generate_surrogate(q=10, mass_scale='M')
        assert isinstance(h, dict)

    def test_valid_mass_scale_m1(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
        t, h = BHPTNRSur2dq1e3.generate_surrogate(q=10, calibrated=False, mass_scale='m1')
        assert isinstance(h, dict)

    def test_physical_waveform_rejects_mass_scale_m1(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
        with pytest.raises(ValueError, match="only valid for geometric waveforms"):
            BHPTNRSur2dq1e3.generate_surrogate(q=10, M_tot=50, dist_mpc=100, calibrated=False, mass_scale='m1')


class TestSpinBranching:
    def test_positive_spin(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
        t, h = BHPTNRSur2dq1e3.generate_surrogate(q=10, spin1=0.5)
        assert isinstance(h, dict)

    def test_negative_spin(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
        t, h = BHPTNRSur2dq1e3.generate_surrogate(q=10, spin1=-0.5)
        assert isinstance(h, dict)

    def test_zero_spin(self):
        from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
        t, h = BHPTNRSur2dq1e3.generate_surrogate(q=10, spin1=0.0)
        assert isinstance(h, dict)
