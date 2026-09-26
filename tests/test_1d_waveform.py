"""Integration tests for BHPTNRSur1dq1e4 waveform generation (requires h5 data)."""

import numpy as np
import pytest


pytestmark = pytest.mark.slow


class TestCalibratedGeneration:
    def test_default_returns_dict_with_22(self, model_1d):
        t, h = model_1d.generate_surrogate(q=10)
        assert isinstance(h, dict)
        assert (2, 2) in h
        assert len(t) > 0
        assert len(h[(2, 2)]) == len(t)

    def test_specific_modes_only(self, model_1d):
        modes = [(2, 2), (3, 3)]
        t, h = model_1d.generate_surrogate(q=10, modes=modes, neg_modes=False)
        assert set(h.keys()) == set(modes)

    def test_neg_modes_true(self, model_1d):
        t, h = model_1d.generate_surrogate(q=10, modes=[(2, 2)], neg_modes=True)
        assert (2, 2) in h
        assert (2, -2) in h

    def test_neg_modes_false(self, model_1d):
        t, h = model_1d.generate_surrogate(q=10, modes=[(2, 2)], neg_modes=False)
        assert (2, 2) in h
        assert (2, -2) not in h


class TestUncalibratedGeneration:
    def test_uncalibrated_warns(self, model_1d):
        with pytest.warns(UserWarning, match="NOT NR calibrated"):
            model_1d.generate_surrogate(q=10, calibrated=False)


class TestOutOfRange:
    def test_q_out_of_range_warns(self, model_1d):
        with pytest.warns(UserWarning, match="outside bounds"):
            model_1d.generate_surrogate(q=1.5)


class TestMassScaleBehavior:
    def test_mass_scale_ratio(self, model_1d):
        """mass_scale='M' vs 'm1' should differ by factor q/(q+1) in time and strain."""
        q = 100.0
        t_M, h_M = model_1d.generate_surrogate(q=q, calibrated=False, mass_scale='M')
        t_m1, h_m1 = model_1d.generate_surrogate(q=q, calibrated=False, mass_scale='m1')

        expected_ratio = q / (q + 1.0)

        # time ratio
        time_ratio = t_M[-1] / t_m1[-1]
        np.testing.assert_allclose(time_ratio, expected_ratio, rtol=1e-10)

        # strain ratio for (2,2) mode
        strain_ratio = h_M[(2, 2)][-1] / h_m1[(2, 2)][-1]
        np.testing.assert_allclose(np.abs(strain_ratio), expected_ratio, rtol=1e-10)
