"""Integration tests for BHPTNRSur2dq1e3 waveform generation (requires h5 data)."""

import numpy as np
import pytest


pytestmark = pytest.mark.slow


class TestCalibratedGeneration:
    def test_positive_spin(self, model_2d):
        t, h = model_2d.generate_surrogate(q=10, spin1=0.5)
        assert isinstance(h, dict)
        assert (2, 2) in h
        assert len(t) > 0

    def test_negative_spin(self, model_2d):
        t, h = model_2d.generate_surrogate(q=10, spin1=-0.5)
        assert isinstance(h, dict)
        assert (2, 2) in h

    def test_zero_spin(self, model_2d):
        t, h = model_2d.generate_surrogate(q=10, spin1=0.0)
        assert isinstance(h, dict)
        assert (2, 2) in h


class TestUncalibratedGeneration:
    def test_uncalibrated_warns(self, model_2d):
        with pytest.warns(UserWarning, match="NOT NR calibrated"):
            model_2d.generate_surrogate(q=10, calibrated=False)


class TestOutOfRange:
    def test_q_out_of_range_warns(self, model_2d):
        with pytest.warns(UserWarning, match="outside bounds"):
            model_2d.generate_surrogate(q=2.0)

    def test_spin_out_of_range_warns(self, model_2d):
        with pytest.warns(UserWarning, match="outside bounds"):
            model_2d.generate_surrogate(q=10, spin1=0.9)
