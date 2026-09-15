"""Regression test: BHPTNRSur1dq1e4 waveforms against stored reference data.

Reference data was generated from gwsurrogate (verified to agree with
BHPTNRSurrogate to machine precision in tutorials/testing_for_pr_no9.ipynb)
and stored at 500 evenly-spaced time samples per waveform.

This test ensures that code changes do not alter waveform output.
"""

import os
import warnings

import numpy as np
import pytest
from scipy.interpolate import interp1d


pytestmark = pytest.mark.gwsurrogate

REFERENCE_FILE = os.path.join(os.path.dirname(__file__), "regression_data_1dq1e4.npz")

MODES_TEST = [
    (2, 2), (2, 1), (3, 1), (3, 2), (3, 3),
    (4, 2), (4, 3), (4, 4), (5, 3), (5, 4), (5, 5),
]

ATOL_STRAIN = 1e-11
ATOL_PHASE = 1e-9


@pytest.fixture(scope="module")
def sur_bhpt():
    """Load BHPTNRSur1dq1e4 from the native package."""
    from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4
    BHPTNRSur1dq1e4._ensure_loaded()
    return BHPTNRSur1dq1e4


@pytest.fixture(scope="module")
def reference_data():
    """Load stored reference waveforms."""
    if not os.path.exists(REFERENCE_FILE):
        pytest.fail("Required regression data file not found: %s" % REFERENCE_FILE)
    return dict(np.load(REFERENCE_FILE))


@pytest.mark.parametrize("q", [3, 10, 100])
def test_regression_against_reference(sur_bhpt, reference_data, q):
    """BHPTNRSur1dq1e4 output must match stored reference data."""
    t_ref = reference_data["t_q%d" % q]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t_bhpt, h_bhpt = sur_bhpt.generate_surrogate(
            q=q, modes=MODES_TEST, neg_modes=False,
        )

    for mode in MODES_TEST:
        ref_key = "h_q%d_l%d_m%d" % (q, mode[0], mode[1])
        h_ref = reference_data[ref_key]

        h_interp = (
            interp1d(t_bhpt, h_bhpt[mode].real, "cubic")(t_ref)
            + 1j * interp1d(t_bhpt, h_bhpt[mode].imag, "cubic")(t_ref)
        )

        max_dh = np.max(np.abs(h_interp - h_ref))
        assert max_dh < ATOL_STRAIN, (
            "q=%d mode=%s: max|dh| = %.4e exceeds %s" % (q, mode, max_dh, ATOL_STRAIN)
        )

        a_ref = np.abs(h_ref)
        a_new = np.abs(h_interp)
        thr = 1e-8 * max(np.max(a_ref), np.max(a_new))
        valid = (a_ref > thr) & (a_new > thr)
        if np.any(valid):
            max_dphi = np.max(np.abs(
                np.unwrap(np.angle(h_interp[valid]))
                - np.unwrap(np.angle(h_ref[valid]))
            ))
            assert max_dphi < ATOL_PHASE, (
                "q=%d mode=%s: max|dphi| = %.4e exceeds %s"
                % (q, mode, max_dphi, ATOL_PHASE)
            )
