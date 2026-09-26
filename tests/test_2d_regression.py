"""Regression test: BHPTNRSur2dq1e3 waveforms against stored reference data.

Reference data was generated from the main branch (04cc0fb) by
Reg_test_2D/generate_reference_2D.py, and stored at 500 evenly-spaced time samples
per waveform. See Reg_test_2D/README.md for how to regenerate it.

This test ensures that code changes do not alter 2D waveform output. It is the 2D
counterpart of test_gwsurrogate_comparison.py, which does the same for the 1D model.

Mass convention: main has no mass_scale argument, so with calibrated=False it returns
waveforms in the primary-mass (m1) convention. The reference is therefore in that
convention, and this test asks for it explicitly.
"""

import os
import warnings

import numpy as np
import pytest
from scipy.interpolate import interp1d


pytestmark = pytest.mark.slow

REFERENCE_FILE = os.path.join(os.path.dirname(__file__), "regression_data_2dq1e3_04cc0fb.npz")

MODES_TEST = [(2, 2), (2, 1), (3, 1), (3, 2), (3, 3), (4, 2), (4, 3), (4, 4)]

# (q, spin1): both spin signs, spanning the q = 3 to 1000 domain
POINTS = [(3.0, 0.55), (12.5, 0.11), (47.0, -0.42), (115.0, 0.0), (501.0, -0.39)]

# Strain is compared relative to the peak of the 22 mode at the same point, so that the
# weaker modes are measured against something meaningful rather than their own noise
# floor. max|h22| spans 0.5 at q=3 down to 0.003 at q=501, so an absolute tolerance
# (as used for the 1D model) would be far stricter at large q than at small q.
#
# The model's own accuracy is ~1e-4, so these leave three orders of headroom while
# still catching any real change in the evaluation path.
RTOL_STRAIN = 1e-6
ATOL_PHASE = 1e-5


def tag_for(q, spin1):
    """Key used in the reference file, e.g. q12p5_chim0p42."""
    fmt = lambda x: ('m' if x < 0 else '') + ('%g' % abs(x)).replace('.', 'p')
    return 'q%s_chi%s' % (fmt(q), fmt(spin1))


@pytest.fixture(scope="module")
def reference_data():
    """Load stored reference waveforms."""
    if not os.path.exists(REFERENCE_FILE):
        pytest.fail("Required regression data file not found: %s" % REFERENCE_FILE)
    return dict(np.load(REFERENCE_FILE))


def _generate(model_2d, q, spin1):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return model_2d.generate_surrogate(
            q=q, spin1=spin1, modes=MODES_TEST, neg_modes=False,
            calibrated=False, mass_scale='m1',
        )


@pytest.mark.parametrize("q, spin1", POINTS, ids=["%s_%s" % p for p in POINTS])
def test_regression_against_reference(model_2d, reference_data, q, spin1):
    """BHPTNRSur2dq1e3 output must match stored reference data."""
    tag = tag_for(q, spin1)
    t_ref = reference_data["t_%s" % tag]

    t_new, h_new = _generate(model_2d, q, spin1)

    # scale for the strain comparison: the peak of the dominant mode at this point
    scale = np.max(np.abs(reference_data["h_%s_l2_m2" % tag]))

    for mode in MODES_TEST:
        h_ref = reference_data["h_%s_l%d_m%d" % (tag, mode[0], mode[1])]

        h_interp = (
            interp1d(t_new, h_new[mode].real, "cubic")(t_ref)
            + 1j * interp1d(t_new, h_new[mode].imag, "cubic")(t_ref)
        )

        max_dh = np.max(np.abs(h_interp - h_ref))
        assert max_dh / scale < RTOL_STRAIN, (
            "q=%g spin1=%g mode=%s: max|dh|/max|h22| = %.4e exceeds %s"
            % (q, spin1, mode, max_dh / scale, RTOL_STRAIN)
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
                "q=%g spin1=%g mode=%s: max|dphi| = %.4e exceeds %s"
                % (q, spin1, mode, max_dphi, ATOL_PHASE)
            )


@pytest.mark.parametrize("q, spin1", POINTS, ids=["%s_%s" % p for p in POINTS])
def test_time_array_unchanged(model_2d, reference_data, q, spin1):
    """The time array does not depend on the fit evaluation, so it must match exactly."""
    tag = tag_for(q, spin1)
    t_ref = reference_data["t_%s" % tag]

    t_new, _ = _generate(model_2d, q, spin1)

    idx = np.linspace(0, len(t_new) - 1, len(t_ref)).astype(int)
    assert np.array_equal(np.asarray(t_new)[idx], t_ref), (
        "q=%g spin1=%g: time array differs from the reference by up to %.4e"
        % (q, spin1, np.max(np.abs(np.asarray(t_new)[idx] - t_ref)))
    )
