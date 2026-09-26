#!/usr/bin/env python3
"""Generate the reference data used by test_2d_regression.py.

Evaluates the installed BHPTNRSur2dq1e3 at a fixed set of (q, spin1) points and
stores the waveforms, so later versions can be checked against them.

    python tests/generate_regression_data_2d.py --out tests/regression_data_2dq1e3.npz

Only run this to establish a new baseline — for a new model release, or when a change
to the waveforms is intended and has been checked. Regenerating it to make a failing
test pass would defeat the point of the test.

The stored waveforms are uncalibrated and in the primary-mass (m1) convention, which
is what the model produced before mass_scale was introduced in 0.2.0.
"""

import argparse
import subprocess
import warnings
from pathlib import Path

import numpy as np

from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3

# (q, spin1): both spin signs, spanning the q = 3 to 1000 domain
POINTS = [(3.0, 0.55), (12.5, 0.11), (47.0, -0.42), (115.0, 0.0), (501.0, -0.39)]

MODES = [(2, 2), (2, 1), (3, 1), (3, 2), (3, 3), (4, 2), (4, 3), (4, 4)]

# the full time array is ~136k samples per mode, which would make the file ~90 MB;
# 500 matches the 1D reference data
N_SAMPLES = 500


def tag_for(q, spin1):
    """Key for a (q, spin1) point, e.g. q12p5_chim0p42. Must match the test."""
    fmt = lambda x: ('m' if x < 0 else '') + ('%g' % abs(x)).replace('.', 'p')
    return 'q%s_chi%s' % (fmt(q), fmt(spin1))


def source_version():
    """Short commit of the checkout being used, for provenance."""
    repo = Path(BHPTNRSur2dq1e3.__file__).resolve().parents[2]
    try:
        out = subprocess.run(['git', '-C', str(repo), 'describe', '--always', '--dirty'],
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return 'unknown'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default="regression_data_2dq1e3.npz",
                        help="output .npz file")
    args = parser.parse_args()

    print("model  :", Path(BHPTNRSur2dq1e3.__file__).resolve())
    print("version:", source_version())

    data = {}
    for q, spin1 in POINTS:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            t, h = BHPTNRSur2dq1e3.generate_surrogate(
                q=q, spin1=spin1, modes=MODES, neg_modes=False,
                calibrated=False, mass_scale='m1',
            )

        idx = np.linspace(0, len(t) - 1, N_SAMPLES).astype(int)
        tag = tag_for(q, spin1)
        data['t_%s' % tag] = np.asarray(t)[idx]
        for mode in MODES:
            data['h_%s_l%d_m%d' % (tag, mode[0], mode[1])] = np.asarray(h[mode])[idx]

        print("  q=%-6g spin1=%-6g  max|h22|=%.6e" % (q, spin1, np.abs(h[(2, 2)]).max()))

    data['points'] = np.array(POINTS, dtype=float)
    data['modes'] = np.array(MODES, dtype=int)
    data['source_version'] = np.array(source_version())

    out = Path(args.out).resolve()
    np.savez_compressed(out, **data)
    print("wrote %s (%.1f KB)" % (out, out.stat().st_size / 1024))


if __name__ == "__main__":
    main()
