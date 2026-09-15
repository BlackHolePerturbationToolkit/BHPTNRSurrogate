# Changelog

## [0.2.0] - 2026-09-14

### Added
- `mass_scale` parameter to `BHPTNRSur1dq1e4` and `BHPTNRSur2dq1e3` surrogate models.
  When `calibrated=False`, allows choosing between total mass `M` (default) and
  primary mass `m1` as the mass convention. The rescaling factor `q/(q+1)` is applied
  to both time and strain when `mass_scale='M'`.
- Local installation instructions (`pip install -e .`) in README.

### Fixed
- H5 data download now prints progress messages to stdout. Previously used
  `logger.info()` with no handler configured, so downloads were silent.
- H5 files are stored in the package's `data/` directory, following the package-local strategy used by `gwsurrogate`.
- Replaced `print()` warnings with `warnings.warn()` throughout surrogate modules.
- Fixed `SyntaxWarning` from invalid escape sequences.
- Removed unused imports across modules.

### Changed
- CI now tests Python 3.8, 3.10, 3.11, and 3.14.
- scikit-learn is now a required dependency because it is needed by the two-dimensional surrogate model.
- Passing `mass_scale='m1'` for a calibrated waveform now raises a `ValueError` instead of being ignored.
- The default mass convention for uncalibrated waveforms changed from the
  historical implicit primary-mass (`m1`) convention to total mass (`M`). Pass
  `mass_scale='m1'` to reproduce the earlier behavior.
- Surrogate data is now lazily loaded on first access instead of at import time.
- Switched from `wget` to `urllib` for h5 file downloads.
- Converted to proper Python package with `pyproject.toml` and relative imports.

## [0.1.0] - 2024-07-01

### Added
- `BHPTNRSur2dq1e3`: spinning surrogate model for mass ratios 3-1000 and
  spins -0.8 to 0.8 on the primary black hole, with modes up to l=4.
- `neg_modes=True` default for 2D surrogate.

## [0.0.1] - 2022-04-01

### Added
- Initial release of `BHPTNRSur1dq1e4`: non-spinning surrogate model for
  mass ratios 2.5-10000 with modes up to l=10, calibrated to NR up to l=5.
- NR calibration framework with alpha/beta scaling parameters.
- Support for geometric and SI unit waveform output.
- Mode summation and sky-location evaluation.
- Automatic h5 data download from Zenodo.
- Tutorial notebooks with NR comparison examples.
