"""Locate downloaded surrogate data files."""

import os
from pathlib import Path


DATA_DIR_ENV = 'BHPTNR_SURROGATE_DATA_DIR'


def get_data_dir():
    """Return the package-local surrogate data directory."""
    override = os.environ.get(DATA_DIR_ENV)
    if override:
        return str(Path(override).expanduser().resolve())

    return str(Path(__file__).resolve().parent / 'data')
