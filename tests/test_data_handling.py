"""Tests for surrogate data paths and downloads."""

import hashlib
from pathlib import Path

import pytest

from BHPTNRSurrogate import _data
from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4, BHPTNRSur2dq1e3
from BHPTNRSurrogate.surrogates.common_utils import filehash


def test_environment_override(monkeypatch, tmp_path):
    custom_dir = tmp_path / "custom data"
    monkeypatch.setenv(_data.DATA_DIR_ENV, str(custom_dir))
    assert Path(_data.get_data_dir()) == custom_dir.resolve()


def test_default_data_dir_is_package_local(monkeypatch, tmp_path):
    fake_package_file = tmp_path / "site-packages" / "BHPTNRSurrogate" / "_data.py"
    monkeypatch.delenv(_data.DATA_DIR_ENV, raising=False)
    monkeypatch.setattr(_data, "__file__", str(fake_package_file))
    assert Path(_data.get_data_dir()) == fake_package_file.parent / "data"


def test_models_use_resolved_data_directory():
    assert BHPTNRSur1dq1e4.h5_data_dir == _data.get_data_dir()
    assert BHPTNRSur2dq1e3.h5_data_dir == _data.get_data_dir()


def test_download_creates_directory_and_uses_atomic_replace(monkeypatch, tmp_path):
    payload = b"surrogate data"
    data_dir = tmp_path / "new cache"
    destination_seen = {}

    def fake_urlretrieve(url, destination):
        destination_seen["path"] = Path(destination)
        Path(destination).write_bytes(payload)
        return destination, None

    monkeypatch.setattr(filehash.urllib.request, "urlretrieve", fake_urlretrieve)
    expected_hash = hashlib.md5(payload).hexdigest()
    digest = filehash.md5("model.h5", data_dir, "123", expected_hash)
    final_path = data_dir / "model.h5"

    assert digest == hashlib.md5(payload).hexdigest()
    assert final_path.read_bytes() == payload
    assert destination_seen["path"].parent == data_dir
    assert destination_seen["path"] != final_path
    assert list(data_dir.iterdir()) == [final_path]


def test_bad_download_checksum_is_not_published(monkeypatch, tmp_path):
    data_dir = tmp_path / "new cache"

    def fake_urlretrieve(url, destination):
        Path(destination).write_bytes(b"corrupt")
        return destination, None

    monkeypatch.setattr(filehash.urllib.request, "urlretrieve", fake_urlretrieve)
    with pytest.raises(RuntimeError, match="checksum mismatch"):
        filehash.md5("model.h5", data_dir, "123", hashlib.md5(b"expected").hexdigest())

    assert not (data_dir / "model.h5").exists()
    assert list(data_dir.iterdir()) == []


def test_failed_download_removes_partial_file(monkeypatch, tmp_path):
    data_dir = tmp_path / "new cache"

    def fail_download(url, destination):
        Path(destination).write_bytes(b"partial")
        raise OSError("network failed")

    monkeypatch.setattr(filehash.urllib.request, "urlretrieve", fail_download)
    with pytest.raises(RuntimeError, match="Failed to download"):
        filehash.md5("model.h5", data_dir, "123")

    assert not (data_dir / "model.h5").exists()
    assert list(data_dir.iterdir()) == []
