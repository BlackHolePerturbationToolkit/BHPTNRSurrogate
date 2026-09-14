##==============================================================================
## BHPTNRSurrogate module
## Description : checks hash of h5 files
## Author : Tousif Islam, Aug 2022 [tislam@umassd.edu / tousifislam24@gmail.com]
##==============================================================================

import os
import hashlib
import tempfile
import urllib.request

#----------------------------------------------------------------------------------------------------
def _md5(file_path):
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

#----------------------------------------------------------------------------------------------------
def md5(fname, h5_data_dir, zenodo_ID, expected_hash=None):
    """ Compute hash from file. code taken from
    https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file"""

    file_path = os.path.join(h5_data_dir, fname)

    # download file if not already there
    if not os.path.isfile(file_path):
        url = 'https://zenodo.org/record/%s/files/%s' % (zenodo_ID, fname)
        print('%s not found in %s' % (fname, h5_data_dir))
        print('Downloading %s from Zenodo ... this may take a few minutes ...' % fname)
        temp_path = None
        try:
            os.makedirs(h5_data_dir, exist_ok=True)
            file_descriptor, temp_path = tempfile.mkstemp(prefix='.%s.' % fname, suffix='.tmp', dir=h5_data_dir)
            os.close(file_descriptor)
            urllib.request.urlretrieve(url, temp_path)
            downloaded_hash = _md5(temp_path)
            if expected_hash is not None and downloaded_hash != expected_hash:
                raise ValueError("checksum mismatch (expected %s, received %s)" % (expected_hash, downloaded_hash))
            os.replace(temp_path, file_path)
        except Exception as e:
            raise RuntimeError("Failed to download %s from %s: %s" % (fname, url, e)) from e
        finally:
            if temp_path is not None:
                try:
                    os.remove(temp_path)
                except FileNotFoundError:
                    pass
        print('Download complete: %s' % fname)
        return downloaded_hash

    return _md5(file_path)

#----------------------------------------------------------------------------------------------------
def check_current_hash(file_hash, zenodo_current_hash, url, fname):
    # chech if the h5file is the most recent one or if it is corrupted
    if file_hash != zenodo_current_hash:
        raise AttributeError("%s out of date.\n \
                             Please download new version from %s"\
                             %(fname,url))
    else:
        pass
