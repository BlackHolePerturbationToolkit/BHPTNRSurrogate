"""Smoke tests: verify that the package and its submodules are importable."""


def test_import_top_level():
    import BHPTNRSurrogate  # noqa: F401


def test_import_1d_model():
    from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4  # noqa: F401


def test_import_2d_model():
    from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3  # noqa: F401


def test_import_check_inputs():
    from BHPTNRSurrogate.surrogates.common_utils import check_inputs  # noqa: F401


def test_generate_surrogate_callable_1d():
    from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4
    assert callable(BHPTNRSur1dq1e4.generate_surrogate)


def test_generate_surrogate_callable_2d():
    from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
    assert callable(BHPTNRSur2dq1e3.generate_surrogate)
