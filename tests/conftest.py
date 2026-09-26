import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="Run slow tests that require h5 data download",
    )
    parser.addoption(
        "--run-gwsurrogate",
        action="store_true",
        default=False,
        help="Run gwsurrogate comparison tests",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: requires h5 data download (deselect with '-m \"not slow\"')")
    config.addinivalue_line("markers", "gwsurrogate: requires gwsurrogate package and h5 data")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--run-slow"):
        skip_slow = pytest.mark.skip(reason="need --run-slow option to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)
    if not config.getoption("--run-gwsurrogate"):
        skip_gws = pytest.mark.skip(reason="need --run-gwsurrogate option to run")
        for item in items:
            if "gwsurrogate" in item.keywords:
                item.add_marker(skip_gws)


@pytest.fixture(scope="module")
def model_1d():
    """Load the 1D surrogate model (triggers h5 download)."""
    from BHPTNRSurrogate.surrogates import BHPTNRSur1dq1e4
    BHPTNRSur1dq1e4._ensure_loaded()
    return BHPTNRSur1dq1e4


@pytest.fixture(scope="module")
def model_2d():
    """Load the 2D surrogate model (triggers h5 download).

    Requires the eval_pysur submodule used for GPR fits.
    """
    from BHPTNRSurrogate.surrogates.common_utils.eval_pysur import evaluate_fit  # noqa: F401
    from BHPTNRSurrogate.surrogates import BHPTNRSur2dq1e3
    BHPTNRSur2dq1e3._ensure_loaded()
    return BHPTNRSur2dq1e3
