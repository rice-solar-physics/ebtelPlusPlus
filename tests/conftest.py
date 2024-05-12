import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--ebtel_idl_path", action="store", default=None, help="Path to EBTEL IDL code"
    )
    parser.addoption(
        '--plot_idl_comparisons', action='store_true', default=False, help='Plot IDL results against ebtel++'
    )


@pytest.fixture
def ebtel_idl_path(request):
    return request.config.getoption("--ebtel_idl_path")


@pytest.fixture
def plot_idl_comparisons(request):
    return request.config.getoption("--plot_idl_comparisons")
