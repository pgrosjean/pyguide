import pytest

def test_import_tiling():
    from pyguide import tiling
    assert hasattr(tiling, 'parse_coordinates')
