""" some math functions
"""
import numpy
from ._core import symbols as _symbols
from ._core import coordinates as _coordinates


def almost_equal(geo1, geo2):
    """ are these geometries almost equal?
    """
    ret = False
    if _symbols(geo1) == _symbols(geo2):
        ret = numpy.allclose(_coordinates(geo1), _coordinates(geo2))
    return ret
