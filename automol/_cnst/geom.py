""" geometry constructor
"""
import numpy
from .. import atom as _atom
from .. import _units


def from_data(symbols, coordinates, angstroms=False):
    """ geometry data structure from symbols and coordinates
    """
    syms = list(map(_atom.standard_case, symbols))
    assert all(sym in _atom.SYMBOLS for sym in syms)
    natms = len(syms)

    xyzs = numpy.array(coordinates, dtype=float)
    assert numpy.ndim(xyzs) == 2 and numpy.shape(xyzs) == (natms, 3)
    xyzs = xyzs if not angstroms else numpy.multiply(xyzs, _units.ANG2BOHR)
    xyzs = list(map(tuple, xyzs))
    geo = tuple(zip(syms, xyzs))
    return geo
