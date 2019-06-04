""" geometry constructor
"""
import numpy
from qcelemental import periodictable as pt
from qcelemental import constants as qcc


def from_data(symbols, coordinates, angstrom=False):
    """ geometry data structure from symbols and coordinates
    """
    syms = list(map(pt.to_E, symbols))
    natms = len(syms)

    xyzs = numpy.array(coordinates, dtype=float)
    assert numpy.ndim(xyzs) == 2 and numpy.shape(xyzs) == (natms, 3)
    xyzs = (xyzs if not angstrom else
            numpy.multiply(xyzs, qcc.conversion_factor('angstrom', 'bohr')))
    xyzs = list(map(tuple, xyzs))
    geo = tuple(zip(syms, xyzs))
    return geo
