""" geometry constructor
"""
import numpy
import phycon.units as pcu
import phycon.elements as pce


def from_data(symbols, coordinates, angstrom=False):
    """ geometry data structure from symbols and coordinates
    """
    syms = list(map(pce.standard_case, symbols))
    assert all(sym in pce.element_keys() for sym in syms)
    natms = len(syms)

    xyzs = numpy.array(coordinates, dtype=float)
    assert numpy.ndim(xyzs) == 2 and numpy.shape(xyzs) == (natms, 3)
    xyzs = xyzs if not angstrom else numpy.multiply(xyzs, pcu.ANG2BOHR)
    xyzs = list(map(tuple, xyzs))
    geo = tuple(zip(syms, xyzs))
    return geo
