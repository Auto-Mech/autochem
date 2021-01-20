""" geometry constructor
"""

import numpy
from phydat import phycon, ptab


def from_data(symbols, coordinates, angstrom=False):
    """ Build a geometry data structure from atomic symbols and coordinates.

        format:
            geo = (((sym1, (xcoord1, ycoord1, zcoord1)),...
                   ((symn, (xcoordn, ycoordn, zcoordn)))

        :param symbols: atomic symbols of the atoms
        :type symbols: tuple(str)
        :param coordinates: xyz coordinates of the atoms
        :type coordinates: tuple(float)
        :param angstrom: convert coordinates to Angstrom units
        :type angstrom: bool
    """

    symbs = list(map(ptab.to_symbol, symbols))
    natms = len(symbs)

    xyzs = numpy.array(coordinates, dtype=float)
    assert numpy.ndim(xyzs) == 2 and numpy.shape(xyzs) == (natms, 3)
    xyzs = (xyzs if not angstrom else
            numpy.multiply(xyzs, phycon.ANG2BOHR))
    xyzs = list(map(tuple, xyzs))
    geo = tuple(zip(symbs, xyzs))

    return geo
