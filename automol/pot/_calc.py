"""
  Various info objects for torsions
"""

import numpy
import automol


def grid(zma, tors_name, scan_increment, sym, from_equil=False):
    """ scan grids
    """

    _grid = _scan_linspace(zma, name, scan_increment, sym)

    if from_equil:
        val_dct = automol.zmat.values(zma)
        ini_val = val_dct[tors_name]
        _grid = tuple(val.item() + ini_val for val in _grid)

    return grid_from_equil


def _scan_linspace(zma, coord_name, scan_increment, sym):
    """ scan grids for torsional dihedrals
    """

    interval = ((2.0 * numpy.pi) / sym) - scan_increment
    npoints = int(interval / scan_increment) + 2
    linspace = numpy.linspace(0.0, interval, npoints)

    return tuple(val.item() for val in linspace)
