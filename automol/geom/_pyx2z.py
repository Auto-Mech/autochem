""" pyx2z interface
"""
import pyx2z
import autoparse.pattern as app
from ..readers.zmatrix import from_string as _zmatrix_reader


def from_geometry(geo):
    """ x2z molecule object from a geometry
    """
    _mg = pyx2z.MolecGeom()
    for sym, xyz in geo:
        _atm = pyx2z.Atom(sym)
        _atm[0], _atm[1], _atm[2] = xyz
        _mg.push_back(_atm)
    _ps = pyx2z.PrimStruct(_mg)
    x2m = pyx2z.MolecStruct(_ps)
    return x2m


def to_zmatrix(x2m):
    """ z-matrix from an x2z molecule object
    """
    zma_str = pyx2z.zmatrix_string(x2m)

    zma = _zmatrix_reader(
        zma_str,
        mat_delim_pattern=',',
        setval_delim_pattern=app.one_of_these([app.LINESPACE, app.NEWLINE]),
        one_indexed=True, angstrom=False, degree=True,
    )

    return zma
