""" pyx2z interface
"""
import pyx2z
import autoparse.pattern as app
import autoparse.find as apf
from ..readers.zmatrix import zmatrix as _zmatrix_reader
from ..readers.zmatrix import (matrix_block_capturing_pattern as
                               _zmatrix_matrix_block_capturing_pattern)
from ..readers.zmatrix import (setval_block_capturing_pattern as
                               _zmatrix_setval_block_capturing_pattern)


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
    mat_block_pattern = _zmatrix_matrix_block_capturing_pattern(
        delim_pattern=',')
    setval_block_pattern = _zmatrix_setval_block_capturing_pattern(
        delim_pattern=app.one_of_these([app.LINESPACE, app.NEWLINE]))

    zma_str = pyx2z.zmatrix_string(x2m)
    mat_str = apf.first_capture(mat_block_pattern, zma_str)
    setval_str = apf.first_capture(setval_block_pattern, zma_str)

    zma = _zmatrix_reader(mat_str, setval_str, mat_delim_pattern=',',
                          one_indexed=True, angstrom=False, degree=True)
    return zma
