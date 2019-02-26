""" I/O for z-matrices
"""
import autoparse.find as apf
from ..readers.zmatrix import zmatrix as _zmatrix_reader
from ..readers.zmatrix import matrix_block_capturing_pattern
from ..readers.zmatrix import setval_block_capturing_pattern
from ._core import symbols as _symbols
from ._core import key_matrix as _key_matrix
from ._core import name_matrix as _name_matrix
from ._core import values as _values


def from_zmat_string(zma_str):
    """ read a z-matrix from a .zmat string
    """
    mat_str = apf.first_capture(matrix_block_capturing_pattern(), zma_str)
    setval_str = apf.first_capture(setval_block_capturing_pattern(), zma_str)
    zma = _zmatrix_reader(mat_str, setval_str)
    return zma


def zmat_string(zma):
    """ write a z-matrix to a .zmat string
    """
    mat_str = matrix_block_string(zma)
    setval_str = setval_block_string(zma)
    zma_str = '\n'.join((mat_str, setval_str))
    return zma_str


def matrix_block_string(zma, delim=' '):
    """ write the .zmat matrix block to a string
    """
    syms = _symbols(zma)
    key_mat = _key_matrix(zma, one_indexed=True)
    name_mat = _name_matrix(zma)

    mat_str = ''
    for row_idx, sym in enumerate(syms):
        mat_str += '{:<2s}'.format(sym)
        for col_idx in range(min(row_idx, 3)):
            mat_str += '{}{:>d}'.format(delim, key_mat[row_idx][col_idx])
            mat_str += '{}{:>5s}'.format(delim, name_mat[row_idx][col_idx])
        mat_str += '\n'

    return mat_str


def setval_block_string(zma, setval_sign='='):
    """ write the .zmat setval block to a string
    """
    val_dct = _values(zma, angstrom=True, degree=True)

    setval_str = ''
    for name, val in val_dct.items():
        setval_str += '{:<5s} {} {:>11.6f}\n'.format(name, setval_sign, val)

    return setval_str
