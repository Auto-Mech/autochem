""" I/O for z-matrices
"""
import autoread as ar
import autowrite as aw
from automol.constructors.zmatrix import from_data as _from_data
from ._core import symbols as _symbols
from ._core import key_matrix as _key_matrix
from ._core import name_matrix as _name_matrix
from ._core import values as _values


def from_zmat_string(zma_str):
    """ read a z-matrix from a .zmat string
    """
    syms, key_mat, name_mat, val_dct = ar.zmatrix.read(zma_str)
    zma = _from_data(syms, key_mat, name_mat, val_dct,
                     one_indexed=True, angstrom=True, degree=True)
    return zma


def zmat_string(zma):
    """ write a z-matrix to a .zmat string
    """
    zma_str = aw.zmatrix.write(
        syms=_symbols(zma),
        key_mat=_key_matrix(zma, one_indexed=True),
        name_mat=_name_matrix(zma),
        val_dct=_values(zma, angstrom=True, degree=True)
    )
    return zma_str


def matrix_block_string(zma):
    """ write a z-matrix to a .zmat string
    """
    mat_str = aw.zmatrix.matrix_block(
        syms=_symbols(zma),
        key_mat=_key_matrix(zma, one_indexed=True),
        name_mat=_name_matrix(zma),
    )
    return mat_str


def setval_block_string(zma):
    """ write a z-matrix to a .zmat string
    """
    zma_str = aw.zmatrix.setval_block(
        val_dct=_values(zma, angstrom=True, degree=True)
    )
    return zma_str
