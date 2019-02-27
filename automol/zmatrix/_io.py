""" I/O for z-matrices
"""
from ..readers import zmatrix as _zmatrix_reader
from ..writers import zmatrix as _zmatrix_writer
from ._core import symbols as _symbols
from ._core import key_matrix as _key_matrix
from ._core import name_matrix as _name_matrix
from ._core import values as _values


def from_zmat_string(zma_str):
    """ read a z-matrix from a .zmat string
    """
    zma = _zmatrix_reader.from_string(zma_str)
    return zma


def zmat_string(zma):
    """ write a z-matrix to a .zmat string
    """
    zma_str = _zmatrix_writer.string(
        syms=_symbols(zma),
        key_mat=_key_matrix(zma, one_indexed=True),
        name_mat=_name_matrix(zma),
        val_dct=_values(zma, angstrom=True, degree=True)
    )
    return zma_str


def matrix_block_string(zma):
    """ write a z-matrix to a .zmat string
    """
    mat_str = _zmatrix_writer.matrix_block_string(
        syms=_symbols(zma),
        key_mat=_key_matrix(zma, one_indexed=True),
        name_mat=_name_matrix(zma),
    )
    return mat_str


def setval_block_string(zma):
    """ write a z-matrix to a .zmat string
    """
    zma_str = _zmatrix_writer.setval_block_string(
        val_dct=_values(zma, angstrom=True, degree=True)
    )
    return zma_str
