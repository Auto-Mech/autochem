""" I/O for variable z-matrices
"""
import autoread as ar
import autowrite as aw
from automol.constructors.vmatrix import from_data as _from_data
from automol.zmatrix.v._core import symbols as _symbols
from automol.zmatrix.v._core import key_matrix as _key_matrix
from automol.zmatrix.v._core import name_matrix as _name_matrix


def from_string(vma_str):
    """ read a z-matrix from a string
    """
    syms, key_mat, name_mat = ar.zmatrix.matrix.read(vma_str)

    vma = _from_data(syms, key_mat, name_mat, one_indexed=True)
    return vma


def string(vma):
    """ write a z-matrix to a string
    """
    vma_str = aw.zmatrix.matrix_block(
        syms=_symbols(vma),
        key_mat=_key_matrix(vma, one_indexed=True),
        name_mat=_name_matrix(vma),
    )
    return vma_str
