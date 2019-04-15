""" I/O for z-matrices
"""
import autoread as ar
import autowrite as aw
from automol.zmatrix import v as _v_
from automol.zmatrix._core import var_ as _var_
from automol.zmatrix._core import values as _values
from automol.constructors.zmatrix import from_data as _from_data


def from_string(zma_str):
    """ read a z-matrix from a string
    """
    syms, key_mat, name_mat, val_dct = ar.zmatrix.read(zma_str)

    zma = _from_data(syms, key_mat, name_mat, val_dct,
                     one_indexed=True, angstrom=True, degree=True)
    return zma


def string(zma):
    """ write a z-matrix to a string
    """
    zma_str = aw.zmatrix.write(
        syms=_v_.symbols(_var_(zma)),
        key_mat=_v_.key_matrix(_var_(zma), shift=1),
        name_mat=_v_.name_matrix(_var_(zma)),
        val_dct=_values(zma, angstrom=True, degree=True)
    )
    return zma_str
