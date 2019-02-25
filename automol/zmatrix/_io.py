""" I/O for z-matrices
"""
import numpy
import phycon.units as pcu
from ._core import from_matrices as _from_matrices
from ._core import symbols as _symbols
from ._core import key_matrix as _key_matrix
from ._core import value_matrix as _value_matrix
from ._core import coordinate_matrix as _coordinate_matrix
from .parse import matrix_block as _matrix_block
from .parse import matrix_symbols as _matrix_symbols
from .parse import variable_block as _variable_block
from .parse import variable_values as _variable_values
from .parse import matrix_keys_and_entries as _matrix_keys_and_entries


def from_zmat_string(zma_str):
    """ read a z-matrix from a .zmat string
    """
    mat_str = _matrix_block(zma_str)
    var_str = _variable_block(zma_str)
    val_dct = _variable_values(var_str)
    syms = _matrix_symbols(mat_str)
    key_mat, val_mat = _matrix_keys_and_entries(mat_str)
    key_mat = numpy.array(key_mat, dtype=numpy.object_)
    val_mat = numpy.array(val_mat, dtype=numpy.object_)
    var_dct = {}
    for col_idx in range(3):
        for row_idx in range(col_idx+1, len(syms)):
            key_mat[row_idx, col_idx] -= 1
            coo = (row_idx,) + tuple(key_mat[row_idx, :col_idx+1])
            val = val_mat[row_idx, col_idx]
            if val in val_dct:
                var_dct[coo] = val
                val = val_dct[val]
            val_mat[row_idx, col_idx] = val * (pcu.ANG2BOHR if col_idx == 0
                                               else pcu.DEG2RAD)

    zma = _from_matrices(syms, key_mat, val_mat)
    return zma, var_dct


def zmat_string(zma, var_dct=None, comment=''):
    """ write z-matrix to a .zmat string

    var_dct maps coordinates such as (0, 2, 4, 5) onto string variable names
    if var_dct is set, the variable block is returned as a separate string
    """
    var_dct = {} if var_dct is None else var_dct
    mat_str = zmat_string_matrix_block(zma, var_dct)
    var_str = zmat_string_variable_block(zma, var_dct)
    zma_str = '\n'.join((comment, mat_str, var_str))
    return zma_str


def zmat_string_matrix_block(zma, var_dct):
    """ write the .zmat matrix block to a string
    """
    syms = _symbols(zma)
    key_mat = _key_matrix(zma, one_indexed=True)
    val_mat = _value_matrix(zma, angstroms=True)
    coo_mat = _coordinate_matrix(zma)

    mat_str = ''
    for row_idx, sym in enumerate(syms):
        mat_str += '{:<2s}'.format(sym)
        for col_idx in range(min(row_idx, 3)):
            key = key_mat[row_idx][col_idx]
            mat_str += ' {:<2d}'.format(key)

            val = val_mat[row_idx][col_idx]
            coo = coo_mat[row_idx][col_idx]
            if coo in var_dct:
                mat_str += ' {:10s}'.format(var_dct[coo])
            else:
                mat_str += ' {:10.6f}'.format(val)
        mat_str += '\n'
    return mat_str


def zmat_string_variable_block(zma, var_dct):
    """ write the .zmat matrix block to a string
    """
    val_mat = _value_matrix(zma, angstroms=True)
    coo_mat = _coordinate_matrix(zma)
    natms = len(val_mat)

    val_dct = {}
    for row_idx in range(natms):
        for col_idx in range(min(row_idx, 3)):
            val = val_mat[row_idx][col_idx]
            coo = coo_mat[row_idx][col_idx]
            if coo in var_dct:
                val_dct[var_dct[coo]] = val

    var_str = ''
    for var_name, val in val_dct.items():
        var_str += '{:10s} = {:10.6f}\n'.format(var_name, val)

    return var_str
