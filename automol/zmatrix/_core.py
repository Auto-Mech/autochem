""" core library defining the zmatrix data structure
"""
import numpy
import phycon.units as pcu
from .._cnst.zmatrix import from_data as _from_data


def from_matrices(syms, key_mat, val_mat):
    """ build a z-matrix from key and value matrices
    """
    assert len(syms) == len(key_mat) == len(val_mat)
    zma = tuple(zip(syms, key_mat, val_mat))
    return _from_data(symbols(zma), distance_column(zma), angle_column(zma),
                      torsion_column(zma))


def symbols(zma):
    """ atomic symbols
    """
    if zma:
        syms, _, _ = zip(*zma)
    else:
        syms = ()
    return syms


def _key_matrix(zma):
    """ atom keys, by row and column
    """
    if zma:
        _, key_mat, _ = zip(*zma)
    else:
        key_mat = ()
    return key_mat


def _value_matrix(zma):
    """ coordinate values, by row and column
    """
    if zma:
        _, _, val_mat = zip(*zma)
    else:
        val_mat = ()
    return val_mat


def distance_column(zma):
    """ the z-matrix distance column, keys and values
    """
    return _column(zma, 0)


def angle_column(zma):
    """ the z-matrix angle column, keys and values
    """
    return _column(zma, 1)


def torsion_column(zma):
    """ the z-matrix torsion column, keys and values
    """
    return _column(zma, 2)


def coordinate_matrix(zma):
    """ matrix of z-matrix coordinates
    """
    key_mat = _object_array(_key_matrix(zma))
    coo_mat = numpy.empty_like(key_mat)
    for row_idx, col_idx in iter_(zma):
        coo = (row_idx,) + tuple(key_mat[row_idx][:col_idx+1])
        coo_mat[row_idx, col_idx] = coo
    return tuple(map(tuple, coo_mat))


def key_matrix(zma, one_indexed=False):
    """ matrix of z-matrix keys
    """
    key_mat = _object_array(_key_matrix(zma))
    if one_indexed:
        key_mat[1:, 0] += 1
        key_mat[2:, 1] += 1
        key_mat[3:, 2] += 1
    return tuple(map(tuple, key_mat))


def value_matrix(zma, angstroms=False):
    """ matrix of z-matrix coordinate values
    """
    val_mat = _object_array(_value_matrix(zma))
    if angstroms:
        val_mat[1:, 0] *= pcu.BOHR2ANG
        val_mat[2:, 1] *= pcu.RAD2DEG
        val_mat[3:, 2] *= pcu.RAD2DEG
    return tuple(map(tuple, val_mat))


def iter_(zma):
    """ iterate over valid z-matrix row and column indices
    """
    natms = len(symbols(zma))
    for row_idx in range(natms):
        for col_idx in range(min(row_idx, 3)):
            yield (row_idx, col_idx)


def iter_column_(zma, col_idx):
    """ iterate over valid z-matrix row indices in a given column
    """
    natms = len(symbols(zma))
    for row_idx in range(col_idx+1, natms):
        yield row_idx


def _column(zma, col_idx):
    key_mat = _object_array(_key_matrix(zma))
    val_mat = _object_array(_value_matrix(zma))
    return tuple((key_mat[row_idx, col_idx], val_mat[row_idx, col_idx])
                 for row_idx in iter_column_(zma, col_idx))


def _object_array(arr):
    shape = numpy.shape(arr) if numpy.size(arr) else (0, 3)
    return numpy.array(arr, dtype=numpy.object_).reshape(shape)
