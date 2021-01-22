""" z-matrix constructor
"""

import numpy
import automol.create.vmat
from phydat import phycon


def from_data(symbols, key_matrix, value_matrix, name_matrix=None,
              one_indexed=False, angstrom=False, degree=False):
    """ Z-Matrix constructor

        :param symbols: atomic symbols
        :type symbols: tuple[str]
        :param key_matrix: key/index columns of the z-matrix, zero-indexed
        :type key_matrix: tuple[tuple[int, int or None, int or None]]
        :param value_matrix: z-matrix coordinate values
        :type value_matrix: tuple[tuple[float, float or None, float or None]]
        :param name_matrix: coordinate name columns of the z-matrix
        :type name_matrix; tuple[tuple[str, str or None, str or None]]
    """
    vma = automol.create.vmat.from_data(
        symbols, key_matrix, name_matrix, one_indexed)
    val_mat = _value_matrix(value_matrix, angstrom, degree)

    syms, key_mat, name_mat = zip(*vma)
    zma = tuple(zip(syms, key_mat, name_mat, val_mat))

    return zma


def _value_matrix(val_mat, angstrom, degree):
    """ Build value matrix of the V-Matrix that contains the
        coordinate keys by row and column.

        :param val_mat: value matrix containing coordinate values
        :type val_mat: tuple(tuple(int))
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
        :rtype: tuple(tuple(str))
    """

    # Check dimensions and ensure proper formatting
    val_mat = [list(row) + [None]*(3-len(row)) for row in val_mat]
    val_mat = numpy.array(val_mat, dtype=numpy.object_)
    natms = val_mat.shape[0]

    assert val_mat.ndim == 2 and val_mat.shape == (natms, 3)
    triu_idxs = numpy.triu_indices(natms, m=3)

    val_mat[1:, 0] *= phycon.ANG2BOHR if angstrom else 1
    val_mat[2:, 1] *= phycon.DEG2RAD if degree else 1
    val_mat[3:, 2] *= phycon.DEG2RAD if degree else 1

    val_mat[triu_idxs] = None

    return tuple(map(tuple, val_mat))
