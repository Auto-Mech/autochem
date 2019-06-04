""" z-matrix constructor
"""
import numpy
from qcelemental import constants as qcc
from automol.create.vmatrix import from_data as _v_from_data


def from_data(symbols, key_matrix, name_matrix, values,
              one_indexed=False, angstrom=False, degree=False):
    """ z-matrix constructor

    :param symbols: atomic symbols
    :type symbols: tuple[str]
    :param key_matrix: key/index columns of the z-matrix, zero-indexed
    :type key_matrix: tuple[tuple[float, float or None, float or None]]
    :param name_matrix: coordinate name columns of the z-matrix
    :type name_matrix; tuple[tuple[str, str or None, str or None]]
    :param values: coordinate values, by coordinate name
    :type values: dict
    """
    vma = _v_from_data(symbols, key_matrix, name_matrix, one_indexed)
    val_dct = _values(values, name_matrix, angstrom, degree)

    zma = (vma, val_dct)
    return zma


def _values(val_dct, name_mat, angstrom, degree):
    ret_val_dct = {}

    name_set = set(numpy.ravel(name_mat)) - {None}
    assert set(val_dct) == name_set

    for name, val in val_dct.items():
        # make sure every entry in the value dictionary corresponds to a
        # coordinate
        assert numpy.any(numpy.equal(name_mat, name))

        # make sure coordinates with the same name are in the same column
        col_idxs = numpy.where(numpy.equal(name_mat, name))[1]
        col_idx_set = set(col_idxs)
        assert len(col_idx_set) == 1

        # get the column index
        col_idx, = col_idx_set

        # convert units according to which column the values come from
        if col_idx == 0:
            if angstrom:
                val *= qcc.conversion_factor('angstrom', 'bohr')
        else:
            if degree:
                val *= qcc.conversion_factor('degree', 'radian')

        ret_val_dct[name] = val

    return ret_val_dct
