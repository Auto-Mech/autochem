""" zmatrix constructor
"""
import numpy
import phycon.units as pcu
import phycon.elements as pce


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
    syms = list(map(pce.standard_case, symbols))
    assert all(sym in pce.element_keys() for sym in syms)
    natms = len(syms)

    key_mat = _key_matrix(key_matrix, natms, one_indexed)
    name_mat = _name_matrix(name_matrix, natms)
    val_dct = _values(values, name_mat, angstrom, degree)

    mat = tuple(zip(syms, key_mat, name_mat))
    zma = (mat, val_dct)
    return zma


def _key_matrix(key_mat, natms, one_indexed):
    key_mat = numpy.array(key_mat, dtype=numpy.object_)
    assert key_mat.ndim == 2 and key_mat.shape == (natms, 3)
    triu_idxs = numpy.triu_indices(natms, m=3)

    # check the key matrix and make it one-indexed
    key_mat[triu_idxs] = -1
    key_mat = key_mat.astype(int)
    key_mat -= 1 if one_indexed else 0
    key_mat = key_mat.astype(numpy.object_)
    key_mat[triu_idxs] = None

    return tuple(map(tuple, key_mat))


def _name_matrix(name_mat, natms):
    name_mat = numpy.array(name_mat, dtype=numpy.object_)
    assert name_mat.ndim == 2 and name_mat.shape == (natms, 3)
    natms = name_mat.shape[0]
    triu_idxs = numpy.triu_indices(natms, m=3)
    tril_idxs = numpy.tril_indices(natms, -1, m=3)

    assert all(isinstance(name, str) for name in name_mat[tril_idxs])
    name_mat[triu_idxs] = None

    return tuple(map(tuple, name_mat))


def _values(val_dct, name_mat, angstrom, degree):
    ret_val_dct = {}

    # makes sure the value dictionary is complete
    assert all(name in val_dct or name is None
               for name in numpy.ravel(name_mat))

    for name, val in val_dct.items():
        # make sure every entry in the value dictionary corresponds to a
        # coordinate
        assert numpy.any(numpy.equal(name_mat, name))

        # make sure coordinates with the same name are in the same column
        _, col_idxs = numpy.where(numpy.equal(name_mat, name))
        col_idx_set = set(col_idxs)
        assert len(col_idx_set) == 1

        # get the column index
        col_idx, = col_idx_set

        # convert units according to which column the values come from
        if col_idx == 0:
            val *= (pcu.ANG2BOHR if angstrom else 1.)
        else:
            val *= (pcu.DEG2RAD if degree else 1.)

        ret_val_dct[name] = val

    return ret_val_dct
