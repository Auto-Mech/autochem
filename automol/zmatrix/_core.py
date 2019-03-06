""" core library defining the z-matrix data structure
"""
import numpy
import phycon.units as pcu
from ..constructors.zmatrix import from_data as _from_data


# value getters
def matrix(zma):
    """ the matrix (atom symbols, atom keys, and coordinate names)
    """
    mat, _ = zma
    return mat


def values(zma, angstrom=False, degree=False):
    """ coordinate values, by coordinate name
    """
    _, val_dct = zma
    return _values(val_dct, zma, angstrom, degree)


def _values(val_dct, zma, angstrom, degree):
    """ values post-processing
    """
    dist_names = distance_names(zma)
    ang_names = angle_names(zma) + dihedral_names(zma)
    orig_val_dct = val_dct

    val_dct = {}
    for name, val in orig_val_dct.items():
        if angstrom and name in dist_names:
            val *= pcu.BOHR2ANG
        if degree and name in ang_names:
            val *= pcu.RAD2DEG
        val_dct[name] = val
    return val_dct


def symbols(zma):
    """ atomic symbols, by z-matrix row
    """
    mat = matrix(zma)
    if mat:
        syms, _, _ = zip(*mat)
    else:
        syms = ()
    return syms


def key_matrix(zma, one_indexed=False):
    """ coordinate atom keys, by z-matrix row and column
    """
    mat = matrix(zma)
    if mat:
        _, key_mat, _ = zip(*mat)
    else:
        key_mat = ()

    return _key_matrix(key_mat, one_indexed)


def _key_matrix(key_mat, one_indexed):
    """ key matrix post-processing
    """
    if one_indexed:
        key_mat = numpy.array(key_mat)
        tril_idxs = numpy.tril_indices(key_mat.shape[0], -1, m=3)
        key_mat[tril_idxs] += 1
    return tuple(map(tuple, key_mat))


def name_matrix(zma):
    """ coordinate names, by z-matrix row and column
    """
    mat = matrix(zma)
    if mat:
        _, _, name_mat = zip(*mat)
    else:
        name_mat = ()
    return name_mat


def value_matrix(zma):
    """ coordinate values, by z-matrix row and column
    """
    val_dct = values(zma)
    name_mat = name_matrix(zma)
    val_mat = tuple(tuple(val_dct[name] if name is not None else None
                          for name in name_mat_row)
                    for name_mat_row in name_mat)
    return val_mat


def coordinate_key_matrix(zma):
    """ coordinate keys, by z-matrix row and column
    """
    key_mat = key_matrix(zma)
    coo_key_mat = tuple(
        tuple((key,) + key_row[:col+1] if key_row[col] is not None else None
              for col in range(3))
        for key, key_row in enumerate(key_mat))
    return coo_key_mat


def distance_names(zma):
    """ distance coordinate names, from top to bottom (no repeats)
    """
    name_mat = numpy.array(name_matrix(zma))
    return tuple(name_mat[1:, 0])


def angle_names(zma):
    """ angle coordinate names, from top to bottom
    """
    name_mat = numpy.array(name_matrix(zma))
    return tuple(name_mat[2:, 1])


def dihedral_names(zma):
    """ dihedral coordinate names, from top to bottom
    """
    name_mat = numpy.array(name_matrix(zma))
    return tuple(name_mat[3:, 2])


def distance_keys(zma):
    """ distance coordinate keys, from top to bottom (no repeats)
    """
    coo_key_mat = numpy.array(coordinate_key_matrix(zma))
    return tuple(coo_key_mat[1:, 0])


def angle_keys(zma):
    """ angle coordinate keys, from top to bottom
    """
    coo_key_mat = numpy.array(coordinate_key_matrix(zma))
    return tuple(coo_key_mat[2:, 1])


def dihedral_keys(zma):
    """ dihedral coordinate keys, from top to bottom
    """
    coo_key_mat = numpy.array(coordinate_key_matrix(zma))
    return tuple(coo_key_mat[3:, 2])


# value setters
def set_names(zma, name_dct):
    """ set coordinate names for the z-matrix
    """
    orig_name_mat = numpy.array(name_matrix(zma))
    tril_idxs = numpy.tril_indices(orig_name_mat.shape[0], -1, m=3)
    orig_names = set(orig_name_mat[tril_idxs])
    assert set(name_dct.keys()) <= orig_names

    name_dct.update({orig_name: orig_name for orig_name in orig_names
                     if orig_name not in name_dct})

    name_mat = numpy.empty(orig_name_mat.shape, dtype=numpy.object_)
    name_mat[tril_idxs] = list(map(name_dct.__getitem__,
                                   orig_name_mat[tril_idxs]))

    orig_val_dct = values(zma)
    val_dct = {name_dct[orig_name]: val
               for orig_name, val in orig_val_dct.items()}

    return _from_data(symbols(zma), key_matrix(zma), name_mat, val_dct)


def set_values(zma, val_dct):
    """ set coordinate values for the z-matrix
    """
    new_val_dct = val_dct
    orig_val_dct = values(zma)
    assert set(new_val_dct.keys()) <= set(orig_val_dct.keys())
    val_dct = {}
    val_dct.update(orig_val_dct)
    val_dct.update(new_val_dct)
    return _from_data(symbols(zma), key_matrix(zma), name_matrix(zma), val_dct)
