""" core library defining the z-matrix data structure
"""
import itertools
import more_itertools
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
    ang_names = central_angle_names(zma) + dihedral_angle_names(zma)
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


def coordinate_key_matrix(zma, one_indexed=False):
    """ coordinate keys, by z-matrix row and column
    """
    key_mat = key_matrix(zma, one_indexed=one_indexed)
    natms = len(key_mat)
    atm_keys = range(natms) if not one_indexed else range(1, natms+1)
    coo_key_mat = [[(atm_key,) + key_row[:col+1]
                    if key_row[col] is not None else None for col in range(3)]
                   for atm_key, key_row in zip(atm_keys, key_mat)]
    return tuple(map(tuple, coo_key_mat))


def coordinates(zma, one_indexed=False):
    """ coordinate keys associated with each coordinate name, as a dictionary

    (the values are sequences of coordinate keys, since there may be multiple)
    """
    _names = numpy.ravel(name_matrix(zma))
    coo_keys = numpy.ravel(coordinate_key_matrix(zma, one_indexed))

    coo_dct = {name: () for name in _names}
    for name, coo_key in zip(_names, coo_keys):
        coo_dct[name] += (coo_key,)

    coo_dct.pop(None)
    return coo_dct


def names(zma):
    """ coordinate names
    """
    _names = filter(lambda x: x is not None,
                    numpy.ravel(numpy.transpose(name_matrix(zma))))
    return tuple(more_itertools.unique_everseen(_names))


def value_names(zma):
    """ value coordinate names
    """
    val_dct = values(zma)
    return tuple(name for name in names(zma) if name in val_dct)


def variable_names(zma):
    """ variable coordinate names
    """
    val_dct = values(zma)
    return tuple(name for name in names(zma) if name not in val_dct)


def distance_names(zma):
    """ distance coordinate names
    """
    name_mat = numpy.array(name_matrix(zma))
    return tuple(more_itertools.unique_everseen(name_mat[1:, 0]))


def distance_value_names(zma):
    """ distance value coordinate names
    """
    val_dct = values(zma)
    return tuple(name for name in distance_names(zma) if name in val_dct)


def distance_variable_names(zma):
    """ distance variable coordinate names
    """
    val_dct = values(zma)
    return tuple(name for name in distance_names(zma) if name not in val_dct)


def central_angle_names(zma):
    """ central angle coordinate names, from top to bottom
    """
    name_mat = numpy.array(name_matrix(zma))
    return tuple(more_itertools.unique_everseen(name_mat[2:, 1]))


def central_angle_value_names(zma):
    """ central angle value coordinate names
    """
    val_dct = values(zma)
    return tuple(name for name in central_angle_names(zma)
                 if name in val_dct)


def central_angle_variable_names(zma):
    """ central angle variable coordinate names
    """
    val_dct = values(zma)
    return tuple(name for name in central_angle_names(zma)
                 if name not in val_dct)


def dihedral_angle_names(zma):
    """ dihedral angle coordinate names, from top to bottom
    """
    name_mat = numpy.array(name_matrix(zma))
    return tuple(more_itertools.unique_everseen(name_mat[3:, 2]))


def dihedral_angle_value_names(zma):
    """ dihedral angle value coordinate names
    """
    val_dct = values(zma)
    return tuple(name for name in dihedral_angle_names(zma)
                 if name in val_dct)


def dihedral_angle_variable_names(zma):
    """ dihedral angle variable coordinate names
    """
    val_dct = values(zma)
    return tuple(name for name in dihedral_angle_names(zma)
                 if name not in val_dct)


def angle_names(zma):
    """ angle coordinate names (dihedral and central)
    """
    return tuple(itertools.chain(central_angle_names(zma),
                                 dihedral_angle_names(zma)))


def angle_value_names(zma):
    """ dihedral angle value coordinate names
    """
    val_dct = values(zma)
    return tuple(name for name in angle_names(zma) if name in val_dct)


def angle_variable_names(zma):
    """ dihedral angle variable coordinate names
    """
    val_dct = values(zma)
    return tuple(name for name in angle_names(zma) if name not in val_dct)


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

    val_dct = {name_dct[orig_name]: val
               for orig_name, val in values(zma).items()}

    return _from_data(symbols(zma), key_matrix(zma), name_mat, val_dct,
                      complete=False)


def set_values(zma, val_dct):
    """ set coordinate values for the z-matrix
    """
    _names = names(zma)
    assert set(val_dct.keys()) <= set(_names)

    new_val_dct = values(zma).copy()
    new_val_dct.update(val_dct)
    return _from_data(symbols(zma), key_matrix(zma), name_matrix(zma),
                      new_val_dct, complete=False)


# misc
def is_valid(zma, require_complete=True):
    """ is this a valid zmatrix?
    """
    ret = hasattr(zma, '__len__') and len(zma) == 2
    if ret:
        mat, val_dct = zma
        ret = all(hasattr(obj, '__len__') and len(obj) == 3 for obj in mat)
        if ret:
            syms, key_mat, name_mat = zip(*mat)
            try:
                _from_data(syms, key_mat, name_mat, val_dct,
                           complete=require_complete)
            except AssertionError:
                ret = False
    return ret


def is_complete(zma):
    """ is this a complete zmatrix? (no variables)
    """
    return value_names(zma) == names(zma) and not variable_names(zma)
