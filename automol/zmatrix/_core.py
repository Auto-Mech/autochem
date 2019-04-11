""" core library defining the z-matrix data structure
"""
import numpy
import phycon.units as pcu
from automol.zmatrix import v as _v_
from automol.constructors.zmatrix import from_data as _from_data


# value getters
def var_(zma):
    """ the variable matrix (atom symbols, atom keys, and coordinate names)
    """
    vma, _ = zma
    return vma


def symbols(zma):
    """ atomic symbols, by z-matrix row
    """
    return _v_.symbols(var_(zma))


def key_matrix(zma, one_indexed=False):
    """ coordinate atom keys, by z-matrix row and column
    """
    return _v_.key_matrix(var_(zma), one_indexed)


def name_matrix(zma):
    """ coordinate names, by z-matrix row and column
    """
    return _v_.name_matrix(var_(zma))


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
    return _v_.coordinate_key_matrix(var_(zma), one_indexed)


def coordinates(zma, one_indexed=False):
    """ coordinate keys associated with each coordinate name, as a dictionary

    (the values are sequences of coordinate keys, since there may be multiple)
    """
    return _v_.coordinates(var_(zma), one_indexed)


def names(zma):
    """ coordinate names
    """
    return _v_.names(var_(zma))


def distance_names(zma):
    """ distance coordinate names
    """
    return _v_.distance_names(var_(zma))


def central_angle_names(zma):
    """ central angle coordinate names
    """
    return _v_.central_angle_names(var_(zma))


def dihedral_angle_names(zma):
    """ dihedral angle coordinate names
    """
    return _v_.dihedral_angle_names(var_(zma))


def angle_names(zma):
    """ angle coordinate names (dihedral and central)
    """
    return _v_.angle_names(var_(zma))


def values(zma, angstrom=False, degree=False):
    """ coordinate values, by coordinate name
    """
    vma, val_dct = zma

    # post-processing for unit convertions
    dist_names = _v_.distance_names(vma)
    ang_names = _v_.angle_names(vma)
    orig_val_dct = val_dct

    val_dct = {}
    for name, val in orig_val_dct.items():
        if angstrom and name in dist_names:
            val *= pcu.BOHR2ANG
        if degree and name in ang_names:
            val *= pcu.RAD2DEG
        val_dct[name] = val
    return val_dct


# value setters
def set_names(zma, name_dct):
    """ set coordinate names for the z-matrix
    """
    orig_vma = var_(zma)
    vma = _v_.set_names(orig_vma, name_dct)
    name_mat = _v_.name_matrix(vma)

    name_dct = dict(zip(numpy.ravel(_v_.name_matrix(orig_vma)),
                        numpy.ravel(name_mat)))
    val_dct = {name_dct[orig_name]: val
               for orig_name, val in values(zma).items()}

    return _from_data(_v_.symbols(vma), _v_.key_matrix(vma), name_mat, val_dct)


def set_values(zma, val_dct):
    """ set coordinate values for the z-matrix
    """
    vma = var_(zma)
    _names = _v_.names(vma)
    assert set(val_dct.keys()) <= set(_names)

    new_val_dct = values(zma).copy()
    new_val_dct.update(val_dct)
    return _from_data(_v_.symbols(vma), _v_.key_matrix(vma),
                      _v_.name_matrix(vma), new_val_dct)


def standard_form(zma):
    """ set standard variable names for the z-matrix
    """
    orig_vma = var_(zma)
    vma = _v_.standard_form(orig_vma)
    name_dct = dict(zip(numpy.ravel(_v_.name_matrix(orig_vma)),
                        numpy.ravel(_v_.name_matrix(vma))))
    name_dct.pop(None)
    return set_names(zma, name_dct)


# misc
def is_valid(zma):
    """ is this a valid zmatrix?
    """
    ret = hasattr(zma, '__len__') and len(zma) == 2
    if ret:
        vma, val_dct = zma
        ret = _v_.is_valid(vma) and set(_v_.names(vma)) == set(val_dct)
        if ret:
            try:
                _from_data(_v_.symbols(vma), _v_.key_matrix(vma),
                           _v_.name_matrix(vma), val_dct)
            except AssertionError:
                ret = False
    return ret


def is_standard_form(zma):
    """ set standard variable names for the z-matrix
    """
    return is_valid(zma) and _v_.is_standard_form(var_(zma))
