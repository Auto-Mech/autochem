""" core library defining the z-matrix data structure
"""
import itertools
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


def count(zma):
    """ the number of z-matrix rows (number of atoms or dummy atoms)
    """
    return _v_.count(var_(zma))


def symbols(zma):
    """ atomic symbols, by z-matrix row
    """
    return _v_.symbols(var_(zma))


def key_matrix(zma, shift=0):
    """ coordinate atom keys, by z-matrix row and column
    """
    return _v_.key_matrix(var_(zma), shift=shift)


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


def coordinate_key_matrix(zma, shift=0):
    """ coordinate keys, by z-matrix row and column
    """
    return _v_.coordinate_key_matrix(var_(zma), shift=shift)


def coordinates(zma, shift=0):
    """ coordinate keys associated with each coordinate name, as a dictionary

    (the values are sequences of coordinate keys, since there may be multiple)
    """
    return _v_.coordinates(var_(zma), shift=shift)


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


def standard_names(zma, shift=0):
    """ standard names for the coordinates, by their current names
    """
    return _v_.standard_names(var_(zma), shift=shift)


def standard_form(zma, shift=0):
    """ set standard variable names for the z-matrix
    """
    return set_names(zma, standard_names(zma, shift=shift))


# operations
def join(zma1, zma2, join_key_mat, join_name_mat, join_val_dct):
    """ join two z-matrices together
    """
    syms1 = symbols(zma1)
    syms2 = symbols(zma2)
    natms1 = count(zma1)
    natms2 = count(zma2)
    key_mat1 = numpy.array(key_matrix(zma1))
    key_mat2 = numpy.array(key_matrix(zma2, shift=natms1))  # note the shift
    name_mat1 = numpy.array(name_matrix(zma1))
    name_mat2 = numpy.array(name_matrix(zma2))
    val_dct1 = values(zma1)
    val_dct2 = values(zma2)

    join_natms = min(natms2, 3)
    assert len(join_key_mat) == len(join_name_mat) == join_natms

    join_key_mat = numpy.array(join_key_mat, dtype=numpy.object_)
    join_name_mat = numpy.array(join_name_mat, dtype=numpy.object_)

    # make sure we aren't overwriting values -- the constructor should take
    # care of the rest of the necessary validation
    assert numpy.all(numpy.equal(join_key_mat, None) ==
                     numpy.equal(join_key_mat, None))

    join_idxs = numpy.not_equal(join_key_mat, None)
    assert numpy.all(numpy.equal(key_mat2[join_idxs], None))
    assert numpy.all(numpy.equal(name_mat2[join_idxs], None))
    key_mat2[join_idxs] = join_key_mat[join_idxs]
    name_mat2[join_idxs] = join_name_mat[join_idxs]

    syms = tuple(itertools.chain(syms1, syms2))
    key_mat = tuple(itertools.chain(key_mat1, key_mat2))
    name_mat = tuple(itertools.chain(name_mat1, name_mat2))

    # Could be made to allow for joins with common zma1 and zma2 names (for
    # symmetry constraints).  Not sure if we really want that.
    val_dct = val_dct1.copy()
    assert not set(val_dct.keys()) & set(val_dct2.keys())
    assert not set(val_dct.keys()) & set(join_val_dct.keys())
    val_dct.update(val_dct2)
    val_dct.update(join_val_dct)

    return _from_data(syms, key_mat, name_mat, val_dct)


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


# helpers
def _coordinate_count(natms):
    """ z-matrix degrees of freedom, by the number of atoms
    """
    assert natms > 0
    if natms == 1:
        ncoo = 0
    elif natms == 2:
        ncoo = 1
    elif natms > 2:
        ncoo = 3 * natms - 6
    return ncoo
