""" z-matrix
"""
import numpy
from qcelemental import constants as qcc
from automol import vmat
import automol.create.zmat
import automol.convert.zmat
import automol.geom
import autoread as ar
import autowrite as aw

BOHR2ANG = qcc.conversion_factor('bohr', 'angstrom')
RAD2DEG = qcc.conversion_factor('radian', 'degree')


# getters
def count(zma):
    """ the number of z-matrix rows (number of atoms or dummy atoms)
    """
    return vmat.count(zma)


def symbols(zma):
    """ atomic symbols, by z-matrix row
    """
    return vmat.symbols(zma)


def key_matrix(zma, shift=0):
    """ coordinate atom keys, by z-matrix row and column
    """
    return vmat.key_matrix(zma, shift=shift)


def name_matrix(zma):
    """ coordinate names, by z-matrix row and column
    """
    return vmat.name_matrix(zma)


def value_matrix(zma, angstrom=False, degree=False):
    """ coordinate values, by z-matrix row and column
    """
    if zma:
        val_mat = tuple(zip(*zma))[3]

        val_mat = [list(row) + [None]*(3-len(row)) for row in val_mat]
        val_mat = numpy.array(val_mat, dtype=numpy.object_)

        val_mat[1:, 0] *= BOHR2ANG if angstrom else 1
        val_mat[2:, 1] *= RAD2DEG if degree else 1
        val_mat[3:, 2] *= RAD2DEG if degree else 1
    else:
        val_mat = ()

    return tuple(map(tuple, val_mat))


def value_dictionary(zma, angstrom=False, degree=False):
    """ coordinate values, by coordinate name (as a dictionary)
    """
    names = numpy.ravel(name_matrix(zma))
    vals = numpy.ravel(value_matrix(zma, angstrom=angstrom, degree=degree))
    val_dct = dict(zip(names, vals))
    val_dct.pop(None)
    return val_dct


def distance(zma, key1, key2, angstrom=False):
    """ measure the distance between atoms
    """
    geo = automol.convert.zmat.geometry(zma)
    return automol.geom.distance(geo, key1, key2, angstrom=angstrom)


def central_angle(zma, key1, key2, key3, degree=False):
    """ measure the angle inscribed by three atoms
    """
    geo = automol.convert.zmat.geometry(zma)
    return automol.geom.central_angle(geo, key1, key2, key3, degree=degree)


def dihedral_angle(zma, key1, key2, key3, key4, degree=False):
    """ measure the angle inscribed by three atoms
    """
    geo = automol.convert.zmat.geometry(zma)
    return automol.geom.dihedral_angle(geo, key1, key2, key3, key4,
                                       degree=degree)


# setters
def set_key_matrix(zma, key_mat):
    """ set the key matrix
    """
    syms = symbols(zma)
    val_mat = value_matrix(zma)
    name_mat = name_matrix(zma)
    zma = automol.create.zmat.from_data(syms, key_mat, val_mat, name_mat)
    return zma


def set_name_matrix(zma, name_mat):
    """ set the name matrix
    """
    syms = symbols(zma)
    val_mat = value_matrix(zma)
    key_mat = key_matrix(zma)
    zma = automol.create.zmat.from_data(syms, key_mat, val_mat, name_mat)
    return zma


def set_value_matrix(zma, val_mat):
    """ set the name matrix
    """
    syms = symbols(zma)
    key_mat = key_matrix(zma)
    name_mat = name_matrix(zma)
    zma = automol.create.zmat.from_data(syms, key_mat, val_mat, name_mat)
    return zma


def standard_name_matrix(zma, shift=0):
    """ standard names for the coordinates (follows x2z format)
    """
    return vmat.standard_name_matrix(zma, shift=shift)


def standard_form(zma, shift=0):
    """ set standard variable names for the z-matrix (x2z format)
    """
    name_mat = standard_name_matrix(zma, shift=shift)
    return set_name_matrix(zma, name_mat)


def rename(zma, name_dct):
    """ set coordinate names for the z-matrix
    """
    syms = symbols(zma)
    key_mat = key_matrix(zma)
    val_mat = value_matrix(zma)

    vma = vmat.rename(zma, name_dct)
    name_mat = vmat.name_matrix(vma)

    zma = automol.create.zmat.from_data(syms, key_mat, val_mat, name_mat)
    return zma


# operations
def add_atom(zma, sym, key_row, val_row, name_row=None,
             one_indexed=False, angstrom=True, degree=True):
    """ add an atom to the z-matrix
    """
    syms = symbols(zma)
    syms += (sym,)

    key_mat = key_matrix(zma, shift=(1 if one_indexed else 0))
    key_mat += (key_row,)

    val_mat = value_matrix(zma, angstrom=angstrom, degree=degree)
    val_mat += (val_row,)

    name_mat = None if name_row is None else name_matrix(zma) + (name_row,)

    zma = automol.create.zmat.from_data(
        syms, key_mat, val_mat, name_mat, one_indexed=one_indexed,
        angstrom=angstrom, degree=degree)
    return zma


def remove_atom(zma, key):
    """ remove an atom from the z-matrix

    complains if attempting to remove an atom that other atoms depend on
    """
    syms = list(symbols(zma))
    syms.pop(key)

    key_mat = list(key_matrix(zma))
    key_mat.pop(key)
    key_mat = numpy.array(key_mat, dtype=numpy.object_)

    if (key_mat == key).any():
        raise ValueError("Other atoms in z-matrix depend on atom {}"
                         .format(key))

    key_map = numpy.vectorize(lambda x: x if (x is None or x < key) else x-1)
    key_mat = key_map(key_mat)

    val_mat = list(value_matrix(zma))
    val_mat.pop(key)

    name_mat = list(name_matrix(zma))
    name_mat.pop(key)

    zma = automol.create.zmat.from_data(syms, key_mat, val_mat, name_mat)
    return zma


def join_replace_one(zma1, zma2, rep1_key, key_mat, val_mat, name_mat=None,
                     degree=True):
    """ join two z-matrices, replacing the first atom in zma2 by an atom in zma1

    :param rep1_key: the key of the atom in zma1 that will replace the first
        atom in zma2
    """
    rep2_key = count(zma1) - 1

    syms = symbols(zma1) + symbols(zma2)[1:]

    key_mat1 = key_matrix(zma1)
    key_mat2 = key_matrix(zma2, shift=rep2_key)[1:]
    key_mat2 = numpy.array(key_mat2, dtype=numpy.object_)
    key_mat2[key_mat2 == rep2_key] = rep1_key
    key_mat2[0, 1:] = key_mat[0][-2:]
    key_mat2[1, 2:] = key_mat[1][-1:]
    key_mat2 = tuple(map(tuple, key_mat2))
    key_mat = key_mat1 + key_mat2

    val_mat1 = value_matrix(zma1, degree=degree)
    val_mat2 = value_matrix(zma2, degree=degree)[1:]
    val_mat2 = numpy.array(val_mat2, dtype=numpy.object_)
    val_mat2[0, 1:] = val_mat[0][-2:]
    val_mat2[1, 2:] = val_mat[1][-1:]
    val_mat2 = tuple(map(tuple, val_mat2))
    val_mat = val_mat1 + val_mat2

    if name_mat is not None:
        name_mat1 = name_matrix(zma1)
        name_mat2 = name_matrix(zma2)[1:]
        name_mat2 = numpy.array(name_mat2, dtype=numpy.object_)
        name_mat2[0, 1:] = name_mat[0][-2:]
        name_mat2[1, 2:] = name_mat[1][-1:]
        name_mat2 = tuple(map(tuple, name_mat2))
        name_mat = name_mat1 + name_mat2

    zma = automol.create.zmat.from_data(syms, key_mat, val_mat, name_mat,
                                        degree=degree)
    return zma


# conversions
def vmatrix(zma):
    """ convert z-matrix to v-matrix
    """
    return vmat.from_data(syms=symbols(zma),
                          key_mat=key_matrix(zma),
                          name_mat=name_matrix(zma))


# comparisons
def almost_equal(zma1, zma2, dist_rtol=2e-5, ang_atol=2e-3, just_dist=False):
    """ are these z-matrices numerically equal?

    :param zma1: The first z-matrix
    :param zma2: The second z-matrix
    :param dist_rtol: Relative tolerance for the distances
    :type dist_rtol: float
    :param ang_atol: Absolute tolerance for the angles
    :type ang_atol: float
    :param just_dist: Only compare distances?
    :type just_dist: bool
    """
    ret = False
    if vmatrix(zma1) == vmatrix(zma2):
        val_mat1 = numpy.array(value_matrix(zma1), dtype=float)
        val_mat2 = numpy.array(value_matrix(zma2), dtype=float)

        # first compare the distances
        dist_vals1 = val_mat1[1:, 0]
        dist_vals2 = val_mat2[1:, 0]
        if numpy.allclose(dist_vals1, dist_vals2, rtol=dist_rtol):
            if just_dist:
                ret = True
            else:
                # now compare the angles
                # see https://gamedev.stackexchange.com/a/4472
                ang_vals1 = numpy.hstack((val_mat1[2:, 1], val_mat1[3:, 2]))
                ang_vals2 = numpy.hstack((val_mat1[2:, 1], val_mat1[3:, 2]))
                ang_vals1 = numpy.mod(ang_vals1, 2*numpy.pi)
                ang_vals2 = numpy.mod(ang_vals2, 2*numpy.pi)

                ang_diffs = numpy.abs(ang_vals1 - ang_vals2)
                ang_diffs = numpy.pi - numpy.abs(ang_diffs - numpy.pi)

                if numpy.allclose(ang_diffs, 0., atol=ang_atol):
                    ret = True
    return ret


# I/O
def from_string(zma_str, one_indexed=True, angstrom=True, degree=True):
    """ read a z-matrix from a string
    """
    syms, key_mat, name_mat, val_dct = ar.zmatrix.read(zma_str)

    val_mat = tuple(tuple(val_dct[name] if name is not None else None
                          for name in name_mat_row)
                    for name_mat_row in name_mat)

    zma = automol.create.zmat.from_data(
        syms, key_mat, val_mat, name_mat, one_indexed=one_indexed,
        angstrom=angstrom, degree=degree)
    return zma


def string(zma, one_indexed=True, angstrom=True, degree=True):
    """ write a z-matrix to a string
    """
    shift = 1 if one_indexed else 0
    zma_str = aw.zmatrix.write(
        syms=symbols(zma),
        key_mat=key_matrix(zma, shift=shift),
        name_mat=name_matrix(zma),
        val_dct=value_dictionary(zma, angstrom=angstrom, degree=degree)
    )
    return zma_str


# validators
def is_valid(zma):
    """ is this a valid vmatrix?
    """
    ret = True

    try:
        assert _is_sequence_of_quadruples(zma)
        syms, key_mat, name_mat, val_mat = zip(*zma)
        automol.create.zmat.from_data(syms, key_mat, val_mat, name_mat)
    except AssertionError:
        ret = False

    return ret


def _is_sequence_of_quadruples(obj):
    ret = hasattr(obj, '__len__')
    if ret:
        ret = all(hasattr(item, '__len__') and len(item) == 4 for item in obj)
    return ret
