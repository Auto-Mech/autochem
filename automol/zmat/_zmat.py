""" z-matrix
"""
import numpy
from qcelemental import constants as qcc
from automol import vmat
from automol.vmat import standard_names
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


# setters
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


def standard_form(zma):
    """ set standard variable names for the z-matrix (x2z format)
    """
    name_dct = standard_names(zma)
    return rename(zma, name_dct)


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
