""" variable z-matrix

(z-matrix without coordinate values)
"""
import itertools
import more_itertools
import numpy
from qcelemental import periodictable as pt
import autoread as ar
import autowrite as aw
import automol.create.vmat
import automol.geom


# constructor
def from_data(syms, key_mat, name_mat=None, one_indexed=False):
    """ v-matrix constructor

    :param syms: atomic symbols
    :type syms: tuple[str]
    :param key_mat: key/index columns of the v-matrix, zero-indexed
    :type key_mat: tuple[tuple[float, float or None, float or None]]
    :param name_mat: coordinate name columns of the v-matrix
    :type name_mat; tuple[tuple[str, str or None, str or None]]
    """
    return automol.create.vmat.from_data(
        symbols=syms, key_matrix=key_mat, name_matrix=name_mat,
        one_indexed=one_indexed)


# getters
def symbols(vma):
    """ atomic symbols, by v-matrix row
    """
    if vma:
        syms = tuple(zip(*vma))[0]
    else:
        syms = ()
    return syms


def count(zma):
    """ the number of v-matrix rows (number of atoms or dummy atoms)
    """
    return len(symbols(zma))


def key_matrix(vma, shift=0):
    """ coordinate atom keys, by v-matrix row and column
    """
    if vma:
        key_mat = tuple(zip(*vma))[1]

        # post-processing for adding the shift
        key_mat = [list(row) + [None]*(3-len(row)) for row in key_mat]
        key_mat = numpy.array(key_mat)
        tril_idxs = numpy.tril_indices(key_mat.shape[0], -1, m=3)
        key_mat[tril_idxs] += shift
    else:
        key_mat = ()

    return tuple(map(tuple, key_mat))


def name_matrix(vma):
    """ coordinate names, by v-matrix row and column
    """
    if vma:
        name_mat = tuple(zip(*vma))[2]
    else:
        name_mat = ()

    name_mat = [list(row) + [None]*(3-len(row)) for row in name_mat]

    return tuple(map(tuple, name_mat))


def coordinate_key_matrix(vma, shift=0):
    """ coordinate keys, by v-matrix row and column
    """
    key_mat = key_matrix(vma, shift=shift)
    natms = len(key_mat)
    atm_keys = range(shift, natms+shift)
    coo_key_mat = [[(atm_key,) + key_row[:col+1]
                    if key_row[col] is not None else None for col in range(3)]
                   for atm_key, key_row in zip(atm_keys, key_mat)]
    return tuple(map(tuple, coo_key_mat))


def coordinates(vma, shift=0, multi=True):
    """ coordinate keys associated with each coordinate name, as a dictionary

    (the values are sequences of coordinate keys, since there may be multiple)
    """
    _names = numpy.ravel(name_matrix(vma))
    coo_keys = numpy.ravel(coordinate_key_matrix(vma, shift))

    if not multi:
        coo_dct = dict(zip(_names, coo_keys))
    else:
        coo_dct = {name: () for name in _names}
        for name, coo_key in zip(_names, coo_keys):
            coo_dct[name] += (coo_key,)

    coo_dct.pop(None)
    return coo_dct


def names(vma):
    """ coordinate names
    """
    name_mat = name_matrix(vma)
    _names = filter(lambda x: x is not None,
                    numpy.ravel(numpy.transpose(name_mat)))
    return tuple(more_itertools.unique_everseen(_names))


def distance_names(vma):
    """ distance coordinate names
    """
    name_mat = numpy.array(name_matrix(vma))
    return tuple(more_itertools.unique_everseen(name_mat[1:, 0]))


def central_angle_names(vma):
    """ central angle coordinate names
    """
    name_mat = numpy.array(name_matrix(vma))
    return tuple(more_itertools.unique_everseen(name_mat[2:, 1]))


def dihedral_angle_names(vma):
    """ dihedral angle coordinate names
    """
    name_mat = numpy.array(name_matrix(vma))
    return tuple(more_itertools.unique_everseen(name_mat[3:, 2]))


def angle_names(vma):
    """ angle coordinate names (dihedral and central)
    """
    return tuple(itertools.chain(central_angle_names(vma),
                                 dihedral_angle_names(vma)))


def dummy_coordinate_names(vma):
    """ names of dummy atom coordinates
    """
    syms = symbols(vma)
    name_mat = numpy.array(name_matrix(vma))
    dummy_keys = [idx for idx, sym in enumerate(syms) if not pt.to_Z(sym)]
    dummy_names = []
    for dummy_key in dummy_keys:
        for col_idx in range(3):
            dummy_name = next(filter(lambda x: x is not None,
                                     name_mat[dummy_key:, col_idx]))
            dummy_names.append(dummy_name)

    dummy_names = tuple(dummy_names)
    return dummy_names


# value setters
def set_key_matrix(zma, key_mat):
    """ set the key matrix
    """
    syms = symbols(zma)
    name_mat = name_matrix(zma)
    zma = automol.create.vmat.from_data(syms, key_mat, name_mat)
    return zma


def set_name_matrix(zma, name_mat):
    """ set the name matrix
    """
    syms = symbols(zma)
    key_mat = key_matrix(zma)
    zma = automol.create.vmat.from_data(syms, key_mat, name_mat)
    return zma


def standard_name_matrix(vma, shift=0):
    """ standard names for the coordinates (follows x2z format)
    """
    natms = count(vma)

    name_mat = numpy.array(name_matrix(vma), dtype=numpy.object_)

    name_mat[1:, 0] = ['R{:d}'.format(num + shift + 1) for num in range(natms)]
    name_mat[2:, 1] = ['A{:d}'.format(num + shift + 2) for num in range(natms)]
    name_mat[3:, 2] = ['D{:d}'.format(num + shift + 3) for num in range(natms)]

    name_mat = tuple(map(tuple, name_mat))
    return name_mat


def standard_form(vma, shift=0):
    """ set standard variable names for the z-matrix (x2z format)
    """
    name_mat = standard_name_matrix(vma, shift=shift)
    return set_name_matrix(vma, name_mat)


def rename(vma, name_dct):
    """ set coordinate names for the variable v-matrix
    """
    orig_name_mat = numpy.array(name_matrix(vma))
    tril_idxs = numpy.tril_indices(orig_name_mat.shape[0], -1, m=3)
    orig_names = set(orig_name_mat[tril_idxs])
    assert set(name_dct.keys()) <= orig_names

    name_dct.update({orig_name: orig_name for orig_name in orig_names
                     if orig_name not in name_dct})

    name_mat = numpy.empty(orig_name_mat.shape, dtype=numpy.object_)
    name_mat[tril_idxs] = list(map(name_dct.__getitem__,
                                   orig_name_mat[tril_idxs]))

    return from_data(symbols(vma), key_matrix(vma), name_mat)


def standard_names(vma, shift=0):
    """ standard names for the coordinates, by their current names

    (follows x2z format)
    """
    dist_names = distance_names(vma)
    cent_ang_names = central_angle_names(vma)
    dih_ang_names = dihedral_angle_names(vma)
    name_dct = {}
    name_dct.update({
        dist_name: 'R{:d}'.format(num + shift + 1)
        for num, dist_name in enumerate(dist_names)})
    name_dct.update({
        cent_ang_name: 'A{:d}'.format(num + shift + 2)
        for num, cent_ang_name in enumerate(cent_ang_names)})
    name_dct.update({
        dih_ang_name: 'D{:d}'.format(num + shift + 3)
        for num, dih_ang_name in enumerate(dih_ang_names)})
    return name_dct


# operations
def add_atom(vma, sym, key_row, name_row=None, one_indexed=False):
    """ add an atom to the z-matrix
    """
    syms = symbols(vma)
    syms += (sym,)

    key_mat = key_matrix(vma, shift=(1 if one_indexed else 0))
    key_mat += (key_row,)

    name_mat = None if name_row is None else name_matrix(vma) + (name_row,)

    vma = automol.create.vmat.from_data(
        syms, key_mat, name_mat, one_indexed=one_indexed)
    return vma


def remove_atom(vma, key):
    """ remove an atom from the z-matrix

    complains if attempting to remove an atom that other atoms depend on
    """
    syms = list(symbols(vma))
    syms.pop(key)

    key_mat = list(key_matrix(vma))
    key_mat.pop(key)
    key_mat = numpy.array(key_mat, dtype=numpy.object_)

    if (key_mat == key).any():
        raise ValueError("Other atoms in z-matrix depend on atom {}"
                         .format(key))

    key_map = numpy.vectorize(lambda x: x if (x is None or x < key) else x-1)
    key_mat = key_map(key_mat)

    name_mat = list(name_matrix(vma))
    name_mat.pop(key)

    vma = automol.create.vmat.from_data(syms, key_mat, name_mat)
    return vma


def is_valid(vma):
    """ is this a valid vmatrix?
    """
    ret = True

    try:
        assert _is_sequence_of_triples(vma)
        syms, key_mat, name_mat = zip(*vma)
        from_data(syms, key_mat, name_mat)
    except AssertionError:
        ret = False

    return ret


def _is_sequence_of_triples(obj):
    ret = hasattr(obj, '__len__')
    if ret:
        ret = all(hasattr(item, '__len__') and len(item) == 3 for item in obj)
    return ret


def is_standard_form(vma):
    """ set standard variable names for the v-matrix

    (follows x2z format)
    """
    return names(vma) == names(standard_form(vma))


# I/O
def from_string(vma_str):
    """ read a v-matrix from a string
    """
    syms, key_mat, name_mat = ar.zmatrix.matrix.read(vma_str)

    vma = from_data(syms, key_mat, name_mat, one_indexed=True)
    return vma


def string(vma, one_indexed=True):
    """ write a v-matrix to a string
    """
    shift = 1 if one_indexed else 0
    vma_str = aw.zmatrix.matrix_block(
        syms=symbols(vma),
        key_mat=key_matrix(vma, shift=shift),
        name_mat=name_matrix(vma),
    )
    return vma_str
