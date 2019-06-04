""" variable z-matrix

(z-matrix without the value dictionary)
"""
import itertools
import more_itertools
import numpy
import autoread as ar
import autowrite as aw
import automol.create.vmatrix


# constructor
def from_data(syms, key_mat, name_mat, one_indexed=False):
    """ z-matrix constructor

    :param syms: atomic symbols
    :type syms: tuple[str]
    :param key_mat: key/index columns of the z-matrix, zero-indexed
    :type key_mat: tuple[tuple[float, float or None, float or None]]
    :param name_mat: coordinate name columns of the z-matrix
    :type name_mat; tuple[tuple[str, str or None, str or None]]
    """
    return automol.create.vmatrix.from_data(
        symbols=syms, key_matrix=key_mat, name_matrix=name_mat,
        one_indexed=one_indexed)


# getters
def count(zma):
    """ the number of z-matrix rows (number of atoms or dummy atoms)
    """
    return len(symbols(zma))


def symbols(vma):
    """ atomic symbols, by z-matrix row
    """
    if vma:
        syms, _, _ = zip(*vma)
    else:
        syms = ()
    return syms


def key_matrix(vma, shift=0):
    """ coordinate atom keys, by z-matrix row and column
    """
    if vma:
        _, key_mat, _ = zip(*vma)
    else:
        key_mat = ()

    # post-processing for adding the shift
    key_mat = numpy.array(key_mat)
    tril_idxs = numpy.tril_indices(key_mat.shape[0], -1, m=3)
    key_mat[tril_idxs] += shift

    return tuple(map(tuple, key_mat))


def name_matrix(vma):
    """ coordinate names, by z-matrix row and column
    """
    if vma:
        _, _, name_mat = zip(*vma)
    else:
        name_mat = ()
    return name_mat


def coordinate_key_matrix(vma, shift=0):
    """ coordinate keys, by z-matrix row and column
    """
    key_mat = key_matrix(vma, shift=shift)
    natms = len(key_mat)
    atm_keys = range(shift, natms+shift)
    coo_key_mat = [[(atm_key,) + key_row[:col+1]
                    if key_row[col] is not None else None for col in range(3)]
                   for atm_key, key_row in zip(atm_keys, key_mat)]
    return tuple(map(tuple, coo_key_mat))


def coordinates(vma, shift=0):
    """ coordinate keys associated with each coordinate name, as a dictionary

    (the values are sequences of coordinate keys, since there may be multiple)
    """
    _names = numpy.ravel(name_matrix(vma))
    coo_keys = numpy.ravel(coordinate_key_matrix(vma, shift))

    coo_dct = {name: () for name in _names}
    for name, coo_key in zip(_names, coo_keys):
        coo_dct[name] += (coo_key,)

    coo_dct.pop(None)
    return coo_dct


def names(vma):
    """ coordinate names
    """
    _names = filter(lambda x: x is not None,
                    numpy.ravel(numpy.transpose(name_matrix(vma))))
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


# value setters
def set_names(vma, name_dct):
    """ set coordinate names for the variable z-matrix
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
    """
    dist_names = distance_names(vma)
    cent_ang_names = central_angle_names(vma)
    dih_ang_names = dihedral_angle_names(vma)
    name_dct = {}
    name_dct.update({
        dist_name: 'r{:d}'.format(num + shift + 1)
        for num, dist_name in enumerate(dist_names)})
    name_dct.update({
        cent_ang_name: 'a{:d}'.format(num + shift + 1)
        for num, cent_ang_name in enumerate(cent_ang_names)})
    name_dct.update({
        dih_ang_name: 'd{:d}'.format(num + shift + 1)
        for num, dih_ang_name in enumerate(dih_ang_names)})
    return name_dct


def standard_form(vma):
    """ set standard variable names for the variable z-matrix
    """
    name_dct = standard_names(vma)
    return set_names(vma, name_dct)


def is_valid(vma):
    """ is this a valid zmatrix?
    """
    ret = _is_sequence_of_triples(vma)
    if ret:
        syms, key_mat, name_mat = zip(*vma)
        ret = (_is_sequence_of_triples(key_mat) and
               _is_sequence_of_triples(name_mat))
        if ret:
            try:
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
    """ set standard variable names for the z-matrix
    """
    return names(vma) == names(standard_form(vma))


# I/O
def from_string(vma_str):
    """ read a z-matrix from a string
    """
    syms, key_mat, name_mat = ar.zmatrix.matrix.read(vma_str)

    vma = from_data(syms, key_mat, name_mat, one_indexed=True)
    return vma


def string(vma):
    """ write a z-matrix to a string
    """
    vma_str = aw.zmatrix.matrix_block(
        syms=symbols(vma),
        key_mat=key_matrix(vma, shift=1),
        name_mat=name_matrix(vma),
    )
    return vma_str
