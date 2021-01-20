""" V-Matrix: Variable Z-Matrix (Z-Matrix without coordinate values)
"""

import itertools
import more_itertools
import numpy
import autoread as ar
import autowrite as aw
import automol.create.vmat
import automol.geom
from phydat import ptab


# constructor
def from_data(symbs, key_mat, name_mat=None, one_indexed=False):
    """ Build a V-Matrix data structure from atomic symbols and coordinates.

        :param symbs: atomic symbols
        :type symbs: tuple[str]
        :param key_mat: key/index columns of the v-matrix, zero-indexed
        :type key_mat: tuple[tuple[float, float or None, float or None]]
        :param name_mat: coordinate name columns of the v-matrix
        :type name_mat; tuple[tuple[str, str or None, str or None]]
        :param one_indexed: parameter to store keys in one-indexing
        :type one_indexed: bool
        :rtype: automol V-Matrix data structure
    """
    return automol.create.vmat.from_data(
        symbols=symbs, key_matrix=key_mat, name_matrix=name_mat,
        one_indexed=one_indexed)


# getters
def symbols(vma):
    """ Obtain the atomic symbols for all atoms defined in the V-Matrix.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :rtype: tuple(str)
    """
    return tuple(zip(*vma))[0] if vma else ()


def count(vma):
    """ Obtain the number of rows of the V-Matrix, which corresponds to
        the number of atoms defined in the V-Matrix. This includes all
        real and dummy atoms.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :rtype: int
    """
    return len(symbols(vma))


def key_matrix(vma, shift=0):
    """ Obtain the key matrix of the V-Matrix that containts the
        coordinate atom keys by row and column.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param shift: value to shift the keys by when obtaining the key matrix
        :type shift: int
        :rtype: tuple(tuple(int))
    """

    if vma:
        key_mat = tuple(zip(*vma))[1]

        # post-processing for adding the shift
        key_mat = [list(row) + [None]*(3-len(row)) for row in key_mat]
        key_mat = numpy.array(key_mat)
        tril_idxs = numpy.tril_indices(key_mat.shape[0], -1, m=3)
        key_mat[tril_idxs] += shift
        key_mat[tril_idxs] = key_mat[tril_idxs].astype(int)
    else:
        key_mat = ()

    return tuple(map(tuple, key_mat))


def name_matrix(vma):
    """ Obtain the name matrix of the V-Matrix that containts the
        coordinate names by row and column.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param shift: value to shift the keys by when obtaining the key matrix
        :type shift: int
        :rtype: tuple(tuple(str))
    """

    if vma:
        name_mat = tuple(zip(*vma))[2]
    else:
        name_mat = ()

    name_mat = [list(row) + [None]*(3-len(row)) for row in name_mat]

    return tuple(map(tuple, name_mat))


def coordinate_key_matrix(vma, shift=0):
    """ Obtain the coordinate key matrix of the V-Matrix that containts the
        coordinate keys by row and column.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param shift: value to shift the keys by when obtaining the key matrix
        :type shift: int
        :rtype: tuple(tuple(int))
    """

    key_mat = key_matrix(vma, shift=shift)
    natms = len(key_mat)
    atm_keys = range(shift, natms+shift)
    coo_key_mat = [[(atm_key,) + key_row[:col+1]
                    if key_row[col] is not None else None for col in range(3)]
                   for atm_key, key_row in zip(atm_keys, key_mat)]

    return tuple(map(tuple, coo_key_mat))


def coordinates(vma, shift=0, multi=True):
    """ Obtain the coordinate keys associated with each coordinate name,
        as a dictionary. Values are sequences of coordinate keys,
        since there may be multiple.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param shift: value to shift the keys by when obtaining the keys
        :type shift: int
        :param multi: parameter to grab multiple coordinate keys
        :type multi: bool
        :rtype: dict[str: tuple(int)]
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
    """ Obtain names of all coordinates defined in the V-Matrix.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :rtype: tuple(str)
    """

    name_mat = name_matrix(vma)
    _names = filter(lambda x: x is not None,
                    numpy.ravel(numpy.transpose(name_mat)))

    return tuple(more_itertools.unique_everseen(_names))


def distance_names(vma):
    """ Obtain names of all distance coordinates defined in the V-Matrix.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :rtype: tuple(str)
    """

    name_mat = numpy.array(name_matrix(vma))

    return tuple(more_itertools.unique_everseen(name_mat[1:, 0]))


def central_angle_names(vma):
    """ Obtain names of all central-angle coordinates defined in the V-Matrix.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :rtype: tuple(str)
    """

    name_mat = numpy.array(name_matrix(vma))

    return tuple(more_itertools.unique_everseen(name_mat[2:, 1]))


def dihedral_angle_names(vma):
    """ Obtain names of all dihedral angle coordinates defined in the V-Matrix.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :rtype: tuple(str)
    """

    name_mat = numpy.array(name_matrix(vma))

    return tuple(more_itertools.unique_everseen(name_mat[3:, 2]))


def angle_names(vma):
    """ Obtain names of all angle coordinates defined in the V-Matrix.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :rtype: tuple(str)
    """
    return tuple(itertools.chain(central_angle_names(vma),
                                 dihedral_angle_names(vma)))


def dummy_coordinate_names(vma):
    """ Obtain names of all coordinates associated with dummy atoms
        defined in the V-Matrix.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :rtype: tuple(str)
    """

    symbs = symbols(vma)
    name_mat = numpy.array(name_matrix(vma))
    dummy_keys = [idx for idx, sym in enumerate(symbs)
                  if not ptab.to_number(sym)]
    dummy_names = []
    for dummy_key in dummy_keys:
        for col_idx in range(3):
            dummy_name = next(filter(lambda x: x is not None,
                                     name_mat[dummy_key:, col_idx]))
            dummy_names.append(dummy_name)

    dummy_names = tuple(dummy_names)

    return dummy_names


# value setters
def set_key_matrix(vma, key_mat):
    """ Re-set the key matrix of a V-Matrix using the input.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param key_mat: key matrix of V-Matrix coordinate keys
        :type key_mat: tuple(tuple(int))
        :rtype: tuple(str)
    """

    symbs = symbols(vma)
    name_mat = name_matrix(vma)
    vma = automol.create.vmat.from_data(symbs, key_mat, name_mat)

    return vma


def set_name_matrix(vma, name_mat):
    """ Re-set the name matrix of a V-Matrix using the input.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param name_mat: name matrix of V-Matrix coordinate names
        :type name_mat: tuple(tuple(int))
        :rtype: tuple(str)
    """

    symbs = symbols(vma)
    key_mat = key_matrix(vma)
    vma = automol.create.vmat.from_data(symbs, key_mat, name_mat)

    return vma


def standard_name_matrix(vma, shift=0):
    """ Builds a name matrix of the V-Matrix where all of the
        coordinate names have been standardized:
            RN: (1<=N<=Ncoords)
            AN: (2<=N<=Ncoords)
            DN: (1<=N<=Ncoords)

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param shift: value to shift the keys by when obtaining the keys
        :type shift: int
        :rtype: tuple(tuple(str))
    """

    natms = count(vma)

    name_mat = numpy.array(name_matrix(vma), dtype=numpy.object_)
    name_mat[1:, 0] = [
        'R{:d}'.format(num + shift) for num in range(1, natms)]
    name_mat[2:, 1] = [
        'A{:d}'.format(num + shift) for num in range(2, natms)]
    name_mat[3:, 2] = [
        'D{:d}'.format(num + shift) for num in range(3, natms)]

    name_mat = tuple(map(tuple, name_mat))

    return name_mat


def standard_form(vma, shift=0):
    """ Build a V-Matrix where all of the coordinate names of an input V-Matrix
        have been put into standard form:
            RN: (1<=N<=Ncoords)
            AN: (2<=N<=Ncoords)
            DN: (1<=N<=Ncoords)

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param shift: value to shift the keys by when obtaining the keys
        :type shift: int
        :rtype: automol V-Matrix data strucutre
    """
    name_mat = standard_name_matrix(vma, shift=shift)
    return set_name_matrix(vma, name_mat)


def rename(vma, name_dct):
    """ Rename a subset of the coordinates of a V-Matrix.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param name_dct: mapping from old coordinate names to new ones
        :type name_dct: dict[str: str]
        :rtype: automol V-Matrix data strucutre
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
    """ Build a dictionary that can mas the coordinate names
        of the input V-Matrix to their name in a standard-form V-Matrix:
            RN: (1<=N<=Ncoords)
            AN: (2<=N<=Ncoords)
            DN: (1<=N<=Ncoords)

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param shift: value to shift the keys by when obtaining the keys
        :type shift: int
        :rtype: dict[str: str]
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
def add_atom(vma, symb, key_row, name_row=None, one_indexed=False):
    """ Add an atom to a V-Matrix.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param symb: symbol of atom to add
        :type symb: str
        :param key_row: row of keys to define new atom added to key matrix
        :type key_row: tuple(int)
        :param name_row: row of names to define new atom added to name matrix
        :type name_row: tuple(str)
        :param one_indexed: parameter to store keys in one-indexing
        :type one_indexed: bool
        :rtype: automol V-Matrix data structure
    """

    symbs = symbols(vma)
    symbs += (symb,)

    key_mat = key_matrix(vma, shift=(1 if one_indexed else 0))
    key_mat += (key_row,)

    name_mat = None if name_row is None else name_matrix(vma) + (name_row,)

    vma = automol.create.vmat.from_data(
        symbs, key_mat, name_mat, one_indexed=one_indexed)

    return vma


def remove_atom(vma, key):
    """ Remove an atom from a V-Matrix. Error raised if attempting
        to remove atom other atoms depend on.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param key: key of atom to remove
        :type key: str
        :rtype: automol V-Matrix data structure
    """

    symbs = list(symbols(vma))
    symbs.pop(key)

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

    vma = automol.create.vmat.from_data(symbs, key_mat, name_mat)

    return vma


def is_valid(vma):
    """ Assess if a V-Matrix has proper structure.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :rtype: bool
    """

    ret = True
    try:
        assert _is_sequence_of_triples(vma)
        symbs, key_mat, name_mat = zip(*vma)
        from_data(symbs, key_mat, name_mat)
    except AssertionError:
        ret = False

    return ret


def _is_sequence_of_triples(obj):
    """ Assess if input object sequence has length of three.

        :param obj: object with __len__ attribute
        :type obj: list, tuple, dict
        :rtype: bool
    """

    ret = hasattr(obj, '__len__')
    if ret:
        ret = all(hasattr(item, '__len__') and len(item) == 3 for item in obj)

    return ret


def is_standard_form(vma):
    """ Assesses if the names of the V-Matrix are in standard form:
            RN: (1<=N<=Ncoords)
            AN: (2<=N<=Ncoords)
            DN: (1<=N<=Ncoords)

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :rtype: bool
    """
    return names(vma) == names(standard_form(vma))


# I/O
def from_string(vma_str):
    """ Parse a V-Matrix object from a string.

        :param vma_str: string containing a V-Matrix
        :type vma_str: str
        :rtype: automol V-Matrix data structure
    """

    symbs, key_mat, name_mat = ar.vmatrix.matrix.read(vma_str)
    vma = from_data(symbs, key_mat, name_mat, one_indexed=True)

    return vma


def string(vma, one_indexed=True):
    """ Write a V-Matrix object to a string.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param one_indexed: parameter to write keys in one-indexing
        :type one_indexed: bool
        :rtype: str
    """

    shift = 1 if one_indexed else 0
    vma_str = aw.vmatrix.matrix_block(
        symbs=symbols(vma),
        key_mat=key_matrix(vma, shift=shift),
        name_mat=name_matrix(vma),
    )

    return vma_str
