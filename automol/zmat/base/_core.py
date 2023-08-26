""" core functionality
"""

import numpy
import pyparsing as pp
from pyparsing import pyparsing_common as ppc
from phydat import phycon
from automol import vmat
from automol import util
import automol.geom.base
from automol.vmat import symbols
from automol.vmat import key_matrix
from automol.vmat import name_matrix
from automol.vmat import count
from automol.vmat import atom_indices
from automol.vmat import coordinates
from automol.vmat import dihedral_angle_names
from automol.vmat import standard_name_matrix


# # constructors
def from_data(symbs, key_mat, val_mat, name_mat=None,
              one_indexed=False, angstrom=False, degree=False):
    """ Z-Matrix constructor

        :param symbs: atomic symbols
        :type symbs: tuple[str]
        :param key_mat: key/index columns of the z-matrix, zero-indexed
        :type key_mat: tuple[tuple[float, float or None, float or None]]
        :param val_mat: z-matrix coordinate values
        :type val_mat: tuple[tuple[float, float or None, float or None]]
        :param name_mat: coordinate name columns of the z-matrix
        :type name_mat; tuple[tuple[str, str or None, str or None]]
        :param one_indexed: are the keys in `key_mat` one-indexed?
        :type one_indexed: bool
        :param angstrom: are distance values in `val_mat` in angstroms?
        :type angstrom: bool
        :param degree: are angle values in `val_mat` in degrees?
        :type degree: bool
        :rtype: automol Z-Matrix data structure
    """
    vma = vmat.from_data(symbs, key_mat, name_mat, one_indexed=one_indexed)
    val_mat = _value_matrix(val_mat, angstrom=angstrom, degree=degree)

    symbs, key_mat, name_mat = zip(*vma)
    zma = tuple(zip(symbs, key_mat, name_mat, val_mat))

    return zma


def from_geometry(vma, geo):
    """  Build a Z-Matrix from a V-Matrix and a molecular geometry.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :rtype: automol Z-Matrix data structure
    """
    assert symbols(vma) == symbols(geo)

    symbs = symbols(vma)
    key_mat = key_matrix(vma)
    name_mat = name_matrix(vma)
    val_mat = numpy.empty(numpy.shape(key_mat), dtype=object)

    for row, key_row in enumerate(key_mat):
        if row > 0:
            val_mat[row, 0] = automol.geom.base.distance(
                    geo, row, *key_row[:1])
        if row > 1:
            val_mat[row, 1] = automol.geom.base.central_angle(
                    geo, row, *key_row[:2])
        if row > 2:
            val_mat[row, 2] = automol.geom.base.dihedral_angle(
                    geo, row, *key_row[:3])

    zma = from_data(symbs, key_mat, val_mat, name_mat)
    return zma


# # getters
def value_matrix(zma, angstrom=False, degree=False):
    """ Obtain the value matrix of the Z-Matrix that contains the
        values of the coordinates by row and column.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
        :rtype: tuple(tuple(str))
    """

    if zma:
        val_mat = tuple(zip(*zma))[3]

        val_mat = [list(row) + [None]*(3-len(row)) for row in val_mat]
        val_mat = numpy.array(val_mat, dtype=object)

        val_mat[1:, 0] *= phycon.BOHR2ANG if angstrom else 1
        val_mat[2:, 1] *= phycon.RAD2DEG if degree else 1
        val_mat[3:, 2] *= phycon.RAD2DEG if degree else 1
    else:
        val_mat = ()

    return tuple(map(tuple, val_mat))


def value_dictionary(zma, angstrom=False, degree=False):
    """ Obtain the values of the coordinates defined in the Z-Matrix
        in the form of a dictionary.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
        :rtype: dict[str: tuple(float)]
    """

    names = numpy.ravel(name_matrix(zma))
    vals = numpy.ravel(value_matrix(zma, angstrom=angstrom, degree=degree))
    val_dct = dict(zip(names, vals))

    # Remove None entries from the None in name mat, convert npflatt to float
    val_dct.pop(None)
    val_dct = {name: float(val) for name, val in val_dct.items()}

    return val_dct


# # setters
def set_key_matrix(zma, key_mat):
    """ Re-set the key matrix of a Z-Matrix using the input key matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param key_mat: matrix of Z-Matrix coordinate keys
        :type key_mat: tuple(tuple(int))
        :rtype: automol Z-Matrix data structure
    """

    symbs = symbols(zma)
    val_mat = value_matrix(zma)
    name_mat = name_matrix(zma)
    zma = from_data(symbs, key_mat, val_mat, name_mat)

    return zma


def set_name_matrix(zma, name_mat):
    """ Re-set the name matrix of a Z-Matrix using the input name matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param name_mat: matrix of Z-Matrix coordinate names
        :type name_mat: tuple(tuple(int))
        :rtype: automol Z-Matrix data structure
    """

    symbs = symbols(zma)
    val_mat = value_matrix(zma)
    key_mat = key_matrix(zma)
    zma = from_data(symbs, key_mat, val_mat, name_mat)

    return zma


def set_value_matrix(zma, val_mat):
    """ Re-set the name matrix of a Z-Matrix using the input value matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param val_mat: matrix of Z-Matrix coordinate values
        :type val_mat: tuple(tuple(int))
        :rtype: automol Z-Matrix data structure
    """

    symbs = symbols(zma)
    key_mat = key_matrix(zma)
    name_mat = name_matrix(zma)
    zma = from_data(symbs, key_mat, val_mat, name_mat)

    return zma


def set_values_by_name(zma, val_dct, angstrom=True, degree=True):
    """ Re-set the name matrix of a Z-Matrix using the input value dictionary.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param val_dct: dictionary of coordinate values
        :type val_dct: dict[str: float]
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
        :rtype: automol Z-Matrix data structure
    """

    val_mat = numpy.array(value_matrix(zma), dtype=object)
    name_mat = numpy.array(name_matrix(zma), dtype=object)

    for (row, col), name in numpy.ndenumerate(name_mat):
        if name in val_dct:
            val = val_dct[name]
            if col == 0:
                val *= phycon.ANG2BOHR if angstrom else 1
            else:
                assert col > 0
                val *= phycon.DEG2RAD if degree else 1
            val_mat[row, col] = val

    return set_value_matrix(zma, val_mat)


# # I/O
def string(zma, one_indexed=True, angstrom=True, degree=True):
    """ Write a Z-Matrix object to a string.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param one_indexed: parameter to store keys in one-indexing
        :type one_indexed: bool
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
        :rtype: str
    """

    # 1. Get the v-matrix string
    vma = vmatrix(zma)
    vma_str = vmat.string(vma, one_indexed=one_indexed)

    # 2. Format the set-value dictionary
    val_dct = value_dictionary(zma, angstrom=angstrom, degree=degree)

    char_dct = {'R': 0, 'A': 1, 'D': 2}

    def _sort_priority(arg):
        """ return a sort priority value for z-matrix variable names
        """
        name, _ = arg
        char, num = name[0], name[1:]
        char_val = char_dct[char] if char in char_dct else 99
        num_val = int(num) if num.isdigit() else numpy.inf
        return (char_val, num_val)

    items = sorted(val_dct.items(), key=_sort_priority)

    setval_str = '\n'.join([
        f'{name:<5s}={val:>11.6f}'
        for name, val in items])

    # 3. Join them together
    zma_str = '\n\n'.join((vma_str, setval_str))

    return zma_str


def from_string(zma_str, one_indexed=True, angstrom=True, degree=True):
    """ Parse a Z-Matrix object from a string.

        :param zma_str: string containing a Z-Matrix
        :type zma_str: str
        :param one_indexed: parameter to store keys in one-indexing
        :type one_indexed: bool
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
        :rtype: automol Z-Matrix data structure
    """
    vma = vmat.from_string(zma_str, one_indexed=one_indexed)
    symbs = vmat.symbols(vma)
    key_mat = vmat.key_matrix(vma)
    name_mat = vmat.name_matrix(vma)

    nrows = len(symbs)
    val_str = '\n'.join(zma_str.splitlines()[nrows:])

    value_line = pp.Group(vmat.VNAME + pp.Suppress('=') + ppc.fnumber)
    value_lines = pp.delimitedList(value_line, delim=pp.lineEnd())

    val_dct = dict(value_lines.parseString(val_str).asList())
    val_dct[None] = None
    val_mat = [list(map(val_dct.__getitem__, nrow)) for nrow in name_mat]

    zma = from_data(
        symbs, key_mat, val_mat, name_mat, one_indexed=False,
        angstrom=angstrom, degree=degree)

    return zma


# # validation
def is_valid(zma):
    """ Assess if a Z-Matrix has proper data structure.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: bool
    """

    ret = True

    try:
        assert _is_sequence_of_quadruples(zma)
        symbs, key_mat, name_mat, val_mat = zip(*zma)
        from_data(symbs, key_mat, val_mat, name_mat)
    except AssertionError:
        ret = False

    return ret


# # conversions
def vmatrix(zma):
    """ Parse and return the V-Matrix component of a Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
    """
    return vmat.from_data(symbs=symbols(zma),
                          key_mat=key_matrix(zma),
                          name_mat=name_matrix(zma))


def formula(zma):
    """ Generate a stoichiometric formula dictionary from a Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :type: dict[str: int]
    """

    syms = symbols(zma)
    fml = util.formula_from_symbols(syms)

    return fml


# # relabelling
def rename(zma, name_dct):
    """ Rename a subset of the coordinates of a Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param name_dct: mapping from old coordinate names to new ones
        :type name_dct: dict[str: str]
        :rtype: automol Z-Matrix data strucutre
    """

    symbs = symbols(zma)
    key_mat = key_matrix(zma)
    val_mat = value_matrix(zma)

    vma = vmat.rename(zma, name_dct)
    name_mat = vmat.name_matrix(vma)

    zma = from_data(symbs, key_mat, val_mat, name_mat)

    return zma


def standard_form(zma, shift=0):
    """ Build a Z-Matrix where all of the coordinate names of an input Z-Matrix
        have been put into standard form:
            RN: (1<=N<=Ncoords)
            AN: (2<=N<=Ncoords)
            DN: (1<=N<=Ncoords)

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param shift: value to shift the keys by when obtaining the keys
        :type shift: int
    """
    name_mat = standard_name_matrix(zma, shift=shift)
    return set_name_matrix(zma, name_mat)


# # add/remove atoms
def add_atom(zma, sym, key_row, val_row, name_row=None,
             one_indexed=False, angstrom=True, degree=True):
    """ Add an atom to a Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param symb: symbol of atom to add
        :type symb: str
        :param key_row: row of keys to define new atom added to key matrix
        :type key_row: tuple(int)
        :param val_row: row of values to define new atom added to name matrix
        :type val_row: tuple(float)
        :param name_row: row of names to define new atom added to name matrix
        :type name_row: tuple(str)
        :param one_indexed: parameter to store keys in one-indexing
        :type one_indexed: bool
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
        :rtype: automol Z-Matrix data structure
    """

    symbs = symbols(zma)
    symbs += (sym,)

    key_mat = key_matrix(zma, shift=(1 if one_indexed else 0))
    key_mat += (key_row,)

    val_mat = value_matrix(zma, angstrom=angstrom, degree=degree)
    val_mat += (val_row,)

    name_mat = None if name_row is None else name_matrix(zma) + (name_row,)

    zma = from_data(
        symbs, key_mat, val_mat, name_mat, one_indexed=one_indexed,
        angstrom=angstrom, degree=degree)

    return zma


def remove_atom(zma, key):
    """ Remove an atom from a Z-Matrix. Error raised if attempting
        to remove atom other atoms depend on.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param key: key of atom to remove
        :type key: str
        :rtype: automol Z-Matrix data structure
    """

    symbs = list(symbols(zma))
    symbs.pop(key)

    key_mat = list(key_matrix(zma))
    key_mat.pop(key)
    key_mat = numpy.array(key_mat, dtype=object)

    if (key_mat == key).any():
        raise ValueError(f"Other atoms in z-matrix depend on atom {key}")

    key_map = numpy.vectorize(lambda x: x if (x is None or x < key) else x-1)
    key_mat = key_map(key_mat)

    val_mat = list(value_matrix(zma))
    val_mat.pop(key)

    name_mat = list(name_matrix(zma))
    name_mat.pop(key)

    zma = from_data(symbs, key_mat, val_mat, name_mat)

    return zma


# # comparisons
def almost_equal(zma1, zma2, dist_rtol=2e-5, ang_atol=2e-3, just_dist=False):
    """ Assss if two Z-Matrices are numerically equal.

        :param zma1: The first z-matrix
        :type zma1: automol Z-Matrix data structure
        :param zma2: The second z-matrix
        :type zma2: automol Z-Matrix data structure
        :param dist_rtol: Relative tolerance for the distances
        :type dist_rtol: float
        :param ang_atol: Absolute tolerance for the angles
        :type ang_atol: float
        :param just_dist: parameter to only compare distances
        :type just_dist: bool
        :rtype: bool
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
                ang_vals2 = numpy.hstack((val_mat2[2:, 1], val_mat2[3:, 2]))
                ang_vals1 = numpy.mod(ang_vals1, 2*numpy.pi)
                ang_vals2 = numpy.mod(ang_vals2, 2*numpy.pi)
                ang_diffs = numpy.abs(ang_vals1 - ang_vals2)
                ang_diffs = numpy.pi - numpy.abs(ang_diffs - numpy.pi)

                if numpy.allclose(ang_diffs, 0., atol=ang_atol):
                    ret = True

    return ret


# # coordinate names
def distance_coordinate_name(zma, key1, key2):
    """ get the name of a distance coordinate for a given bond

        :param zma: the z-matrix
        :type zma: automol Z-Matrix data structure
        :param key1: the first key in the torsion axis (rotational bond)
        :type key1: int
        :param key2: the second key in the torsion axis (rotational bond)
        :type key2: int
        :rtype: str
    """

    key1, key2 = sorted([key1, key2])
    name_mat = name_matrix(zma)
    key_mat = key_matrix(zma)
    assert key_mat[key2][0] == key1, (
        f"{key1:d}-{key2:d} is not a distance coordinate in this zmatrix:"
        f"\n{string(zma, one_indexed=False)}")
    name = name_mat[key2][0]

    return name


def central_angle_coordinate_name(zma, key1, key2, key3):
    """ get the name of angle coordinate for a set of 3 atoms

        :param zma: the z-matrix
        :type zma: automol Z-Matrix data structure
        :param key1: the first key
        :type key1: int
        :param key2: the second key (central atom)
        :type key2: int
        :param key3: the third key
        :type key3: int
        :rtype: str
    """

    key1, key3 = sorted([key1, key3])
    name_mat = name_matrix(zma)
    key_mat = key_matrix(zma)
    assert key_mat[key3][0] == key2 and key_mat[key3][1] == key1, (
        f"{key1:d}-{key2:d}-{key3:d} "
        "is not a angle coordinate in this zmatrix:"
        f"\n{string(zma, one_indexed=False)}")
    name = name_mat[key3][1]

    return name


def dihedral_angle_coordinate_name(zma, key1, key2, key3, key4):
    """ get the name of dihedral coordinate for a set of 4 atoms

        :param zma: the z-matrix
        :type zma: automol Z-Matrix data structure
        :param key1: the first key
        :type key1: int
        :param key2: the second key
        :type key2: int
        :param key3: the third key
        :type key3: int
        :param key4: the fourth key
        :type key4: int
        :rtype: str
    """

    if key1 > key4:
        key1, key2, key3, key4 = key4, key3, key2, key1

    name_mat = name_matrix(zma)
    key_mat = key_matrix(zma)
    assert (
        key_mat[key4][0] == key3 and key_mat[key4][1] == key2 and
        key_mat[key4][2] == key1
    ), (
        f"{key1:d}-{key2:d}-{key3:d}-{key4:d} "
        "is not a dihedral coordinate in this zmatrix:"
        f"\n{string(zma, one_indexed=False)}")

    name = name_mat[key4][2]

    return name


def dihedral_axis_name(zma, axis):
    """ gives this name of a dihedral angle that has
        the given axis atoms

        currently fails if indices of axis not in zma, i.e.
        if 5,1 in zma, code fails if 1,5 given
    """
    coords = coordinates(zma)
    angles = dihedral_angle_names(zma)
    name = None
    for ang in angles:
        _coord_idxs = coords[ang]
        if (tuple(list(_coord_idxs[0])[1:3]) == axis
                or tuple(list(_coord_idxs[0])[3:1:-1]) == axis):
            name = ang
            break
    return name


# # dummy atom functions
def dummy_keys(zma):
    """ Obtain keys to dummy atoms in the Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: tuple[int]
    """
    keys = tuple(key for key, sym in enumerate(symbols(zma)) if sym == 'X')
    return keys


def dummy_key_dictionary(zma):
    """ Obtain keys to linear atoms in the Z-matrix along with their associated
        dummy atoms.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :returns: a dictionary with the linear atoms as keys and their dummy
            atoms as values
        :rtype: dict[int: int]
    """
    dummy_key_dct = dict(map(reversed, dummy_neighbor_keys(zma).items()))
    return dummy_key_dct


def dummy_neighbor_keys(zma):
    """ Obtain keys to dummy atoms in the Z-Matrix, along with their
        neighboring atoms.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :returns: a dictionary with the dummy atoms as keys and their neighbors
            as values
        :rtype: dict[int: int]
    """
    key_mat = key_matrix(zma)
    dum_keys = dummy_keys(zma)
    key_dct = {}
    for dum_key in dum_keys:
        ngb_key = key_mat[dum_key][0]
        if ngb_key is None:
            ngb_key = next(row for row, (k, _, _) in enumerate(key_mat)
                           if k == dum_key)

        key_dct[dum_key] = ngb_key
    return key_dct


def linear_atom_keys(zma, geom_indexing=False):
    """ Obtain keys to linear atoms in the Z-matrix. Any atom neighboring a
        dummy atom is considered to be linear.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param geom_indexing: use geometry indexing?
        :type geom_indexing: bool
        :returns: the linear atom keys
        :rtype: tuple[int]
    """
    lin_key_dct = dummy_neighbor_keys(zma)
    lin_keys = tuple(sorted(lin_key_dct.values()))
    if geom_indexing:
        dum_keys = numpy.array(sorted(lin_key_dct.keys()))
        lin_keys = tuple(k-sum(k > dum_keys) for k in lin_keys)
    return lin_keys


def shift_down(zma, vals):
    """
    Build a remdummy list that tells how to shift the groups
    """

    dummy_idxs = sorted(atom_indices(zma, 'X', match=True))
    if dummy_idxs:
        remdummy = [0 for _ in range(count(zma))]
        for dummy in dummy_idxs:
            for idx, _ in enumerate(remdummy):
                if dummy < idx:
                    remdummy[idx] += 1

        vals1 = tuple(val+1 for val in vals)
        vals2 = tuple(val-remdummy[val-1] for val in vals1)
        shift_vals = tuple(val-1 for val in vals2)

    else:
        shift_vals = vals

    return shift_vals


def shift_up(zma, idxs):
    """ shift up from the dummy idxs
    """

    dummy_idxs = sorted(atom_indices(zma, 'X', match=True))

    shift_idxs = []
    for idx in idxs:
        new_idx = idx
        for dummy_idx in dummy_idxs:
            if idx >= dummy_idx:
                new_idx += 1
        shift_idxs.append(new_idx)

    return tuple(shift_idxs)


# # helpers
def _value_matrix(val_mat, angstrom, degree):
    """ Format value matrix of the V-Matrix that contains the
        coordinate values by row and column.

        :param val_mat: value matrix containing coordinate values
        :type val_mat: tuple(tuple(int))
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
        :rtype: tuple(tuple(str))
    """

    # Check dimensions and ensure proper formatting
    val_mat = [list(row) + [None]*(3-len(row)) for row in val_mat]
    val_mat = numpy.array(val_mat, dtype=object)
    natms = val_mat.shape[0]

    assert val_mat.ndim == 2 and val_mat.shape == (natms, 3)
    triu_idxs = numpy.triu_indices(natms, m=3)

    val_mat[1:, 0] *= phycon.ANG2BOHR if angstrom else 1
    val_mat[2:, 1] *= phycon.DEG2RAD if degree else 1
    val_mat[3:, 2] *= phycon.DEG2RAD if degree else 1

    val_mat[triu_idxs] = None

    return tuple(map(tuple, val_mat))


def _is_sequence_of_quadruples(obj):
    """ Assess if input object sequence has length of four.

        :param obj: object with __len__ attribute
        :type obj: list, tuple, dict
        :rtype: bool
    """

    ret = hasattr(obj, '__len__')
    if ret:
        ret = all(hasattr(item, '__len__') and len(item) == 4 for item in obj)

    return ret
