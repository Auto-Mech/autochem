""" core functionality
"""

import itertools

import numpy
import pyparsing as pp
from pyparsing import pyparsing_common as ppc

from phydat import phycon

from ... import util, vmat
from ...geom import base as geom_base
from ...graph import base as graph_base


# # constructors
def from_data(
    symbs,
    key_mat,
    val_mat,
    name_mat=None,
    one_indexed=False,
    angstrom=False,
    degree=False,
):
    """Z-Matrix constructor

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
    """Build a Z-Matrix from a V-Matrix and a molecular geometry.

    :param vma: V-Matrix
    :type vma: automol V-Matrix data structure
    :param geo: molecular geometry
    :type geo: automol molecular geometry data structure
    :rtype: automol Z-Matrix data structure
    """
    assert vmat.symbols(vma) == vmat.symbols(geo), (
        f"\n{vmat.symbols(vma)} !=\n{vmat.symbols(geo)}"
    )

    symbs = vmat.symbols(vma)
    key_mat = vmat.key_matrix(vma)
    name_mat = vmat.name_matrix(vma)
    val_mat = numpy.empty(numpy.shape(key_mat), dtype=object)

    for row, key_row in enumerate(key_mat):
        if row > 0:
            val_mat[row, 0] = geom_base.distance(geo, row, *key_row[:1])
        if row > 1:
            val_mat[row, 1] = geom_base.central_angle(geo, row, *key_row[:2])
        if row > 2:
            val_mat[row, 2] = geom_base.dihedral_angle(geo, row, *key_row[:3])

    zma = from_data(symbs, key_mat, val_mat, name_mat)
    return zma


# # getters
def value_matrix(zma, angstrom=False, degree=False):
    """Obtain the value matrix of the Z-Matrix that contains the
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

        val_mat = [list(row) + [None] * (3 - len(row)) for row in val_mat]
        val_mat = numpy.array(val_mat, dtype=object)

        val_mat[1:, 0] *= phycon.BOHR2ANG if angstrom else 1
        val_mat[2:, 1] *= phycon.RAD2DEG if degree else 1
        val_mat[3:, 2] *= phycon.RAD2DEG if degree else 1

        tril_idxs = numpy.tril_indices(val_mat.shape[0], -1, m=3)
        val_mat[tril_idxs] = val_mat[tril_idxs].astype(float)
    else:
        val_mat = ()

    return tuple(map(tuple, val_mat))


def value_dictionary(zma, angstrom=False, degree=False):
    """Obtain the values of the coordinates defined in the Z-Matrix
    in the form of a dictionary.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param angstrom: parameter to control Bohr->Angstrom conversion
    :type angstrom: bool
    :param degree: parameter to control radian->degree conversion
    :type degree: bool
    :rtype: dict[str: tuple(float)]
    """

    names = numpy.ravel(vmat.name_matrix(zma))
    vals = numpy.ravel(value_matrix(zma, angstrom=angstrom, degree=degree))
    val_dct = dict(zip(names, vals))

    # Remove None entries from the None in name mat, convert npflatt to float
    val_dct.pop(None)
    val_dct = {name: float(val) for name, val in val_dct.items()}

    return val_dct


def value(zma, name: str, angstrom: bool = False, degree: bool = False) -> float:
    """Obtain the value of a Z-Matrix coordinate by name

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param name: The coordinate name
    :type name: str
    :param angstrom: parameter to control Bohr->Angstrom conversion
    :type angstrom: bool
    :param degree: parameter to control radian->degree conversion
    :type degree: bool
    :rtype: dict[str: tuple(float)]
    """
    val_dct = value_dictionary(zma, angstrom=angstrom, degree=degree)
    return val_dct[name]


# # setters
def set_key_matrix(zma, key_mat):
    """Re-set the key matrix of a Z-Matrix using the input key matrix.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param key_mat: matrix of Z-Matrix coordinate keys
    :type key_mat: tuple(tuple(int))
    :rtype: automol Z-Matrix data structure
    """

    symbs = vmat.symbols(zma)
    val_mat = value_matrix(zma)
    name_mat = vmat.name_matrix(zma)
    zma = from_data(symbs, key_mat, val_mat, name_mat)

    return zma


def set_name_matrix(zma, name_mat):
    """Re-set the name matrix of a Z-Matrix using the input name matrix.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param name_mat: matrix of Z-Matrix coordinate names
    :type name_mat: tuple(tuple(int))
    :rtype: automol Z-Matrix data structure
    """

    symbs = vmat.symbols(zma)
    val_mat = value_matrix(zma)
    key_mat = vmat.key_matrix(zma)
    zma = from_data(symbs, key_mat, val_mat, name_mat)

    return zma


def set_value_matrix(zma, val_mat):
    """Re-set the name matrix of a Z-Matrix using the input value matrix.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param val_mat: matrix of Z-Matrix coordinate values
    :type val_mat: tuple(tuple(int))
    :rtype: automol Z-Matrix data structure
    """

    symbs = vmat.symbols(zma)
    key_mat = vmat.key_matrix(zma)
    name_mat = vmat.name_matrix(zma)
    zma = from_data(symbs, key_mat, val_mat, name_mat)

    return zma


def set_values_by_name(zma, val_dct, angstrom=True, degree=True):
    """Re-set the name matrix of a Z-Matrix using the input value dictionary.

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
    name_mat = numpy.array(vmat.name_matrix(zma), dtype=object)

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
def string(zma, one_indexed=False, angstrom=True, degree=True):
    """Write a Z-Matrix object to a string.

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

    char_dct = {"R": 0, "A": 1, "D": 2}

    def _sort_priority(arg):
        """return a sort priority value for z-matrix variable names"""
        name, _ = arg
        char, num = name[0], name[1:]
        char_val = char_dct[char] if char in char_dct else 99
        num_val = int(num) if num.isdigit() else numpy.inf
        return (char_val, num_val)

    items = sorted(val_dct.items(), key=_sort_priority)

    setval_str = "\n".join([f"{name:<5s}={val:>11.6f}" for name, val in items])

    # 3. Join them together
    zma_str = "\n\n".join((vma_str, setval_str))

    return zma_str


def from_string(zma_str, one_indexed=None, angstrom=True, degree=True):
    """Parse a Z-Matrix object from a string.

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
    val_str = "\n".join(zma_str.strip().splitlines()[nrows:])

    value_line = pp.Group(vmat.VNAME + pp.Suppress("=") + ppc.fnumber)
    value_lines = pp.delimitedList(value_line, delim=pp.lineEnd())

    if nrows > 1:
        val_dct = dict(value_lines.parseString(val_str).asList())
        val_dct[None] = None
        val_mat = [list(map(val_dct.__getitem__, nrow)) for nrow in name_mat]
    else:
        assert nrows == 1, f"Failed to parse z-matrix string: {zma_str}"
        val_mat = [(None, None, None)]

    zma = from_data(
        symbs,
        key_mat,
        val_mat,
        name_mat,
        one_indexed=False,
        angstrom=angstrom,
        degree=degree,
    )

    return zma


def yaml_data(zma) -> list:
    """A yaml-friendly data format for the z-matrix

    :param zma: molecular z-matrix
    :type zma: automol molecular z-matrix data structure
    :returns: A yaml-formatted molecular z-matrix
    :rtype: list
    """
    zma = round_(zma)
    symbs = vmat.symbols(zma)
    key_mat = vmat.key_matrix(zma)
    val_mat = value_matrix(zma)
    zma_yml = [
        [s, *itertools.chain(*zip(k, v))] for s, k, v in zip(symbs, key_mat, val_mat)
    ]
    return zma_yml


def from_yaml_data(zma_yml):
    """Put a yaml-formatted z-matrix back into standard format

    :param zma_yml: A yaml-formatted molecular z-matrix
    :type zma_yml: list
    :returns: molecular z-matrix
    :rtype: automol molecular z-matrix data structure
    """
    symbs = [row[0] for row in zma_yml]
    key_mat = [row[1::2] for row in zma_yml]
    val_mat = [row[2::2] for row in zma_yml]
    return from_data(symbs, key_mat=key_mat, val_mat=val_mat)


# # validation
def is_valid(zma):
    """Assess if a Z-Matrix has proper data structure.

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


# # properties
def torsion_coordinates(
    zma, gra, with_h_rotors: bool = True, with_ch_rotors: bool = True
):
    """Get the names and indices of torsional coordinates, as a dictionary

    :param zma: A z-matrix
    :type zma: automol zmat data structure
    :param gra: A graph, specifying connectivity for the z-matrix
    :type gra: automol graph data structure
    :param with_h_rotors: Include H rotors?
    :type with_h_rotors: bool
    :param with_ch_rotors: Include CH rotors?
    :type with_ch_rotors: bool
    """
    # 1. Identify rotational bond keys from graph
    dum_src_dct = vmat.dummy_source_dict(zma, dir_=False)
    lin_keys = list(dum_src_dct.values())
    rot_bkeys = graph_base.rotational_bond_keys(
        gra,
        lin_keys=lin_keys,
        with_h_rotors=with_h_rotors,
        with_ch_rotors=with_ch_rotors,
    )

    # 2. Identify the first dihedral coordinate for each bond
    dih_names = vmat.dihedral_angle_names(zma)
    ckey_dct = vmat.coordinates(zma, multi=False)
    bkey_dct = {n: frozenset(k[1:3]) for n, k in ckey_dct.items()}

    tors_dct = {}
    for bkey in rot_bkeys:
        # The names are ordered, so this is guaranteed to find the first one
        name = next(n for n in dih_names if bkey_dct[n] == bkey)
        tors_dct[name] = ckey_dct[name]

    return tors_dct


# # conversions
def vmatrix(zma):
    """Parse and return the V-Matrix component of a Z-Matrix.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    """
    return vmat.from_data(
        symbs=vmat.symbols(zma),
        key_mat=vmat.key_matrix(zma),
        name_mat=vmat.name_matrix(zma),
    )


def formula(zma):
    """Generate a stoichiometric formula dictionary from a Z-Matrix.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :type: dict[str: int]
    """

    syms = vmat.symbols(zma)
    fml = util.formula_from_symbols(syms)

    return fml


# # relabelling
def rename(zma, name_dct):
    """Rename a subset of the coordinates of a Z-Matrix.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param name_dct: mapping from old coordinate names to new ones
    :type name_dct: dict[str: str]
    :rtype: automol Z-Matrix data strucutre
    """

    symbs = vmat.symbols(zma)
    key_mat = vmat.key_matrix(zma)
    val_mat = value_matrix(zma)

    vma = vmat.rename(zma, name_dct)
    name_mat = vmat.name_matrix(vma)

    zma = from_data(symbs, key_mat, val_mat, name_mat)

    return zma


def standard_form(zma, shift=0):
    """Build a Z-Matrix where all of the coordinate names of an input Z-Matrix
    have been put into standard form:
        RN: (1<=N<=Ncoords)
        AN: (2<=N<=Ncoords)
        DN: (1<=N<=Ncoords)

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param shift: value to shift the keys by when obtaining the keys
    :type shift: int
    """
    name_mat = vmat.standard_name_matrix(zma, shift=shift)
    return set_name_matrix(zma, name_mat)


def round_(zma, decimals=6):
    """Round the coordinates of this z-matrix to a certain number of decimal
    places

    :param zma: molecular z-matrix
    :type zma: automol molecular z-matrix data structure
    :param decimals: the number of decimals to round to
    :type decimals: int
    :rtype: automol molecular z-matrix data structure
    """
    val_mat = numpy.array(value_matrix(zma), dtype=object)
    r_ = numpy.vectorize(
        lambda x: x if x is None else numpy.round(x, decimals=decimals)
    )
    return from_data(
        symbs=vmat.symbols(zma),
        key_mat=vmat.key_matrix(zma),
        val_mat=r_(val_mat),
        name_mat=vmat.name_matrix(zma),
    )


# # add/remove atoms
def add_atom(
    zma,
    sym,
    key_row,
    val_row,
    name_row=None,
    one_indexed=False,
    angstrom=True,
    degree=True,
):
    """Add an atom to a Z-Matrix.

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

    symbs = vmat.symbols(zma)
    symbs += (sym,)

    key_mat = vmat.key_matrix(zma, shift=(1 if one_indexed else 0))
    key_mat += (key_row,)

    val_mat = value_matrix(zma, angstrom=angstrom, degree=degree)
    val_mat += (val_row,)

    name_mat = None if name_row is None else vmat.name_matrix(zma) + (name_row,)

    zma = from_data(
        symbs,
        key_mat,
        val_mat,
        name_mat,
        one_indexed=one_indexed,
        angstrom=angstrom,
        degree=degree,
    )

    return zma


def remove_atom(zma, key):
    """Remove an atom from a Z-Matrix. Error raised if attempting
    to remove atom other atoms depend on.

    :param zma: Z-Matrix
    :type zma: automol Z-Matrix data structure
    :param key: key of atom to remove
    :type key: str
    :rtype: automol Z-Matrix data structure
    """

    symbs = list(vmat.symbols(zma))
    symbs.pop(key)

    key_mat = list(vmat.key_matrix(zma))
    key_mat.pop(key)
    key_mat = numpy.array(key_mat, dtype=object)

    if (key_mat == key).any():
        raise ValueError(f"Other atoms in z-matrix depend on atom {key}")

    key_map = numpy.vectorize(lambda x: x if (x is None or x < key) else x - 1)
    key_mat = key_map(key_mat)

    val_mat = list(value_matrix(zma))
    val_mat.pop(key)

    name_mat = list(vmat.name_matrix(zma))
    name_mat.pop(key)

    zma = from_data(symbs, key_mat, val_mat, name_mat)

    return zma


# # comparisons
def almost_equal(zma1, zma2, dist_rtol=2e-5, ang_atol=2e-3, just_dist=False):
    """Assss if two Z-Matrices are numerically equal.

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
                ang_vals1 = numpy.mod(ang_vals1, 2 * numpy.pi)
                ang_vals2 = numpy.mod(ang_vals2, 2 * numpy.pi)
                ang_diffs = numpy.abs(ang_vals1 - ang_vals2)
                ang_diffs = numpy.pi - numpy.abs(ang_diffs - numpy.pi)

                if numpy.allclose(ang_diffs, 0.0, atol=ang_atol):
                    ret = True

    return ret


# # coordinate names
# # helpers
def _value_matrix(val_mat, angstrom, degree):
    """Format value matrix of the V-Matrix that contains the
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
    val_mat = [list(row) + [None] * (3 - len(row)) for row in val_mat]
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
    """Assess if input object sequence has length of four.

    :param obj: object with __len__ attribute
    :type obj: list, tuple, dict
    :rtype: bool
    """

    ret = hasattr(obj, "__len__")
    if ret:
        ret = all(hasattr(item, "__len__") and len(item) == 4 for item in obj)

    return ret
