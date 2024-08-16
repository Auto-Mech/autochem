"""V-Matrix: Variable V-Matrix (V-Matrix without coordinate values)."""

import itertools
from collections import defaultdict

import more_itertools
import numpy
import pyparsing as pp
from pyparsing import pyparsing_common as ppc

from phydat import ptab

from .util import ZmatConv, dict_, zmat_conv

Symbol = str
Key = int | None
Name = str | None
KeyRow = tuple[Key, Key, Key]
NameRow = tuple[Name, Name, Name]
KeyMatrix = tuple[KeyRow, ...]
NameMatrix = tuple[NameRow, ...]
VMatrixRow = tuple[Symbol, KeyRow, NameRow] | tuple[Symbol, KeyRow, NameRow, object]
VMatrix = tuple[VMatrixRow, ...]
CoordinateKey = tuple[int, ...] | None
CoordinateKeyRow = tuple[CoordinateKey, CoordinateKey, CoordinateKey]
CoordinateKeyMatrix = tuple[CoordinateKeyRow, ...]


# Build the v-xmatrix parser
CHAR = pp.Char(pp.alphas)
SYMBOL = pp.Combine(CHAR + pp.Opt(CHAR))
VNAME = pp.Combine(pp.Word(pp.alphas) + pp.Opt(pp.Word(pp.nums)))
LINE_END = pp.Suppress(pp.lineEnd())

LINE0 = pp.Group(SYMBOL)
LINE1 = pp.Group(SYMBOL + ppc.integer + VNAME)
LINE2 = pp.Group(SYMBOL + (ppc.integer + VNAME) * 2)
LINE3 = pp.Group(SYMBOL + (ppc.integer + VNAME) * 3)

LINES0 = LINE0
LINES1 = LINE0 + LINE_END + LINE1
LINES2 = LINE0 + LINE_END + LINE1 + LINE_END + LINE2
LINES3 = (
    LINE0
    + LINE_END
    + LINE1
    + LINE_END
    + LINE2
    + LINE_END
    + pp.delimitedList(LINE3, delim=pp.lineEnd())
)

VMAT_LINES = LINES0 ^ LINES1 ^ LINES2 ^ LINES3


# # constructors
def from_data(
    symbs: tuple[str, ...],
    key_mat: KeyMatrix,
    name_mat: NameMatrix = None,
    one_indexed: bool | None = None,
) -> VMatrix:
    """V-Matrix constructor (V-Matrix without numerical coordinate values).

    :param symbs: Atomic symbols
    :param key_mat: Key/index columns of the v-matrix, zero-indexed
    :param name_mat: Coordinate name columns of the v-matrix
    :param one_indexed: Parameter to store keys in one-indexing
    :return: Automol V-Matrix data structure
    """
    symbs = list(map(ptab.to_symbol, symbs))
    natms = len(symbs)

    key_mat = _key_matrix(key_mat, natms, one_indexed)
    name_mat = _name_matrix(name_mat, natms)

    vma = tuple(zip(symbs, key_mat, name_mat, strict=False))

    return vma


# # V-Matrix/V-Matrix common functions (document these as z-matrix functions)
# # # getters
def symbols(vma: VMatrix, idxs: list[int] | None = None) -> tuple[Symbol, ...]:
    """Obtain the atomic symbols for all atoms defined in the V-Matrix.

    :param vma: V-Matrix
    :param idxs: Indices of atoms to obtain information for
    :return: The list of atomic symbols
    """
    if vma:
        symbs, *_ = tuple(zip(*vma, strict=False))
    else:
        symbs = ()

    return symbs if idxs is None else tuple(map(symbs.__getitem__, idxs))


def key_matrix(vma: VMatrix, shift: int = 0) -> KeyMatrix:
    """Obtain the key matrix of the V-Matrix that contains the
    coordinate atom keys by row and column.

    :param vma: V-Matrix
    :param shift: Value to shift the keys by when obtaining the key matrix
    :return: Key matrix
    """
    if vma:
        key_mat = tuple(zip(*vma, strict=True))[1]

        # post-processing for adding the shift
        key_mat = [list(row) + [None] * (3 - len(row)) for row in key_mat]
        key_mat = numpy.array(key_mat)
        tril_idxs = numpy.tril_indices(key_mat.shape[0], -1, m=3)
        key_mat[tril_idxs] += shift
        key_mat[tril_idxs] = key_mat[tril_idxs].astype(int)
    else:
        key_mat = ()

    return tuple(map(tuple, key_mat))


def name_matrix(vma: VMatrix) -> NameMatrix:
    """Obtain the name matrix of the V-Matrix that contains the
    coordinate names by row and column.

    :param vma: V-Matrix
    :type vma: Automol V-Matrix data structure
    :return: Name matrix
    """
    if vma:
        name_mat = tuple(zip(*vma, strict=True))[2]
    else:
        name_mat = ()

    name_mat = [list(row) + [None] * (3 - len(row)) for row in name_mat]

    return tuple(map(tuple, name_mat))


# # # properties
def count(vma: VMatrix) -> int:
    """Obtain the number of rows of the V-Matrix, which corresponds to
    the number of atoms defined in the V-Matrix. This includes all
    real and dummy atoms.

    :param vma: V-Matrix
    :type vma: Automol V-Matrix data structure
    :return: Number of rows
    """
    return len(symbols(vma))


def atom_indices(vma: VMatrix, symb: Symbol, match: bool = True) -> tuple[int]:
    """Obtain the indices of a atoms of a particular type in the v-matrix.

    :param vma: V-Matrix
    :param match: Grab idxs that match given atom type
    :param symb: Atomic symbol
    :param match: Obtain indices of symbols that match the type?
    :return: Indices
    """
    symbs = symbols(vma)
    idxs = ()
    for idx, symb_ in enumerate(symbs):
        if symb_ == symb and match:
            idxs += (idx,)
        elif symb_ != symb and not match:
            idxs += (idx,)

    return idxs


def coordinate_key_matrix(vma: VMatrix, shift: int = 0) -> CoordinateKeyMatrix:
    """Obtain the coordinate key matrix of the V-Matrix that contains the
    coordinate keys by row and column.

    :param vma: V-Matrix
    :param shift: value to shift the keys by when obtaining the key matrix
    :return: Coordinate key matrix
    """
    key_mat = key_matrix(vma, shift=shift)
    natms = len(key_mat)
    atm_keys = range(shift, natms + shift)
    coo_key_mat = [
        [
            (atm_key,) + key_row[: col + 1] if key_row[col] is not None else None
            for col in range(3)
        ]
        for atm_key, key_row in zip(atm_keys, key_mat, strict=True)
    ]

    return tuple(map(tuple, coo_key_mat))


def coordinates(
    vma: VMatrix, shift: int = 0, multi: bool = True
) -> dict[Name, CoordinateKey]:
    """Obtain the coordinate keys associated with each coordinate name,
    as a dictionary. Values are sequences of coordinate keys,
    since there may be multiple.

    :param vma: V-Matrix
    :param shift: Value to shift the keys by when obtaining the keys
    :param multi: Parameter to grab multiple coordinate keys
    :return: Dictionary of coordinate keys
    """
    _names = numpy.ravel(name_matrix(vma))
    coo_keys = numpy.ravel(numpy.array(coordinate_key_matrix(vma, shift), dtype=object))

    if not multi:
        coo_dct = dict(zip(_names, coo_keys, strict=True))
    else:
        coo_dct = {name: () for name in _names}
        for name, coo_key in zip(_names, coo_keys, strict=False):
            coo_dct[name] += (coo_key,)

    coo_dct.pop(None)

    return coo_dct


def distance_coordinates(vma: VMatrix) -> dict[Name, CoordinateKey]:
    """Get the distance coordinates by coordinate name.

    :param vma: V-Matrix
    :return: The coordinates, by coordinate name
    """
    return dict_.by_key(coordinates(vma, multi=False), distance_names(vma))


def central_angle_coordinates(vma: VMatrix) -> dict[Name, CoordinateKey]:
    """Get the central angle coordinates by coordinate name.

    :param vma: V-Matrix
    :return: The coordinates, by coordinate name
    """
    return dict_.by_key(coordinates(vma, multi=False), central_angle_names(vma))


def dihedral_angle_coordinates(vma: VMatrix) -> dict[Name, CoordinateKey]:
    """Get the dihedral angle coordinates by coordinate name.

    :param vma: V-Matrix
    :return: The coordinates, by coordinate name
    """
    return dict_.by_key(coordinates(vma, multi=False), dihedral_angle_names(vma))


def coordinate(vma: VMatrix, name: Name) -> CoordinateKey:
    """Get the atom keys defining a coordinate by name.

    :param vma: A v-matrix or z-matrix
    :param name: The name of the coordinate, e.g. "R5"
    :return: The atom keys defining the coordinate
    """
    coo, *_ = coordinates(vma)[name]
    return coo


def torsion_axis(vma: VMatrix, dih_name: Name) -> tuple[int, int]:
    """Get the rotational axis of a torsion from the dihedral angle name.

    :param vma: A v-matrix or z-matrix
    :param dih_name: The dihedral angle name of the torsion
    :return: The axis, i.e. the central two atoms in the coordinate
    """
    dih_coo = coordinate(vma, dih_name)
    assert len(dih_coo) == 4, f"{dih_name} is not a dihedral angle:\n{vma}"
    _, ax_key1, ax_key2, _ = dih_coo
    return ax_key1, ax_key2


# # # names and standard naming
def names(vma: VMatrix) -> tuple[Name, ...]:
    """Obtain names of all coordinates defined in the V-Matrix.

    :param vma: V-Matrix
    :return: Names of coordinates
    """
    name_mat = name_matrix(vma)
    _names = filter(lambda x: x is not None, numpy.ravel(numpy.transpose(name_mat)))

    return tuple(more_itertools.unique_everseen(_names))


def distance_names(vma: VMatrix) -> tuple[str]:
    """Obtain names of all distance coordinates defined in the V-Matrix.

    :param vma: V-Matrix
    :return: Names of defined coordinates
    """
    name_mat = numpy.array(name_matrix(vma))

    return tuple(more_itertools.unique_everseen(name_mat[1:, 0]))


def central_angle_names(vma: VMatrix) -> tuple[Name, ...]:
    """Obtain names of all central-angle coordinates defined in the V-Matrix.

    :param vma: V-Matrix
    :return: Names of central angle coordinates
    """
    name_mat = numpy.array(name_matrix(vma))

    return tuple(more_itertools.unique_everseen(name_mat[2:, 1]))


def dihedral_angle_names(vma: VMatrix) -> tuple[Name, ...]:
    """Obtain names of all dihedral angle coordinates defined in the V-Matrix.

    :param vma: V-Matrix
    :return: Names of dihedral angles
    """
    name_mat = numpy.array(name_matrix(vma))

    return tuple(more_itertools.unique_everseen(name_mat[3:, 2]))


def angle_names(vma: VMatrix) -> tuple[Name, ...]:
    """Obtain names of all angle coordinates defined in the V-Matrix.

    :param vma: V-Matrix
    :return: Name of angle coordinates
    """
    return tuple(itertools.chain(central_angle_names(vma), dihedral_angle_names(vma)))


def standard_names(vma: VMatrix, shift: int = 0) -> dict[Name, Name]:
    """Build a dictionary that can mas the coordinate names
    of the input V-Matrix to their name in a standard-form V-Matrix:
        RN: (1<=N<=Ncoords)
        AN: (2<=N<=Ncoords)
        DN: (1<=N<=Ncoords).

    :param vma: V-Matrix
    :param shift: Value to shift the keys by when obtaining the keys
    :return: Dictionary with names
    """
    dist_names = distance_names(vma)
    cent_ang_names = central_angle_names(vma)
    dih_ang_names = dihedral_angle_names(vma)
    name_dct = {}
    name_dct.update(
        {dist_name: f"R{num+shift+1:d}" for num, dist_name in enumerate(dist_names)}
    )
    name_dct.update(
        {
            cent_ang_name: f"A{num+shift+2:d}"
            for num, cent_ang_name in enumerate(cent_ang_names)
        }
    )
    name_dct.update(
        {
            dih_ang_name: f"D{num+shift+3:d}"
            for num, dih_ang_name in enumerate(dih_ang_names)
        }
    )

    return name_dct


def standard_name_matrix(vma: VMatrix, shift: int = 0) -> NameMatrix:
    """Build a name matrix of the V-Matrix where all of the
    coordinate names have been standardized:
        RN: (1<=N<=Ncoords)
        AN: (2<=N<=Ncoords)
        DN: (1<=N<=Ncoords).

    :param vma: V-Matrix
    :param shift: value to shift the keys by when obtaining the keys
    :return: Name matrix
    """
    natms = count(vma)

    name_mat = numpy.array(name_matrix(vma), dtype=object)
    name_mat[1:, 0] = [f"R{num+shift:d}" for num in range(1, natms)]
    name_mat[2:, 1] = [f"A{num+shift:d}" for num in range(2, natms)]
    name_mat[3:, 2] = [f"D{num+shift:d}" for num in range(3, natms)]

    name_mat = tuple(map(tuple, name_mat))

    return name_mat


def distance_coordinate_name(zma: VMatrix, key1: int, key2: int) -> Name:
    """Get the name of a distance coordinate for a given bond.

    :param zma: the z-matrix
    :param key1: the first key
    :param key2: the second key
    :return: Name of coordinate
    """
    key1, key2 = sorted([key1, key2])
    name_mat = name_matrix(zma)
    key_mat = key_matrix(zma)
    assert (
        key_mat[key2][0] == key1
    ), f"{key1}-{key2} is not a coordinate in this zmatrix:\n{string(zma)}"
    name = name_mat[key2][0]

    return name


def central_angle_coordinate_name(
    zma: VMatrix, key1: int, key2: int, key3: int
) -> Name:
    """Get the name of angle coordinate for a set of 3 atoms.

    :param zma: The z-matrix
    :param key1: The first key
    :param key2: The second key (central atom)
    :param key3: The third key
    :return: Name of angle coordinates
    """
    key1, key3 = sorted([key1, key3])
    name_mat = name_matrix(zma)
    key_mat = key_matrix(zma)
    assert (
        key_mat[key3][0] == key2 and key_mat[key3][1] == key1
    ), f"{key1}-{key2}-{key3} is not a coordinate in this zmatrix:\n{string(zma)}"
    name = name_mat[key3][1]

    return name


def dihedral_angle_coordinate_name(
    zma: VMatrix, key1: int, key2: int, key3: int, key4: int
) -> Name:
    """Get the name of dihedral coordinate for a set of 4 atoms.

    :param zma: The z-matrix
    :param key1: The first key
    :param key2: The second key
    :param key3: The third key
    :param key4: The fourth key
    :return:Name for dihedral coordinate
    """
    if key1 > key4:
        key1, key2, key3, key4 = key4, key3, key2, key1

    name_mat = name_matrix(zma)
    key_mat = key_matrix(zma)
    assert (
        key_mat[key4][0] == key3
        and key_mat[key4][1] == key2
        and key_mat[key4][2] == key1
    ), f"{key1}-{key2}-{key3}-{key4} is not a coordinate in this zmat:\n{string(zma)}"

    name = name_mat[key4][2]

    return name


# # dummy atom functions
def dummy_keys(zma: VMatrix) -> tuple[Key, ...]:
    """Obtain keys to dummy atoms in the Z-Matrix.

    :param zma: Z-Matrix
    :return: Key to dummy atoms
    """
    keys = tuple(key for key, sym in enumerate(symbols(zma)) if sym == "X")
    return keys


def dummy_coordinate_names(vma: VMatrix) -> tuple[Name, ...]:
    """Obtain names of all coordinates associated with dummy atoms
    defined in the V-Matrix.

    :param vma: V-Matrix
    :return: Name of coordinates
    """
    symbs = symbols(vma)
    name_mat = numpy.array(name_matrix(vma))
    dummy_keys = [idx for idx, sym in enumerate(symbs) if not ptab.to_number(sym)]
    dummy_names = []
    for dummy_key in dummy_keys:
        for col_idx in range(3):
            dummy_name = next(
                filter(lambda x: x is not None, name_mat[dummy_key:, col_idx])
            )
            dummy_names.append(dummy_name)

    dummy_names = tuple(dummy_names)

    return dummy_names


def dummy_source_dict(
    zma: VMatrix, dir_: bool = True
) -> dict[int, int | tuple[int, int]]:
    """Obtain keys to dummy atoms in the Z-Matrix, along with their
    parent atoms.

    :param zma: Z-Matrix
    :param dir_: Include linear direction atoms? defaults to True
    :returns: A dictionary mapping dummy atoms onto their parent atoms
    """
    key_mat = key_matrix(zma)
    dum_keys = dummy_keys(zma)
    src_dct = {}
    for dum_key in dum_keys:
        lin_key, dir_key, _ = key_mat[dum_key]
        if lin_key is None:
            lin_key = next(lk for lk, (k, _, _) in enumerate(key_mat) if k == dum_key)
        if dir_key is None:
            dir_key = next(dk for dk, (_, k, _) in enumerate(key_mat) if k == dum_key)

        if not dir_:
            src_dct[dum_key] = lin_key
        else:
            src_dct[dum_key] = (lin_key, dir_key)

    return src_dct


def conversion_info(zma: VMatrix) -> ZmatConv:
    """Get the conversion information for this z-matrix, relative to geometry following
    the same atom order.

    :param zma: Z-Matrix
    :return: The z-matrix conversion
    """
    zcount = count(zma)
    src_zkeys_dct = dummy_source_dict(zma)
    return zmat_conv.from_zmat_data(zcount, src_zkeys_dct)


def neighbor_keys(vma: VMatrix) -> dict[Key, frozenset[Key]]:
    """Identify which atoms are explicit neighbors in the V-Matrix.

    :param vma: V-Matrix
    :return: A dictionary mapping atoms onto their neighbors
    """
    dist_coos = list(distance_coordinates(vma).values())
    nkeys_dct = defaultdict(set)
    for key1, key2 in dist_coos:
        nkeys_dct[key1].add(key2)
        nkeys_dct[key2].add(key1)
    return dict_.transform_values(nkeys_dct, frozenset)


# # V-Matrix-specific functions
# # # setters
def set_key_matrix(vma: VMatrix, key_mat: KeyMatrix) -> VMatrix:
    """Re-set the key matrix of a V-Matrix using the input key matrix.

    :param vma: V-Matrix
    :param key_mat: Key matrix of V-Matrix coordinate keys
    :return: VMatrix with input keys
    """
    symbs = symbols(vma)
    name_mat = name_matrix(vma)
    vma = from_data(symbs, key_mat, name_mat)

    return vma


def set_name_matrix(vma: VMatrix, name_mat: NameMatrix) -> VMatrix:
    """Re-set the name matrix of a V-Matrix using the input name matrix.

    :param vma: V-Matrix
    :param name_mat: name matrix of V-Matrix coordinate names
    :return: V-Matrix reset
    """
    symbs = symbols(vma)
    key_mat = key_matrix(vma)
    vma = from_data(symbs, key_mat, name_mat)

    return vma


# # # names and naming
def rename(vma: VMatrix, name_dct: dict[Name, Name]) -> VMatrix:
    """Rename a subset of the coordinates of a V-Matrix.

    :param vma: V-Matrix
    :param name_dct: Mapping from old coordinate names to new ones
    :return: VMatrix with new names
    """
    orig_name_mat = numpy.array(name_matrix(vma))
    tril_idxs = numpy.tril_indices(orig_name_mat.shape[0], -1, m=3)
    orig_names = set(orig_name_mat[tril_idxs])
    assert set(name_dct.keys()) <= orig_names

    name_dct.update(
        {orig_name: orig_name for orig_name in orig_names if orig_name not in name_dct}
    )

    name_mat = numpy.empty(orig_name_mat.shape, dtype=object)
    name_mat[tril_idxs] = list(map(name_dct.__getitem__, orig_name_mat[tril_idxs]))

    return from_data(symbols(vma), key_matrix(vma), name_mat)


def standard_form(vma: VMatrix, shift: int = 0) -> VMatrix:
    """Build a V-Matrix where all of the coordinate names of an input V-Matrix
    have been put into standard form:
        RN: (1<=N<=Ncoords)
        AN: (2<=N<=Ncoords)
        DN: (1<=N<=Ncoords).

    :param vma: V-Matrix
    :param shift: value to shift the keys by when obtaining the keys
    :return: V-Matrix of coordinate names
    """
    name_mat = standard_name_matrix(vma, shift=shift)
    return set_name_matrix(vma, name_mat)


# # # add/remove atoms
def add_atom(
    vma: VMatrix,
    symb: Symbol,
    key_row: KeyRow,
    name_row: NameRow = None,
    one_indexed: bool = False,
) -> VMatrix:
    """Add an atom to a V-Matrix.

    :param vma: V-Matrix
    :param symb: Symbol of atom to add
    :param key_row: Row of keys to define new atom added to key matrix
    :param name_row: Row of names to define new atom added to name matrix
    :param one_indexed: Parameter to store keys in one-indexing
    :return: V-Matrix with new atom
    """
    symbs = symbols(vma)
    symbs += (symb,)

    key_mat = key_matrix(vma, shift=(1 if one_indexed else 0))
    key_mat += (key_row,)

    name_mat = None if name_row is None else (*name_matrix(vma), name_row)

    vma = from_data(symbs, key_mat, name_mat, one_indexed=one_indexed)

    return vma


def remove_atom(vma: VMatrix, key: Key) -> VMatrix:
    """Remove an atom from a V-Matrix. Error raised if attempting
    to remove atom other atoms depend on.

    :param vma: V-Matrix
    :param key: key of atom to remove
    :return: V-Matrix without atom
    """
    symbs = list(symbols(vma))
    symbs.pop(key)

    key_mat = list(key_matrix(vma))
    key_mat.pop(key)
    key_mat = numpy.array(key_mat, dtype=object)

    if (key_mat == key).any():
        raise ValueError(f"Other atoms in z-matrix depend on atom {key}")

    key_map = numpy.vectorize(lambda x: x if (x is None or x < key) else x - 1)
    key_mat = key_map(key_mat)

    name_mat = list(name_matrix(vma))
    name_mat.pop(key)

    vma = from_data(symbs, key_mat, name_mat)

    return vma


# # # validation
def is_valid(vma: VMatrix) -> bool:
    """Assess if a V-Matrix has proper structure.

    :param vma: V-Matrix
    :return: True, if V-Matrix has proper structure; False otherwise
    """
    ret = True
    try:
        assert _is_sequence_of_triples(vma)
        symbs, key_mat, name_mat = zip(*vma, strict=True)
        from_data(symbs, key_mat, name_mat)
    except AssertionError:
        ret = False

    return ret


def is_standard_form(vma: VMatrix) -> bool:
    """Assesses if the names of the V-Matrix are in standard form:
        RN: (1<=N<=Ncoords)
        AN: (2<=N<=Ncoords)
        DN: (1<=N<=Ncoords).

    :param vma: V-Matrix
    :return: True if V-Matrix is in standard form, False if not
    """
    return names(vma) == names(standard_form(vma))


# # # I/O
def string(vma: VMatrix, one_indexed: bool = False) -> str:
    """Write a V-Matrix object to a string.

    :param vma: V-Matrix
    :param one_indexed: Parameter to write keys in one-indexing
    :return: String
    """
    shift = 1 if one_indexed else 0
    symbs = symbols(vma)
    key_mat = key_matrix(vma, shift=shift)
    name_mat = name_matrix(vma)

    def _line_string(row_idx):
        line_str = f"{symbs[row_idx]:<2s} "
        keys = key_mat[row_idx]
        names_ = name_mat[row_idx]
        line_str += " ".join(
            [
                f"{keys[col_idx]:>d} {names_[col_idx]:>5s} "
                for col_idx in range(min(row_idx, 3))
            ]
        )
        return line_str

    natms = len(symbs)
    vma_str = "\n".join([_line_string(row_idx) for row_idx in range(natms)])

    return vma_str


def from_string(vma_str: str, one_indexed: bool | None = None) -> VMatrix:
    """Parse a V-Matrix object from a string.

    :param vma_str: string containing a V-Matrix
    :param one_indexed: Read a one-indexed string?
    :return: V-Matrix of string
    """
    rows = VMAT_LINES.parseString(vma_str).asList()
    symbs = [r.pop(0) for r in rows]
    key_mat = [r[::2] for r in rows]
    name_mat = [r[1::2] for r in rows]

    vma = from_data(symbs, key_mat, name_mat, one_indexed=one_indexed)

    return vma


# # helpers
def _key_matrix(
    key_mat: KeyMatrix, natms: int, one_indexed: bool | None = None
) -> KeyMatrix:
    """Build name matrix of the V-Matrix that contains the
    coordinate keys by row and column.

    :param key_mat: key matrix of V-Matrix coordinate keys
    :param natms: number of atoms
    :param one_indexed: parameter to write keys in one-indexing
    :return: Name matrix of V-Matrix
    """
    if natms == 1:
        return ((None, None, None),)

    # Check dimensions and ensure proper formatting
    key_mat = [list(row) + [None] * (3 - len(row)) for row in key_mat]
    key_mat = numpy.array(key_mat, dtype=object)

    assert key_mat.ndim == 2 and key_mat.shape == (natms, 3)
    triu_idxs = numpy.triu_indices(natms, m=3)

    one_indexed = bool(min(key_mat[1:, 0])) if one_indexed is None else one_indexed
    key_mat[1:, 0] -= 1 if one_indexed else 0
    key_mat[2:, 1] -= 1 if one_indexed else 0
    key_mat[3:, 2] -= 1 if one_indexed else 0

    key_mat[triu_idxs] = None

    return tuple(map(tuple, key_mat))


def _name_matrix(name_mat: NameMatrix, natms: int) -> NameMatrix:
    """Build name matrix of the V-Matrix that contains the
    coordinate names by row and column.

    :param name_mat: key matrix of V-Matrix coordinate keys
    :param natms: number of atoms
    :return: Name matrix of Vmatrix
    """
    if name_mat is None:
        name_mat = numpy.empty((natms, 3), dtype=object)
        for row in range(0, natms):
            if row > 0:
                name_mat[row, 0] = f"R{row:d}"
            if row > 1:
                name_mat[row, 1] = f"A{row:d}"
            if row > 2:
                name_mat[row, 2] = f"D{row:d}"

    # Check dimensions and make sure there are Nones in the right places
    name_mat = [list(row) + [None] * (3 - len(row)) for row in name_mat]
    name_mat = numpy.array(name_mat, dtype=object)

    assert name_mat.ndim == 2 and name_mat.shape == (natms, 3)
    natms = name_mat.shape[0]
    triu_idxs = numpy.triu_indices(natms, m=3)
    tril_idxs = numpy.tril_indices(natms, -1, m=3)

    assert all(isinstance(name, str) for name in name_mat[tril_idxs])
    name_mat[triu_idxs] = None

    return tuple(map(tuple, name_mat))


def _is_sequence_of_triples(obj) -> bool:
    """Assess if input object sequence has length of three.

    :param obj: Object with __len__ attribute
    :return: True if object is less than len=3, False if not
    """
    ret = hasattr(obj, "__len__")
    if ret:
        ret = all(hasattr(item, "__len__") and len(item) == 3 for item in obj)

    return ret
