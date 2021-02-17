""" z-matrix
"""

import numpy
from automol import vmat
import automol.create.zmat
import automol.convert.zmat
import automol.geom
import autoread as ar
import autowrite as aw
from phydat import phycon


# constructors
def from_geometry(vma, geo):
    """  Build a Z-Matrix from a V-Matrix and a molecular geometry.

        :param vma: V-Matrix
        :type vma: automol V-Matrix data structure
        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :rtype: automol Z-Matrix data structure
    """
    assert symbols(vma) == symbols(geo)

    syms = symbols(vma)
    key_mat = key_matrix(vma)
    name_mat = name_matrix(vma)
    val_mat = numpy.empty(numpy.shape(key_mat), dtype=numpy.object_)

    for row, key_row in enumerate(key_mat):
        if row > 0:
            val_mat[row, 0] = automol.geom.distance(geo, row,
                                                    *key_row[:1])
        if row > 1:
            val_mat[row, 1] = automol.geom.central_angle(geo, row,
                                                         *key_row[:2])
        if row > 2:
            val_mat[row, 2] = automol.geom.dihedral_angle(geo, row,
                                                          *key_row[:3])

    zma = automol.create.zmat.from_data(syms, key_mat, val_mat, name_mat)
    return zma


# converters
def geometry(zma, dummy=False):
    """ Convert a Z-Matrix to a molecular geometry.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: automol molecular geometry data structure
    """

    if dummy:
        geo = automol.convert.zmat.geometry_with_dummy_atoms(zma)
    else:
        geo, _ = automol.convert.zmat.geometry(zma)

    return geo


# getters
def count(zma):
    """ Obtain the number of rows of the Z-Matrix, which corresponds to
        the number of atoms defined in the Z-Matrix. This includes all
        real and dummy atoms.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: int
    """
    return vmat.count(zma)


def symbols(zma):
    """ Obtain the atomic symbols for all atoms defined in the Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: tuple(str)
    """
    return vmat.symbols(zma)


def atom_indices(zma, symb, match=True):
    """ Obtain the indices of a atoms of a particular type in the geometry.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param match: grab idxs that match given atom type
        :param symb: atomic symbol
        :type symb: str
        :param match: obtain indices of symbols that match the type?
        :type match: bool
    """
    return vmat.atom_indices(zma, symb, match=match)




def key_matrix(zma, shift=0):
    """ Obtain the key matrix of the Z-Matrix that contains the
        coordinate atom keys by row and column.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param shift: value to shift the keys by when obtaining the key matrix
        :type shift: int
        :rtype: tuple(tuple(int))
    """
    return vmat.key_matrix(zma, shift=shift)


def name_matrix(zma):
    """ Obtain the name matrix of the Z-Matrix that contains the
        coordinate names by row and column.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: tuple(tuple(str))
    """
    return vmat.name_matrix(zma)


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
        val_mat = numpy.array(val_mat, dtype=numpy.object_)

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
        :rtype: dict[str: float]
    """

    names = numpy.ravel(name_matrix(zma))
    vals = numpy.ravel(value_matrix(zma, angstrom=angstrom, degree=degree))
    val_dct = dict(zip(names, vals))
    val_dct.pop(None)

    return val_dct


def dummy_keys(zma):
    """ Obtain keys to dummy atoms in the Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: tuple[int]
    """
    keys = tuple(key for key, sym in enumerate(symbols(zma)) if sym == 'X')
    return keys


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
    ngb_keys = [key_mat[k][0] for k in dum_keys]
    key_dct = dict(zip(dum_keys, ngb_keys))
    return key_dct


def distance(zma, key1, key2, angstrom=False):
    """ Measure the distance between two atoms defined in a Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param key1: key of atom 1 in the pair to be measured
        :type key1: int
        :param key2: key of atom 2 in the pair to be measured
        :type key2: int
        :param angstrom: parameter to control Bohr->Angstrom conversion
        :type angstrom: bool
    """
    geo, _ = automol.convert.zmat.geometry(zma)
    return automol.geom.distance(geo, key1, key2, angstrom=angstrom)


def central_angle(zma, key1, key2, key3, degree=False):
    """ Measure the angle inscribed by three atoms in a Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param key1: key of atom 1 in the triplet to be measured
        :type key1: int
        :param key2: key of atom 2 in the triplet to be measured
        :type key2: int
        :param key3: key of atom 3 in the triplet to be measured
        :type key3: int
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
    """
    geo = automol.convert.zmat.geometry(zma)
    return automol.geom.central_angle(geo, key1, key2, key3, degree=degree)


def dihedral_angle(zma, key1, key2, key3, key4, degree=False):
    """ Measure the angle inscribed by three atoms in a molecular geometry.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param key1: key of atom 1 in the quartet to be measured
        :type key1: int
        :param key2: key of atom 2 in the quartet to be measured
        :type key2: int
        :param key3: key of atom 3 in the quartet to be measured
        :type key3: int
        :param key4: key of atom 4 in the quartet to be measured
        :type key4: int
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
    """
    geo = automol.convert.zmat.geometry(zma)
    return automol.geom.dihedral_angle(geo, key1, key2, key3, key4,
                                       degree=degree)


# setters
def set_key_matrix(zma, key_mat):
    """ Re-set the key matrix of a Z-Matrix using the input key matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param key_mat: matrix of Z-Matrix coordinate keys
        :type key_mat: tuple(tuple(int))
        :rtype: automol Z-Matrix data structure
    """

    syms = symbols(zma)
    val_mat = value_matrix(zma)
    name_mat = name_matrix(zma)
    zma = automol.create.zmat.from_data(syms, key_mat, val_mat, name_mat)

    return zma


def set_name_matrix(zma, name_mat):
    """ Re-set the name matrix of a Z-Matrix using the input name matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param name_mat: matrix of Z-Matrix coordinate names
        :type name_mat: tuple(tuple(int))
        :rtype: automol Z-Matrix data structure
    """

    syms = symbols(zma)
    val_mat = value_matrix(zma)
    key_mat = key_matrix(zma)
    zma = automol.create.zmat.from_data(syms, key_mat, val_mat, name_mat)

    return zma


def set_value_matrix(zma, val_mat):
    """ Re-set the name matrix of a Z-Matrix using the input value matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param val_mat: matrix of Z-Matrix coordinate values
        :type val_mat: tuple(tuple(int))
        :rtype: automol Z-Matrix data structure
    """

    syms = symbols(zma)
    key_mat = key_matrix(zma)
    name_mat = name_matrix(zma)
    zma = automol.create.zmat.from_data(syms, key_mat, val_mat, name_mat)

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

    val_mat = numpy.array(value_matrix(zma), dtype=numpy.object_)
    name_mat = numpy.array(name_matrix(zma), dtype=numpy.object_)

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


def standard_name_matrix(zma, shift=0):
    """ Builds a name matrix of the Z-Matrix where all of the
        coordinate names have been standardized:
            RN: (1<=N<=Ncoords)
            AN: (2<=N<=Ncoords)
            DN: (1<=N<=Ncoords)

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param shift: value to shift the keys by when obtaining the keys
        :type shift: int
        :rtype: tuple(tuple(str))
    """
    return vmat.standard_name_matrix(zma, shift=shift)


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


def rename(zma, name_dct):
    """ Rename a subset of the coordinates of a Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param name_dct: mapping from old coordinate names to new ones
        :type name_dct: dict[str: str]
        :rtype: automol Z-Matrix data strucutre
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
    """ Remove an atom from a Z-Matrix. Error raised if attempting
        to remove atom other atoms depend on.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :param key: key of atom to remove
        :type key: str
        :rtype: automol Z-Matrix data structure
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
    """ Join two z-matrices together replacing the first atom in zma2
        by an atom in zma1.

        :param zma1: Z-Matrix 1
        :type zma: automol Z-Matrix data structure
        :param zma2: Z-Matrix 2
        :type zma: automol Z-Matrix data structure
        :param rep1_key: the key of the atom in zma1 that will replace
            the first atom in zma2
        :type rep1_key: int
        :param key_mat: matrix of Z-Matrix coordinate keys
        :type key_mat: tuple(tuple(int))
        :param val_mat: matrix of Z-Matrix coordinate values
        :type val_mat: tuple(tuple(int))
        :param name_mat: matrix of Z-Matrix coordinate names
        :type name_mat: tuple(tuple(int))
        :param degree: parameter to control radian->degree conversion
        :type degree: bool
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

    zma = automol.create.zmat.from_data(
        syms, key_mat, val_mat, name_mat, degree=degree)

    return zma


# conversions
def vmatrix(zma):
    """ Parse and return the V-Matrix component of a Z-Matrix.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
    """
    return vmat.from_data(symbs=symbols(zma),
                          key_mat=key_matrix(zma),
                          name_mat=name_matrix(zma))


# comparisons
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


# I/O
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

    syms, key_mat, name_mat, val_mat = ar.zmat.read(zma_str)
    zma = automol.create.zmat.from_data(
        syms, key_mat, val_mat, name_mat, one_indexed=one_indexed,
        angstrom=angstrom, degree=degree)

    return zma


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

    shift = 1 if one_indexed else 0
    zma_str = aw.zmat.write(
        symbs=symbols(zma),
        key_mat=key_matrix(zma, shift=shift),
        name_mat=name_matrix(zma),
        val_dct=value_dictionary(zma, angstrom=angstrom, degree=degree)
    )

    return zma_str


# validators
def is_valid(zma):
    """ Assess if a Z-Matrix has proper data structure.

        :param zma: Z-Matrix
        :type zma: automol Z-Matrix data structure
        :rtype: bool
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
    """ Assess if input object sequence has length of four.

        :param obj: object with __len__ attribute
        :type obj: list, tuple, dict
        :rtype: bool
    """

    ret = hasattr(obj, '__len__')
    if ret:
        ret = all(hasattr(item, '__len__') and len(item) == 4 for item in obj)

    return ret


def samples(zma, nsamp, range_dct):
    """ randomly sample over torsional dihedrals
    """
    _names = tuple(range_dct.keys())
    ranges = tuple(range_dct.values())
    vals_lst = _sample_over_ranges(ranges, nsamp)

    zmas = tuple(set_values_by_name(zma, dict(zip(_names, vals)))
                 for vals in vals_lst)
    return zmas


def _sample_over_ranges(rngs, nsamp):
    """ randomly sample over several ranges
    """
    nrng = len(rngs)
    samp_mat = numpy.random.rand(nsamp, nrng)
    for i, (start, stop) in enumerate(rngs):
        samp_mat[:, i] = samp_mat[:, i] * (stop - start) + start
    return tuple(map(tuple, samp_mat))


# Functions that need to be here
def coord_idxs(zma, key):
    """give a bond length key, return the indices of involved bonded atoms
    """
    coords = coordinates(zma)
    idxs = coords.get(key, [None])
    return idxs[0]


def bond_key_from_idxs(zma, idxs):
    """given indices of involved bonded atoms, return bond name
    """
    idxs = list(idxs)
    idxs.sort(reverse=True)
    idxs = tuple(idxs)
    bond_key = None
    coords = coordinates(zma)
    for key in coords:
        for coord in coords.get(key, [None]):
            if idxs == coord:
                bond_key = key
    return bond_key


def get_babs1(zma, dist_name):
    """ get name of torsional coordinate associated with babs1 pre-reformatting
    """
    idxs = coord_idxs(zma, dist_name)
    idx = max(idxs)
    babs1 = 'D{:g}'.format(idx)
    return babs1


def get_babs2(zma, dist_name):
    """ get name of torsional coordinate associated with babs2 pre-reformatting
    """
    idxs = coord_idxs(zma, dist_name)
    idx = max(idxs)
    babs2 = 'D{:g}'.format(idx+1)
    return babs2


def atom_indices(zma, sym, match=True):
    """ indices for a particular atom type
        :param match: grab idxs that match given atom type
    """

    syms = symbols(zma)
    idxs = tuple()
    for idx, sym_ in enumerate(syms):
        if sym_ == sym and match:
            idxs += (idx,)
        elif sym_ != sym and not match:
            idxs += (idx,)

    return idxs


# wrappers to vmat
def dummy_coordinate_names(zma):
    """ names of dummy atom coordinates
    """
    return vmat.dummy_coordinate_names(vmatrix(zma))


def coordinates(zma, shift=0, multi=True):
    """ coordinate keys associated with each coordinate name,
        as a dictionary
    """
    return vmat.coordinates(vmatrix(zma), shift=shift, multi=multi)


def dihedral_angle_names(zma):
    """ dihedral angle coordinate names
    """
    return vmat.dihedral_angle_names(vmatrix(zma))


# z-matrix torsional degrees of freedom
def torsional_symmetry_numbers(zma, tors_names,
                               frm_bnd_key=None, brk_bnd_key=None):
    """ symmetry numbers for torsional dihedrals
    """
    dih_edg_key_dct = _dihedral_edge_keys(zma)
    assert set(tors_names) <= set(dih_edg_key_dct.keys())
    edg_keys = tuple(map(dih_edg_key_dct.__getitem__, tors_names))

    gra = automol.convert.zmat.graph(zma, stereo=False)
    bnd_sym_num_dct = automol.graph.bond_symmetry_numbers(
        gra, frm_bnd_key, brk_bnd_key)
    tors_sym_nums = []
    for edg_key in edg_keys:
        if edg_key in bnd_sym_num_dct.keys():
            sym_num = bnd_sym_num_dct[edg_key]
        else:
            sym_num = 1
        tors_sym_nums.append(sym_num)

    tors_sym_nums = tuple(tors_sym_nums)
    return tors_sym_nums


def _dihedral_edge_keys(zma):
    """ dihedral bonds, by name
    """
    coo_dct = coordinates(zma)
    dih_names = dihedral_angle_names(zma)
    dih_keys_lst = tuple(map(coo_dct.__getitem__, dih_names))
    dih_edg_key_dct = {dih_name: frozenset(dih_key[1:3])
                       for dih_name, dih_keys in zip(dih_names, dih_keys_lst)
                       for dih_key in dih_keys}
    return dih_edg_key_dct


def torsional_sampling_ranges(tors_names):
    """ sampling ranges for torsional dihedrals
    """
    # sym_nums = torsional_symmetry_numbers(zma, tors_names,
    # frm_bnd_key=None, brk_bnd_key=None)
    # return tuple((0, 2*numpy.pi/sym_num) for sym_num in sym_nums)
    # originally restricted range by sym_num.
    # But after all it appears that using the
    # full range is best.
    # after all it appears that using the full sampling range is most effective
    sym_nums = 1.
    return tuple((0, 2*numpy.pi/sym_nums) for tors_name in tors_names)


def torsional_scan_linspaces(zma, tors_names, increment=0.5,
                             frm_bnd_key=None, brk_bnd_key=None):
    """ scan grids for torsional dihedrals
    """
    sym_nums = torsional_symmetry_numbers(
        zma, tors_names, frm_bnd_key=frm_bnd_key, brk_bnd_key=brk_bnd_key)
    intervals = tuple(2*numpy.pi/sym_num - increment for sym_num in sym_nums)
    npoints_lst = tuple(
        (int(interval / increment)+1) for interval in intervals)
    return tuple((0, interval, npoints)
                 for interval, npoints in zip(intervals, npoints_lst))


def constraint_dct(zma, const_names, var_names=()):
    """ Build a dictionary of constraints

        has to round the values for the filesystem
    """

    # Get the list names sorted for dictionary
    rnames = (name for name in const_names if 'R' in name)
    anames = (name for name in const_names if 'A' in name)
    dnames = (name for name in const_names if 'D' in name)
    rnames = tuple(sorted(rnames, key=lambda x: int(x.split('R')[1])))
    anames = tuple(sorted(anames, key=lambda x: int(x.split('A')[1])))
    dnames = tuple(sorted(dnames, key=lambda x: int(x.split('D')[1])))
    constraint_names = rnames + anames + dnames

    # Remove the scan coordinates so they are not placed in the dict
    constraint_names = tuple(name for name in constraint_names
                             if name not in var_names)

    # Build dictionary
    if constraint_names:
        zma_vals = automol.zmat.value_dictionary(zma)
        zma_coords = automol.zmat.coordinates(zma)
        assert set(constraint_names) <= set(zma_coords.keys()), (
            'Attempting to constrain coordinates not in zma:\n{}\n{}'.format(
                constraint_names, zma_coords)
        )
        constraint_dct = dict(zip(
            constraint_names,
            (round(zma_vals[name], 2) for name in constraint_names)
        ))
    else:
        constraint_dct = None

    return constraint_dct
