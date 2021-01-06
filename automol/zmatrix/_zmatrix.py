""" z-matrix
"""

import itertools
import numpy
from qcelemental import periodictable as pt
from qcelemental import constants as qcc
import autoread as ar
import autowrite as aw
import automol.graph
import automol.geom
import automol.create.zmatrix
import automol.convert.zmatrix
import automol.convert.geom
from automol import vmatrix as _v_
from automol import util


# constructors
def from_geometry(vma, geo):
    """ determine z-matrix from v-matrix and geometry
    """
    assert _v_.symbols(vma) == automol.geom.symbols(geo)
    val_dct = {}
    coo_dct = _v_.coordinates(vma, multi=False)
    dist_names = _v_.distance_names(vma)
    cent_names = _v_.central_angle_names(vma)
    dih_names = _v_.dihedral_angle_names(vma)
    for name, coo in coo_dct.items():
        if name in dist_names:
            val_dct[name] = automol.geom.distance(geo, *coo)
        elif name in cent_names:
            val_dct[name] = automol.geom.central_angle(geo, *coo)
        elif name in dih_names:
            val_dct[name] = automol.geom.dihedral_angle(geo, *coo)

    zma = automol.create.zmatrix.from_data(
        symbols=_v_.symbols(vma), key_matrix=_v_.key_matrix(vma),
        name_matrix=_v_.name_matrix(vma), values=val_dct)
    return zma


# getters
def var_(zma):
    """ the variable matrix (atom symbols, atom keys, and coordinate names)
    """
    vma, _ = zma
    return vma


def count(zma):
    """ the number of z-matrix rows (number of atoms or dummy atoms)
    """
    return _v_.count(var_(zma))


def atom_count(zma, sym, match=True):
    """ the number of entries for desired atom types (put in vma code)
    """
    return len(atom_count(zma, sym, match=match))


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
            val *= qcc.conversion_factor('bohr', 'angstrom')
        if degree and name in ang_names:
            val *= qcc.conversion_factor('radian', 'degree')
        val_dct[name] = val
    return val_dct


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


def coordinates(zma, shift=0, multi=True):
    """ coordinate keys associated with each coordinate name, as a dictionary

    (the values are sequences of coordinate keys, since there may be multiple)
    """
    return _v_.coordinates(var_(zma), shift=shift, multi=multi)


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


def new_distance_name(zma):
    """ a distance coordinate name that hasn't been used

    (if standard-form, returns the next standard-form name
    """
    names_ = names(zma)
    dist_name_iter = iter(map('R{:1}'.format, itertools.count(1)))
    dist_name = next(filter(lambda x: x not in names_, dist_name_iter))
    return dist_name


def new_central_angle_name(zma):
    """ a central angle coordinate name that hasn't been used

    (if standard-form, returns the next standard-form name
    """
    names_ = central_angle_names(zma)
    cang_name_iter = iter(map('A{:1}'.format, itertools.count(2)))
    cang_name = next(filter(lambda x: x not in names_, cang_name_iter))
    return cang_name


def new_dihedral_angle_name(zma):
    """ a dihedral angle coordinate name that hasn't been used

    (if standard-form, returns the next standard-form name
    """
    names_ = dihedral_angle_names(zma)
    dih_name_iter = iter(map('D{:1}'.format, itertools.count(3)))
    dih_name = next(filter(lambda x: x not in names_, dih_name_iter))
    return dih_name


def angle_names(zma):
    """ angle coordinate names (dihedral and central)
    """
    return _v_.angle_names(var_(zma))


def dummy_coordinate_names(zma):
    """ names of dummy atom coordinates
    """
    return _v_.dummy_coordinate_names(var_(zma))


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


def atom_specifiers(zma, without_dummies=False):
    """ For each atom, this gives the indices of atoms specifying its position
    """
    syms = symbols(zma)
    key_mat = key_matrix(zma)

    def _specifiers(idx):
        key_mat_row = key_mat[idx]
        idxs = tuple(filter(lambda x: x is not None and syms, key_mat_row))
        if without_dummies:
            idxs = tuple(idx for idx in idxs if pt.to_Z(syms[idx]))
        return idxs

    idxs = range(count(zma))
    atm_dep_dct = dict(zip(idxs, map(_specifiers, idxs)))
    return atm_dep_dct


def atom_dependents(zma, without_dummies=False):
    """ For each atom, this gives the indices of atoms depending on its position
    """
    syms = symbols(zma)
    key_mat = key_matrix(zma)

    def _dependents(idx):
        idxs = [idx_ for idx_, row in enumerate(key_mat) if idx in row]
        if without_dummies:
            idxs = tuple(idx for idx in idxs if pt.to_Z(syms[idx]))
        return tuple(idxs)

    idxs = range(count(zma))
    atm_dep_dct = dict(zip(idxs, map(_dependents, idxs)))
    return atm_dep_dct


def dummy_atom_anchors(zma):
    """ Three atoms specifying each dummy atom's position

    If the dummy atom is below the third row of the z-matrix, the anchoring
    atoms will be the three atoms used to specify its position in `zma`
    (see :meth:`autofile.zmatrix.atom_specifiers`).
    If the dummy atom is in the first three rows of the z-matrix, atoms from
    the rest of `zma` will be included to arrive at a total of three.
    Where possible, the third atom will be chosen to be non-collinear with the
    first two.

    (Could be generalized to non-dummy atoms, but I can't see a use case.)
    """
    geo = automol.convert.zmatrix.geometry(zma)
    xyzs = automol.geom.coordinates(geo)

    atm_spe_dct = atom_specifiers(zma)
    atm_dep_dct = atom_dependents(zma, without_dummies=True)

    dummy_idxs = atom_indices(zma, 'X', match=True)

    atm_ach_dct = {}
    for dummy_idx in dummy_idxs:
        pool = atm_spe_dct[dummy_idx] + atm_dep_dct[dummy_idx]

        # find a non-collinear atom for the third index
        atm_achs = list(pool[:2])
        for atm3_idx in pool[2:]:
            atm1_xyz, atm2_xyz = map(xyzs.__getitem__, atm_achs)
            atm3_xyz = xyzs[atm3_idx]
            if not util.vec.are_parallel(atm1_xyz, atm2_xyz,
                                         orig_xyz=atm3_xyz):
                atm_achs.append(atm3_idx)
                break

        atm_ach_dct[dummy_idx] = tuple(atm_achs)

    return atm_ach_dct


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


# validation
def is_valid(zma):
    """ is this a valid zmatrix?
    """
    ret = hasattr(zma, '__len__') and len(zma) == 2
    if ret:
        vma, val_dct = zma
        ret = _v_.is_valid(vma) and set(_v_.names(vma)) == set(val_dct)
        if ret:
            try:
                automol.create.zmatrix.from_data(
                    _v_.symbols(vma), _v_.key_matrix(vma),
                    _v_.name_matrix(vma), val_dct)
            except AssertionError:
                ret = False
    return ret


# setters
def set_keys(zma, key_dct):
    """ set keys for the z-matrix
    """
    natoms = count(zma)
    for key, idxs in key_dct.items():
        assert key <= (natoms - 1)
        assert all(idx <= (natoms - 1) for idx in idxs)
        assert all(idx != key for idx in idxs)
        assert len(idxs) == len(set(idxs))

    orig_key_matrix = key_matrix(zma)
    new_key_matrix = tuple()
    for i, row in enumerate(orig_key_matrix):
        if i in key_dct.keys():
            new_key_matrix += (key_dct[i],)
        else:
            new_key_matrix += (row,)

    return automol.create.zmatrix.from_data(
        symbols(zma), new_key_matrix,
        name_matrix(zma), values(zma))


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

    return automol.create.zmatrix.from_data(
        _v_.symbols(vma), _v_.key_matrix(vma), name_mat, val_dct)


def set_values(zma, val_dct):
    """ set coordinate values for the z-matrix
    """
    vma = var_(zma)
    _names = _v_.names(vma)
    assert set(val_dct.keys()) <= set(_names)

    new_val_dct = values(zma).copy()
    new_val_dct.update(val_dct)
    return automol.create.zmatrix.from_data(
        _v_.symbols(vma), _v_.key_matrix(vma), _v_.name_matrix(vma),
        new_val_dct)


def shift_row_to_end(zma, row_idx, frm_bnd_key, brk_bnd_key):
    """ move a single row of the zmatrix to the end

        set to only act if keys in row_idx are larger than row_idx
        i.e., if the shift is clearly needed
        function only works if atom row that needs to be moved is not
        used to define other atoms
    """

    # set the original
    orig_symbols = symbols(zma)
    orig_keys = key_matrix(zma)
    orig_names = name_matrix(zma)
    orig_values = values(zma)
    h_idx = row_idx

    # Check if shift will break things
    will_break = False
    for i, key_row in enumerate(orig_keys):
        if i <= row_idx:
            continue
        for key in key_row:
            if key == row_idx:
                will_break = True

    # Check if a shift is needed
    shift_needed = False
    for key in orig_keys[row_idx]:
        if key > row_idx:
            shift_needed = True

    # Build new tuples with the row corresponding to the row_idx removed
    if shift_needed and not will_break:
        new_symbols = tuple()
        new_key_matrix = tuple()
        new_name_matrix = tuple()
        zma_data = zip(orig_symbols, orig_keys, orig_names)
        for i, (sym, key_row, name_row) in enumerate(zma_data):

            if i >= row_idx:
                # Decrement the keys if they are greater than row idx
                new_key_row = []
                for key in key_row:
                    if key >= row_idx:
                        new_key = key - 1
                    else:
                        new_key = key
                    new_key_row.append(new_key)
                new_key_row = tuple(new_key_row)

                # Set key row for row_idx row, else add to overal lsts
                if i == row_idx:
                    new_keys_at_row = (new_key_row,)
                else:
                    new_symbols += (sym,)
                    new_key_matrix += (new_key_row,)
                    new_name_matrix += (name_row,)
            else:
                # Simply add the zmatrix pieces to new parts
                new_symbols += (sym,)
                new_key_matrix += (key_row,)
                new_name_matrix += (name_row,)

        # Append the row corresponding to row_idx
        new_symbols += ((orig_symbols[row_idx]),)
        new_key_matrix += (new_keys_at_row)
        new_name_matrix += ((orig_names[row_idx]),)

        zma = automol.create.zmatrix.from_data(
            new_symbols, new_key_matrix, new_name_matrix, orig_values)
        h_idx = automol.zmatrix.count(zma) - 1

        frm_bnd_key = list(frm_bnd_key)
        brk_bnd_key = list(brk_bnd_key)

        if shift_needed:
            if frm_bnd_key[0] > row_idx:
                frm_bnd_key[0] = frm_bnd_key[0] - 1
            elif frm_bnd_key[0] == row_idx:
                frm_bnd_key[0] = h_idx
            if frm_bnd_key[1] > row_idx:
                frm_bnd_key[1] = frm_bnd_key[1] - 1
            elif frm_bnd_key[1] == row_idx:
                frm_bnd_key[1] = h_idx
            if brk_bnd_key[0] > row_idx:
                brk_bnd_key[0] = brk_bnd_key[0] - 1
            elif brk_bnd_key[0] == row_idx:
                brk_bnd_key[0] = h_idx
            if brk_bnd_key[1] > row_idx:
                brk_bnd_key[1] = brk_bnd_key[1] - 1
            elif brk_bnd_key[1] == row_idx:
                brk_bnd_key[1] = h_idx

        frm_bnd_key = frozenset(frm_bnd_key)
        brk_bnd_key = frozenset(brk_bnd_key)

    return zma, h_idx, frm_bnd_key, brk_bnd_key


def standard_names(zma, shift=0):
    """ standard names for the coordinates, by their current names

    (follows x2z format)
    """
    return _v_.standard_names(var_(zma), shift=shift)


def standard_form(zma, shift=0):
    """ set standard variable names for the z-matrix

    (follows x2z format)
    """
    return set_names(zma, standard_names(zma, shift=shift))


# operations
def append(zma, sym, key_row, name_row, val_row):
    """ append an atom to the end of a z-matrix
    """
    syms = list(symbols(zma))
    key_mat = list(key_matrix(zma))
    nam_mat = list(name_matrix(zma))
    val_dct = values(zma)

    syms.append(sym)
    key_mat.append(key_row)
    nam_mat.append(name_row)
    val_dct.update(dict(zip(name_row, val_row)))

    if None in val_dct:
        val_dct.pop(None)

    return automol.create.zmatrix.from_data(syms, key_mat, nam_mat, val_dct)


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
    assert numpy.all(numpy.equal(key_mat2[:3][join_idxs], None))
    assert numpy.all(numpy.equal(name_mat2[:3][join_idxs], None))
    key_mat2[:3][join_idxs] = join_key_mat[join_idxs]
    name_mat2[:3][join_idxs] = join_name_mat[join_idxs]

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

    return automol.create.zmatrix.from_data(syms, key_mat, name_mat, val_dct)


def insert_dummy_atom(zma, x_key, x_key_mat, x_name_mat, x_val_dct):
    """ insert a dummy atom at a given position in the z-matrix
    """
    syms = symbols(zma)
    key_mat = numpy.array(key_matrix(zma))
    name_mat = numpy.array(name_matrix(zma))

    # check whether x_name_mat overlaps with name_mat

    natms = len(syms)
    keys = numpy.arange(natms)
    keys[x_key:] += 1
    key_dct = dict(enumerate(keys))
    key_dct[None] = None
    key_mat = [list(map(key_dct.__getitem__, row)) for row in key_mat]

    syms = list(syms)
    key_mat = list(key_mat)
    name_mat = list(name_mat)
    syms.insert(x_key, 'X')
    key_mat.insert(x_key, [None, None, None])
    name_mat.insert(x_key, [None, None, None])

    key_mat = numpy.array(key_mat)
    name_mat = numpy.array(name_mat)

    x_key_mat = numpy.array(x_key_mat, dtype=numpy.object_)
    x_name_mat = numpy.array(x_name_mat, dtype=numpy.object_)

    x_idxs = numpy.not_equal(x_key_mat, None)

    offset = min(3, len(key_mat)-x_key)
    key_mat[x_key:x_key+offset][x_idxs] = x_key_mat[:offset][x_idxs]
    name_mat[x_key:x_key+offset][x_idxs] = x_name_mat[:offset][x_idxs]

    val_dct = values(zma)
    assert not set(val_dct.keys()) & set(x_val_dct.keys())
    val_dct.update(x_val_dct)

    return automol.create.zmatrix.from_data(syms, key_mat, name_mat, val_dct)


def convert(zma1, zma2, frm_bnd_keys=(), brk_bnd_keys=()):
    """ Convert the geometry of zma1 into the z-matrix form of zma2

    Each dummy atom is placed based on its position in `zma2`, relative to
    three "anchoring atoms" (see :meth:`autofile.zmatrix.dummy_atom_anchors`).

    :param zma1: A z-matrix specifying the desired geometry
    :param zma2: A z-matrix giving the desired z-matrix form. If no
        formed/broken bond keys are provided, its connectivity graph must be
        isomorphic to the connectivity graph of `zma1`.
    :param frm_bnd_keys: Bonds that must be formed in order to make `zma1`
        isomorphic to `zma2`.
    :type frm_bnd_keys: list, tuple, set
    :param frm_bnd_keys: Bonds that must be broken in order to make `zma1`
        isomorphic to `zma2`.
    :type brk_bnd_keys: list, tuple, set
    """

    # Remove dummy atoms and reorder the atoms in geo1 to match geo2
    geo1 = automol.convert.zmatrix.geometry(zma1)
    geo2 = automol.convert.zmatrix.geometry(zma2)

    geo1 = automol.geom.without_dummy_atoms(geo1)
    geo2_no_dummies = automol.geom.without_dummy_atoms(geo2)

    gra1 = automol.convert.geom.graph(geo1)
    gra2 = automol.convert.geom.graph(geo2_no_dummies)

    # If bonds are broken and formed, apply the transformation so that there is
    # an isomorphism between gra1 and gra2
    if frm_bnd_keys or brk_bnd_keys:
        tra = automol.graph.trans.from_data(frm_bnd_keys, brk_bnd_keys)
        gra1 = automol.graph.trans.apply(tra, gra1)

    iso_dct = automol.graph.full_isomorphism(gra1, gra2)

    assert iso_dct is not None

    geo1 = automol.geom.reorder(geo1, iso_dct)

    # Insert the dummy atoms according to their relative positions in geo2
    atm_ach_dct2 = dummy_atom_anchors(zma2)

    # First, insert the dummies without coordinates to fully line up the atom
    # indices in geo1 and geo2
    for dummy_idx in sorted(atm_ach_dct2):
        geo1 = automol.geom.insert(geo1, 'X', [0., 0., 0.], idx=dummy_idx)

    for dummy_idx, (idx1, idx2, idx3) in sorted(atm_ach_dct2.items()):
        dist = automol.geom.distance(geo2, dummy_idx, idx1)
        ang = automol.geom.central_angle(geo2, dummy_idx, idx1, idx2)
        dih = automol.geom.dihedral_angle(geo2, dummy_idx, idx1, idx2, idx3)

        xyzs = automol.geom.coordinates(geo1)
        xyz1, xyz2, xyz3 = map(xyzs.__getitem__, (idx1, idx2, idx3))

        dummy_xyz = util.vec.from_internals(dist=dist, xyz1=xyz1,
                                            ang=ang, xyz2=xyz2,
                                            dih=dih, xyz3=xyz3)
        geo1 = automol.geom.set_coordinates(geo1, {dummy_idx: dummy_xyz})

    # At this point, the two geometries should be fully lined up and we can
    # simply call the constructor
    vma2 = var_(zma2)
    zma = from_geometry(vma2, geo1)
    return zma


# misc
def is_standard_form(zma):
    """ set standard variable names for the z-matrix

    (follows x2z format)
    """
    return is_valid(zma) and _v_.is_standard_form(var_(zma))


# I/O
def from_string(zma_str, one_indexed=True, angstrom=True, degree=True):
    """ read a z-matrix from a string
    """
    syms, key_mat, name_mat, val_dct = ar.zmatrix.read(zma_str)

    zma = automol.create.zmatrix.from_data(
        syms, key_mat, name_mat, val_dct, one_indexed=one_indexed,
        angstrom=angstrom, degree=degree)
    return zma


def string(zma, one_indexed=True, angstrom=True, degree=True):
    """ write a z-matrix to a string
    """
    shift = 1 if one_indexed else 0
    zma_str = aw.zmatrix.write(
        syms=_v_.symbols(var_(zma)),
        key_mat=_v_.key_matrix(var_(zma), shift=shift),
        name_mat=_v_.name_matrix(var_(zma)),
        val_dct=values(zma, angstrom=angstrom, degree=degree)
    )
    return zma_str


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
    if var_(zma1) == var_(zma2):
        # first compare the distances
        val_dct1 = values(zma1)
        val_dct2 = values(zma2)
        dist_names = distance_names(zma1)
        dist_vals1 = tuple(map(val_dct1.__getitem__, dist_names))
        dist_vals2 = tuple(map(val_dct2.__getitem__, dist_names))
        if numpy.allclose(dist_vals1, dist_vals2, rtol=dist_rtol):
            if just_dist:
                ret = True
            else:
                # now compare the angles
                # see https://gamedev.stackexchange.com/a/4472
                ang_names = angle_names(zma1)

                ang_vals1 = tuple(map(val_dct1.__getitem__, ang_names))
                ang_vals2 = tuple(map(val_dct2.__getitem__, ang_names))

                ang_vals1 = numpy.mod(ang_vals1, 2*numpy.pi)
                ang_vals2 = numpy.mod(ang_vals2, 2*numpy.pi)

                ang_diffs = numpy.abs(ang_vals1 - ang_vals2)
                ang_diffs = numpy.pi - numpy.abs(ang_diffs - numpy.pi)

                if numpy.allclose(ang_diffs, 0., atol=ang_atol):
                    ret = True
    return ret


# random sampling
def samples(zma, nsamp, range_dct):
    """ randomly sample over torsional dihedrals
    """
    _names = tuple(range_dct.keys())
    ranges = tuple(range_dct.values())
    vals_lst = _sample_over_ranges(ranges, nsamp)

    zmas = tuple(set_values(zma, dict(zip(_names, vals)))
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


# z-matrix torsional degrees of freedom
def torsional_symmetry_numbers(zma, tors_names,
                               frm_bnd_key=None, brk_bnd_key=None):
    """ symmetry numbers for torsional dihedrals
    """
    dih_edg_key_dct = _dihedral_edge_keys(zma)
    assert set(tors_names) <= set(dih_edg_key_dct.keys())
    edg_keys = tuple(map(dih_edg_key_dct.__getitem__, tors_names))

    gra = automol.convert.zmatrix.graph(zma, remove_stereo=True)
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
