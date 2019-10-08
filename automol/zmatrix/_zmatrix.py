""" z-matrix
"""
import itertools
import numpy
from qcelemental import constants as qcc
import autoread as ar
import autowrite as aw
import automol.graph
import automol.create.zmatrix
import automol.convert.zmatrix
import automol.convert.geom
from automol import vmatrix as _v_


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


def angle_names(zma):
    """ angle coordinate names (dihedral and central)
    """
    return _v_.angle_names(var_(zma))


def dummy_coordinate_names(zma):
    """ names of dummy atom coordinates
    """
    return _v_.dummy_coordinate_names(var_(zma))


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


def bond_idxs(zma, key):
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
    coords = coordinates(zma)
    for key in coords:
        for coord in coords.get(key, [None]):
            if idxs == coord:
                bond_key = key

    return bond_key


def get_babs1(zma, dist_name):
    """ get name of torsional coordinate associated with babs1 pre-reformatting
    """
    idxs = bond_idxs(zma, dist_name)
    idx = max(idxs)
    babs1 = 'D{:g}'.format(idx)
    return babs1


def get_babs2(zma, dist_name):
    """ get name of torsional coordinate associated with babs2 pre-reformatting
    """
    idxs = bond_idxs(zma, dist_name)
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


def shift_row_to_end(zma, row_idx):
    """ move a single row of the zmatrix to the end

        function only works if atom row that needs to be moved is not
        used to define other atoms
    """

    # set the original
    orig_symbols = symbols(zma)
    orig_keys = key_matrix(zma)
    orig_names = name_matrix(zma)
    orig_values = values(zma)

    # Build new tuples with the row corresponding to the row_idx removed
    new_symbols, new_key_matrix, new_name_matrix = tuple(), tuple(), tuple()
    zma_data = zip(orig_symbols, orig_keys, orig_names)
    for i, (sym, key_row, name_row) in enumerate(zma_data):
        if i == row_idx:
            continue
        else:
            new_symbols += (sym,)
            new_key_matrix += (key_row,)
            new_name_matrix += (name_row,)

    # Append the row corresponding to row_idx
    new_symbols += ((orig_symbols[row_idx]),)
    new_key_matrix += ((orig_keys[row_idx]),)
    new_name_matrix += ((orig_names[row_idx]),)

    return automol.create.zmatrix.from_data(
        new_symbols, new_key_matrix, new_name_matrix, orig_values)


# def _shift_method2(zma, row_idx, syms, key_mat, name_mat, value_dct):
#     """ attempt at generic shift method - not functioning
#
#         method needs some post processing
#         since atom to be moved is used to define number of rows
#     """
#
#     # Build new tuples with the row corresponding to the row_idx removed
#     new_symbols, new_key_matrix, new_name_matrix = tuple(), tuple(), tuple()
#     zmat_data = zip(syms, key_mat, name_mat)
#     bad_names = []
#     for i, (sym, key_row, name_row) in enumerate(zmat_data):
#         new_symbols += (sym,)
#         if i == 0:
#             new_key_matrix += ((None, None, None),)
#             new_name_matrix += ((None, None, None),)
#             bad_names.append([name_row[i] for i in (0, 1, 2)
#                               if name_row[i] is not None])
#         elif i == 1:
#             new_key_matrix += ((key_row[0], None, None),)
#             new_name_matrix += ((name_row[0], None, None),)
#             bad_names.append([name_row[i] for i in (1, 2)
#                               if name_row[i] is not None])
#         elif i == 2:
#             new_key_matrix += ((key_row[0], key_row[1], None),)
#             new_name_matrix += ((name_row[0], name_row[1], None),)
#             bad_names.append([name_row[i] for i in (2,)
#                               if name_row[i] is not None])
#         elif i == count(zma)-1:
#             new_key_matrix += ((0, 1, 2),)
#             new_name_matrix += (('R1000', 'A1000', 'D1000'),)
#             bad_names.append([name_row[i] for i in (0, 1, 2)
#                               if name_row[i] is not None])
#         else:
#             new_key_matrix += (key_row,)
#             new_name_matrix += (name_row,)
#
#     # lower the index of keys if needd
#     # for i, keys in enumerate(new_key_matrix):
#     #     if i in keys:
#     #         new_keys = ()
#     #         for key in keys:
#     #         new_key_matrix[i] =
#     #     else:
#     #         continue
#
#     # build new dict: delete bad, rename
#     flat_bad = [item for sublist in bad_names for item in sublist]
#     new_values = {}
#     for key, val in value_dct.items():
#         if key in flat_bad:
#             continue
#         else:
#             new_values[key] = val
#     new_values['R1000'] = 1.000
#     new_values['A1000'] = 1.000
#     new_values['D1000'] = 1.000
#
#     # determine the coordinates in the name matrix needed to recalculate
#     names_to_recalc = []
#     zma_info = zip(new_key_matrix, new_name_matrix)
#     for i, (key_row, name_row) in enumerate(zma_info):
#         # figure out what element of the keys corresponds to moving atom idx
#         if i != count(zma)-1:
#             if row_idx in key_row:
#                 idx_ref = key_row.index(row_idx)
#                 name = name_row[idx_ref]
#                 names_to_recalc.append([name, idx_ref])
#             else:
#                 names_to_recalc.append([None, None])
#
#     # form the initial new z-matrix
#     new_zma = automol.create.zmatrix.from_data(
#         new_symbols, new_key_matrix, new_name_matrix, new_values)
#
#     # get a reordered geometry
#     geo = automol.zmatrix.geometry(zma)
#     new_geo = (*geo[:row_idx],) + (*geo[row_idx+1:],) + tuple([geo[row_idx]])
#
#     # loop through recalc list and recalc stuff as necessary
#     for atom_idx, (name, idx) in enumerate(names_to_recalc):
#         # Get the keys needed to recalculate the coordinate
#         keys = new_key_matrix[i]
#         if idx in keys:
#             # Calculate the coordinate
#             if idx == 0:
#                 new_val = automol.geom.distance(
#                     new_geo, atom_idx, keys[0])
#             elif idx == 1:
#                 new_val = automol.geom.central_angle(
#                     new_geo, atom_idx, keys[0], keys[1])
#             elif idx == 2:
#                 new_val = automol.geom.dihedral_angle(
#                     new_geo, atom_idx, keys[0], keys[1], keys[2])
#             # Reset the value in the value dictionary
#             new_zma = set_values(new_zma, {name: new_val})
#
#     # set final values
#     atom_idx = count(zma) - 1
#     keys = new_key_matrix[-1]
#     print(keys)
#     new_zma = set_values(
#         new_zma, {'R1000': automol.geom.distance(
#             new_geo, atom_idx, keys[0])})
#     new_zma = set_values(
#         new_zma, {'A1000': automol.geom.central_angle(
#             new_geo, atom_idx, keys[0], keys[1])})
#     new_zma = set_values(
#         new_zma, {'D1000': automol.geom.dihedral_angle(
#             new_geo, atom_idx, keys[0], keys[1], keys[2])})
#
#     print(string(new_zma))
#     print('\n')
#     print(automol.geom.string(new_geo))
#     sys.exit()
#     return new_zma


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


# misc
def is_standard_form(zma):
    """ set standard variable names for the z-matrix

    (follows x2z format)
    """
    return is_valid(zma) and _v_.is_standard_form(var_(zma))


# I/O
def from_string(zma_str):
    """ read a z-matrix from a string
    """
    syms, key_mat, name_mat, val_dct = ar.zmatrix.read(zma_str)

    zma = automol.create.zmatrix.from_data(
        syms, key_mat, name_mat, val_dct, one_indexed=True, angstrom=True,
        degree=True)
    return zma


def string(zma):
    """ write a z-matrix to a string
    """
    zma_str = aw.zmatrix.write(
        syms=_v_.symbols(var_(zma)),
        key_mat=_v_.key_matrix(var_(zma), shift=1),
        name_mat=_v_.name_matrix(var_(zma)),
        val_dct=values(zma, angstrom=True, degree=True)
    )
    return zma_str


# comparisons
def almost_equal(zma1, zma2, rtol=2e-5, just_dist=False):
    """ are these z-matrices numerically equal?
    """
    ret = False
    if var_(zma1) == var_(zma2):
        val_dct1 = values(zma1)
        val_dct2 = values(zma2)
        dist_names = distance_names(zma1)
        dist_vals1 = tuple(map(val_dct1.__getitem__, dist_names))
        dist_vals2 = tuple(map(val_dct2.__getitem__, dist_names))
        if numpy.allclose(dist_vals1, dist_vals2, rtol=rtol):
            if just_dist:
                ret = True
            else:
                ang_names = angle_names(zma1)
                ang_vals1 = tuple(map(val_dct1.__getitem__, ang_names))
                ang_vals2 = tuple(map(val_dct2.__getitem__, ang_names))
                for shift in (0., numpy.pi/10.):
                    ang_vals1 = numpy.mod(numpy.add(ang_vals1, shift), 2*numpy.pi)
                    ang_vals2 = numpy.mod(numpy.add(ang_vals2, shift), 2*numpy.pi)
                    if numpy.allclose(ang_vals1, ang_vals2, rtol=rtol):
                        ret = True
                        break
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
def torsional_symmetry_numbers(zma, tors_names):
    """ symmetry numbers for torsional dihedrals
    """
    dih_edg_key_dct = _dihedral_edge_keys(zma)
    assert set(tors_names) <= set(dih_edg_key_dct.keys())
    edg_keys = tuple(map(dih_edg_key_dct.__getitem__, tors_names))

    gra = automol.convert.zmatrix.graph(zma, remove_stereo=True)
    bnd_sym_num_dct = automol.graph.bond_symmetry_numbers(gra)
    tors_sym_nums = []
    for edg_key in edg_keys:
        if edg_key in bnd_sym_num_dct.keys():
            sym_num = bnd_sym_num_dct[edg_key]
        else:
            sym_num = 1.
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


def torsional_sampling_ranges(zma, tors_names):
    """ sampling ranges for torsional dihedrals
    """
    sym_nums = torsional_symmetry_numbers(zma, tors_names)
    return tuple((0, 2*numpy.pi/sym_num) for sym_num in sym_nums)


def torsional_scan_linspaces(zma, tors_names, increment=0.5):
    """ scan grids for torsional dihedrals
    """
    sym_nums = torsional_symmetry_numbers(zma, tors_names)
    intervals = tuple(2*numpy.pi/sym_num - increment for sym_num in sym_nums)
    npoints_lst = tuple(
        (int(interval / increment)+1) for interval in intervals)
    return tuple((0, interval, npoints)
                 for interval, npoints in zip(intervals, npoints_lst))
