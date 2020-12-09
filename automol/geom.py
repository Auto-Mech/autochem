""" cartesian geometries
"""

import itertools
import functools
import more_itertools as mit
import numpy
from qcelemental import periodictable as pt
from qcelemental import constants as qcc
import autoread as ar
import autowrite as aw
import automol.create.geom
import automol.convert.geom
import automol.convert.inchi
from automol import cart
from automol.convert._pyx2z import to_oriented_geometry


# constructor
def from_data(syms, xyzs, angstrom=False):
    """ geometry data structure from symbols and coordinates
    """
    return automol.create.geom.from_data(
        symbols=syms, coordinates=xyzs, angstrom=angstrom)


def from_subset(geo, idxs):
    """ generate a new geometry from a subset of the atoms
    """
    syms = symbols(geo)
    xyzs = coordinates(geo)

    syms = list(map(syms.__getitem__, idxs))
    xyzs = list(map(xyzs.__getitem__, idxs))
    return from_data(syms, xyzs)


# getters
def symbols(geo, idxs=None):
    """ atomic symbols
    """
    idxs = list(range(count(geo))) if idxs is None else idxs

    if geo:
        syms, _ = zip(*geo)
    else:
        syms = ()

    syms = tuple(sym for idx, sym in enumerate(syms) if idx in idxs)
    return syms


def coordinates(geo, idxs=None, angstrom=False):
    """ atomic coordinates
    """
    idxs = list(range(count(geo))) if idxs is None else idxs
    if geo:
        _, xyzs = zip(*geo)
    else:
        xyzs = ()
    xyzs = xyzs if not angstrom else numpy.multiply(
        xyzs, qcc.conversion_factor('bohr', 'angstrom'))
    xyzs = tuple(xyz for idx, xyz in enumerate(xyzs) if idx in idxs)
    return xyzs


def count(geo):
    """ count the number of atoms in the geometry
    """
    return len(geo)


def atom_count(geo, sym, match=True):
    """ count the number of some atom type in the geometry
    """
    return len(atom_indices(geo, sym, match=match))


def atom_indices(geo, sym, match=True):
    """ indices for a particular atom type
        :param match: grab idxs that match given atom type
    """

    syms = symbols(geo)
    idxs = tuple()
    for idx, sym_ in enumerate(syms):
        if sym_ == sym and match:
            idxs += (idx,)
        elif sym_ != sym and not match:
            idxs += (idx,)

    return idxs


def dummy_atom_indices(geo):
    """ indices of dummy atoms in this geometry (Replace w/ above at some pt)
    """
    syms = symbols(geo)
    dummy_idxs = [idx for idx, sym in enumerate(syms) if not pt.to_Z(sym)]
    return tuple(dummy_idxs)


# validation
def is_valid(geo):
    """ is this a valid geometry?
    """
    ret = hasattr(geo, '__iter__')
    if ret:
        ret = all(hasattr(obj, '__len__') and len(obj) == 2 for obj in geo)
        if ret:
            syms, xyzs = zip(*geo)
            try:
                from_data(syms, xyzs)
            except AssertionError:
                ret = False
    return ret


def connected(geo, remove_stereo=False):
    """ Determine if all atoms in geometry are completely connected
    """
    return bool(len(components_graph(geo, remove_stereo=remove_stereo)) == 1)


# setters
def set_coordinates(geo, xyz_dct):
    """ set coordinate values for the geometry, using a dictionary by index
    """
    syms = symbols(geo)
    xyzs = coordinates(geo)

    natms = len(syms)
    assert all(idx in range(natms) for idx in xyz_dct)

    xyzs = [xyz_dct[idx] if idx in xyz_dct else xyz
            for idx, xyz in enumerate(xyzs)]
    return from_data(syms, xyzs)


def without_dummy_atoms(geo):
    """ return a copy of the geometry without dummy atoms
    """
    syms = symbols(geo)

    non_dummy_idxs = [idx for idx, sym in enumerate(syms) if pt.to_Z(sym)]
    return from_subset(geo, non_dummy_idxs)


# conversions
def zmatrix(geo, ts_bnds=()):
    """ geometry => z-matrix
    """
    return automol.convert.geom.zmatrix(geo, ts_bnds)


def zmatrix_torsion_coordinate_names(geo, ts_bnds=()):
    """ z-matrix torsional coordinate names
    """
    return automol.convert.geom.zmatrix_torsion_coordinate_names(geo, ts_bnds)


def zmatrix_atom_ordering(geo, ts_bnds=()):
    """ z-matrix atom ordering
    """
    return automol.convert.geom.zmatrix_atom_ordering(geo, ts_bnds)


def graph(geo, remove_stereo=False):
    """ geometry => graph
    """
    return automol.convert.geom.graph(
        geo, remove_stereo=remove_stereo)


def connectivity_graph(geo,
                       rqq_bond_max=3.5, rqh_bond_max=2.6, rhh_bond_max=1.9):
    """ geometry => connectivity graph
    """
    gra = automol.convert.geom.connectivity_graph(
        geo, rqq_bond_max=rqq_bond_max, rqh_bond_max=rqh_bond_max,
        rhh_bond_max=rhh_bond_max)
    return gra


def components_graph(geo, remove_stereo=False):
    """ geometry => connected components graphs
    """
    return automol.graph.connected_components(
        automol.convert.geom.graph(geo, remove_stereo=remove_stereo))


def inchi(geo, remove_stereo=False):
    """ geometry => inchi
    """
    return automol.convert.geom.inchi(geo, remove_stereo=remove_stereo)


def smiles(geo, remove_stereo=False):
    """ geometry => inchi
    """
    ich = inchi(geo, remove_stereo=remove_stereo)
    return automol.convert.inchi.smiles(ich)


def formula(geo):
    """ geometry => formula
    """
    return automol.convert.geom.formula(geo)


def remove(geo, idxs=()):
    """ Remove idxs from a geometry
    """
    new_geo = tuple(row for i, row in enumerate(geo) if i not in idxs)
    return new_geo


def end_group_sym_factor(geo, frm_bnd_keys=(), brk_bnd_keys=()):
    """ Determine sym factor for terminal groups in a geometry
    """

    # Set saddle based on frm and brk keys existing
    saddle = bool(frm_bnd_keys or brk_bnd_keys)

    gra = graph(geo, remove_stereo=True)
    term_atms = {}
    all_hyds = []
    neighbor_dct = automol.graph.atom_neighbor_keys(gra)

    # determine if atom is a part of a double bond
    unsat_atms = automol.graph.unsaturated_atom_keys(gra)
    if not saddle:
        rad_atms = automol.graph.sing_res_dom_radical_atom_keys(gra)
        res_rad_atms = automol.graph.resonance_dominant_radical_atom_keys(gra)
        rad_atms = [atm for atm in rad_atms if atm not in res_rad_atms]
    else:
        rad_atms = []

    gra = gra[0]
    for atm in gra:
        if gra[atm][0] == 'H':
            all_hyds.append(atm)
    for atm in gra:
        if atm in unsat_atms and atm not in rad_atms:
            pass
        else:
            if atm not in frm_bnd_keys and atm not in brk_bnd_keys:
                nonh_neighs = []
                h_neighs = []
                neighs = neighbor_dct[atm]
                for nei in neighs:
                    if nei in all_hyds:
                        h_neighs.append(nei)
                    else:
                        nonh_neighs.append(nei)
                if len(nonh_neighs) < 2 and len(h_neighs) > 1:
                    term_atms[atm] = h_neighs
    factor = 1.
    remove_atms = []
    for atm in term_atms:
        hyds = term_atms[atm]
        if len(hyds) > 1:
            factor *= len(hyds)
            remove_atms.extend(hyds)
    geo = remove(geo, remove_atms)
    return geo, factor


# operations
def join(geo1, geo2,
         dist_cutoff=3.*qcc.conversion_factor('angstrom', 'bohr'),
         theta=0.*qcc.conversion_factor('degree', 'radian'),
         phi=0.*qcc.conversion_factor('degree', 'radian')):
    """ join two geometries together
    """
    orient_vec = numpy.array([numpy.sin(theta) * numpy.cos(phi),
                              numpy.sin(theta) * numpy.sin(phi),
                              numpy.cos(theta)])
    neg_orient_vec = -1.0 * orient_vec

    # get the correct distance apart
    geo1 = mass_centered(geo1)
    geo2 = mass_centered(geo2)
    ext1 = max(numpy.vdot(orient_vec, xyz) for xyz in coordinates(geo1))
    ext2 = max(numpy.vdot(neg_orient_vec, xyz) for xyz in coordinates(geo2))

    cm_dist = ext1 + dist_cutoff + ext2
    dist_grid = numpy.arange(cm_dist, 0., -0.1)
    for dist in dist_grid:
        trans_geo2 = translate(geo2, orient_vec * dist)
        min_dist = minimum_distance(geo1, trans_geo2)
        if numpy.abs(min_dist - dist_cutoff) < 0.1:
            break

    geo2 = trans_geo2

    # now, join them together
    syms = symbols(geo1) + symbols(geo2)
    xyzs = coordinates(geo1) + coordinates(geo2)
    return from_data(syms, xyzs)


# I/O
def from_string(geo_str, angstrom=True):
    """ read a cartesian geometry from a string
    """
    syms, xyzs = ar.geom.read(geo_str)
    geo = from_data(syms, xyzs, angstrom=angstrom)
    return geo


def from_xyz_string(xyz_str):
    """ read a cartesian geometry from a .xyz string
    """
    syms, xyzs = ar.geom.read_xyz(xyz_str)
    geo = from_data(syms, xyzs, angstrom=True)
    return geo


def string(geo, angstrom=True):
    """ write the cartesian geometry as a string
    """
    syms = symbols(geo)
    xyzs = coordinates(geo, angstrom=angstrom)
    geo_str = aw.geom.write(syms=syms, xyzs=xyzs)
    return geo_str


def xyz_string(geo, comment=None):
    """ write the cartesian geometry to a .xyz string
    """
    syms = symbols(geo)
    xyzs = coordinates(geo, angstrom=True)
    geo_str = aw.geom.write_xyz(syms=syms, xyzs=xyzs, comment=comment)
    return geo_str


def xyz_trajectory_string(geo_lst, comments=None):
    """ write a series of cartesian geometries to a .xyz string
    """
    syms_lst = [symbols(geo) for geo in geo_lst]
    xyzs_lst = [coordinates(geo, angstrom=True) for geo in geo_lst]
    assert len(set(syms_lst)) == 1
    syms = syms_lst[0]
    xyz_traj_str = aw.geom.write_xyz_trajectory(syms, xyzs_lst,
                                                comments=comments)
    return xyz_traj_str


def from_xyz_trajectory_string(geo_str):
    """ read a series of cartesian geometries from a .xyz string
    """

    def _blocks(lst, size):
        """ Split list into parts of size n"""
        split_lst = []
        for i in range(0, len(lst), size):
            split_lst.append(lst[i:i + size])
        return split_lst

    # Split the lines for iteration
    geo_lines = [line for line in geo_str.splitlines()
                 if line != '']

    # Get the number of atoms used to partition the trajectory file
    line1 = geo_lines[0].strip()
    natoms = int(line1)

    geoms, comments = tuple(), tuple()
    for block in _blocks(geo_lines, natoms+2):
        comments += (block[1],)
        geoms += (from_string('\n'.join(block[2:])),)

    return (geoms, comments)


# representations
def coulomb_spectrum(geo):
    """ (sorted) coulomb matrix eigenvalue spectrum
    """
    mat = _coulomb_matrix(geo)
    vals = tuple(sorted(numpy.linalg.eigvalsh(mat)))
    return vals


def _coulomb_matrix(geo):
    nums = numpy.array(list(map(pt.to_Z, symbols(geo))))
    xyzs = numpy.array(coordinates(geo))

    _ = numpy.newaxis
    natms = len(nums)
    diag_idxs = numpy.diag_indices(natms)
    tril_idxs = numpy.tril_indices(natms, -1)
    triu_idxs = numpy.triu_indices(natms, 1)

    zxz = numpy.outer(nums, nums)
    rmr = numpy.linalg.norm(xyzs[:, _, :] - xyzs[_, :, :], axis=2)

    mat = numpy.zeros((natms, natms))
    mat[diag_idxs] = nums ** 2.4 / 2.
    mat[tril_idxs] = zxz[tril_idxs] / rmr[tril_idxs]
    mat[triu_idxs] = zxz[triu_idxs] / rmr[triu_idxs]
    return mat


# comparisons
def almost_equal(geo1, geo2, rtol=2e-3):
    """ are these geometries numerically equal?
    """
    ret = False
    if symbols(geo1) == symbols(geo2):
        ret = numpy.allclose(coordinates(geo1), coordinates(geo2), rtol=rtol)
    return ret


def minimum_distance(geo1, geo2):
    """ get the minimum distance between atoms in geo1 and those in geo2
    """
    xyzs1 = coordinates(geo1)
    xyzs2 = coordinates(geo2)
    return min(cart.vec.distance(xyz1, xyz2)
               for xyz1, xyz2 in itertools.product(xyzs1, xyzs2))


def almost_equal_coulomb_spectrum(geo1, geo2, rtol=1e-2):
    """ do these geometries have similar coulomb spectrums?
    """
    ret = numpy.allclose(coulomb_spectrum(geo1), coulomb_spectrum(geo2),
                         rtol=rtol)
    return ret


def argunique_coulomb_spectrum(geos, seen_geos=(), rtol=1e-2):
    """ get indices of unique geometries, by coulomb spectrum
    """
    comp_ = functools.partial(almost_equal_coulomb_spectrum, rtol=rtol)
    idxs = _argunique(geos, comp_, seen_items=seen_geos)
    return idxs


def _argunique(items, comparison, seen_items=()):
    """ get the indices of unique items using some comparison function
    """
    idxs = []
    seen_items = list(seen_items)
    for idx, item in enumerate(items):
        if not any(comparison(item, seen_item) for seen_item in seen_items):
            idxs.append(idx)
            seen_items.append(item)
    idxs = tuple(idxs)
    return idxs


# transformations
def reorder(geo, idx_dct):
    """ Reorder the atoms in this geometry

    :param geo: The geometry
    :param idx_dct: The new order of the atoms, by index
    :type idx_dct: dict
    """
    syms = symbols(geo)
    xyzs = coordinates(geo)

    idxs = [idx for idx, _ in sorted(idx_dct.items(), key=lambda x: x[1])]
    assert len(syms) == len(xyzs) == len(idxs)

    syms = [syms[idx] for idx in idxs]
    xyzs = [xyzs[idx] for idx in idxs]
    return from_data(syms, xyzs)


def swap_coordinates(geo, idx1, idx2):
    """ swap the order of the coordinates of the two atoms
    """
    geo = [list(x) for x in geo]
    geo[idx1], geo[idx2] = geo[idx2], geo[idx1]
    geo_swp = tuple(tuple(x) for x in geo)
    return geo_swp


def insert(geo, sym, xyz, idx=None):
    """ Insert an atom into this geometry
    """
    syms = list(symbols(geo))
    xyzs = list(coordinates(geo))

    idx = idx if idx is not None else len(syms)

    syms.insert(idx, sym)
    xyzs.insert(idx, xyz)
    return from_data(syms, xyzs)


def displace(geo, xyzs):
    """ displacement of the geometry
    """
    syms = symbols(geo)
    orig_xyzs = coordinates(geo)
    xyzs = numpy.add(orig_xyzs, xyzs)
    return from_data(syms, xyzs)


def translate(geo, xyz):
    """ translation of the geometry
    """
    syms = symbols(geo)
    xyzs = coordinates(geo)
    xyzs = numpy.add(xyzs, xyz)
    return from_data(syms, xyzs)


def invert(geo):
    """ inversion of the geometry
    """
    syms = symbols(geo)
    xyzs = numpy.array(coordinates(geo))
    xyzs *= -1.
    xyzs = list(map(list, xyzs))
    return from_data(syms, xyzs)


def transform(geo, func, idxs=None):
    """ transform the coordinates of a geometry by a function

    if transforming a subset, specify the atoms indices with `idxs`
    """
    idxs = list(range(count(geo))) if idxs is None else idxs
    syms = symbols(geo)
    xyzs = coordinates(geo)
    xyzs = [func(xyz) if idx in idxs else xyz for idx, xyz in enumerate(xyzs)]
    return from_data(syms, xyzs)


def transform_by_matrix(geo, mat):
    """ transform the coordinates of a geometry by a matrix
    """
    syms = symbols(geo)
    xyzs = coordinates(geo)
    xyzs = numpy.dot(xyzs, numpy.transpose(mat))
    return from_data(syms, xyzs)


def rotate(geo, axis, angle, orig_xyz=None, idxs=None):
    """ axis-angle rotation of the geometry

    if rotating a subset, specify the atoms indices with `idxs`
    """
    func = cart.vec.rotater(axis, angle, orig_xyz=orig_xyz)
    return transform(geo, func, idxs=idxs)


def euler_rotate(geo, theta, phi, psi):
    """ axis-angle rotation of the geometry
    """
    mat = cart.mat.euler_rotation(theta, phi, psi)
    return transform_by_matrix(geo, mat)


def move_coordinates(geo, idx1, idx2):
    """ move the atom at position idx1 to idx2, shifting all other atoms
    """

    # Get the coordinates at idx1 that are to be moved
    geo = [list(x) for x in geo]
    moving_coords = geo[idx1]

    # move the coordinates to idx2
    geo.remove(moving_coords)
    geo.insert(idx2, moving_coords)
    geo_move = tuple(tuple(x) for x in geo)
    return geo_move


def reflect_coordinates(geo, idxs, axes):
    """ reflect each coordinate about the requested axes
    """

    # check input
    assert all(idx < len(geo) for idx in idxs)
    assert all(axis in ('x', 'y', 'z') for axis in axes)

    # get coords
    coords = coordinates(geo)

    # convert x,y,z to nums
    axis_dct = {'x': 0, 'y': 1, 'z': 2}
    axes = [axis_dct[axis] for axis in axes]

    # build set atom dct with relected coords
    reflect_dct = {}
    for idx in idxs:
        coord_lst = list(coords[idx])
        for axis in axes:
            coord_lst[axis] *= -1.0
        reflect_dct[idx] = coord_lst

    # Reflect coords with dct
    geo_reflected = set_coordinates(geo, reflect_dct)

    return geo_reflected


def rot_permutated_geoms(geo, frm_bnd_keys=(), brk_bnd_keys=()):
    """ convert an input geometry to a list of geometries
        corresponding to the rotational permuations of all the terminal groups
    """

    # Set saddle based on frm and brk keys existing
    saddle = bool(frm_bnd_keys or brk_bnd_keys)

    gra = graph(geo, remove_stereo=True)
    term_atms = {}
    all_hyds = []
    neighbor_dct = automol.graph.atom_neighbor_keys(gra)

    # determine if atom is a part of a double bond
    unsat_atms = automol.graph.unsaturated_atom_keys(gra)
    if not saddle:
        rad_atms = automol.graph.sing_res_dom_radical_atom_keys(gra)
        res_rad_atms = automol.graph.resonance_dominant_radical_atom_keys(gra)
        rad_atms = [atm for atm in rad_atms if atm not in res_rad_atms]
    else:
        rad_atms = []

    gra = gra[0]
    for atm in gra:
        if gra[atm][0] == 'H':
            all_hyds.append(atm)
    for atm in gra:
        if atm in unsat_atms and atm not in rad_atms:
            pass
        else:
            if atm not in frm_bnd_keys and atm not in brk_bnd_keys:
                nonh_neighs = []
                h_neighs = []
                neighs = neighbor_dct[atm]
                for nei in neighs:
                    if nei in all_hyds:
                        h_neighs.append(nei)
                    else:
                        nonh_neighs.append(nei)
                if len(nonh_neighs) < 2 and len(h_neighs) > 1:
                    term_atms[atm] = h_neighs
    geo_final_lst = [geo]
    for atm in term_atms:
        hyds = term_atms[atm]
        geo_lst = []
        for geom in geo_final_lst:
            geo_lst.extend(_swap_for_one(geom, hyds))
        geo_final_lst = geo_lst
    return geo_final_lst


def _swap_for_one(geo, hyds):
    """ rotational permuation for one rotational group
    """
    geo_lst = []
    if len(hyds) > 1:
        new_geo = geo
        if len(hyds) > 2:
            geo_lst.append(new_geo)
            new_geo = swap_coordinates(new_geo, hyds[0], hyds[1])
            new_geo = swap_coordinates(new_geo, hyds[0], hyds[2])
            geo_lst.append(new_geo)
            new_geo = swap_coordinates(new_geo, hyds[0], hyds[1])
            new_geo = swap_coordinates(new_geo, hyds[0], hyds[2])
            geo_lst.append(new_geo)
        else:
            geo_lst.append(new_geo)
            new_geo = swap_coordinates(new_geo, hyds[0], hyds[1])
            geo_lst.append(new_geo)
    return geo_lst


# geometric properties
def distance(geo, idx1, idx2):
    """ measure the distance between atoms
    """
    xyzs = coordinates(geo)
    xyz1 = xyzs[idx1]
    xyz2 = xyzs[idx2]
    return cart.vec.distance(xyz1, xyz2)


def central_angle(geo, idx1, idx2, idx3):
    """ measure the angle inscribed by three atoms
    """
    xyzs = coordinates(geo)
    xyz1 = xyzs[idx1]
    xyz2 = xyzs[idx2]
    xyz3 = xyzs[idx3]
    return cart.vec.central_angle(xyz1, xyz2, xyz3)


def dihedral_angle(geo, idx1, idx2, idx3, idx4):
    """ measure the dihedral angle defined by four atoms
    """
    xyzs = coordinates(geo)
    xyz1 = xyzs[idx1]
    xyz2 = xyzs[idx2]
    xyz3 = xyzs[idx3]
    xyz4 = xyzs[idx4]
    return cart.vec.dihedral_angle(xyz1, xyz2, xyz3, xyz4)


def distance_matrix(geo):
    """form distance matrix for a set of xyz coordinates
    """
    mat = numpy.zeros((len(geo), len(geo)))
    for i in range(len(geo)):
        for j in range(len(geo)):
            mat[i][j] = distance(geo, i, j)
    return mat


def almost_equal_dist_matrix(geo1, geo2, thresh=0.1):
    """form distance matrix for a set of xyz coordinates
    """
    for i in range(len(geo1)):
        for j in range(len(geo1)):
            dist_mat1_ij = distance(geo1, i, j)
            dist_mat2_ij = distance(geo2, i, j)
            if abs(dist_mat1_ij - dist_mat2_ij) > thresh:
                return False
    return True


def external_symmetry_factor(geo):
    """ obtain external symmetry factor for a geometry using x2z
    """
    # Get initial external symmetry number
    if automol.geom.is_atom(geo):
        ext_sym_fac = 1.
    else:
        oriented_geom = to_oriented_geometry(geo)
        ext_sym_fac = oriented_geom.sym_num()
        # Divide symmetry number by enantiomeric factor
        if oriented_geom.is_enantiomer():
            ext_sym_fac *= 0.5
    return ext_sym_fac


def find_xyzp_using_internals(xyz1, xyz2, xyz3, pdist, pangle, pdihed):
    """ geometric approach for calculating the xyz coordinates of atom A
        when the xyz coordinates of the A B and C are known and
        the position is defined w/r to A B C with internal coordinates

    TODO: This already exists: cart.vec.from_internals does exactly the same
    thing; find where this is used and replace it with that
    """

    # Set to numpy arrays
    xyz1 = numpy.array(xyz1)
    xyz2 = numpy.array(xyz2)
    xyz3 = numpy.array(xyz3)

    # Set the coordinates of Point P in the RT system
    xyzp_rt = numpy.array([pdist * numpy.sin(pangle) * numpy.cos(pdihed),
                           pdist * numpy.cos(pangle),
                           -(pdist * numpy.sin(pangle) * numpy.sin(pdihed))
                           ])

    # Set the coordinates of the Point 2 and 3 in the RT system
    dist12 = numpy.linalg.norm(xyz1 - xyz2)
    dist13 = numpy.linalg.norm(xyz1 - xyz3)
    dist23 = numpy.linalg.norm(xyz2 - xyz3)
    xyz2_rt = numpy.array([0.0, dist12, 0.0])

    val = ((dist12**2 + dist13**2 - dist23**2) / 2.0 / dist12)
    valx3 = numpy.sqrt(dist13**2 - val**2)
    valy3 = ((dist12**2 + dist13**2 - dist23**2) / 2.0 / dist12)
    xyz3_rt = numpy.array([valx3, valy3, 0.0])

    # Translate original frame of ref coors so that xyz1 is at (0, 0, 0)
    xyz2_t = xyz2 - xyz1
    xyz3_t = xyz3 - xyz1

    # Rotation matrix to rotate back to the original ref system
    r12 = (xyz2[0] - xyz1[0]) / xyz2_rt[1]
    r22 = (xyz2[1] - xyz1[1]) / xyz2_rt[1]
    r32 = (xyz2[2] - xyz1[2]) / xyz2_rt[1]

    r11 = (xyz3[0] - xyz1[0] - xyz3_rt[1]*r12) / xyz3_rt[0]
    r21 = (xyz3[1] - xyz1[1] - xyz3_rt[1]*r22) / xyz3_rt[0]
    r31 = (xyz3[2] - xyz1[2] - xyz3_rt[1]*r32) / xyz3_rt[0]

    anum_aconst = xyz2_t[1] - (xyz3_t[1] / xyz3_t[0]) * xyz2_t[0]
    den_aconst = xyz2_t[2] - (xyz3_t[2] / xyz3_t[0]) * xyz2_t[0]

    if abs(anum_aconst) < 1.0e-6 and abs(den_aconst) < 1.0e-6:
        if anum_aconst < 0.0:
            aconst = -1.0e20
        else:
            aconst = 1.0e20
    elif abs(den_aconst) < 1.0e-6:
        if anum_aconst < 0.0:
            aconst = -1.0e20
        else:
            aconst = 1.0e20
    else:
        anum = xyz2_t[1] - (xyz3_t[1] / xyz3_t[0]) * xyz2_t[0]
        aden = xyz2_t[2] - (xyz3_t[2] / xyz3_t[0]) * xyz2_t[0]
        aconst = anum / aden

    den1 = (xyz3_t[1] / xyz3_t[0]) - aconst * (xyz3_t[2] / xyz3_t[0])
    if den1 == 0.0:
        den1 = 1.0e-20
    bconst = 1.0 / den1

    # Set vals for another point
    valx = -(1.0 / numpy.sqrt(1.0 + (bconst**2) * (1.0 + aconst**2)))
    valy = -(valx * bconst)
    xyz4_t = numpy.array([valx, valy, -(valy * aconst)])

    r13 = xyz4_t[0]
    r23 = xyz4_t[1]
    r33 = xyz4_t[2]
    r13n = -r13
    r23n = -r23
    r33n = -r33

    # Now rotate and translate back
    # Here I check  the (001) vector direction to decide whether
    # To take the positive of negative results of square root taken above
    xap = (xyz1[0] + (r11 * xyzp_rt[0]) +
           (r12 * xyzp_rt[1]) + (r13 * xyzp_rt[2]))
    yap = (xyz1[1] + (r21 * xyzp_rt[0]) +
           (r22 * xyzp_rt[1]) + (r33 * xyzp_rt[2]))
    zap = (xyz1[2] + (r31 * xyzp_rt[0]) +
           (r32 * xyzp_rt[1]) + (r33 * xyzp_rt[2]))

    xan = (xyz1[0] + (r11 * xyzp_rt[0]) +
           (r12 * xyzp_rt[1]) + (r13n * xyzp_rt[2]))
    yan = (xyz1[1] + (r21 * xyzp_rt[0]) +
           (r22 * xyzp_rt[1]) + (r23n * xyzp_rt[2]))
    zan = (xyz1[2] + (r31 * xyzp_rt[0]) +
           (r32 * xyzp_rt[1]) + (r33n * xyzp_rt[2]))

    bvec = xyz1 - xyz2
    cvec = xyz2 - xyz3
    vec1 = (bvec[1] * cvec[2]) - (bvec[2] * cvec[1])
    vec2 = (bvec[2] * cvec[0]) - (bvec[0] * cvec[2])
    vec3 = (bvec[0] * cvec[1]) - (bvec[1] * cvec[0])

    if abs(xyz4_t[0]) > 1.0e-5:
        checkv = vec1 / xyz4_t[0]
    elif abs(xyz4_t[1]) > 1.0e-5:
        checkv = vec2 / xyz4_t[1]
    else:
        checkv = vec3 / xyz4_t[2]

    if checkv >= 0.0:
        xyzp = numpy.array([xap, yap, zap])
    else:
        xyzp = numpy.array([xan, yan, zan])

    return xyzp[0], xyzp[1], xyzp[2]


def closest_unbonded_atoms(geo, gra):
    """ determine the closest unbonded pair of atoms in this geometry
    """
    atm_keys = automol.graph.atom_keys(gra)
    bnd_keys = automol.graph.bond_keys(gra)
    poss_bnd_keys = set(map(frozenset, itertools.combinations(atm_keys, r=2)))

    # The set of candidates includes all unbonded pairs of atoms
    cand_bnd_keys = poss_bnd_keys - bnd_keys

    min_bnd_key = None
    min_dist_val = 1000.
    for bnd_key in cand_bnd_keys:
        dist_val = distance(geo, *bnd_key)
        if dist_val < min_dist_val:
            min_dist_val = dist_val
            min_bnd_key = bnd_key

    return min_bnd_key, min_dist_val


# chemical properties
def is_atom(geo):
    """ return return the atomic masses
    """
    syms = symbols(geo)
    ret = False
    if len(syms) == 1:
        ret = True
    return ret


def masses(geo, amu=True):
    """ return the atomic masses
    """
    syms = symbols(geo)
    amas = list(map(pt.to_mass, syms))

    if not amu:
        conv = qcc.conversion_factor("atomic_mass_unit", "electron_mass")
        amas = numpy.multiply(amas, conv)

    amas = tuple(amas)
    return amas


def total_mass(geo):
    """ Calculate the total mass
    """
    return sum(masses(geo))


def center_of_mass(geo):
    """ center of mass
    """
    xyzs = coordinates(geo)
    amas = masses(geo)
    cm_xyz = tuple(
        sum(numpy.multiply(xyz, ama) for xyz, ama in zip(xyzs, amas)) /
        sum(amas))

    return cm_xyz


def reduced_mass(geo1, geo2):
    """ Calculate the reduced mass for two species.

        :param geo1: geometry of species 1 (Bohr)
        :type geo1: list(float)
        :param geo2: geometry of species 2 (Bohr)
        :type geo2: list(float)
        :return: reduced mass (amu)
        :rtype: float
    """

    mass1 = automol.geom.total_mass(geo1)
    mass2 = automol.geom.total_mass(geo2)

    return (mass1 * mass2) / (mass1 + mass2)


def mass_centered(geo):
    """ mass-centered geometry
    """
    geo = translate(geo, numpy.negative(center_of_mass(geo)))
    return geo


def inertia_tensor(geo, amu=True):
    """ molecula# r inertia tensor (atomic units if amu=False)
    """
    geo = mass_centered(geo)
    amas = masses(geo, amu=amu)
    xyzs = coordinates(geo)
    ine = tuple(map(tuple, sum(
        ama * (numpy.vdot(xyz, xyz) * numpy.eye(3) - numpy.outer(xyz, xyz))
        for ama, xyz in zip(amas, xyzs))))

    return ine


def principal_axes(geo, amu=True):
    """ principal inertial axes (atomic units if amu=False)
    """
    ine = inertia_tensor(geo, amu=amu)
    _, paxs = numpy.linalg.eigh(ine)
    paxs = tuple(map(tuple, paxs))
    return paxs


def moments_of_inertia(geo, amu=True):
    """ principal inertial axes (atomic units if amu=False)
    """
    ine = inertia_tensor(geo, amu=amu)
    moms, _ = numpy.linalg.eigh(ine)
    moms = tuple(moms)
    return moms


def rotational_constants(geo, amu=True):
    """ rotational constants (atomic units if amu=False)
    """
    moms = moments_of_inertia(geo, amu=amu)
    sol = (qcc.get('speed of light in vacuum') *
           qcc.conversion_factor('meter / second', 'bohr hartree / h'))
    cons = numpy.divide(1., moms) / 4. / numpy.pi / sol
    cons = tuple(cons)
    return cons


def is_linear(geo, tol=2.*qcc.conversion_factor('degree', 'radian')):
    """ is this geometry linear?
    """
    ret = True

    if len(geo) == 1:
        ret = False
    elif len(geo) == 2:
        ret = True
    else:
        keys = range(len(symbols(geo)))
        for key1, key2, key3 in mit.windowed(keys, 3):
            cangle = numpy.abs(central_angle(geo, key1, key2, key3))
            if not (numpy.abs(cangle) < tol or
                    numpy.abs(cangle - numpy.pi) < tol):
                ret = False
    return ret


def permutation(geo, ref_geo, thresh=1e-4):
    """ determine the permutation of one geometry that reproduces another

    (If there isn't one -- the geometries are not aligned, return None)
    """
    natms = count(geo)
    syms = symbols(geo)
    xyzs = coordinates(geo)

    perm_idxs = [None] * natms
    for idx, (sym, xyz) in enumerate(zip(syms, xyzs)):
        # Loop over atoms in the reference geometry with the same symbol
        ref_idxs = atom_indices(ref_geo, sym=sym)
        ref_xyzs = coordinates(ref_geo, idxs=ref_idxs)
        perm_idx = next(
            (ref_idx for ref_idx, ref_xyz in zip(ref_idxs, ref_xyzs)
             if cart.vec.distance(xyz, ref_xyz) < thresh), None)
        perm_idxs[idx] = perm_idx

    perm_idxs = tuple(perm_idxs)

    if any(perm_idx is None for perm_idx in perm_idxs):
        perm_idxs = None

    return perm_idxs


if __name__ == '__main__':
    GEO_STR = """
C    0.000000   0.000000   0.000000
C    0.000000   0.000000   1.528500
H    0.000000   1.021271  -0.395008
H   -0.884960  -0.505333  -0.402273
H    0.884431  -0.514168  -0.394620
C    0.131848  -1.396208   2.150831
H    0.830707   0.618974   1.889645
H   -0.915682   0.483661   1.891358
C   -0.935918  -2.407074   1.716757
H    0.105716  -1.297371   3.244463
H    1.121314  -1.802704   1.901432
C   -2.368078  -1.949468   2.003020
H   -0.754135  -3.357275   2.235849
H   -0.822109  -2.622070   0.646658
C   -3.429891  -2.998316   1.719404
H   -2.443822  -1.645020   3.056185
H   -2.587464  -1.045228   1.414070
C   -3.441270  -3.623096   0.332366
H   -3.196344  -3.971965   2.538467
H   -4.422939  -2.656165   2.030354
C   -3.487708  -2.607740  -0.819639
H   -2.548685  -4.250408   0.207484
H   -4.303285  -4.298844   0.256996
C   -4.697641  -1.673392  -0.777069
H   -2.564651  -2.015142  -0.805424
H   -3.484882  -3.150427  -1.774953
C   -4.756162  -0.678230  -1.939836
H   -4.695662  -1.111648   0.167050
H   -5.614566  -2.278478  -0.770618
C   -3.584557   0.306860  -1.978679
H   -5.693215  -0.109020  -1.871759
H   -4.798868  -1.225952  -2.892242
C   -3.760205   1.390162  -3.041302
H   -3.472805   0.773991  -0.989972
H   -2.648859  -0.232888  -2.169286
H   -4.662764   1.983230  -2.853383
H   -2.905957   2.074798  -3.059709
H   -3.856588   0.947583  -4.039544
X   -3.239107  -3.334837   3.308037
C   -2.887701  -5.027299   3.409914
H   -3.850477  -5.459083   3.681271
H   -2.368177  -4.557577   4.244852
H   -2.257867  -5.697614   2.825405
"""
    GEO = from_string(GEO_STR)
