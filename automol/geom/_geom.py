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
import automol.graph
import automol.create.geom
import automol.convert.geom
import automol.convert.inchi
from automol import util

BOHR2ANG = qcc.conversion_factor('bohr', 'angstrom')
RAD2DEG = qcc.conversion_factor('radian', 'degree')


# constructors
def from_subset(geo, idxs):
    """ generate a new geometry from a subset of the atoms

    (Rename this and put it under operations?)
    """
    symbs = symbols(geo)
    xyzs = coordinates(geo)

    symbs = list(map(symbs.__getitem__, idxs))
    xyzs = list(map(xyzs.__getitem__, idxs))
    return automol.create.geom.from_data(symbs, xyzs)


# getters
def symbols(geo, idxs=None):
    """ atomic symbols
    """
    idxs = list(range(count(geo))) if idxs is None else idxs

    if geo:
        symbs, _ = zip(*geo)
    else:
        symbs = ()

    symbs = tuple(symb for idx, symb in enumerate(symbs) if idx in idxs)
    return symbs


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


def atom_count(geo, symb, match=True):
    """ count the number of some atom type in the geometry
    """
    return len(atom_indices(geo, symb, match=match))


def atom_indices(geo, symb, match=True):
    """ indices for a particular atom type
        :param match: grab idxs that match given atom type
    """

    symbs = symbols(geo)
    idxs = tuple()
    for idx, symb_ in enumerate(symbs):
        if symb_ == symb and match:
            idxs += (idx,)
        elif symb_ != symb and not match:
            idxs += (idx,)

    return idxs


def dummy_atom_indices(geo):
    """ indices of dummy atoms in this geometry (Replace w/ above at some pt)
    """
    symbs = symbols(geo)
    dummy_idxs = [idx for idx, symb in enumerate(symbs) if not pt.to_Z(symb)]
    return tuple(dummy_idxs)


def components_graph(geo, remove_stereo=False):
    """ geometry => connected components graphs
    """
    return automol.graph.connected_components(
        automol.convert.geom.graph(geo, remove_stereo=remove_stereo))


def connected(geo, remove_stereo=False):
    """ Determine if all atoms in geometry are completely connected
    """
    comps = components_graph(geo, remove_stereo=remove_stereo)
    return bool(len(comps) == 1)


# validation
def is_valid(geo):
    """ is this a valid geometry?
    """
    ret = hasattr(geo, '__iter__')
    if ret:
        ret = all(hasattr(obj, '__len__') and len(obj) == 2 for obj in geo)
        if ret:
            symbs, xyzs = zip(*geo)
            try:
                automol.create.geom.from_data(symbs, xyzs)
            except AssertionError:
                ret = False
    return ret


# setters
def set_coordinates(geo, xyz_dct):
    """ set coordinate values for the geometry, using a dictionary by index
    """
    symbs = symbols(geo)
    xyzs = coordinates(geo)

    natms = len(symbs)
    assert all(idx in range(natms) for idx in xyz_dct)

    xyzs = [xyz_dct[idx] if idx in xyz_dct else xyz
            for idx, xyz in enumerate(xyzs)]
    return automol.create.geom.from_data(symbs, xyzs)


def without_dummy_atoms(geo):
    """ return a copy of the geometry without dummy atoms
    """
    symbs = symbols(geo)

    non_dummy_idxs = [idx for idx, symb in enumerate(symbs) if pt.to_Z(symb)]
    return from_subset(geo, non_dummy_idxs)


# I/O
def from_string(geo_str, angstrom=True):
    """ read a cartesian geometry from a string
    """
    symbs, xyzs = ar.geom.read(geo_str)
    geo = automol.create.geom.from_data(symbs, xyzs, angstrom=angstrom)
    return geo


def from_xyz_string(xyz_str):
    """ read a cartesian geometry from a .xyz string
    """
    symbs, xyzs = ar.geom.read_xyz(xyz_str)
    geo = automol.create.geom.from_data(symbs, xyzs, angstrom=True)
    return geo


def string(geo, angstrom=True):
    """ write the cartesian geometry as a string
    """
    symbs = symbols(geo)
    xyzs = coordinates(geo, angstrom=angstrom)
    geo_str = aw.geom.write(symbs=symbs, xyzs=xyzs)
    return geo_str


def xyz_string(geo, comment=None):
    """ write the cartesian geometry to a .xyz string
    """
    symbs = symbols(geo)
    xyzs = coordinates(geo, angstrom=True)
    geo_str = aw.geom.write_xyz(symbs=symbs, xyzs=xyzs, comment=comment)
    return geo_str


def xyz_trajectory_string(geo_lst, comments=None):
    """ write a series of cartesian geometries to a .xyz string
    """
    symbs_lst = [symbols(geo) for geo in geo_lst]
    xyzs_lst = [coordinates(geo, angstrom=True) for geo in geo_lst]
    assert len(set(symbs_lst)) == 1
    symbs = symbs_lst[0]
    xyz_traj_str = aw.geom.write_xyz_trajectory(symbs, xyzs_lst,
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
    return min(util.vec.distance(xyz1, xyz2)
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


# operations
def remove(geo, idxs=()):
    """ Remove idxs from a geometry
    """
    new_geo = tuple(row for i, row in enumerate(geo) if i not in idxs)
    return new_geo


def join(geo1, geo2,
         dist_cutoff=3.*qcc.conversion_factor('angstrom', 'bohr'),
         theta=0.*qcc.conversion_factor('degree', 'radian'),
         phi=0.*qcc.conversion_factor('degree', 'radian')):
    """ join two geometries together
    """
    if not geo1:
        symbs = symbols(geo2)
        xyzs = coordinates(geo2)
    elif not geo2:
        symbs = symbols(geo1)
        xyzs = coordinates(geo1)
    else:
        orient_vec = numpy.array([numpy.sin(theta) * numpy.cos(phi),
                                  numpy.sin(theta) * numpy.sin(phi),
                                  numpy.cos(theta)])
        neg_orient_vec = -1.0 * orient_vec

        # get the correct distance apart
        geo1 = mass_centered(geo1)
        geo2 = mass_centered(geo2)
        ext1 = max(numpy.vdot(orient_vec, xyz) for xyz in coordinates(geo1))
        ext2 = max(numpy.vdot(neg_orient_vec, xyz)
                   for xyz in coordinates(geo2))

        cm_dist = ext1 + dist_cutoff + ext2
        dist_grid = numpy.arange(cm_dist, 0., -0.1)
        for dist in dist_grid:
            trans_geo2 = translate(geo2, orient_vec * dist)
            min_dist = minimum_distance(geo1, trans_geo2)
            if numpy.abs(min_dist - dist_cutoff) < 0.1:
                break

        geo2 = trans_geo2

        # now, join them together
        symbs = symbols(geo1) + symbols(geo2)
        xyzs = coordinates(geo1) + coordinates(geo2)

    return automol.create.geom.from_data(symbs, xyzs)


def reorder(geo, idx_dct):
    """ Reorder the atoms in this geometry

    :param geo: The geometry
    :param idx_dct: The new order of the atoms, by index
    :type idx_dct: dict
    """
    symbs = symbols(geo)
    xyzs = coordinates(geo)

    idxs = [idx for idx, _ in sorted(idx_dct.items(), key=lambda x: x[1])]
    assert len(symbs) == len(xyzs) == len(idxs)

    symbs = [symbs[idx] for idx in idxs]
    xyzs = [xyzs[idx] for idx in idxs]
    return automol.create.geom.from_data(symbs, xyzs)


def swap_coordinates(geo, idx1, idx2):
    """ swap the order of the coordinates of the two atoms
    """
    geo = [list(x) for x in geo]
    geo[idx1], geo[idx2] = geo[idx2], geo[idx1]
    geo_swp = tuple(tuple(x) for x in geo)
    return geo_swp


def insert(geo, symb, xyz, idx=None, angstrom=False):
    """ Insert an atom into this geometry
    """
    symbs = list(symbols(geo))
    xyzs = list(coordinates(geo, angstrom=angstrom))

    idx = idx if idx is not None else len(symbs)

    symbs.insert(idx, symb)
    xyzs.insert(idx, xyz)
    return automol.create.geom.from_data(symbs, xyzs, angstrom=angstrom)


def insert_dummies_on_linear_atoms(geo, lin_idxs=None, gra=None, dist=1.,
                                   tol=5.):
    """ Insert dummy atoms over linear atoms in the geometry

    :param geo: the geometry
    :param lin_idxs: the indices of the linear atoms; if None, indices will be
        automatically determined from the geometry based on the graph
    :param gra: the graph describing connectivity; if None, a connectivity
        graph will be generated using default distance thresholds
    :param dist: distance of the dummy atom from the linear atom, in angstroms
    :param tol: the tolerance threshold for linearity, in degrees
    :returns: the geometry with dummy atoms inserted, along with a dictionary
        mapping the linear atoms onto their associated dummy atoms
    """
    lin_idxs = linear_atoms(geo) if lin_idxs is None else lin_idxs
    gra = automol.convert.geom.connectivity_graph(geo) if gra is None else gra

    dummy_ngb_idxs = automol.graph.dummy_atom_neighbor_keys(gra)
    assert not dummy_ngb_idxs & set(lin_idxs), (
        "Attempting to add dummy atoms on atoms that already have them: {}"
        .format(dummy_ngb_idxs & set(lin_idxs)))

    ngb_idxs_dct = automol.graph.sorted_atom_neighbor_keys(gra)

    xyzs = coordinates(geo, angstrom=True)

    def _perpendicular_direction(idxs):
        """ find a nice perpendicular direction for a series of linear atoms
        """
        triplets = []
        for idx in idxs:
            for n1idx in ngb_idxs_dct[idx]:
                for n2idx in ngb_idxs_dct[n1idx]:
                    if n2idx != idx:
                        ang = central_angle(geo, idx, n1idx, n2idx,
                                            degree=True)
                        if numpy.abs(ang - 180.) > tol:
                            triplets.append((idx, n1idx, n2idx))

        if triplets:
            idx1, idx2, idx3 = min(triplets, key=lambda x: x[1:])
            xyz1, xyz2, xyz3 = map(xyzs.__getitem__, (idx1, idx2, idx3))
            r12 = util.vec.unit_direction(xyz1, xyz2)
            r23 = util.vec.unit_direction(xyz2, xyz3)
            direc = util.vec.orthogonalize(r12, r23, normalize=True)
        else:
            assert len(idxs) >= 2, "This should never happen."
            idx1, idx2 = idxs[:2]
            xyz1, xyz2 = map(xyzs.__getitem__, (idx1, idx2))
            r12 = util.vec.unit_direction(xyz1, xyz2)
            for i in range(3):
                disp = numpy.zeros((3,))
                disp[i] = -1.
                alt = numpy.add(r12, disp)
                direc = util.vec.unit_perpendicular(r12, alt)
                if numpy.linalg.norm(direc) > 1e-2:
                    break
        return direc

    # partition the linear atoms into adjacent groups, to be handled together
    lin_idxs_lst = sorted(map(sorted, util.equivalence_partition(
        lin_idxs, lambda x, y: x in ngb_idxs_dct[y])))

    dummy_idx_dct = {}

    for idxs in lin_idxs_lst:
        direc = _perpendicular_direction(idxs)
        for idx in idxs:
            xyz = numpy.add(xyzs[idx], numpy.multiply(dist, direc))
            dummy_idx_dct[idx] = count(geo)

            geo = insert(geo, 'X', xyz, angstrom=True)

    return geo, dummy_idx_dct


def displace(geo, xyzs):
    """ displacement of the geometry
    """
    symbs = symbols(geo)
    orig_xyzs = coordinates(geo)
    xyzs = numpy.add(orig_xyzs, xyzs)
    return automol.create.geom.from_data(symbs, xyzs)


def translate(geo, xyz):
    """ translation of the geometry
    """
    symbs = symbols(geo)
    xyzs = coordinates(geo)
    xyzs = numpy.add(xyzs, xyz)
    return automol.create.geom.from_data(symbs, xyzs)


def invert(geo):
    """ inversion of the geometry
    """
    symbs = symbols(geo)
    xyzs = numpy.array(coordinates(geo))
    xyzs *= -1.
    xyzs = list(map(list, xyzs))
    return automol.create.geom.from_data(symbs, xyzs)


def transform(geo, func, idxs=None):
    """ transform the coordinates of a geometry by a function

    if transforming a subset, specify the atoms indices with `idxs`
    """
    idxs = list(range(count(geo))) if idxs is None else idxs
    symbs = symbols(geo)
    xyzs = coordinates(geo)
    xyzs = [func(xyz) if idx in idxs else xyz for idx, xyz in enumerate(xyzs)]
    return automol.create.geom.from_data(symbs, xyzs)


def transform_by_matrix(geo, mat):
    """ transform the coordinates of a geometry by a matrix
    """
    symbs = symbols(geo)
    xyzs = coordinates(geo)
    xyzs = numpy.dot(xyzs, numpy.transpose(mat))
    return automol.create.geom.from_data(symbs, xyzs)


def rotate(geo, axis, angle, orig_xyz=None, idxs=None):
    """ axis-angle rotation of the geometry

    if rotating a subset, specify the atoms indices with `idxs`
    """
    func = util.vec.rotater(axis, angle, orig_xyz=orig_xyz)
    return transform(geo, func, idxs=idxs)


def euler_rotate(geo, theta, phi, psi):
    """ axis-angle rotation of the geometry
    """
    mat = util.mat.euler_rotation(theta, phi, psi)
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


# geometric properties
def distance(geo, idx1, idx2, angstrom=False):
    """ measure the distance between atoms
    """
    xyzs = coordinates(geo)
    xyz1 = xyzs[idx1]
    xyz2 = xyzs[idx2]
    dist = util.vec.distance(xyz1, xyz2)
    dist *= BOHR2ANG if angstrom else 1
    return dist


def central_angle(geo, idx1, idx2, idx3, degree=False):
    """ measure the angle inscribed by three atoms
    """
    xyzs = coordinates(geo)
    xyz1 = xyzs[idx1]
    xyz2 = xyzs[idx2]
    xyz3 = xyzs[idx3]
    ang = util.vec.central_angle(xyz1, xyz2, xyz3)
    ang *= RAD2DEG if degree else 1
    return ang


def dihedral_angle(geo, idx1, idx2, idx3, idx4, degree=False):
    """ measure the dihedral angle defined by four atoms
    """
    xyzs = coordinates(geo)
    xyz1 = xyzs[idx1]
    xyz2 = xyzs[idx2]
    xyz3 = xyzs[idx3]
    xyz4 = xyzs[idx4]
    dih = util.vec.dihedral_angle(xyz1, xyz2, xyz3, xyz4)
    dih *= RAD2DEG if degree else 1
    return dih


def zmatrix_row_values(geo, idx, idx1=None, idx2=None, idx3=None,
                       angstrom=True, degree=True):
    """ get coordinate values for a row in a z-matrix
    """
    val_row = [None, None, None]

    if idx1 is not None:
        val_row[0] = distance(geo, idx, idx1, angstrom=angstrom)

    if idx2 is not None:
        assert idx1 is not None
        val_row[1] = central_angle(geo, idx, idx1, idx2, degree=degree)

    if idx3 is not None:
        assert idx2 is not None
        val_row[2] = dihedral_angle(geo, idx, idx1, idx2, idx3, degree=degree)

    val_row = tuple(val_row)
    return val_row


def linear_atoms(geo, gra=None, tol=5.):
    """ find linear atoms in a geometry (atoms with 180 degree bond angle)

    :param geo: the geometry
    :param gra: the graph describing connectivity; if None, a connectivity
        graph will be generated using default distance thresholds
    :param tol: the tolerance threshold for linearity, in degrees
    """
    gra = automol.convert.geom.connectivity_graph(geo) if gra is None else gra
    ngb_idxs_dct = automol.graph.atom_neighbor_keys(gra)

    lin_idxs = []
    for idx in range(count(geo)):
        nidxs = ngb_idxs_dct[idx]
        if len(nidxs) >= 2:
            for nidx1, nidx2 in itertools.combinations(nidxs, 2):
                ang = central_angle(geo, nidx1, idx, nidx2, degree=True)
                if numpy.abs(ang - 180.) < tol:
                    lin_idxs.append(idx)

    lin_idxs = tuple(lin_idxs)
    return lin_idxs


def distance_matrix(geo):
    """form distance matrix for a set of xyz coordinates
    """
    mat = numpy.zeros((len(geo), len(geo)))
    for i in range(len(geo)):
        for j in range(len(geo)):
            mat[i][j] = distance(geo, i, j)
    return mat


def closest_unbonded_atoms(geo, gra=None):
    """ determine the closest unbonded pair of atoms in this geometry
    """
    gra = automol.convert.geom.connectivity_graph(geo) if gra is None else gra
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
    """ does this geometry describe a single atom?
    """
    return len(symbols(geo)) == 1


def masses(geo, amu=True):
    """ return the atomic masses
    """
    symbs = symbols(geo)
    amas = list(map(pt.to_mass, symbs))

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


def is_linear(geo, tol=5.):
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
            ang = central_angle(geo, key1, key2, key3, degree=True)
            if not (numpy.abs(ang) < tol or
                    numpy.abs(ang - 180.) < tol):
                ret = False
                break
    return ret


def permutation(geo, ref_geo, thresh=1e-4):
    """ determine the permutation of one geometry that reproduces another

    (If there isn't one -- the geometries are not aligned, return None)
    """
    natms = count(geo)
    symbs = symbols(geo)
    xyzs = coordinates(geo)

    perm_idxs = [None] * natms
    for idx, (symb, xyz) in enumerate(zip(symbs, xyzs)):
        # Loop over atoms in the reference geometry with the same symbol
        ref_idxs = atom_indices(ref_geo, symb=symb)
        ref_xyzs = coordinates(ref_geo, idxs=ref_idxs)
        perm_idx = next(
            (ref_idx for ref_idx, ref_xyz in zip(ref_idxs, ref_xyzs)
             if util.vec.distance(xyz, ref_xyz) < thresh), None)
        perm_idxs[idx] = perm_idx

    perm_idxs = tuple(perm_idxs)

    if any(perm_idx is None for perm_idx in perm_idxs):
        perm_idxs = None

    return perm_idxs


if __name__ == '__main__':
    # ICH = automol.smiles.inchi('C=C=C')
    ICH = automol.smiles.inchi('C#CCCCC#C')
    # ICH = automol.smiles.inchi('C#C')
    GEO = automol.inchi.geometry(ICH)
    ZMA = automol.geom.zmatrix(GEO)
    GEO = automol.zmat.geometry(ZMA)
    GEO = automol.geom.without_dummy_atoms(GEO)
    LIN_IDXS = linear_atoms(GEO)
    print(LIN_IDXS)
    GEO = insert_dummies_on_linear_atoms(GEO)
    print(automol.geom.string(GEO))
