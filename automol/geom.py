""" cartesian geometries
"""
import itertools
import functools
import more_itertools as mit
import numpy
import qcelemental as qcel
from qcelemental import periodictable as pt
from qcelemental import constants as qcc
import autoread as ar
import autowrite as aw
from automol import cart
import automol.create.geom
import automol.convert.geom
import automol.convert.inchi


# constructor
def from_data(syms, xyzs, angstrom=False):
    """ geometry data structure from symbols and coordinates
    """
    return automol.create.geom.from_data(
        symbols=syms, coordinates=xyzs, angstrom=angstrom)


# getters
def symbols(geo):
    """ atomic symbols
    """
    if geo:
        syms, _ = zip(*geo)
    else:
        syms = ()
    return syms


def coordinates(geo, angstrom=False):
    """ atomic coordinates
    """
    if geo:
        _, xyzs = zip(*geo)
    else:
        xyzs = ()
    xyzs = xyzs if not angstrom else numpy.multiply(
        xyzs, qcc.conversion_factor('bohr', 'angstrom'))
    return xyzs


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
    xyzs = coordinates(geo)

    non_dummy_keys = [idx for idx, sym in enumerate(syms) if pt.to_Z(sym)]
    syms = list(map(syms.__getitem__, non_dummy_keys))
    xyzs = list(map(xyzs.__getitem__, non_dummy_keys))
    return from_data(syms, xyzs)


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

    # get the correct distance apart
    geo1 = mass_centered(geo1)
    geo2 = mass_centered(geo2)
    ext1 = max(numpy.vdot(orient_vec, xyz) for xyz in coordinates(geo1))
    ext2 = max(numpy.vdot(-orient_vec, xyz) for xyz in coordinates(geo2))

    cm_dist = ext1 + dist_cutoff + ext2
    dist_grid = numpy.arange(cm_dist, 0., -0.1)
    for dist in dist_grid:
        trans_geo2 = translated(geo2, orient_vec * dist)
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


def view(geo):
    """ view the geometry using nglview
    """
    import nglview

    xyz_str = xyz_string(geo)
    qcm = qcel.models.Molecule.from_data(xyz_str)
    ngv = nglview.show_qcelemental(qcm)
    return ngv


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
def displaced(geo, xyzs):
    """ displacement of the geometry
    """
    syms = symbols(geo)
    orig_xyzs = coordinates(geo)
    xyzs = numpy.add(orig_xyzs, xyzs)
    return from_data(syms, xyzs)


def translated(geo, xyz):
    """ translation of the geometry
    """
    syms = symbols(geo)
    xyzs = coordinates(geo)
    xyzs = numpy.add(xyzs, xyz)
    return from_data(syms, xyzs)


def rotated(geo, axis, angle):
    """ axis-angle rotation of the geometry
    """
    syms = symbols(geo)
    xyzs = coordinates(geo)
    rot_mat = cart.mat.rotation(axis, angle)
    xyzs = numpy.dot(xyzs, numpy.transpose(rot_mat))
    return from_data(syms, xyzs)


def euler_rotated(geo, theta, phi, psi):
    """ axis-angle rotation of the geometry
    """
    syms = symbols(geo)
    xyzs = coordinates(geo)
    rot_mat = cart.mat.euler_rotation(theta, phi, psi)
    xyzs = numpy.dot(xyzs, numpy.transpose(rot_mat))
    return from_data(syms, xyzs)


def swap_coordinates(geo, idx1, idx2):
    """ swap the order of the coordinates of the two atoms
    """

    # convert geo (tuple-of-tuples) to a list-of-lists
    geo2 = [list(x) for x in geo]

    # swap the two coordinates
    geo2[idx1], geo2[idx2] = geo2[idx2], geo2[idx1]

    # convert geo back to a tuple-of-tuples
    geo_swp = tuple(tuple(x) for x in geo2)

    return geo_swp


def move_coordinates(geo, idx1, idx2):
    """ move the atom at position idx1 to idx2, shifting all other atoms
    """

    # convert geo (tuple-of-tuples) to a list-of-lists
    geo2 = [list(x) for x in geo]

    # Get the coordinates at idx1 that are to be moved
    moving_coords = geo2[idx1]

    # move the coordinates to idx2
    geo2.remove(moving_coords)
    geo2.insert(idx2, moving_coords)

    # convert geo back to a tuple-of-tuples
    geo_move = tuple(tuple(x) for x in geo2)

    return geo_move


# geometric properties
def distance(geo, key1, key2):
    """ measure the distance between atoms
    """
    xyzs = coordinates(geo)
    xyz1 = xyzs[key1]
    xyz2 = xyzs[key2]
    return cart.vec.distance(xyz1, xyz2)


def central_angle(geo, key1, key2, key3):
    """ measure the angle inscribed by three atoms
    """
    xyzs = coordinates(geo)
    xyz1 = xyzs[key1]
    xyz2 = xyzs[key2]
    xyz3 = xyzs[key3]
    return cart.vec.central_angle(xyz1, xyz2, xyz3)


def dihedral_angle(geo, key1, key2, key3, key4):
    """ measure the dihedral angle defined by four atoms
    """
    xyzs = coordinates(geo)
    xyz1 = xyzs[key1]
    xyz2 = xyzs[key2]
    xyz3 = xyzs[key3]
    xyz4 = xyzs[key4]
    return cart.vec.dihedral_angle(xyz1, xyz2, xyz3, xyz4)


def dist_mat(geo):
    """form distance matrix for a set of xyz coordinates
    """
    mat = numpy.zeros((len(geo),len(geo)))
    for i in range(len(geo)):
        for j  in range(len(geo)):
            mat[i][j] = distance(geo, i, j)
    return mat

def almost_equal_dist_mat(geo1, geo2, thresh = 0.1):
    """form distance matrix for a set of xyz coordinates
    """
    dist_mat1 = dist_mat(geo1)
    dist_mat2 = dist_mat(geo2)
    diff_mat  = numpy.zeros((len(geo),len(geo)))
    almost_equal = True
    for i in range(len(dist_mat1)):
        for j in range(len(dist_mat1)):
            diff_mat[i][j] = abs(dist_mat1[i][j] - dist_mat2[i][j])
    if numpy.amax(diff_mat) > thresh:
        almost_equal = False
    return almost_equal


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


def center_of_mass(geo):
    """ center of mass
    """
    xyzs = coordinates(geo)
    amas = masses(geo)
    cm_xyz = tuple(
        sum(numpy.multiply(xyz, ama) for xyz, ama in zip(xyzs, amas)) /
        sum(amas))

    return cm_xyz


def mass_centered(geo):
    """ mass-centered geometry
    """
    geo = translated(geo, numpy.negative(center_of_mass(geo)))
    return geo


def inertia_tensor(geo, amu=True):
    """ molecular inertia tensor (atomic units if amu=False)
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
            if numpy.abs(central_angle(geo, key1, key2, key3) % numpy.pi) > tol:
                ret = False
    return ret


# conversions
def zmatrix(geo):
    """ geometry => z-matrix
    """
    return automol.convert.geom.zmatrix(geo)


def zmatrix_torsion_coordinate_names(geo):
    """ z-matrix torsional coordinate names
    """
    return automol.convert.geom.zmatrix_torsion_coordinate_names(geo)


def zmatrix_atom_ordering(geo):
    """ z-matrix atom ordering
    """
    return automol.convert.geom.zmatrix_atom_ordering(geo)


def graph(geo, remove_stereo=False):
    """ geometry => graph
    """
    return automol.convert.geom.graph(geo, remove_stereo=remove_stereo)


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
