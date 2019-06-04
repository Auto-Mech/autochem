""" cartesian geometries
"""
import functools
import numpy
import more_itertools as mit
from qcelemental import periodictable as pt
from qcelemental import constants as qcc
import autoread as ar
import autowrite as aw
from automol import cart
import automol.create.geom
import automol.convert.geom


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


def without_ghost_atoms(geo):
    """ return a copy of the geometry without ghost atoms
    """
    syms = symbols(geo)
    xyzs = coordinates(geo)

    non_ghost_keys = [idx for idx, sym in enumerate(syms) if pt.to_Z(sym)]
    syms = list(map(syms.__getitem__, non_ghost_keys))
    xyzs = list(map(xyzs.__getitem__, non_ghost_keys))
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


# representations
def formula(geo):
    """ molecular formula, as a dictionary
    """
    syms = list(map(pt.to_E, symbols(geo)))
    frm_dct = {sym: syms.count(sym) for sym in mit.unique_everseen(syms)}
    return frm_dct


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
def almost_equal(geo1, geo2, rtol=2e-5):
    """ are these geometries numerically equal?
    """
    ret = False
    if symbols(geo1) == symbols(geo2):
        ret = numpy.allclose(coordinates(geo1), coordinates(geo2), rtol=rtol)
    return ret


def almost_equal_coulomb_spectrum(geo1, geo2, rtol=2e-5):
    """ do these geometries have similar coulomb spectrums?
    """
    ret = numpy.allclose(coulomb_spectrum(geo1), coulomb_spectrum(geo2),
                         rtol=rtol)
    return ret


def argunique_coulomb_spectrum(geos, seen_geos=(), rtol=2e-5):
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


# chemical properties
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


# conversions
def zmatrix(geo):
    """ geometry => z-matrix
    """
    return automol.convert.geom.zmatrix(geo)


def zmatrix_torsion_coordinate_names(geo):
    """ z-matrix torsional coordinate names
    """
    return automol.convert.geom.zmatrix_torsion_coordinate_names(geo)


def graph(geo):
    """ geometry => graph
    """
    return automol.convert.geom.graph(geo)


def inchi(geo, remove_stereo=False):
    """ geometry => inchi
    """
    return automol.convert.geom.inchi(geo, remove_stereo=remove_stereo)
