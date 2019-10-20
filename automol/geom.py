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

    geo = [list(x) for x in geo]
    geo[idx1], geo[idx2] = geo[idx2], geo[idx1]
    geo_swp = tuple(tuple(x) for x in geo)
    return geo_swp


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



def rot_permutated_geoms(geo, saddle=False, frm_bnd_key=[], brk_bnd_key=[], form_coords=[]):
    """ convert an input geometry to a list of geometries
        corresponding to the rotational permuations of all the terminal groups
    """
    gra = graph(geo, remove_stereo=True)
    term_atms = {}
    all_hyds = []
    neighbor_dct = automol.graph.atom_neighbor_keys(gra)

    # determine if atom is a part of a double bond
    unsat_atms = automol.graph.unsaturated_atom_keys(gra)
    if not saddle:
        rad_atms = automol.graph.sing_res_dom_radical_atom_keys(gra)
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
            if atm not in frm_bnd_key and atm not in brk_bnd_key:
            #if atm not in form_coords:
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
    mat = numpy.zeros((len(geo), len(geo)))
    for i in range(len(geo)):
        for j in range(len(geo)):
            mat[i][j] = distance(geo, i, j)
    return mat


def almost_equal_dist_mat(geo1, geo2, thresh=0.1):
    """form distance matrix for a set of xyz coordinates
    """
    dist_mat1 = dist_mat(geo1)
    dist_mat2 = dist_mat(geo2)
    diff_mat = numpy.zeros((len(geo1), len(geo2)))
    almost_equal_dm = True
    for i, _ in enumerate(dist_mat1):
        for j, _ in enumerate(dist_mat1):
            diff_mat[i][j] = abs(dist_mat1[i][j] - dist_mat2[i][j])
    if numpy.amax(diff_mat) > thresh:
        almost_equal_dm = False
    return almost_equal_dm


def external_symmetry_factor(geo):
    """ obtain external symmetry number for a geometry using x2z
    """
    # Get initial external symmetry number
    if automol.geom.is_atom(geo):
        ext_sym_fac = 1.
    else:
        oriented_geom = to_oriented_geometry(geo)
        ext_sym_fac = oriented_geom.sym_num()
        # Change symmetry number if geometry has enantiomers
        if oriented_geom.is_enantiomer():
            ext_sym_fac *= 0.5
    return ext_sym_fac


def find_xyzp_using_internals(xyz1, xyz2, xyz3, pdist, pangle, pdihed):
    """ geometric approach for calculating the xyz coordinates of atom A
        when the xyz coordinates of the A B and C are known and
        the position is defined w/r to A B C with internal coordinates
    """
    # Set to numpy arrays
    xyz1 = numpy.array(xyz1)
    xyz2 = numpy.array(xyz2)
    xyz3 = numpy.array(xyz3)

    print('coords')
    print(pdist)
    print(pangle)
    print(pdihed)

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

    print('dists')
    print(dist12)
    print(dist13)
    print(dist23)

    val = ((dist12**2 + dist13**2 - dist23**2) / 2.0 / dist12)
    valx3 = numpy.sqrt(dist13**2 - val**2)
    valy3 = ((dist12**2 + dist13**2 - dist23**2) / 2.0 / dist12)
    xyz3_rt = numpy.array([valx3, valy3, 0.0])

    print('pointd rt')
    print(xyz3_rt)

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

    print('r11 test')
    print(xyz3[0])
    print(xyz1[0])
    print(xyz3_rt[1]*r12)
    print(xyz3_rt[0])

    anum_aconst = xyz2_t[1] - (xyz3_t[1] / xyz3_t[0]) * xyz2_t[0]
    den_aconst = xyz2_t[2] - (xyz3_t[2] / xyz3_t[0]) * xyz2_t[0]

    print('rvals')
    print(r12)
    print(r22)
    print(r32)
    print(r11)
    print(r21)
    print(r31)
    print(anum_aconst)
    print(den_aconst)

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

    print('aconst')
    print(aconst)

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

    print('rn3 vals')
    print(r13)
    print(r23)
    print(r33)

    print('xyzp vals')
    print(xyzp_rt[0])
    print(xyzp_rt[1])
    print(xyzp_rt[2])

    # Now rotate and translate back
    # Here I check  the (001) vector direction to decide whether
    # To take the positive of negative results of the
    # Square root taken above
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

    print('xp')
    print(xap)
    print(yap)
    print(zap)
    print('xn')
    print(xan)
    print(yan)
    print(zan)

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
            if cangle % numpy.pi > tol:
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
    return automol.convert.geom.graph(
        geo, remove_stereo=remove_stereo)


def weakly_connected_graph(geo, remove_stereo=False):
    """ geometry => graph
    """
    return automol.convert.geom.weakly_connected_graph(
        geo, remove_stereo=remove_stereo)


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
