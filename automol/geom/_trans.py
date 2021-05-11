"""
  Performs several operations and manipulations of geometries
"""

import numpy
from phydat import phycon
import automol.create.geom
from automol import util
from automol.geom import _base as geom_base
from automol.graph.geom import center_of_mass
from automol.graph.geom import translate as _translate
from automol.graph.geom import geometry_join as _geometry_join

AXIS_DCT = {'x': 0, 'y': 1, 'z': 2}


# General transformation functions
def transform(geo, func, idxs=None):
    """ Transform the coordinates of a geometry by a function.
        A set of `idxs` can be supplied to transform a subset of coordinates.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param func: transformation function
        :type func: function object
        :param idxs: indices representing the subset of atoms
        :type idxs: tuple(int)
    """

    idxs = list(range(geom_base.count(geo))) if idxs is None else idxs
    symbs = geom_base.symbols(geo)
    xyzs = geom_base.coordinates(geo)
    xyzs = [func(xyz) if idx in idxs else xyz for idx, xyz in enumerate(xyzs)]

    return automol.create.geom.from_data(symbs, xyzs)


def transform_by_matrix(geo, mat):
    """ Transform the coordinates of a molecular geometry by multiplying
        it by some input transfomration matrix.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param mat: transformation matrix
        :type mat: tuple(tuple(float))
        :rtype: automol moleculer geometry data structure
    """

    symbs = geom_base.symbols(geo)
    xyzs = geom_base.coordinates(geo)
    xyzs = numpy.dot(xyzs, numpy.transpose(mat))

    return automol.create.geom.from_data(symbs, xyzs)


# transformations
def remove_coordinates(geo, idxs=()):
    """ Remove atoms by from the molecular geometry by their indices.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param idxs: indices of atoms to remove
        :type idxs: tuple(int)
        :rtype: automol molecular geometry data structure
    """
    return tuple(row for i, row in enumerate(geo) if i not in idxs)


def join(geo1, geo2,
         dist_cutoff=3.0*phycon.ANG2BOHR, theta=0.0, phi=0.0):
    """ Join two molecular geometries together where the intermolecular
        separation and orientation can be specified.

        :param geo1: molecular geometry 1
        :type geo1: automol molecular geometry data structure
        :param geo2: molecular geometry 2
        :type geo2: automol molecular geometry data structure
        :param dist_cutoff: threshhold for center-of-mass distance
        :type: dist_cutoff: float
        :param theta: theta angle for intermolecular orientation
        :type theta: float
        :param phi: phi angle for intermolecular orientation
        :type phi: float
        :rtype: automol molecular geometry data structure
    """
    return _geometry_join(
        geo1, geo2, dist_cutoff=dist_cutoff, theta=theta, phi=phi)


def reorder_coordinates(geo, idx_dct):
    """ Reorder the atoms of a molecular geometry using
        the mapping of an input dictionary.

        :param geo: The geometry
        :param idx_dct: The new order of the atoms, by index
        :type idx_dct: dict
        :rtype: automol geometry data structure
    """

    symbs = geom_base.symbols(geo)
    xyzs = geom_base.coordinates(geo)

    idxs = [idx for idx, _ in sorted(idx_dct.items(), key=lambda x: x[1])]
    assert len(symbs) == len(xyzs) == len(idxs)

    symbs = [symbs[idx] for idx in idxs]
    xyzs = [xyzs[idx] for idx in idxs]

    return automol.create.geom.from_data(symbs, xyzs)


def swap_coordinates(geo, idx1, idx2):
    """ Swap the order of the coordinates of two atoms in a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param idx1: index for one atom to swap coordinates
        :type idx1: int
        :param idx2: index for one atom to swap coordinates
        :type idx2: int
        :rtype: molecular geometry
    """

    geo = [list(x) for x in geo]
    geo[idx1], geo[idx2] = geo[idx2], geo[idx1]
    geo_swp = tuple(tuple(x) for x in geo)

    return geo_swp


def insert(geo, symb, xyz, idx=None, angstrom=False):
    """ Insert an atom into a molecular geometry.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param symb: symbol of atom to add
        :type symb: str
        :param xyz: xyz coordinates of atom to add
        :type xyz: tuple(float)
        :param idx: index of geometry to place atom
        :type idx: int
        :rtype: automol geometry date structure
    """

    return automol.convert.geom.insert(
        geo, symb, xyz, idx=idx, angstrom=angstrom)


def insert_dummies_on_linear_atoms(geo, lin_idxs=None, gra=None, dist=1.,
                                   tol=5.):
    """ Insert dummy atoms over linear atoms in the geometry.

        :param geo: the geometry
        :type geo: automol molecular geometry data structure
        :param lin_idxs: the indices of the linear atoms; if None, indices are
            automatically determined from the geometry based on the graph
        :type lin_idxs: tuple(int)
        :param gra: the graph describing connectivity; if None, a connectivity
            graph will be generated using default distance thresholds
        :type gra: automol molecular graph data structure
        :param dist: distance of dummy atom from the linear atom, in angstroms
        :type dist: float
        :param tol: the tolerance threshold for linearity, in degrees
        :type tol: float
        :returns: geometry with dummy atoms inserted, along with a dictionary
            mapping the linear atoms onto their associated dummy atoms
        :rtype: automol molecular geometry data structure
    """
    return automol.convert.geom.insert_dummies_on_linear_atoms(
        geo, lin_idxs=lin_idxs, gra=gra, dist=dist, tol=tol)


def displace(geo, xyzs):
    """ Displace the coordinates of a geometry along a vector.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param xyzs: vector to displace along
        :type xyzs: tuple(float)
        :rtype: automol molecular geometry data structure
    """

    symbs = geom_base.symbols(geo)
    orig_xyzs = geom_base.coordinates(geo)
    xyzs = numpy.add(orig_xyzs, xyzs)

    return automol.create.geom.from_data(symbs, xyzs)


# redundant with above
def translate(geo, xyz):
    """ Translate the coordinates of a molecular geometry along
        a three-dimensiona vector.

        :param geo: molecular geometry
        :type geo: automol molecular geometry data structure
        :param xyz: vector to translate along
        :type xyz: tuple(float)
        :rtype: automol molecular geometry data structure
    """
    return _translate(geo, xyz)


def perturb(geo, atm_idx, pert_xyz):
    """ Perturb the position of one atom by
        changing the value of an xyz coord by some amount.
    """

    # Get the xyz coordinates of the atom to perturb
    atm_coords = list(geom_base.coordinates(geo)[atm_idx])

    # Get the perturbed set of atomic coordinates
    for idx, val in enumerate(pert_xyz):
        atm_coords[idx] += val
    pert_dct = {atm_idx: atm_coords}

    # Perturb the coordinates of the atom
    pert_geo = geom_base.set_coordinates(geo, pert_dct)

    return pert_geo


def mass_centered(geo):
    """ Generate a new geometry where the coordinates of the input geometry
        have been translated to the center-of-mass.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :rtype: tuple(float)
    """
    return translate(geo, numpy.negative(center_of_mass(geo)))


def rotate(geo, axis, angle, orig_xyz=None, idxs=None):
    """ Rotate the coordinates of a molecular geometry about
        an axis by a specified angle. A set of `idxs` can be supplied
        to transform a subset of coordinates.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param axis: axis to rotate about
        :type axis: tuple(float)
        :param angle: angle of rotation
        :type angle: float
        :param orig_xyz: xyz coordinates of the origin
        :type orig_xyz: tuple(float)
        :param idxs: indices of atoms whose coordinates are to be rotated
        :type idxs: tuple(int)
    """

    func = util.vec.rotater(axis, angle, orig_xyz=orig_xyz)

    return transform(geo, func, idxs=idxs)


def euler_rotate(geo, theta, phi, psi):
    """ Rotate the coordinates of a molecular geometry about
        the three Euler angles.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param theta: angle to rotate about z-axis
        :type theta: float
        :param phi: angle to rotate about x'-axis
        :type phi: float
        :param psi: angle to rotate about z'-axis
        :type psi: float
        :rtype: automol geometry data structure
    """

    mat = util.mat.euler_rotation_matrix(theta, phi, psi)

    return transform_by_matrix(geo, mat)


def shift_atom_position(geo, idx1, idx2):
    """ Move the atom at position idx1 to idx2, shifting all other atoms.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param idx1: index of atom 1 in the pair to be measured
        :type idx1: int
        :param idx2: index of atom 2 in the pair to be measured
        :type idx2: int
        :rtype: automol geometry data structure
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
    """ Reflect a specified set of coordinates of a molecular geometry
        about some each of the requested axes.

        A set of `idxs` can be supplied to transform a subset of coordinates.

        :param geo: molecular geometry
        :type geo: automol geometry data structure
        :param idxs: indices of atoms whose coordinates are to be reflected
        :type idxs: tuple(int)
        :param axes: axes to reflect about
        :type axes: tuple(str)
        :rtype: automol geometry data structure
    """

    # check input
    assert all(idx < len(geo) for idx in idxs)
    assert all(axis in ('x', 'y', 'z') for axis in axes)

    # get coords
    coords = geom_base.coordinates(geo)

    # convert x,y,z to nums
    axes = [AXIS_DCT[axis] for axis in axes]

    # build set atom dct with relected coords
    reflect_dct = {}
    for idx in idxs:
        coord_lst = list(coords[idx])
        for axis in axes:
            coord_lst[axis] *= -1.0
        reflect_dct[idx] = coord_lst

    # Reflect coords with dct
    geo_reflected = geom_base.set_coordinates(geo, reflect_dct)

    return geo_reflected
