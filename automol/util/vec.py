""" vector functions
"""
import numbers
import numpy
import transformations as tf


def unit_norm(xyz):
    """ Normalize a vector (xyz) to 1.0.

        :param xyz: vector
        :type xyz: tuple, list, or numpy nd.array
        :rtype: float
    """
    norm = numpy.linalg.norm(xyz)
    uxyz = numpy.divide(xyz, norm)
    assert numpy.allclose(numpy.linalg.norm(uxyz), 1.0)
    return uxyz


def unit_direction(xyz1, xyz2):
    """ Calculate a unit direction vector from `xyz1` to `xyz2`.

        :param xyz1: 3D vector
        :type xyz1: tuple, list, or numpy nd.array
        :param xyz2: 3D vector
        :type xyz2: tuple, list, or numpy nd.array
        :rtype: float
    """
    dxyz12 = numpy.subtract(xyz2, xyz1)
    uxyz12 = unit_norm(dxyz12)
    return uxyz12


def are_parallel(xyz1, xyz2, orig_xyz=(0., 0., 0.), tol=1e-7):
    """ Assess if two vectors are parallel to each other.

        :param xyz1: 3D vector
        :type xyz1: tuple, list, or numpy nd.array
        :param xyz2: 3D vector
        :type xyz2: tuple, list, or numpy nd.array
        :param orig_xyz: origin of coordinate system `xyz1` and `xyz2` are in
        :type orig_xyz: tuple, list, or numpy nd.array
        :param tol: tolerance for checking determinant
        :type tol: float
        :rtype: bool
    """

    det = numpy.linalg.det([list(orig_xyz), list(xyz1), list(xyz2)])
    return det > tol


def orthogonalize(xyz1, xyz2, normalize=False):
    """ orthogonalize `xyz2` against `xyz1`
    """
    overlap = numpy.dot(xyz1, xyz2)
    norm = numpy.dot(xyz1, xyz1)
    oxyz2 = numpy.subtract(xyz2, numpy.multiply(overlap/norm, xyz1))
    if normalize:
        oxyz2 = unit_norm(oxyz2)
    return oxyz2


def arbitrary_unit_perpendicular(xyz, orig_xyz=(0., 0., 0.)):
    """ determine an arbitrary perpendicular vector
    """
    for xyz2 in ([1., 0., 0.], [0., 1., 0.], [0., 0., 1]):
        uxyz = unit_perpendicular(xyz, xyz2, orig_xyz=orig_xyz)
        if numpy.linalg.norm(uxyz) > 1e-7:
            break

    return uxyz


def unit_perpendicular(xyz1, xyz2, orig_xyz=(0., 0., 0.), allow_parallel=True):
    """ Calculate a unit perpendicular on `xyz1` and `xyz2`.

        :param xyz1: 3D vector
        :type xyz1: tuple, list, or numpy nd.array
        :param xyz2: 3D vector
        :type xyz2: tuple, list, or numpy nd.array
        :param orig_xyz: origin of coordinate system `xyz1` and `xyz2` are in
        :type orig_xyz: tuple, list, or numpy nd.array
        :param allow_parallel: parameter to allow if vector can be parallel
        :type allow_parallel: bool
        :rtype: numpy.ndarray
    """

    xyz1 = numpy.subtract(xyz1, orig_xyz)
    xyz2 = numpy.subtract(xyz2, orig_xyz)
    xyz3 = numpy.cross(xyz1, xyz2)

    if numpy.linalg.norm(xyz3) > 1e-7:
        uxyz3 = unit_norm(xyz3)
    elif allow_parallel:
        uxyz3 = numpy.zeros((3,))
    else:
        raise ValueError

    return uxyz3


def unit_bisector(xyz1, xyz2, orig_xyz):
    """ Calculate a unit bisector.

        :param xyz1: 3D vector
        :type xyz1: tuple, list, or numpy nd.array
        :param xyz2: 3D vector
        :type xyz2: tuple, list, or numpy nd.array
        :param orig_xyz: origin of coordinate system `xyz1` and `xyz2` are in
        :type orig_xyz: tuple, list, or numpy nd.array
        :rtype: numpy.ndarray
    """

    ang = central_angle(xyz1, orig_xyz, xyz2)
    rot_ = rotater(
        axis=unit_perpendicular(xyz1, xyz2, orig_xyz), angle=ang/2.0,
        orig_xyz=orig_xyz)
    xyz = unit_norm(numpy.subtract(rot_(xyz1), orig_xyz))
    return xyz


def from_internals(dist=0., xyz1=(0., 0., 0.), ang=0., xyz2=(0., 0., 1.),
                   dih=0., xyz3=(0., 1., 0.)):
    """ Determine the position of a point (xyz4) in Cartesian coordinates whose
        position is related to three other points (xyz1, xyz2, xyz3) via a set
        of internal coordinates.

        :param dist: distance between `xyz1` and `xyz2` (in bohr)
        :type dist: float
        :param xyz1: 3D vector to point 1
        :type xyz1: tuple, list, or numpy nd.array
        :param angle: angle between `xyz1`, `xyz2`, `xyz3` (in radians)
        :type angle: float
        :param xyz2: 3D vector to point 2
        :type xyz2: tuple, list, or numpy nd.array
        :param dih: dihedral from `xyz1`, `xyz2`, `xyz3` to `xyz4` (in radians)
        :type dih: float
        :param xyz3: 3D vector to point 2
        :type xyz3: tuple, list, or numpy nd.array
        :rtyp: tuple(float)
    """

    local_xyz = _local_position(dist=dist, ang=ang, dih=dih)
    local_basis = _local_axes(xyz1=xyz1, xyz2=xyz2, xyz3=xyz3)
    xyz4 = tuple(xyz1 + numpy.dot(local_xyz, local_basis))

    return xyz4


def _local_position(dist=0., ang=0., dih=0.):
    """ Determine the xyz coordinates of a point in the local axis frame
        defined by a set of internal coordinates.

        :param dist: distance between `xyz1` and `xyz2` (in bohr)
        :type dist: float
        :param angle: angle between `xyz1`, `xyz2`, `xyz3` (in radians)
        :type angle: float
        :param dih: dihedral from `xyz1`, `xyz2`, `xyz3` to `xyz4` (in radians)
        :type dih: float
        :rtyp: tuple(float)
    """

    x_comp = dist * numpy.sin(ang) * numpy.sin(dih)
    y_comp = dist * numpy.sin(ang) * numpy.cos(dih)
    z_comp = dist * numpy.cos(ang)

    return (x_comp, y_comp, z_comp)


def _local_axes(xyz1=(0., 0., 0.), xyz2=(0., 0., 1.), xyz3=(0., 1., 0.)):
    """ Determine the  local axes for defining bond, angle, dihedral from
        the Cartesian coordinates of three support atoms.

        :param xyz1: 3D vector to point 1
        :type xyz1: tuple, list, or numpy nd.array
        :param xyz2: 3D vector to point 2
        :type xyz2: tuple, list, or numpy nd.array
        :param xyz3: 3D vector to point 2
        :type xyz3: tuple, list, or numpy nd.array
        :rtyp: tuple(float)
    """

    uxyz12 = unit_direction(xyz1, xyz2)
    uxyz23 = unit_direction(xyz2, xyz3)
    uxyz123_perp = unit_perpendicular(uxyz23, uxyz12)

    z_ax = tuple(uxyz12)
    y_ax = tuple(unit_perpendicular(uxyz12, uxyz123_perp))
    x_ax = tuple(unit_perpendicular(y_ax, z_ax))
    return (x_ax, y_ax, z_ax)


def distance(xyz1, xyz2):
    """ Measure the distance between points.

        :param xyz1: 3D vector to point 1
        :type xyz1: tuple, list, or numpy nd.array
        :param xyz2: 3D vector to point 2
        :type xyz2: tuple, list, or numpy nd.array
        :rtype: float
    """

    dist = numpy.linalg.norm(numpy.subtract(xyz1, xyz2))

    return dist


def central_angle(xyz1, xyz2, xyz3):
    """ Measure the angle inscribed by three atoms.

        :param xyz1: 3D vector to point 1
        :type xyz1: tuple, list, or numpy nd.array
        :param xyz2: 3D vector to point 2
        :type xyz2: tuple, list, or numpy nd.array
        :param xyz3: 3D vector to point 3
        :type xyz3: tuple, list, or numpy nd.array
        :rtype: float
    """
    uxyz21 = unit_direction(xyz2, xyz1)
    uxyz23 = unit_direction(xyz2, xyz3)
    cos = numpy.dot(uxyz21, uxyz23)
    if cos < -1.0:
        assert numpy.allclose(cos, -1.0)
        cos = -1.0
    elif cos > 1.0:
        assert numpy.allclose(cos, 1.0)
        cos = 1.0
    ang = numpy.arccos(cos)
    return ang


def projected_central_angle(xyz1, xyz2, xyz3):
    """ Measure the angle inscribed by three atoms,
        projected onto the normal plane of the central atom.

        :param xyz1: 3D vector to point 1
        :type xyz1: tuple, list, or numpy nd.array
        :param xyz2: 3D vector to point 2
        :type xyz2: tuple, list, or numpy nd.array
        :param xyz3: 3D vector to point 3
        :type xyz3: tuple, list, or numpy nd.array
        :rtype: float
    """

    uxyz21 = unit_perpendicular(xyz2, xyz1)
    uxyz23 = unit_perpendicular(xyz2, xyz3)
    cos = numpy.dot(uxyz21, uxyz23)
    if cos < -1.0:
        assert numpy.allclose(cos, -1.0)
        cos = -1.0
    elif cos > 1.0:
        assert numpy.allclose(cos, 1.0)
        cos = 1.0
    ang = numpy.arccos(cos)
    return ang


def dihedral_angle(xyz1, xyz2, xyz3, xyz4):
    """ Measure the dihedral angle defined by four atoms.

        :param xyz1: 3D vector to point 1
        :type xyz1: tuple, list, or numpy nd.array
        :param xyz2: 3D vector to point 2
        :type xyz2: tuple, list, or numpy nd.array
        :param xyz3: 3D vector to point 3
        :type xyz3: tuple, list, or numpy nd.array
        :param xyz4: 3D vector to point 4
        :type xyz4: tuple, list, or numpy nd.array
        :rtype: float
    """

    # Get the cosine of the angle
    uxyz21 = unit_direction(xyz2, xyz1)
    uxyz23 = unit_direction(xyz2, xyz3)
    uxyz32 = unit_direction(xyz3, xyz2)
    uxyz34 = unit_direction(xyz3, xyz4)
    uxyz123_perp = unit_perpendicular(uxyz21, uxyz23)
    uxyz234_perp = unit_perpendicular(uxyz32, uxyz34)
    cos = numpy.dot(uxyz123_perp, uxyz234_perp)

    # Get the sign of the angle
    val = numpy.dot(uxyz123_perp, uxyz34)
    val = max(min(val, 1.), -1.)
    sign = 2 * (val < 0) - 1

    # Before plugging it into the arccos function, make sure we haven't
    # Slightly run out of bounds
    if cos < -1.0:
        assert numpy.allclose(cos, -1.0)
        cos = -1.0
    elif cos > 1.0:
        assert numpy.allclose(cos, 1.0)
        cos = 1.0

    dih = sign * numpy.arccos(cos)
    dih = numpy.mod(dih, 2*numpy.pi)
    return dih


# transformations
def rotater(axis, angle, orig_xyz=None):
    """ A function to rotate vectors about an axis at a particular point.

        :param axis: axis to rotate about
        :type axis: list, tuple, or nd.array
        :param angle: angle by which to rotate the vectors
        :type angle: float
        :param orig_xyz: origin of system
        :type: list, tuple, or nd.array
        :rtype: tuple(tuple(float))
    """
    aug_rot_mat = tf.rotation_matrix(angle, axis, point=orig_xyz)
    return _transformer(aug_rot_mat)


def aligner(xyz1a, xyz1b, xyz2a, xyz2b):
    """ a function to translate and rotate a system, bringing two points into
    alignment

    Takes 1a-1b into 2a-2b if they are equidistant; otherwise, takes 1a into 2a
    and takes 1b onto the 2a-2b line
    """
    xyz1a, xyz1b, xyz2a, xyz2b = map(numpy.array, (xyz1a, xyz1b, xyz2a, xyz2b))

    trans = xyz2a - xyz1a

    trans_mat = tf.translation_matrix(trans)

    xyz1b = xyz1b + trans

    rot_vec = unit_perpendicular(xyz1b, xyz2b, orig_xyz=xyz2a)
    rot_ang = central_angle(xyz1b, xyz2a, xyz2b)
    rot_mat = tf.rotation_matrix(rot_ang, rot_vec, point=xyz2a)

    mat = numpy.dot(rot_mat, trans_mat)

    return _transformer(mat)


def _transformer(aug_mat):
    def _transform(xyz):
        aug_xyz = _augmented(xyz)
        rot_aug_xyz = numpy.dot(aug_mat, aug_xyz)
        rot_xyz = tuple(rot_aug_xyz[:3])
        return rot_xyz
    return _transform


def _augmented(xyz):
    assert len(xyz) == 3
    assert all(isinstance(val, numbers.Real) for val in xyz)
    xyz_aug = tuple(xyz) + (1.,)
    return xyz_aug


# I/O
def string(vec, num_per_row=None, val_format='{0:>8.3f}'):
    """ Write a vector to a string.

        :param vec: vector to form string with
        :type vec: list, tuple, or nd.array
        :param num_per_row: number of vector elements to write to a row
        :type num_per_row: int
        :rtype: str
    """

    if num_per_row is None:
        num_per_row = len(vec)

    assert isinstance(num_per_row, int), 'num_per_row must be an integer'

    nvals = len(vec)
    vec_str = ''
    for i, val in enumerate(vec):
        vec_str += val_format.format(val)
        if ((i+1) % num_per_row) == 0 and (i+1) != nvals:
            vec_str += '\n'

    return vec_str
