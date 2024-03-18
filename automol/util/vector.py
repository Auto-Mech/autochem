""" vector functions
"""
from typing import Callable, Optional, Tuple

import numpy
from scipy.spatial.transform import Rotation
from phydat import phycon

Vector = Tuple[float, float, float]


def unit_norm(xyz):
    """Normalize a vector (xyz) to 1.0.

    Returns 0. for null vectors

    :param xyz: vector
    :type xyz: tuple, list, or numpy nd.array
    :rtype: float
    """
    uxyz = xyz
    norm = numpy.linalg.norm(xyz)
    if not numpy.isclose(norm, 0.):
        uxyz = numpy.divide(xyz, norm)
    return uxyz


def unit_direction(xyz1, xyz2):
    """Calculate a unit direction vector from `xyz1` to `xyz2`.

    :param xyz1: 3D vector
    :type xyz1: tuple, list, or numpy nd.array
    :param xyz2: 3D vector
    :type xyz2: tuple, list, or numpy nd.array
    :rtype: float
    """
    dxyz12 = numpy.subtract(xyz2, xyz1)
    uxyz12 = unit_norm(dxyz12)
    return uxyz12


def are_parallel(xyz1, xyz2, orig_xyz=(0.0, 0.0, 0.0), anti=False, tol=1e-3):
    """Assess if two vectors are parallel to each other.

    :param xyz1: 3D vector
    :type xyz1: tuple, list, or numpy nd.array
    :param xyz2: 3D vector
    :type xyz2: tuple, list, or numpy nd.array
    :param orig_xyz: origin of coordinate system `xyz1` and `xyz2` are in
    :type orig_xyz: tuple, list, or numpy nd.array
    :param anti: Test for anti-parallelity too?
    :type anti: bool
    :param tol: tolerance for checking determinant
    :type tol: float
    :rtype: bool
    """
    vec1 = numpy.subtract(xyz1, orig_xyz)
    vec2 = numpy.subtract(xyz2, orig_xyz)
    len1 = numpy.linalg.norm(vec1)
    len2 = numpy.linalg.norm(vec2)

    ratio = numpy.dot(vec1, vec2) / (len1 * len2)

    ret = abs(ratio - 1) < tol
    if anti:
        ret |= abs(ratio + 1) < tol

    return ret


def orthogonalize(xyz1, xyz2, normalize=False):
    """orthogonalize `xyz2` against `xyz1`"""
    overlap = numpy.dot(xyz1, xyz2)
    norm = numpy.dot(xyz1, xyz1)
    oxyz2 = numpy.subtract(xyz2, numpy.multiply(overlap / norm, xyz1))
    if normalize:
        oxyz2 = unit_norm(oxyz2)
    return oxyz2


def flip_if_left_handed(xvec, yvec, zvec) -> Vector:
    """Given three vectors, flips the third one, if necessary, to make a right-handed
    coordinate system

    The vectors need not be orthogonal, but they must be linearly independent for this
    to work.  Does not change direction -- only flips the sign, if necessary.

    :param xvec: A vector defining the x direction
    :type xvec: Vector
    :param yvec: A vector defining the y direction
    :type yvec: Vector
    :param zvec: A vector defining the z direction
    :type zvec: Vector
    :returns: `zvec` or, if they made a left-handed system, `-zvec`
    :rtype: Vector
    """
    zdir = numpy.cross(xvec, yvec)
    proj = numpy.dot(zvec, zdir)
    if proj < 0.0:
        zvec = numpy.negative(zvec)
    return tuple(map(float, zvec))


def best_unit_perpendicular(xyzs):
    """Find a vector that is perpendicular to a series of points as much as possible

    For 1 point, this is an arbitrary vector.
    For 2 points, this is an arbitrary perpendicular vector.
    For 3 points, this is the exact normal to the plane.
    For 4+ points, this is the normal to the best-fit plane.

    Solution to the last problem was found here: https://stackoverflow.com/a/51132260

    :param xyzs: The points
    :type xyzs: List[Vector]
    """
    if len(xyzs) <= 1:
        nvec = [1, 0, 0]
    elif len(xyzs) == 2:
        nvec = arbitrary_unit_perpendicular(*xyzs)
    elif len(xyzs) == 3:
        nvec = unit_perpendicular(*xyzs)
    else:
        xyz0 = numpy.sum(xyzs, axis=0) / len(xyzs)
        _, _, vmat = numpy.linalg.svd(numpy.subtract(xyzs, xyz0))
        nvec = vmat[2, :]

    nvec = tuple(map(float, nvec))
    return nvec


def arbitrary_unit_perpendicular(xyz, orig_xyz=(0.0, 0.0, 0.0)):
    """determine an arbitrary perpendicular vector"""
    for xyz2 in ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1]):
        uxyz = unit_perpendicular(xyz, xyz2, orig_xyz=orig_xyz)
        if numpy.linalg.norm(uxyz) > 1e-7:
            break

    return uxyz


def unit_perpendicular(
    xyz1, xyz2, orig_xyz=(0.0, 0.0, 0.0), value_if_parallel: Vector = (0.0, 0.0, 0.0)
) -> Vector:
    """Calculate a unit perpendicular on `xyz1` and `xyz2`.

    :param xyz1: 3D vector
    :type xyz1: tuple, list, or numpy nd.array
    :param xyz2: 3D vector
    :type xyz2: tuple, list, or numpy nd.array
    :param orig_xyz: origin of coordinate system `xyz1` and `xyz2` are in
    :type orig_xyz: tuple, list, or numpy nd.array
    :param value_if_parallel: What to return if the vectors are parallel
    :type value_if_parallel: Vector
    :rtype: numpy.ndarray
    """

    xyz1 = numpy.subtract(xyz1, orig_xyz)
    xyz2 = numpy.subtract(xyz2, orig_xyz)
    xyz3 = numpy.cross(xyz1, xyz2)

    if numpy.linalg.norm(xyz3) > 1e-7:
        uxyz3 = unit_norm(xyz3)
    else:
        uxyz3 = value_if_parallel

    return uxyz3


def unit_bisector(xyz1, xyz2, orig_xyz, outer: bool = False) -> Vector:
    """Calculate a unit bisector.

    :param xyz1: 3D vector
    :type xyz1: tuple, list, or numpy nd.array
    :param xyz2: 3D vector
    :type xyz2: tuple, list, or numpy nd.array
    :param orig_xyz: origin of coordinate system `xyz1` and `xyz2` are in
    :type orig_xyz: tuple, list, or numpy nd.array
    :param outer: Get the outer, instead of the inner bisector?, defaults to False
    :type outer: bool, optional
    :rtype: numpy.ndarray
    """

    ang = central_angle(xyz1, orig_xyz, xyz2)
    rot_ = rotator(
        axis=unit_perpendicular(xyz1, xyz2, orig_xyz),
        ang=ang / 2.0,
        orig_xyz=orig_xyz,
    )
    xyz = unit_norm(numpy.subtract(rot_(xyz1), orig_xyz))
    if outer:
        xyz = numpy.negative(xyz)
    return xyz


def from_internals(
    dist=0.0,
    xyz1=(0.0, 0.0, 0.0),
    ang=0.0,
    xyz2=(0.0, 0.0, 1.0),
    dih=0.0,
    xyz3=(0.0, 1.0, 0.0),
):
    """Determine the position of a point (xyz4) in Cartesian coordinates whose
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


def _local_position(dist=0.0, ang=0.0, dih=0.0):
    """Determine the xyz coordinates of a point in the local axis frame
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


def _local_axes(xyz1=(0.0, 0.0, 0.0), xyz2=(0.0, 0.0, 1.0), xyz3=(0.0, 1.0, 0.0)):
    """Determine the  local axes for defining bond, angle, dihedral from
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
    """Measure the distance between points.

    :param xyz1: 3D vector to point 1
    :type xyz1: tuple, list, or numpy nd.array
    :param xyz2: 3D vector to point 2
    :type xyz2: tuple, list, or numpy nd.array
    :rtype: float
    """

    dist = numpy.linalg.norm(numpy.subtract(xyz1, xyz2))

    return dist


def angle(xyz1, xyz2, orig_xyz=(0.0, 0.0, 0.0)):
    """Measure the angle inscribed by three atoms.

    :param xyz1: 3D vector to point 1
    :type xyz1: tuple, list, or numpy nd.array
    :param xyz2: 3D vector to point 2
    :type xyz2: tuple, list, or numpy nd.array
    :param xyz3: 3D vector to point 3
    :type xyz3: tuple, list, or numpy nd.array
    :rtype: float
    """
    uxyz1 = unit_direction(orig_xyz, xyz1)
    uxyz2 = unit_direction(orig_xyz, xyz2)
    cos = numpy.dot(uxyz1, uxyz2)
    if cos < -1.0:
        assert numpy.allclose(cos, -1.0)
        cos = -1.0
    elif cos > 1.0:
        assert numpy.allclose(cos, 1.0)
        cos = 1.0
    ang = numpy.arccos(cos)
    return ang


def central_angle(xyz1, xyz2, xyz3):
    """Measure the angle inscribed by three atoms.

    :param xyz1: 3D vector to point 1
    :type xyz1: tuple, list, or numpy nd.array
    :param xyz2: 3D vector to point 2
    :type xyz2: tuple, list, or numpy nd.array
    :param xyz3: 3D vector to point 3
    :type xyz3: tuple, list, or numpy nd.array
    :rtype: float
    """
    return angle(xyz1=xyz1, xyz2=xyz3, orig_xyz=xyz2)


def projected_central_angle(xyz1, xyz2, xyz3):
    """Measure the angle inscribed by three atoms,
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
    """Measure the dihedral angle defined by four atoms.

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
    val = max(min(val, 1.0), -1.0)
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
    dih = numpy.mod(dih, 2 * numpy.pi)
    return dih


# transformations
def rotator(
    axis: Vector, ang: float, degree: bool = False, orig_xyz: Optional[Vector] = None
) -> Callable[[Vector], Vector]:
    """Get a function for axis-angle rotations, optionally specifying the origin

    :param axis: Rotational axis (norm is ignored)
    :type axis: Vector
    :param ang: Rotational angle
    :type ang: float
    :param degree: _description_, defaults to False
    :type degree: bool, optional
    :param orig_xyz: _description_, defaults to None
    :type orig_xyz: Optional[Vector], optional
    :return: _description_
    :rtype: Callable[[Vector], Vector]
    """
    orig_xyz = numpy.array([0.0, 0.0, 0.0] if orig_xyz is None else orig_xyz)
    ang = ang * phycon.DEG2RAD if degree else ang

    axis = unit_norm(axis)

    rot_mat = Rotation.from_rotvec(numpy.multiply(axis, ang)).as_matrix()

    def rotate_(xyz: Vector) -> Vector:
        xyz = numpy.array(xyz)
        return numpy.dot(rot_mat, xyz - orig_xyz) + orig_xyz

    return rotate_


# I/O
def string(vec, num_per_row=None, val_format="{0:>8.3f}"):
    """Write a vector to a string.

    :param vec: vector to form string with
    :type vec: list, tuple, or nd.array
    :param num_per_row: number of vector elements to write to a row
    :type num_per_row: int
    :rtype: str
    """

    if num_per_row is None:
        num_per_row = len(vec)

    assert isinstance(num_per_row, int), "num_per_row must be an integer"

    nvals = len(vec)
    vec_str = ""
    for i, val in enumerate(vec):
        vec_str += val_format.format(val)
        if ((i + 1) % num_per_row) == 0 and (i + 1) != nvals:
            vec_str += "\n"

    return vec_str
