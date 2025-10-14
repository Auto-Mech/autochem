"""vector functions."""
from collections.abc import Callable, Sequence

import numpy
from scipy.spatial.transform import Rotation
from phydat import phycon

Vector = numpy.ndarray


def is_vector_like(obj: Sequence | Vector) -> bool:
    """Check if an object has the right shape and dtype to be a vector.

    :param obj: The object
    :return: `True` if it does, `False` if it doesn't
    """
    arr = numpy.array(obj)
    return arr.shape == (3,) and numpy.issubdtype(arr.dtype, numpy.number)


def require_vector_like(obj: Sequence | Vector):
    """Require an object to have the right shape and dtype to be a vector.

    :param obj: The object
    :raises: Raises an AssertionError if the object is not vector-like
    """
    assert is_vector_like(obj), f"{obj} is not vector-like"


def unit_norm(xyz: Sequence | Vector) -> Vector:
    """Normalize a vector (xyz) to 1.0.

    Does nothing to null vectors

    :param xyz:3D Vector
    :return: Numbers in the vector
    """
    require_vector_like(xyz)
    uxyz = xyz
    norm = numpy.linalg.norm(xyz)
    if not numpy.isclose(norm, 0.0):
        uxyz = numpy.divide(xyz, norm)
    return uxyz


def unit_direction(xyz1: Sequence | Vector, xyz2: Sequence | Vector) -> Vector:
    """Calculate a unit direction vector from `xyz1` to `xyz2`.

    :param xyz1: 3D vector
    :param xyz2: 3D vector
    :return: Unit direction vector from two vectors
    """
    require_vector_like(xyz1)
    require_vector_like(xyz2)
    dxyz12 = numpy.subtract(xyz2, xyz1)
    uxyz12 = unit_norm(dxyz12)
    return uxyz12


def are_parallel(
    xyz1: Sequence | Vector,
    xyz2: Sequence | Vector,
    orig_xyz: Sequence | Vector = (0.0, 0.0, 0.0),
    anti: bool = False,
    tol: float = 1e-3,
) -> bool:
    """Assess if two vectors are parallel to each other.

    :param xyz1: 3D vector
    :param xyz2: 3D vector
    :param orig_xyz: Origin of coordinate system `xyz1` and `xyz2` are in
    :param anti: Test for anti-parallelity too?

    :param tol: Tolerance for checking determinant
    :return: True if vectors are parallel, false if not
    """
    require_vector_like(xyz1)
    require_vector_like(xyz2)
    require_vector_like(orig_xyz)
    vec1 = numpy.subtract(xyz1, orig_xyz)
    vec2 = numpy.subtract(xyz2, orig_xyz)
    len1 = numpy.linalg.norm(vec1)
    len2 = numpy.linalg.norm(vec2)

    ratio = numpy.dot(vec1, vec2) / (len1 * len2)

    ret = abs(ratio - 1) < tol
    if anti:
        ret |= abs(ratio + 1) < tol

    return ret


def orthogonalize(
    xyz1: Sequence | Vector, xyz2: Sequence | Vector, normalize: bool = False
) -> Vector:
    """Orthogonalize `xyz2` against `xyz1`.
    :param xyz1: 3D vector
    :param xyz2: 3D vector
    :param normalize: Whether to normalize the result
    :return: The orthogonal component.
    """
    require_vector_like(xyz1)
    require_vector_like(xyz2)
    overlap = numpy.dot(xyz1, xyz2)
    norm = numpy.dot(xyz1, xyz1)
    oxyz2 = numpy.subtract(xyz2, numpy.multiply(overlap / norm, xyz1))
    if normalize:
        oxyz2 = unit_norm(oxyz2)
    return oxyz2


def flip_if_left_handed(
    xvec: Sequence | Vector, yvec: Sequence | Vector, zvec: Sequence | Vector
) -> Vector:
    """Given three vectors, flips the third one, if necessary, to make a right-handed
    coordinate system.

    The vectors need not be orthogonal, but they must be linearly independent for this
    to work.  Does not change direction -- only flips the sign, if necessary.

    :param xvec: A vector defining the x direction
    :param yvec: A vector defining the y direction
    :param zvec: A vector defining the z direction
    :returns: `zvec` or, if they made a left-handed system, `-zvec`
    """
    require_vector_like(xvec)
    require_vector_like(yvec)
    require_vector_like(zvec)
    zdir = numpy.cross(xvec, yvec)
    proj = numpy.dot(zvec, zdir)
    if proj < 0.0:
        return numpy.negative(zvec)
    return zvec


def best_unit_perpendicular(xyzs: Sequence[Sequence | Vector]) -> Vector:
    """Find a vector that is perpendicular to a series of points as much as possible.

    For 1 point, this is an arbitrary vector.
    For 2 points, this is an arbitrary perpendicular vector.
    For 3 points, this is the exact normal to the plane.
    For 4+ points, this is the normal to the best-fit plane.

    Solution to the last problem was found here: https://stackoverflow.com/a/51132260

    :param xyzs: The points
    :return: New vector
    """
    assert all(
        is_vector_like(xyz) for xyz in xyzs
    ), f"One of these is not vector like: {xyzs}"
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

    return nvec


def arbitrary_unit_perpendicular(
    xyz: Sequence | Vector, orig_xyz: Sequence | Vector = (0.0, 0.0, 0.0)
) -> Vector:
    """Determine an arbitrary perpendicular vector."""
    require_vector_like(xyz)
    require_vector_like(orig_xyz)
    for xyz2 in ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1]):
        uxyz = unit_perpendicular(xyz, xyz2, orig_xyz=orig_xyz)
        if numpy.linalg.norm(uxyz) > 1e-7:
            break

    return uxyz


def unit_perpendicular(
    xyz1: Sequence | Vector,
    xyz2: Sequence | Vector,
    orig_xyz: Sequence | Vector = (0.0, 0.0, 0.0),
    value_if_parallel: Vector = (0.0, 0.0, 0.0),
) -> Vector:
    """Calculate a unit perpendicular on `xyz1` and `xyz2`.

    :param xyz1: 3D vector
    :param xyz2: 3D vector
    :param orig_xyz: origin of coordinate system `xyz1` and `xyz2` are in
    :param value_if_parallel: What to return if the vectors are parallel
    :return: Unit vector perpendicular to xyz1 and xyz2.
    """
    require_vector_like(xyz1)
    require_vector_like(xyz2)
    require_vector_like(orig_xyz)
    require_vector_like(value_if_parallel)
    xyz1 = numpy.subtract(xyz1, orig_xyz)
    xyz2 = numpy.subtract(xyz2, orig_xyz)
    xyz3 = numpy.cross(xyz1, xyz2)

    if numpy.linalg.norm(xyz3) > 1e-7:
        uxyz3 = unit_norm(xyz3)
    else:
        uxyz3 = numpy.array(value_if_parallel)

    return uxyz3


def unit_bisector(
    xyz1: Sequence | Vector,
    xyz2: Sequence | Vector,
    orig_xyz: Sequence | Vector = (0.0, 0.0, 0.0),
    outer: bool = False,
) -> Vector:
    """Calculate a unit bisector.

    :param xyz1: 3D vector
    :param xyz2: 3D vector
    :param orig_xyz: Origin of coordinate system `xyz1` and `xyz2` are in
    :param outer: Get the outer, instead of the inner bisector?, defaults to False
    :return: Unit bisector
    """
    require_vector_like(xyz1)
    require_vector_like(xyz2)
    require_vector_like(orig_xyz)
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
    dist: float = 0.0,
    xyz1: Sequence | Vector = (0.0, 0.0, 0.0),
    ang: float = 0.0,
    xyz2: Sequence | Vector = (0.0, 0.0, 1.0),
    dih: float = 0.0,
    xyz3: Sequence | Vector = (0.0, 1.0, 0.0),
) -> Vector:
    """Determine the position of a point (xyz4) in Cartesian coordinates whose
    position is related to three other points (xyz1, xyz2, xyz3) via a set
    of internal coordinates.

    :param dist: Distance between `xyz1` and `xyz2` (in bohr)
    :param angle: Angle between `xyz1`, `xyz2`, `xyz3` (in radians)
    :param xyz2: 3D vector to point 2
    :param dih: Dihedral from `xyz1`, `xyz2`, `xyz3` to `xyz4` (in radians)
    :param xyz3: 3D vector to point 2
    :return: New point 'xyz4'
    """
    require_vector_like(xyz1)
    require_vector_like(xyz2)
    require_vector_like(xyz3)
    local_xyz = _local_position(dist=dist, ang=ang, dih=dih)
    local_basis = _local_axes(xyz1=xyz1, xyz2=xyz2, xyz3=xyz3)
    xyz4 = xyz1 + numpy.dot(local_xyz, local_basis)

    return xyz4


def _local_position(dist: float = 0.0, ang: float = 0.0, dih: float = 0.0) -> Vector:
    """Determine the xyz coordinates of a point in the local axis frame
    defined by a set of internal coordinates.

    :param dist: distance between `xyz1` and `xyz2` (in bohr)
    :param angle: angle between `xyz1`, `xyz2`, `xyz3` (in radians)
    :param dih: dihedral from `xyz1`, `xyz2`, `xyz3` to `xyz4` (in radians)
    :return: New coordinates
    """
    x_comp = dist * numpy.sin(ang) * numpy.sin(dih)
    y_comp = dist * numpy.sin(ang) * numpy.cos(dih)
    z_comp = dist * numpy.cos(ang)

    return numpy.array([x_comp, y_comp, z_comp])


def _local_axes(
    xyz1: Sequence | Vector = (0.0, 0.0, 0.0),
    xyz2: Sequence | Vector = (0.0, 0.0, 1.0),
    xyz3: Sequence | Vector = (0.0, 1.0, 0.0),
) -> Vector:
    """Determine the  local axes for defining bond, angle, dihedral from
    the Cartesian coordinates of three support atoms.

    :param xyz1: 3D vector to point 1
    :param xyz2: 3D vector to point 2
    :param xyz3: 3D vector to point 2
    :return: Local axes
    """
    require_vector_like(xyz1)
    require_vector_like(xyz2)
    require_vector_like(xyz3)
    uxyz12 = unit_direction(xyz1, xyz2)
    uxyz23 = unit_direction(xyz2, xyz3)
    uxyz123_perp = unit_perpendicular(uxyz23, uxyz12)

    z_ax = uxyz12
    y_ax = unit_perpendicular(uxyz12, uxyz123_perp)
    x_ax = unit_perpendicular(y_ax, z_ax)
    return numpy.array([x_ax, y_ax, z_ax])


def distance(xyz1: Sequence | Vector, xyz2: Sequence | Vector) -> float:
    """Measure the distance between points.

    :param xyz1: 3D vector to point 1
    :param xyz2: 3D vector to point 2
    :return: Distance
    """
    require_vector_like(xyz1)
    require_vector_like(xyz2)
    dist = numpy.linalg.norm(numpy.subtract(xyz1, xyz2))

    return dist


def angle(
    xyz1: Sequence | Vector,
    xyz2: Sequence | Vector,
    orig_xyz: Sequence | Vector = (0.0, 0.0, 0.0),
) -> float:
    """Measure the angle inscribed by three atoms.

    :param xyz1: 3D vector to point 1
    :param xyz2: 3D vector to point 2
    :param orig_xyz: Origin to measure the vector angle from
    :return: Angle between three atoms
    """
    require_vector_like(xyz1)
    require_vector_like(xyz2)
    require_vector_like(orig_xyz)
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


def central_angle(
    xyz1: Sequence | Vector,
    xyz2: Sequence | Vector,
    xyz3: Sequence | Vector,
) -> float:
    """Measure the angle inscribed by three atoms.

    :param xyz1: 3D vector to point 1
    :param xyz2: 3D vector to point 2
    :param xyz3: 3D vector to point 3
    :return: Angle inscribed by three atoms
    """
    require_vector_like(xyz1)
    require_vector_like(xyz2)
    require_vector_like(xyz3)
    return angle(xyz1=xyz1, xyz2=xyz3, orig_xyz=xyz2)


def projected_central_angle(
    xyz1: Sequence | Vector,
    xyz2: Sequence | Vector,
    xyz3: Sequence | Vector,
) -> float:
    """Measure the angle inscribed by three atoms,
    projected onto the normal plane of the central atom.

    :param xyz1: 3D vector to point 1
    :param xyz2: 3D vector to point 2
    :param xyz3: 3D vector to point 3
    :return: Angle projected on normal plane
    """
    require_vector_like(xyz1)
    require_vector_like(xyz2)
    require_vector_like(xyz3)
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


def dihedral_angle(
    xyz1: Sequence | Vector,
    xyz2: Sequence | Vector,
    xyz3: Sequence | Vector,
    xyz4: Sequence | Vector,
) -> float:
    """Measure the dihedral angle defined by four atoms.

    :param xyz1: 3D vector to point 1
    :param xyz2: 3D vector to point 2
    :param xyz3: 3D vector to point 3
    :param xyz4: 3D vector to point 4
    :return: Dihedral angle of four atoms
    """
    require_vector_like(xyz1)
    require_vector_like(xyz2)
    require_vector_like(xyz3)
    require_vector_like(xyz4)
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
# Note: does not have updated format because that causes issues with geometry tests
def rotator(
    axis: Vector,
    ang: float,
    degree: bool = False,
    orig_xyz: Vector | Sequence = (0.0, 0.0, 0.0),
) -> Callable[[Vector], Vector]:
    """Get a function for axis-angle rotations, optionally specifying the origin.

    :param axis: Rotational axis (norm is ignored)
    :param ang: Rotational angle
    :param degree: _description_, defaults to False
    :return: _description_
    """
    orig_xyz = numpy.array(orig_xyz)
    ang = ang * phycon.DEG2RAD if degree else ang

    axis = unit_norm(axis)

    rot_mat = Rotation.from_rotvec(numpy.multiply(axis, ang)).as_matrix()

    def rotate_(xyz: Vector) -> Vector:
        xyz = numpy.array(xyz)
        return numpy.dot(rot_mat, xyz - orig_xyz) + orig_xyz

    return rotate_


# I/O
def string(
    vec: Sequence | Vector,
    num_per_row: int | None = None,
    val_format: str = "{0:>8.3f}",
) -> str:
    """Write a vector to a string.

    :param vec: vector to form string with
    :param num_per_row: number of vector elements to write to a row
    :return: String
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
