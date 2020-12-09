""" vector functions
"""
import numbers
import numpy
import transformations as tf


def unit_norm(xyz):
    """ vector normalized to 1
    """
    norm = numpy.linalg.norm(xyz)
    uxyz = numpy.divide(xyz, norm)
    assert numpy.allclose(numpy.linalg.norm(uxyz), 1.)
    return uxyz


def unit_direction(xyz1, xyz2):
    """ calculate a unit direction vector from `xyz1` to `xyz2`
    """
    dxyz12 = numpy.subtract(xyz2, xyz1)
    uxyz12 = unit_norm(dxyz12)
    return uxyz12


def are_parallel(xyz1, xyz2, orig_xyz=(0., 0., 0.), tol=1e-7):
    """ Are these vectors parallel to each other?
    """
    det = numpy.linalg.det([list(orig_xyz), list(xyz1), list(xyz2)])
    return det > tol


def unit_perpendicular(xyz1, xyz2, orig_xyz=(0., 0., 0.), allow_parallel=True):
    """ calculate a unit perpendicular on `xyz1` and `xyz2`
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
    """ calculate a unit bisector
    """
    ang = central_angle(xyz1, orig_xyz, xyz2)
    rot_ = rotater(
        axis=unit_perpendicular(xyz1, xyz2, orig_xyz), angle=ang/2.,
        orig_xyz=orig_xyz)
    xyz = unit_norm(numpy.subtract(rot_(xyz1), orig_xyz))
    return xyz


def from_internals(dist=0., xyz1=(0., 0., 0.), ang=0., xyz2=(0., 0., 1.),
                   dih=0., xyz3=(0., 1., 0.)):
    """ cartesian position from internal coordinates
    """
    local_xyz = _local_position(dist=dist, ang=ang, dih=dih)
    local_basis = _local_axes(xyz1=xyz1, xyz2=xyz2, xyz3=xyz3)
    xyz = tuple(xyz1 + numpy.dot(local_xyz, local_basis))
    return xyz


def _local_position(dist=0., ang=0., dih=0.):
    """ position by internal coordinates in the local axis frame
    """
    x_comp = dist * numpy.sin(ang) * numpy.sin(dih)
    y_comp = dist * numpy.sin(ang) * numpy.cos(dih)
    z_comp = dist * numpy.cos(ang)
    return (x_comp, y_comp, z_comp)


def _local_axes(xyz1=(0., 0., 0.), xyz2=(0., 0., 1.), xyz3=(0., 1., 0.)):
    """ local axes for defining bond, angle, dihedral from support atoms
    """
    uxyz12 = unit_direction(xyz1, xyz2)
    uxyz23 = unit_direction(xyz2, xyz3)
    uxyz123_perp = unit_perpendicular(uxyz23, uxyz12)

    z_ax = tuple(uxyz12)
    y_ax = tuple(unit_perpendicular(uxyz12, uxyz123_perp))
    x_ax = tuple(unit_perpendicular(y_ax, z_ax))
    return (x_ax, y_ax, z_ax)


def distance(xyz1, xyz2):
    """ measure the distance between points
    """
    dist = numpy.linalg.norm(numpy.subtract(xyz1, xyz2))
    return dist


def central_angle(xyz1, xyz2, xyz3):
    """ measure the angle inscribed by three atoms
    """
    uxyz21 = unit_direction(xyz2, xyz1)
    uxyz23 = unit_direction(xyz2, xyz3)
    ang = numpy.arccos(numpy.dot(uxyz21, uxyz23))
    return ang


def projected_central_angle(xyz1, xyz2, xyz3):
    """ measure the angle inscribed by three atoms, projected onto the normal
    plane of the central atom
    """
    uxyz21 = unit_perpendicular(xyz2, xyz1)
    uxyz23 = unit_perpendicular(xyz2, xyz3)
    ang = numpy.arccos(numpy.dot(uxyz21, uxyz23))
    return ang


def dihedral_angle(xyz1, xyz2, xyz3, xyz4):
    """ measure the dihedral angle defined by four atoms
    """
    # get the cosine of the angle
    uxyz21 = unit_direction(xyz2, xyz1)
    uxyz23 = unit_direction(xyz2, xyz3)
    uxyz32 = unit_direction(xyz3, xyz2)
    uxyz34 = unit_direction(xyz3, xyz4)
    uxyz123_perp = unit_perpendicular(uxyz21, uxyz23)
    uxyz234_perp = unit_perpendicular(uxyz32, uxyz34)
    cos = numpy.dot(uxyz123_perp, uxyz234_perp)

    # get the sign of the angle
    val = numpy.dot(uxyz123_perp, uxyz34)
    val = max(min(val, 1.), -1.)
    sign = 2 * (val < 0) - 1

    # before plugging it into the arccos function, make sure we haven't
    # slightly run out of bounds
    if cos < -1.:
        assert numpy.allclose(cos, -1.)
        cos = -1.
    elif cos > 1.:
        assert numpy.allclose(cos, 1.)
        cos = 1.

    dih = sign * numpy.arccos(cos)
    dih = numpy.mod(dih, 2*numpy.pi)
    return dih


def distance_matrix(xyzs):
    """ determine the distance matrix for a series of points
    """
    dist_mat = tuple(tuple(distance(xyz1, xyz2) if idx1 != idx2 else 0.
                           for idx2, xyz2 in enumerate(xyzs))
                     for idx1, xyz1 in enumerate(xyzs))
    return dist_mat


# transformations
def rotater(axis, angle, orig_xyz=None):
    """ a function to rotate vectors about an axis at a particular point
    """
    aug_rot_mat = tf.rotation_matrix(angle, axis, point=orig_xyz)
    return _transformer(aug_rot_mat)


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
def string(vec, num_per_row=3):
    """ Iterate over a 1D vector and print the values over some number of rows
    """
    nvals = len(vec)
    vec_str = ''
    for i, val in enumerate(vec):
        if ((i+1) % num_per_row) == 0 and (i+1) != nvals:
            vec_str += '{0:>8.3f}\n'.format(val)
        else:
            vec_str += '{0:>8.3f}'.format(val)

    return vec_str
