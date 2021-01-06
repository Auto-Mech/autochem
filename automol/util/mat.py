""" transformation matrices
"""
import numpy
import automol.util.vec
import transformations as tf


def random_rotation():
    """ random rotation matrix
    """
    rot_mat = _slice_affine(tf.random_rotation_matrix())
    return tuple(map(tuple, rot_mat))


def axis_alignment(xyz1, xyz2):
    """ a rotation matrix for aligning two points onto the x-axis and xy-plane

    (xyz1 will point along the x-axis and xyz2 will point in the xy-plane
    twoards the y axis)
    """
    assert numpy.linalg.norm(xyz1) > 1e-3, 'norm({}) != 0'.format(xyz1)
    assert numpy.linalg.norm(xyz2) > 1e-3, 'norm({}) != 0'.format(xyz2)
    mat = numpy.zeros((3, 3))
    mat[0] = automol.util.vec.unit_norm(xyz1)
    mat[2] = automol.util.vec.unit_perpendicular(xyz1, xyz2)
    mat[1] = numpy.cross(mat[2], mat[1])
    return tuple(map(tuple, mat))


def superimposition(points1, points2, keep_origin=False):
    """ superimposition matrix
    """
    if keep_origin:
        points1 = numpy.concatenate(([[0., 0., 0.]], points1), axis=0)
        points2 = numpy.concatenate(([[0., 0., 0.]], points2), axis=0)

    points1 = numpy.transpose(points1)
    points2 = numpy.transpose(points2)
    rot_mat = _slice_affine(tf.superimposition_matrix(points1, points2))
    return tuple(map(tuple, rot_mat))


def rotation(axis, angle):
    """ axis-angle rotation matrix
    """
    rot_mat = _slice_affine(tf.rotation_matrix(angle, axis))
    return tuple(map(tuple, rot_mat))


def euler_rotation(theta, phi, psi):
    """ euler rotation matrix
    """
    rot_mat = _slice_affine(tf.euler_matrix(theta, phi, psi))
    return tuple(map(tuple, rot_mat))


def _slice_affine(mat):
    return mat[:3, :3]


# I/O
def string(mat):
    """ Convert the Cartesian matrix
    """
    if mat is not None:
        mat_str = ''
        for row in mat:
            mat_str += '  '.join('{0:>8.3f}'.format(val) for val in row)
            mat_str += '\n'
        mat_str.rstrip()
    else:
        mat_str = ''

    return mat_str
