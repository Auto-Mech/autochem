""" transformation matrices
"""
import numpy
import automol.util.vec
import transformations as tf


def random_rotation_matrix():
    """ Build a 3x3 rotation matrix to rotate vectors about some random angle.

        :rtype: tuple(tuple(float))
    """

    rot_mat = _slice_affine(tf.random_rotation_matrix())

    return tuple(map(tuple, rot_mat))


def axis_alignment_matrix(xyz1, xyz2):
    """ Build a 3x3 rotation matrix that aligns two points
        onto the x-axis and xy-plane (xyz1 will point along the
        x-axis and xyz2 will point in the xy-plane towards the y axis).

        :param xyz1: 3D vector to point 1
        :type xyz1: tuple, list, or numpy nd.array
        :param xyz2: 3D vector to point 2
        :type xyz2: tuple, list, or numpy nd.array
        :rtype: tuple(tuple(float))
    """

    assert numpy.linalg.norm(xyz1) > 1e-3, f'norm({xyz1}) != 0'
    assert numpy.linalg.norm(xyz2) > 1e-3, f'norm({xyz2}) != 0'
    mat = numpy.zeros((3, 3))
    mat[0] = automol.util.vec.unit_norm(xyz1)
    mat[2] = automol.util.vec.unit_perpendicular(xyz1, xyz2)
    mat[1] = numpy.cross(mat[2], mat[1])

    return tuple(map(tuple, mat))


def superimposition_matrix(points1, points2, keep_origin=False):
    """ superimposition matrix
    """
    if keep_origin:
        points1 = numpy.concatenate(([[0., 0., 0.]], points1), axis=0)
        points2 = numpy.concatenate(([[0., 0., 0.]], points2), axis=0)

    points1 = numpy.transpose(points1)
    points2 = numpy.transpose(points2)
    rot_mat = _slice_affine(tf.superimposition_matrix(points1, points2))
    return tuple(map(tuple, rot_mat))


def rotation_matrix(axis, angle):
    """ Build a 3x3 rotation matrix to rotate vectors about a given
        axis by a given angle.

        :param axis: axis to rotate about
        :type axis: list, tuple, or nd.array
        :param angle: angle by which to rotate the vectors
        :type angle: float
        :rtype: tuple(tuple(float))
    """
    rot_mat = _slice_affine(tf.rotation_matrix(angle, axis))
    return tuple(map(tuple, rot_mat))


def euler_rotation_matrix(theta, phi, psi):
    """ Build a 3x3 Euler rotation matrix to rotate vectors about
        three given Euler angles.

        :param theta: angle to rotate about z-axis
        :type theta: float
        :param phi: angle to rotate about x'-axis
        :type phi: float
        :param psi: angle to rotate about z'-axis
        :type psi: float
        :rtype: tuple(tuple(float))
    """
    rot_mat = _slice_affine(tf.euler_matrix(theta, phi, psi))
    return tuple(map(tuple, rot_mat))


def _slice_affine(mat):
    """ Return slice of a matrix up to 3 rows and 3 columns.

        :param mat: matrix to slice
        :type mat: tuple(tuple(float))
        :rtype: tuple(tuple(float))
    """
    return mat[:3, :3]


# I/O
def string(mat, val_format='{0:>8.3f}'):
    """ Write a matrix to a string.

        :param mat: matrix to form string with
        :type mat: tuple(tuple(float))
        :param precision: number of integers past decimal
        :type precision: int
        :rtype: str
    """

    mat_str = ''
    for row in mat:
        mat_str += ''.join(val_format.format(val) for val in row)
        mat_str += '\n'
    mat_str = mat_str.rstrip()

    return mat_str
