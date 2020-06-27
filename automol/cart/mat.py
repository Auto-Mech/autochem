""" transformation matrices
"""
import numpy
import transformations as tf


def random_rotation():
    """ random rotation matrix
    """
    rot_mat = _slice_affine(tf.random_rotation_matrix())
    return tuple(map(tuple, rot_mat))


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
