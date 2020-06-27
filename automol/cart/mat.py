""" transformation matrices
"""
import numpy
import transformations as tf


def random_rotation():
    """ random rotation matrix
    """
    rot_mat = _slice_affine(tf.random_rotation_matrix())
    return tuple(map(tuple, rot_mat))


def superimposition(points1, points2):
    """ superimposition matrix
    """
    pts1 = numpy.transpose(points1)
    pts2 = numpy.transpose(points2)
    rot_mat = _slice_affine(tf.superimposition_matrix(pts1, pts2))
    return rot_mat


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
