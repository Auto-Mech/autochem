""" transformation matrices
"""

from transformations import rotation_matrix
from transformations import euler_matrix


def rotation(axis, angle):
    """ axis-angle rotation matrix
    """
    rot_mat = _slice_affine(rotation_matrix(angle, axis))
    return tuple(map(tuple, rot_mat))


def euler_rotation(theta, phi, psi):
    """ euler rotation matrix
    """
    rot_mat = _slice_affine(euler_matrix(theta, phi, psi))
    return tuple(map(tuple, rot_mat))


def _slice_affine(mat):
    return mat[:3, :3]
