""" transformation matrices
"""
import transformations


def rotation(axis, angle):
    """ axis-angle rotation matrix
    """
    rot_mat = _slice_affine(transformations.rotation_matrix(angle, axis))
    return tuple(map(tuple, rot_mat))


def euler_rotation(theta, phi, psi):
    """ euler rotation matrix
    """
    rot_mat = _slice_affine(transformations.euler_matrix(theta, phi, psi))
    return tuple(map(tuple, rot_mat))


def _slice_affine(mat):
    return mat[:3, :3]
