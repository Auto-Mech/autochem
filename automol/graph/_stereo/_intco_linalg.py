""" integer 3-vector library
"""
from numbers import Real as _Real
from numbers import Integral as _Integral
import numpy


def unit_direction(int_xyz1, int_xyz2):
    """ unit direction vector pointing from `int_xyz1` to `int_xyz2`
    """
    int_xyz = numpy.subtract(int_xyz2, int_xyz1)
    uint_xyz = numpy.divide(int_xyz, numpy.linalg.norm(int_xyz))
    assert is_unit_integer_triple(uint_xyz)
    return tuple(numpy.array(uint_xyz, dtype=int))


def local_coordinate_interpreter(trans, rot):
    """ interpret local coordinates with a given origin and rotation matrix

    :param trans: translation vector
    :type trans: (int, int, int)
    :param rot: rotation matrix
    :type rot: ((int, int, int), (int, int, int), (int, int, int))

    :returns: integer triple local-> global coordinate transformation
    :rtype: callable
    """

    def _global_coordinates(loc_xyz):
        return tuple(numpy.add(trans, numpy.dot(rot, loc_xyz)))

    return _global_coordinates


def aligning_rotation_matrix(uint_xyz, uint_xyz_target):
    """ rotation matrix to align `uint_xyz` with `uint_xyz_target`

    :param uint_xyz: unit integer vector (+/- a standard basis vector)
    :type uint_xyz: (int, int, int)
    :param uint_xyz_target: unit integer vector (+/- a standard basis vector)
    :type uint_xyz_target: (int, int, int)

    :returns: integer rotation mapping
    :rtype: callable
    """
    assert is_unit_integer_triple(uint_xyz)
    assert is_unit_integer_triple(uint_xyz_target)
    dot = numpy.vdot(uint_xyz, uint_xyz_target)
    cross = numpy.cross(uint_xyz, uint_xyz_target)
    if numpy.any(cross):
        rot = rotation_matrix(cross, clicks=1)
    else:
        assert dot in (-1, 1)
        clicks = 0 if dot == 1 else 2
        perp = unit_perpendicular(uint_xyz)
        rot = rotation_matrix(perp, clicks=clicks)
    return rot


def rotation_matrix(uint_xyz, clicks):
    """ 3d integer rotation matrix

    :param uint_xyz: unit integer vector (+/- a standard basis vector)
    :type uint_xyz: (int, int, int)
    :param clicks: angle of roation, in units of pi/2
    :type clicks: int

    :returns: matrix
    :rtype: ((int, int, int), (int, int, int), (int, int, int))
    """
    assert is_unit_integer_triple(uint_xyz)
    assert isinstance(clicks, _Integral)
    angle = clicks * numpy.pi / 2.
    cos = int(numpy.cos(angle))
    sin = int(numpy.sin(angle))
    uix, uiy, uiz = uint_xyz
    eye = numpy.eye(3, dtype=int)
    ux_ = numpy.array([[0, -uiz, uiy],
                       [uiz, 0, -uix],
                       [-uiy, uix, 0]], dtype=int)
    uxu = numpy.outer(uint_xyz, uint_xyz)
    rot = cos * eye + sin * ux_ + (1 - cos) * uxu
    # assert that this is actually a rotation matrix
    assert numpy.all(numpy.dot(rot, rot.T) == eye)
    assert numpy.all(numpy.dot(rot.T, rot) == eye)
    assert numpy.linalg.det(rot) == 1
    return rot


def unit_perpendicular(uint_xyz):
    """ perpendicular unit vector to this one

    :param uint_xyz: unit integer vector (+/- a standard basis vector)
    :type uint_xyz: (int, int, int)
    """
    assert is_unit_integer_triple(uint_xyz)
    uix, uiy, uiz = uint_xyz
    if uix:
        uint_xyz_perp = numpy.cross(uint_xyz, [0, 0, -1])
    elif uiy:
        uint_xyz_perp = numpy.cross(uint_xyz, [0, 0, +1])
    else:
        assert numpy.abs(uiz) == 1
        uint_xyz_perp = numpy.cross(uint_xyz, [0, -1, 0])
    assert is_unit_integer_triple(uint_xyz_perp)
    return uint_xyz_perp


def is_unit_integer_triple(vec):
    """ is this vector a unit integer triple?
    """
    return is_integer_triple(vec) and is_unit_vector(vec)


def is_unit_vector(vec):
    """ is this a unit vector
    """
    return numpy.linalg.norm(vec) == 1


def is_integer_triple(vec):
    """ is this an integer triple?
    """
    return len(vec) == 3 and all(map(is_integer, vec))


def is_integer(num):
    """ is this number an integer?
    """
    assert isinstance(num, _Real)
    return isinstance(num, _Integral) or float.is_integer(num)
