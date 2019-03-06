""" z-matrix conversion

borrows from https://github.com/tmpchem/computational_chemistry
"""
import numpy
from ..cart import unit_direction as _unit_direction
from ..cart import unit_perpendicular as _unit_perpendicular
from ._core import coordinates as _coordinates
from ._pyx2z import from_geometry as _x2m_from_geometry
from ._pyx2z import to_zmatrix as _x2m_to_zmatrix
from ._pyx2z import (zmatrix_rotational_coordinate_names as
                     _x2m_zmatrix_rotational_coordinate_names)


def zmatrix(geo):
    """ z-matrix
    """
    x2m = _x2m_from_geometry(geo)
    zma = _x2m_to_zmatrix(x2m)
    return zma


def zmatrix_rotational_coordinate_names(geo):
    """ z-matrix rotational coordinate names
    """
    x2m = _x2m_from_geometry(geo)
    names = _x2m_zmatrix_rotational_coordinate_names(x2m)
    return names


def distance(geo, key1, key2):
    """ measure the distance between atoms
    """
    xyzs = _coordinates(geo)
    xyz1 = xyzs[key1]
    xyz2 = xyzs[key2]
    dist = numpy.linalg.norm(numpy.subtract(xyz1, xyz2))
    return dist


def angle(geo, key1, key2, key3):
    """ measure the angle inscribed by three atoms
    """
    xyzs = _coordinates(geo)
    xyz1 = xyzs[key1]
    xyz2 = xyzs[key2]
    xyz3 = xyzs[key3]
    uxyz21 = _unit_direction(xyz2, xyz1)
    uxyz23 = _unit_direction(xyz2, xyz3)
    ang = numpy.arccos(numpy.dot(uxyz21, uxyz23))
    return ang


def torsion(geo, key1, key2, key3, key4):
    """ measure the torsion angle defined by four atoms
    """
    xyzs = _coordinates(geo)
    xyz1 = xyzs[key1]
    xyz2 = xyzs[key2]
    xyz3 = xyzs[key3]
    xyz4 = xyzs[key4]

    # get the cosine of the angle
    uxyz21 = _unit_direction(xyz2, xyz1)
    uxyz23 = _unit_direction(xyz2, xyz3)
    uxyz32 = _unit_direction(xyz3, xyz2)
    uxyz34 = _unit_direction(xyz3, xyz4)
    uxyz123_perp = _unit_perpendicular(uxyz21, uxyz23)
    uxyz234_perp = _unit_perpendicular(uxyz32, uxyz34)
    cos = numpy.dot(uxyz123_perp, uxyz234_perp)

    # get the sign of the angle
    val = numpy.dot(uxyz123_perp, uxyz34)
    val = max(min(val, 1.), -1.)
    sign = 2 * (val < 0) - 1

    tors = sign * numpy.arccos(cos)
    return tors
