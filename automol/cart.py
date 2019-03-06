""" cartesian coordinate linear algbera helpers
"""
import numpy


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


def unit_perpendicular(xyz1, xyz2):
    """ calculate a unit perpendicular on `xyz1` and `xyz2`
    """
    xyz3 = numpy.cross(xyz1, xyz2)
    uxyz3 = unit_norm(xyz3)
    return uxyz3
