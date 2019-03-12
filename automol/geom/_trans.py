""" molecular geometry transformations
"""
import numpy
from ..cart.mat import rotation as _rotation_matrix
from ..constructors.geom import from_data
from ._core import symbols as _symbols
from ._core import coordinates as _coordinates


def rotate(geo, axis, angle):
    """ axis-angle rotation of the geometry
    """
    syms = _symbols(geo)
    xyzs = _coordinates(geo)
    rot_mat = _rotation_matrix(axis, angle)
    xyzs = numpy.dot(xyzs, numpy.transpose(rot_mat))
    return from_data(syms, xyzs)
