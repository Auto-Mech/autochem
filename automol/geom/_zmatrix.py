""" z-matrix conversion
"""
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
