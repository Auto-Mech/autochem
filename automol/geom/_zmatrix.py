""" z-matrix conversion
"""
from ._pyx2z import from_geometry as _x2m_from_geometry
from ._pyx2z import to_zmatrix as _x2m_to_zmatrix


def zmatrix(geo):
    """ z-matrix
    """
    x2m = _x2m_from_geometry(geo)
    zma = _x2m_to_zmatrix(x2m)
    return zma
