""" I/O for cartesian geometries
"""
import numpy
import phycon.units as pcu
from ..readers import geom as _geom_reader
from ._core import symbols as _symbols
from ._core import coordinates as _coordinates


def from_xyz_string(xyz_str):
    """ read a cartesian geometry from a .xyz string
    """
    geo = _geom_reader.from_xyz_string(xyz_str)
    return geo


def xyz_string(geo, comment=''):
    """ write the cartesian geometry to a .xyz string
    """
    natms = len(_symbols(geo))
    geo_str = string(geo)
    xyz_str = '{:d}\n{:s}\n{:s}'.format(natms, comment, geo_str)
    return xyz_str


def from_string(geo_str, angstrom=True, strict=True):
    """ read a cartesian geometry from a string
    """
    geo = _geom_reader.from_string(geo_str, angstrom=angstrom, strict=strict)
    return geo


def string(geo, to_angstrom=True):
    """ write the cartesian geometry as a string
    """
    syms = _symbols(geo)
    xyzs = _coordinates(geo)
    xyzs = xyzs if not to_angstrom else numpy.multiply(xyzs, pcu.BOHR2ANG)

    geo_str = '\n'.join('{:2s} {:10.6f} {:10.6f} {:10.6f}'.format(sym, *xyz)
                        for sym, xyz in zip(syms, xyzs))
    return geo_str
