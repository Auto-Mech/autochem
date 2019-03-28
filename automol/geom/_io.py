""" I/O for cartesian geometries
"""
import numpy
import phycon.units as pcu
import autoparser as apr
from automol.constructors.geom import from_data as _from_data
from ._core import symbols as _symbols
from ._core import coordinates as _coordinates


def from_xyz_string(xyz_str):
    """ read a cartesian geometry from a .xyz string
    """
    lines = xyz_str.splitlines()
    try:
        natms = int(lines[0])
    except ValueError:
        raise ValueError('Invalid xyz string')

    geo_str = '\n'.join(lines[2:natms+2])
    geo = from_string(geo_str, angstrom=True)
    return geo


def xyz_string(geo, comment=''):
    """ write the cartesian geometry to a .xyz string
    """
    natms = len(_symbols(geo))
    geo_str = string(geo)
    xyz_str = '{:d}\n{:s}\n{:s}'.format(natms, comment, geo_str)
    return xyz_str


def from_string(geo_str, angstrom=True):
    """ read a cartesian geometry from a string
    """
    syms, xyzs = apr.geom.read(geo_str)
    geo = _from_data(syms, xyzs, angstrom=angstrom)
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
