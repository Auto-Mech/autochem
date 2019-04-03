""" I/O for cartesian geometries
"""
import autoread as ar
import autowrite as aw
from automol.constructors.geom import from_data as _from_data
from ._core import symbols as _symbols
from ._core import coordinates as _coordinates


def from_xyz_string(xyz_str):
    """ read a cartesian geometry from a .xyz string
    """
    syms, xyzs = ar.geom.read_xyz(xyz_str)
    geo = _from_data(syms, xyzs, angstrom=True)
    return geo


def from_string(geo_str, angstrom=True):
    """ read a cartesian geometry from a string
    """
    syms, xyzs = ar.geom.read(geo_str)
    geo = _from_data(syms, xyzs, angstrom=angstrom)
    return geo


def string(geo, angstrom=True):
    """ write the cartesian geometry as a string
    """
    syms = _symbols(geo)
    xyzs = _coordinates(geo, angstrom=angstrom)
    geo_str = aw.geom.write(syms=syms, xyzs=xyzs)
    return geo_str


def xyz_string(geo):
    """ write the cartesian geometry to a .xyz string
    """
    syms = _symbols(geo)
    xyzs = _coordinates(geo, angstrom=True)
    geo_str = aw.geom.write_xyz(syms=syms, xyzs=xyzs)
    return geo_str
