""" I/O for cartesian geometries
"""
import numpy
import autoparse.pattern as app
import autoparse.find as apf
import autoparse.conv as apc
import phycon.units as pcu
from .._cnst.geom import from_data
from ._core import symbols as _symbols
from ._core import coordinates as _coordinates

ATOM_SYMBOL_PATTERN = app.LETTER + app.maybe(app.LETTER)


def from_xyz_string(xyz_str):
    """ read a cartesian geometry from a .xyz string
    """
    lines = xyz_str.splitlines()
    assert apf.has_match(app.UNSIGNED_INTEGER, lines[0])
    natms = int(lines[0])
    # comment_line = lines[1]
    geo_str = '\n'.join(lines[2:natms+2])
    geo = from_string(geo_str, angstroms=True, strict=True)
    return geo


def xyz_string(geo, comment=''):
    """ write the cartesian geometry to a .xyz string
    """
    natms = len(_symbols(geo))
    assert not apf.has_match(app.NEWLINE, comment)
    geo_str = string(geo)
    xyz_str = '{:d}\n{:s}\n{:s}'.format(natms, comment, geo_str)
    return xyz_str


def from_string(geo_str, angstroms=True, strict=True):
    """ read a cartesian geometry from a string
    """
    pattern = app.LINESPACES.join([
        app.capturing(ATOM_SYMBOL_PATTERN),
        app.capturing(app.FLOAT),
        app.capturing(app.FLOAT),
        app.capturing(app.FLOAT),
    ])

    if strict:
        # first check the string
        line_pattern = app.maybe(app.LINESPACES).join(
            [app.LINE_START, pattern, app.LINE_END])
        lines = apf.strip_spaces(geo_str).splitlines()
        assert all(apf.has_match(line_pattern, line) for line in lines)

    mcaps = apf.all_captures(pattern, geo_str)
    mvals = apc.multis(mcaps, dtypes=(str, float, float, float))
    syms = tuple(mval[0] for mval in mvals)
    xyzs = tuple(mval[1:] for mval in mvals)
    geo = from_data(syms, xyzs, angstroms=angstroms)
    return geo


def string(geo, to_angstroms=True):
    """ write the cartesian geometry as a string
    """
    syms = _symbols(geo)
    xyzs = _coordinates(geo)
    xyzs = xyzs if not to_angstroms else numpy.multiply(xyzs, pcu.BOHR2ANG)

    geo_str = '\n'.join('{:2s} {:10.6f} {:10.6f} {:10.6f}'.format(sym, *xyz)
                        for sym, xyz in zip(syms, xyzs))
    return geo_str
