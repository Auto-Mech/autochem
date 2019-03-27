""" cartesian geometry readers
"""
from autoparse import cast as _cast
import autoparse.pattern as app
import autoparse.find as apf
from ..constructors.geom import from_data


SYM_PATTERN = app.LETTER + app.maybe(app.LETTER)
VAL_PATTERN = app.one_of_these([app.FLOAT, app.INTEGER])


def from_xyz_string(xyz_str):
    """ read a cartesian geometry from a .xyz string
    """
    lines = xyz_str.splitlines()
    assert apf.has_match(app.UNSIGNED_INTEGER, lines[0])
    natms = int(lines[0])
    # comment_line = lines[1]
    geo_str = '\n'.join(lines[2:natms+2])
    geo = from_string(geo_str, angstrom=True, strict=True)
    return geo


def from_string(geo_str,
                sym_pattern=SYM_PATTERN,
                val_pattern=VAL_PATTERN,
                angstrom=True,
                strict=False):
    """ read a geometry from a string
    """
    _block_pattern = block_pattern(sym_pattern=sym_pattern,
                                   val_pattern=val_pattern)
    _line_pattern_ = line_pattern(app.capturing(sym_pattern),
                                  app.capturing(val_pattern))

    if strict:
        assert apf.full_match(_block_pattern, geo_str)

    block_str = apf.first_capture(app.capturing(_block_pattern), geo_str)
    mcaps = apf.all_captures(_line_pattern_, block_str)
    mvals = _cast(mcaps)
    syms = tuple(mval[0] for mval in mvals)
    xyzs = tuple(mval[1:] for mval in mvals)
    geo = from_data(syms, xyzs, angstrom=angstrom)
    return geo


def block_pattern(sym_pattern=SYM_PATTERN, val_pattern=VAL_PATTERN):
    """ pattern for the whole geometry block
    """
    return app.series(
        line_pattern(sym_pattern=sym_pattern, val_pattern=val_pattern),
        app.NEWLINE)


def line_pattern(sym_pattern=SYM_PATTERN, val_pattern=VAL_PATTERN):
    """ pattern to capture values from a line of the geometry block
    """
    pattern = app.LINE_START + app.padded(
        app.LINESPACES.join([sym_pattern] + 3 * [val_pattern]))
    return pattern
