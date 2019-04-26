""" geometry parsers
"""
from autoparse import cast as _cast
import autoparse.find as apf
import autoparse.pattern as app
from autoread import par


def read(string,
         sym_ptt=par.Pattern.ATOM_SYMBOL,
         val_ptt=par.Pattern.NUMERIC_VALUE,
         start_ptt=None,
         line_sep_ptt=None,
         line_start_ptt=None,
         last=True,
         case=False):
    """ read geometry from a string
    """
    line_ptt_ = line_pattern(
        sym_ptt=app.capturing(sym_ptt), val_ptt=app.capturing(val_ptt),
        sep_ptt=line_sep_ptt, start_ptt=line_start_ptt)
    block_ptt_ = app.capturing(block_pattern(
        sym_ptt=sym_ptt, val_ptt=val_ptt, line_sep_ptt=line_sep_ptt,
        line_start_ptt=line_start_ptt))

    block_ptt_ = block_ptt_ if start_ptt is None else start_ptt + block_ptt_

    block_str = (apf.last_capture(block_ptt_, string, case=case) if last else
                 apf.first_capture(block_ptt_, string, case=case))

    caps = apf.all_captures(line_ptt_, block_str)
    syms, xcomps, ycomps, zcomps = zip(*_cast(caps))
    xyzs = tuple(zip(xcomps, ycomps, zcomps))
    return syms, xyzs


def read_xyz(string, sym_ptt=par.Pattern.ATOM_SYMBOL,
             val_ptt=par.Pattern.NUMERIC_VALUE):
    """ read geometry from a .xyz string
    """
    lines = string.splitlines()
    try:
        natms = int(lines[0])
    except ValueError:
        raise ValueError('Invalid xyz string')

    geo_str = '\n'.join(lines[2:natms+2])
    syms, xyzs = read(geo_str, sym_ptt=sym_ptt, val_ptt=val_ptt)
    return syms, xyzs


def block_pattern(sym_ptt=par.Pattern.ATOM_SYMBOL,
                  val_ptt=par.Pattern.NUMERIC_VALUE,
                  line_sep_ptt=None,
                  line_start_ptt=None):
    """ geometry block pattern
    """
    line_ptt = line_pattern(
        sym_ptt=sym_ptt, val_ptt=val_ptt, sep_ptt=line_sep_ptt,
        start_ptt=line_start_ptt)
    block_ptt = app.series(line_ptt, app.padded(app.NEWLINE))
    return block_ptt


def line_pattern(sym_ptt=par.Pattern.ATOM_SYMBOL,
                 val_ptt=par.Pattern.NUMERIC_VALUE,
                 sep_ptt=None,
                 start_ptt=None):
    """ geometry line pattern
    """
    parts = (
        ([] if start_ptt is None else [start_ptt]) +
        [sym_ptt] +
        ([] if sep_ptt is None else [sep_ptt]) +
        3 * [val_ptt]
    )
    ptt = app.LINE_START + app.padded(app.LINESPACES.join(parts))
    return ptt
