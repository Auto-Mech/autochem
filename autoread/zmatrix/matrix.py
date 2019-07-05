""" z-matrix matrix block parsers
"""
import numpy
from autoparse import cast as _cast
import autoparse.find as apf
import autoparse.pattern as app
from autoread import par

KEY_PATTERN = app.UNSIGNED_INTEGER
NAME_PATTERN = app.VARIABLE_NAME
ENTRY_SEP_PATTERN = app.LINESPACE


def read(string,
         start_ptt=None,
         sym_ptt=par.Pattern.ATOM_SYMBOL,
         key_ptt=KEY_PATTERN,
         name_ptt=NAME_PATTERN,
         entry_start_ptt=None,
         entry_sep_ptt=ENTRY_SEP_PATTERN,
         entry_end_ptt=None,
         line_start_ptt=None,
         line_end_ptt=None,
         last=True,
         case=False):
    """ read matrix from a string
    """
    line_ptts_ = [
        line_pattern(
            num,
            sym_ptt=app.capturing(sym_ptt),
            key_ptt=app.capturing(key_ptt),
            name_ptt=app.capturing(name_ptt),
            entry_start_ptt=entry_start_ptt,
            entry_sep_ptt=entry_sep_ptt,
            entry_end_ptt=entry_end_ptt,
            start_ptt=line_start_ptt,
            end_ptt=line_end_ptt,
        )
        for num in range(4)]

    block_ptt_ = app.capturing(block_pattern(
        sym_ptt=sym_ptt,
        key_ptt=key_ptt,
        name_ptt=name_ptt,
        entry_start_ptt=entry_start_ptt,
        entry_sep_ptt=entry_sep_ptt,
        entry_end_ptt=entry_end_ptt,
        line_start_ptt=line_start_ptt,
        line_end_ptt=line_end_ptt))

    block_ptt_ = block_ptt_ if start_ptt is None else start_ptt + block_ptt_

    block_str = (apf.last_capture(block_ptt_, string, case=case) if last else
                 apf.first_capture(block_ptt_, string, case=case))

    lines = block_str.splitlines()
    nrows = len(lines)
    syms = []
    key_mat = numpy.empty((nrows, 3), dtype=numpy.object_)
    name_mat = numpy.empty((nrows, 3), dtype=numpy.object_)
    for row_idx, line in enumerate(lines):
        ncols = min(row_idx, 3)
        caps = _cast(apf.first_capture(line_ptts_[ncols], line, case=case))
        sym = caps if ncols == 0 else caps[0]
        keys = caps[1::2]
        names = caps[2::2]

        syms.append(sym)
        key_mat[row_idx, :ncols] = keys
        name_mat[row_idx, :ncols] = names

    syms = tuple(syms)
    key_mat = tuple(map(tuple, key_mat))
    name_mat = tuple(map(tuple, name_mat))
    return syms, key_mat, name_mat


def block_pattern(sym_ptt=par.Pattern.ATOM_SYMBOL,
                  key_ptt=KEY_PATTERN,
                  name_ptt=NAME_PATTERN,
                  entry_start_ptt=None,
                  entry_sep_ptt=ENTRY_SEP_PATTERN,
                  entry_end_ptt=None,
                  line_start_ptt=None,
                  line_end_ptt=None):
    """ matrix pattern (assumes more than one atom)
    """
    line_ptts = [
        line_pattern(
            num,
            sym_ptt=sym_ptt,
            key_ptt=key_ptt,
            name_ptt=name_ptt,
            entry_start_ptt=entry_start_ptt,
            entry_sep_ptt=entry_sep_ptt,
            entry_end_ptt=entry_end_ptt,
            start_ptt=line_start_ptt,
            end_ptt=line_end_ptt,
        )
        for num in range(4)]

    block_end_ptt = app.series(line_ptts[3], app.padded(app.NEWLINE))

    block_ptt = app.one_of_these([
        app.padded(app.NEWLINE).join(line_ptts[:3] + [block_end_ptt]),
        app.padded(app.NEWLINE).join(line_ptts[:3]),
        app.padded(app.NEWLINE).join(line_ptts[:2]),
        app.padded(app.NEWLINE).join(line_ptts[:1]),
    ])
    return block_ptt


def line_pattern(num,
                 sym_ptt=par.Pattern.ATOM_SYMBOL,
                 key_ptt=KEY_PATTERN,
                 name_ptt=NAME_PATTERN,
                 entry_start_ptt=None,
                 entry_sep_ptt=ENTRY_SEP_PATTERN,
                 entry_end_ptt=None,
                 start_ptt=None,
                 end_ptt=None):
    """ matrix line pattern
    """
    assert num in range(0, 4)
    entry_ptt = entry_pattern(
        key_ptt=key_ptt,
        name_ptt=name_ptt,
        start_ptt=entry_start_ptt,
        sep_ptt=entry_sep_ptt,
        end_ptt=entry_end_ptt)

    parts = (
        [app.LINE_START] +
        ([] if start_ptt is None else [start_ptt]) +
        [sym_ptt] + num * [entry_ptt] +
        ([] if end_ptt is None else [end_ptt])
    )
    ptt = app.PADDING.join(parts)
    return ptt


def entry_pattern(key_ptt=KEY_PATTERN,
                  name_ptt=NAME_PATTERN,
                  start_ptt=None,
                  sep_ptt=ENTRY_SEP_PATTERN,
                  end_ptt=None):
    """ matrix entry pattern
    """
    parts = (
        ([] if start_ptt is None else [start_ptt]) +
        [key_ptt, sep_ptt, name_ptt] +
        ([] if end_ptt is None else [end_ptt])
    )
    ptt = app.PADDING.join(parts)
    return ptt
