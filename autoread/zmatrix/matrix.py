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

    captures symbols (column 1), keys (columns 2, 4, 6), and variable
    names/values (columns 3, 5, 7) from the z-matrix and returns them

    :param start_ptt: pattern before the start of the z-matrix
    :type start_ptt: str
    :param sym_ptt: matches atom symbol in the first column of the z-matrix
    :type sym_ptt: str
    :param key_ptt: matches key/index in columns 2, 4, 6 of the z-matrix
    :type key_ptt: str
    :param name_ptt: matches z-matrix variable names in columns 3, 5, 7; can
        also be used to match numbers at these positions
    :type name_ptt: str
    :param entry_start_ptt: matches before key_ptt
    :type entry_start_ptt: str
    :param entry_sep_ptt: matches between key_ptt and name_ptt
    :type entry_sep_ptt: str
    :param entry_end_ptt: matches after name_ptt
    :type entry_end_ptt: str
    :param line_start_ptt: matches at the start of a z-matrix line
    :type line_start_ptt: str
    :param line_end_ptt: matches at the end of a z-matrix line
    :type line_end_ptt: str
    :param last: capture the last match, instead of the first?
    :type last: bool
    :param case: make the match case-sensitive?
    :type case: bool
    :return: symbols, key matrix, variable name matrix
    :rtype: tuple
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
