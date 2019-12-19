""" matrix parsers
"""
import numpy
from autoparse import cast as _cast
import autoparse.find as apf
import autoparse.pattern as app

VALUE_PATTERN = app.one_of_these([app.FLOAT])


def read(string,
         val_ptt=VALUE_PATTERN,
         start_ptt=None,
         block_start_ptt=None,
         line_start_ptt=None,
         last=True,
         tril=False,
         case=False):
    """ read matrix from a string

    captures values from a matrix, which can be rectangular or
    lower-triangular, and may or may not be broken into multiple blocks

    :param val_ptt: matches numeric matrix entry
    :type val_ptt: str
    :param start_ptt: matches before the start of the matrix
    :type start_ptt: str
    :param block_start_ptt: matches at the start of each block within the
        matrix
    :type block_start_ptt: str
    :param line_start_ptt: matches at the start of each line in the matrix or
        matrix block
    :param line_start_ptt: str
    :param last: capture the last match, instead of the first?
    :type last: bool
    :param case: make the match case-sensitive?
    :type case: bool
    """
    line_ptt_ = line_pattern(val_ptt=val_ptt, start_ptt=line_start_ptt,
                             capture_values=True)
    block_ptt_ = block_pattern(val_ptt=val_ptt, start_ptt=block_start_ptt,
                               line_start_ptt=line_start_ptt,
                               capture_block=True)
    blocks_ptt_ = blocks_pattern(val_ptt=val_ptt, start_ptt=start_ptt,
                                 block_start_ptt=block_start_ptt,
                                 line_start_ptt=line_start_ptt,
                                 capture_blocks=True)

    blocks_str = (apf.last_capture(blocks_ptt_, string, case=case) if last else
                  apf.first_capture(blocks_ptt_, string, case=case))

    block_strs = apf.all_captures(block_ptt_, blocks_str, case=case)

    if not tril:
        rows = numpy.concatenate(
            [_block_rows(block_str, val_ptt, line_ptt_, case=case)
             for block_str in block_strs], axis=1)
        mat = _matrix(rows)
    else:
        rows = list(_block_rows(block_strs[0], val_ptt, line_ptt_, case=case))
        nrows = len(rows)
        for block_str in block_strs[1:]:
            block_rows = _block_rows(block_str, val_ptt, line_ptt_, case=case)
            nblock_rows = len(block_rows)
            for block_row_idx, row_idx in enumerate(
                    range(nrows-nblock_rows, nrows)):
                rows[row_idx] += block_rows[block_row_idx]

        mat = _symmetric_matrix_from_lower_triangle(rows)

    return mat


def _matrix(rows):
    mat = numpy.array(rows)
    assert mat.ndim == 2
    return tuple(map(tuple, mat))


def _symmetric_matrix_from_lower_triangle(tril_rows):
    nrows = len(tril_rows)
    mat = numpy.zeros((nrows, nrows), dtype=numpy.object_)
    for row_idx, tril_row in enumerate(tril_rows):
        mat[row_idx, :row_idx+1] = tril_row
        mat[:row_idx+1, row_idx] = tril_row
    # make sure we ended up with a square symmetric matrix
    assert mat.ndim == 2 and mat.shape[0] == mat.shape[1]
    return tuple(map(tuple, mat))


def _block_rows(block_str, val_ptt, line_ptt_, case=False):
    rows = []
    val_ptt_ = app.capturing(val_ptt)
    for line_str in apf.all_captures(line_ptt_, block_str, case=case):
        row = _cast(apf.all_captures(val_ptt_, line_str))
        rows.append(row)
    return rows


def blocks_pattern(val_ptt=VALUE_PATTERN,
                   start_ptt=None,
                   block_start_ptt=None,
                   line_start_ptt=None,
                   capture_blocks=False):
    """ multi-block matrix pattern
    """
    block_ptt = block_pattern(val_ptt=val_ptt, start_ptt=block_start_ptt,
                              line_start_ptt=line_start_ptt)

    blocks_ptt_ = app.series(block_ptt, app.padded(app.NEWLINE))

    if capture_blocks:
        blocks_ptt_ = app.capturing(blocks_ptt_)

    blocks_ptt_ = blocks_ptt_ if start_ptt is None else start_ptt + blocks_ptt_

    return blocks_ptt_


def block_pattern(val_ptt=VALUE_PATTERN,
                  start_ptt=None,
                  line_start_ptt=None,
                  capture_block=False):
    """ matrix block pattern
    """
    line_ptt = line_pattern(val_ptt=val_ptt, start_ptt=line_start_ptt)
    block_ptt_ = app.series(line_ptt, app.padded(app.NEWLINE))

    if capture_block:
        block_ptt_ = app.capturing(block_ptt_)

    block_ptt_ = block_ptt_ if start_ptt is None else start_ptt + block_ptt_

    return app.padded(block_ptt_)


def line_pattern(val_ptt=VALUE_PATTERN,
                 start_ptt=None,
                 capture_values=False):
    """ matrix line pattern
    """
    vals_ptt = app.series(val_ptt, app.LINESPACES)

    if capture_values:
        vals_ptt = app.capturing(vals_ptt)

    parts = (
        ([] if start_ptt is None else [start_ptt]) + [vals_ptt])

    ptt = app.LINE_START + app.padded(app.LINESPACES.join(parts))
    return ptt
