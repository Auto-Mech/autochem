""" z-matrix readers

(this can be cleaned up with new autoparse functionality)
"""
import numpy
from autoparse import cast as _cast
import autoparse.pattern as app
import autoparse.find as apf
from ..constructors.zmatrix import from_data as _from_data


SYM_PATTERN = app.LETTER + app.maybe(app.LETTER)
KEY_PATTERN = app.UNSIGNED_INTEGER
VAL_PATTERN = app.one_of_these([app.FLOAT, app.INTEGER])
NAME_PATTERN = app.VARIABLE_NAME


def from_string(zma_str,
                sym_ptt=SYM_PATTERN,
                key_ptt=KEY_PATTERN,
                val_ptt=VAL_PATTERN,
                name_ptt=NAME_PATTERN,
                mat_delim_ptt=app.LINESPACE,
                setv_ptt=app.escape('='),
                setval_delim_ptt=app.NEWLINE,
                block_delim_ptt=app.one_or_more(app.NEWLINE),
                one_indexed=True,
                angstrom=True,
                degree=True):
    """ read a z-matrix from a single string
    """
    zma_ptt = block_pattern(
        sym_ptt, key_ptt, val_ptt, name_ptt, mat_delim_ptt, setv_ptt,
        setval_delim_ptt, block_delim_ptt,
        capture_matrix_block=True,
        capture_setval_block=True)

    mat_str, setv_str = apf.first_capture(zma_ptt, zma_str)

    zma = from_matrix_and_setv_strings(
        mat_str, setv_str,
        sym_ptt=sym_ptt,
        key_ptt=key_ptt,
        val_ptt=val_ptt,
        name_ptt=name_ptt,
        mat_delim_ptt=mat_delim_ptt,
        setv_ptt=setv_ptt,
        one_indexed=one_indexed,
        angstrom=angstrom,
        degree=degree
    )
    return zma


def from_matrix_and_setv_strings(mat_str, setv_str,
                                 sym_ptt=SYM_PATTERN,
                                 key_ptt=KEY_PATTERN,
                                 val_ptt=VAL_PATTERN,
                                 name_ptt=NAME_PATTERN,
                                 mat_delim_ptt=app.LINESPACE,
                                 setv_ptt=app.escape('='),
                                 one_indexed=True,
                                 angstrom=True,
                                 degree=True):
    """ read a z-matrix from matrix and setval strings

    (use: first capture the matrix and setval blocks using the functions below,
    then call this function)
    """
    # parse the matrix string
    syms, key_mat, name_mat = matrix_block_information(
        mat_str, sym_ptt=sym_ptt, key_ptt=key_ptt, name_ptt=name_ptt,
        delim_ptt=mat_delim_ptt)

    # parse the setval string
    val_dct = setval_block_information(
        setv_str, name_ptt=name_ptt, setv_ptt=setv_ptt, val_ptt=val_ptt)

    return _from_data(syms, key_mat, name_mat, val_dct,
                      one_indexed=one_indexed,
                      angstrom=angstrom,
                      degree=degree)


def matrix_block_information(mat_str,
                             sym_ptt=SYM_PATTERN,
                             key_ptt=KEY_PATTERN,
                             name_ptt=NAME_PATTERN,
                             delim_ptt=app.LINESPACE):
    """ read the atomic symbols, the key matrix, and the name matrix from a
    z-matrix matrix block string
    """
    mat_str = mat_str.strip()
    lines = mat_str.splitlines()
    nrows = len(lines)
    syms = []
    key_mat = numpy.empty((nrows, 3), dtype=numpy.object_)
    name_mat = numpy.empty((nrows, 3), dtype=numpy.object_)
    for row_idx, line in enumerate(lines):
        ncols = min(row_idx, 3)
        line_ptt = matrix_line_pattern(
            num=ncols,
            sym_ptt=app.capturing(sym_ptt),
            key_ptt=app.capturing(key_ptt),
            name_ptt=app.capturing(name_ptt),
            delim_ptt=delim_ptt,
        )
        caps = apf.first_capture(line_ptt, line)
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


def setval_block_information(setv_str,
                             name_ptt=NAME_PATTERN,
                             setv_ptt=app.escape('='),
                             val_ptt=VAL_PATTERN):
    """ read z-matrix coordinate values from the z-matrix setval block string
    """
    setv_str = setv_str.strip()
    setv_ptt_ = setval_term_pattern(
        name_ptt=app.capturing(name_ptt),
        setv_ptt=setv_ptt,
        val_ptt=app.capturing(val_ptt),
    )
    caps = apf.all_captures(setv_ptt_, setv_str)
    val_dct = dict(_cast(caps))
    return val_dct


# patterns
def block_pattern(sym_ptt=SYM_PATTERN,
                  key_ptt=KEY_PATTERN,
                  val_ptt=VAL_PATTERN,
                  name_ptt=NAME_PATTERN,
                  mat_delim_ptt=app.LINESPACE,
                  setv_ptt=app.escape('='),
                  setval_delim_ptt=app.NEWLINE,
                  block_delim_ptt=app.one_or_more(app.NEWLINE),
                  capture_matrix_block=False,
                  capture_setval_block=False):
    """ pattern for a z-matrix string with subblocks separated by newlines
    """
    mat_ptt = matrix_block_pattern(
        sym_ptt=sym_ptt,
        key_ptt=key_ptt,
        name_ptt=name_ptt,
        delim_ptt=mat_delim_ptt,
    )
    setv_ptt = setval_block_pattern(
        name_ptt=name_ptt,
        setv_ptt=setv_ptt,
        val_ptt=val_ptt,
        delim_ptt=setval_delim_ptt,
    )
    mat_ptt = app.capturing(mat_ptt) if capture_matrix_block else mat_ptt
    setv_ptt = app.capturing(setv_ptt) if capture_setval_block else setv_ptt

    return block_delim_ptt.join([mat_ptt, setv_ptt])


def matrix_block_pattern(sym_ptt=SYM_PATTERN,
                         key_ptt=KEY_PATTERN,
                         name_ptt=NAME_PATTERN,
                         delim_ptt=app.LINESPACE):
    """ the matrix block of a z-matrix string
    """
    line_ptts = [
        matrix_line_pattern(num, sym_ptt, key_ptt, name_ptt, delim_ptt)
        for num in range(4)
    ]

    return app.one_of_these([
        app.NEWLINE.join(
            line_ptts[:3] + [app.series(line_ptts[3], app.NEWLINE)]),
        app.NEWLINE.join(line_ptts[:3]),
        app.NEWLINE.join(line_ptts[:2]),
        app.NEWLINE.join(line_ptts[:1]),
    ])


def matrix_line_pattern(num,
                        sym_ptt=SYM_PATTERN,
                        key_ptt=KEY_PATTERN,
                        name_ptt=NAME_PATTERN,
                        delim_ptt=app.LINESPACE):
    """ pattern to capture a line of the zmatrix
    """
    assert num in range(0, 4)
    return app.LINE_START + app.padded(
        app.padded(delim_ptt).join([sym_ptt] + num * [key_ptt, name_ptt]))


def setval_block_pattern(name_ptt=NAME_PATTERN,
                         setv_ptt=app.escape('='),
                         val_ptt=VAL_PATTERN,
                         delim_ptt=app.NEWLINE):
    """ captures the whole setval block of a z-matrix
    """
    term_ptt = setval_term_pattern(name_ptt=name_ptt, setv_ptt=setv_ptt,
                                   val_ptt=val_ptt)
    return app.LINE_START + app.padded(
        app.series(term_ptt, app.padded(delim_ptt)))


def setval_term_pattern(name_ptt=NAME_PATTERN, setv_ptt=app.escape('='),
                        val_ptt=VAL_PATTERN):
    """ pattern to capture names and values from the setval block
    """
    return app.LINESPACES.join([name_ptt, setv_ptt, val_ptt])
