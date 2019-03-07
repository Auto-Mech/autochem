""" z-matrix readers

(this can be cleaned up with new autoparse functionality)
"""
import numpy
import autoparse.pattern as app
import autoparse.find as apf
import autoparse.conv as apc
from ..constructors.zmatrix import from_data as _from_data


_LSTART = app.LINE_START + app.maybe(app.LINESPACES)
_LEND = app.maybe(app.LINESPACES) + app.LINE_END

SYM_PATTERN = app.LETTER + app.maybe(app.LETTER)
KEY_PATTERN = app.UNSIGNED_INTEGER
VAL_PATTERN = app.one_of_these([app.FLOAT, app.INTEGER])
NAME_PATTERN = app.VARIABLE_NAME


def from_string(zma_str,
                sym_pattern=SYM_PATTERN,
                key_pattern=KEY_PATTERN,
                val_pattern=VAL_PATTERN,
                name_pattern=NAME_PATTERN,
                mat_delim_pattern=app.LINESPACE,
                setval_pattern=app.escape('='),
                setval_delim_pattern=app.NEWLINE,
                one_indexed=True,
                angstrom=True,
                degree=True):
    """ read a z-matrix from a single string

    (grabs the first matrix and setval blocks it finds)
    """
    mat_block_pattern = matrix_block_capturing_pattern(
        sym_pattern=sym_pattern,
        key_pattern=key_pattern,
        val_pattern=val_pattern,
        name_pattern=name_pattern,
        delim_pattern=mat_delim_pattern,
    )
    setval_block_pattern = setval_block_capturing_pattern(
        name_pattern=name_pattern,
        setval_pattern=setval_pattern,
        val_pattern=val_pattern,
        delim_pattern=setval_delim_pattern,
    )

    mat_str = apf.first_capture(mat_block_pattern, zma_str)
    setval_str = apf.first_capture(setval_block_pattern, zma_str)

    zma = from_matrix_and_setval_strings(
        mat_str, setval_str,
        sym_pattern=sym_pattern,
        key_pattern=key_pattern,
        val_pattern=val_pattern,
        name_pattern=name_pattern,
        mat_delim_pattern=mat_delim_pattern,
        setval_pattern=setval_pattern,
        one_indexed=one_indexed,
        angstrom=angstrom,
        degree=degree
    )
    return zma


def from_matrix_and_setval_strings(mat_str, setval_str,
                                   sym_pattern=SYM_PATTERN,
                                   key_pattern=KEY_PATTERN,
                                   val_pattern=VAL_PATTERN,
                                   name_pattern=NAME_PATTERN,
                                   mat_delim_pattern=app.LINESPACE,
                                   setval_pattern=app.escape('='),
                                   one_indexed=True,
                                   angstrom=True,
                                   degree=True):
    """ read a z-matrix from matrix and setval strings

    (use: first capture the matrix and setval blocks using the functions below,
    then call this function)
    """
    # parse the matrix string
    col_patterns = matrix_column_capturing_patterns(
        sym_pattern=sym_pattern, key_pattern=key_pattern,
        val_pattern=val_pattern, name_pattern=name_pattern,
        delim_pattern=mat_delim_pattern)
    syms = apf.all_captures(col_patterns[0], mat_str)
    natms = len(syms)

    key_mat = numpy.empty((natms, 3), dtype=numpy.object_)
    name_mat = numpy.empty((natms, 3), dtype=numpy.object_)

    for col, col_pattern in enumerate(col_patterns[1:]):
        caps = apf.all_captures(col_pattern, mat_str)
        if caps:
            vals = apc.multis(caps, (int, str))
            keys, names = zip(*vals)
            key_mat[col+1:, col] = keys
            name_mat[col+1:, col] = names

    # parse the setval string
    pattern = setval_capturing_pattern(
        name_pattern=name_pattern, setval_pattern=setval_pattern,
        val_pattern=val_pattern)
    caps = apf.all_captures(pattern, setval_str)
    val_dct = dict(apc.multis(caps, (str, float)))

    return _from_data(syms, key_mat, name_mat, val_dct,
                      one_indexed=one_indexed,
                      angstrom=angstrom,
                      degree=degree)


# matrix block parsers
def matrix_block_capturing_pattern(sym_pattern=SYM_PATTERN,
                                   key_pattern=KEY_PATTERN,
                                   val_pattern=VAL_PATTERN,
                                   name_pattern=NAME_PATTERN,
                                   delim_pattern=app.LINESPACE):
    """ the matrix block of a z-matrix string
    """
    _entry_pattern = app.one_of_these([val_pattern, name_pattern])
    _delim_pattern = (app.maybe(app.LINESPACES) + delim_pattern +
                      app.maybe(app.LINESPACES))

    def _parts(nentries):
        return [sym_pattern] + nentries * [key_pattern, _entry_pattern]

    line_patterns = [_LSTART + _delim_pattern.join(_parts(n)) + _LEND
                     for n in range(0, 4)]
    pattern = app.one_of_these([
        app.NEWLINE.join(line_patterns[:3])
        + app.zero_or_more(app.NEWLINE + line_patterns[3]),
        app.NEWLINE.join(line_patterns[:2]),
        app.NEWLINE.join(line_patterns[:1]),
    ])
    pattern = app.capturing(pattern)
    return pattern


def matrix_column_capturing_patterns(sym_pattern=SYM_PATTERN,
                                     key_pattern=KEY_PATTERN,
                                     val_pattern=VAL_PATTERN,
                                     name_pattern=NAME_PATTERN,
                                     delim_pattern=app.LINESPACE):
    """ patterns to capture each of the four columns of the matrix block
    """
    _entry_pattern = app.one_of_these([val_pattern, name_pattern])
    _delim_pattern = (app.maybe(app.LINESPACES) + delim_pattern +
                      app.maybe(app.LINESPACES))

    sym_capturing_pattern = _LSTART + app.capturing(sym_pattern)

    column_patterns = [sym_capturing_pattern]
    for col in range(3):
        parts = ([sym_pattern] + col * [key_pattern, _entry_pattern] +
                 [app.capturing(key_pattern), app.capturing(_entry_pattern)])
        pattern = _LSTART + _delim_pattern.join(parts)
        column_patterns.append(pattern)

    return column_patterns


# setval block parsers
def setval_block_capturing_pattern(name_pattern=NAME_PATTERN,
                                   setval_pattern=app.escape('='),
                                   val_pattern=VAL_PATTERN,
                                   delim_pattern=app.NEWLINE):
    """ captures the whole setval block of a z-matrix
    """
    setval_pattern = app.LINESPACES.join([name_pattern, setval_pattern,
                                          val_pattern])
    _delim_pattern = (app.maybe(app.LINESPACES) + delim_pattern +
                      app.maybe(app.LINESPACES))
    pattern = (
        _LSTART + app.capturing(
            setval_pattern +
            app.zero_or_more(_delim_pattern + setval_pattern)) + _LEND)
    return pattern


def setval_capturing_pattern(name_pattern=NAME_PATTERN,
                             setval_pattern=app.escape('='),
                             val_pattern=VAL_PATTERN):
    """ pattern to capture names and values from the setval block
    """
    pattern = app.LINESPACES.join([app.capturing(name_pattern), setval_pattern,
                                   app.capturing(val_pattern)])
    return pattern
