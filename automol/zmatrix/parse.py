""" z-matrix parsing functions
"""
import autoparse.pattern as app
import autoparse.find as apf


_LSTART = app.LINE_START + app.maybe(app.LINESPACES)
_LEND = app.maybe(app.LINESPACES) + app.LINE_END

SYMBOL_PATTERN = app.LETTER + app.maybe(app.LETTER)
KEY_PATTERN = app.UNSIGNED_INTEGER
VARIABLE_PATTERN = app.VARIABLE_NAME
DELIM_PATTERN = app.escape('=')


def comment_line(zmat_str,
                 symbol_pattern=SYMBOL_PATTERN):
    """ get the comment line (the header, if there is one)
    """
    head_line = apf.strip_spaces(zmat_str.splitlines()[0])
    return '' if apf.full_match(symbol_pattern, head_line) else head_line


# matrix block
def matrix_symbols(mat_str,
                   symbol_pattern=SYMBOL_PATTERN):
    """ the atom symbols of a matrix block
    """
    pattern = matrix_line_symbol_pattern(symbol_pattern)
    syms = apf.all_captures(pattern, mat_str)
    return syms


def matrix_keys_and_entries(mat_str,
                            symbol_pattern=SYMBOL_PATTERN,
                            key_pattern=KEY_PATTERN,
                            variable_pattern=VARIABLE_PATTERN):
    """ the keys and entries of a matrix block, as matrices

    (float entries are converted to floats)
    """
    column_patterns = matrix_line_column_key_entry_patterns(symbol_pattern,
                                                            key_pattern,
                                                            variable_pattern)
    key_cols = []
    entry_cols = []
    for col_idx, column_pattern in enumerate(column_patterns):
        caps = apf.all_captures(column_pattern, mat_str)
        keys, entries = zip(*caps) if caps else ([], [])
        keys = list(map(int, keys))
        entries = [float(entry) if apf.is_number(entry) else entry
                   for entry in entries]
        key_cols.append([None] * (col_idx + 1) + keys)
        entry_cols.append([None] * (col_idx + 1) + entries)

    key_mat = tuple(zip(*key_cols))
    entry_mat = tuple(zip(*entry_cols))
    return key_mat, entry_mat


def matrix_block(zma_str,
                 symbol_pattern=SYMBOL_PATTERN,
                 key_pattern=KEY_PATTERN,
                 variable_pattern=VARIABLE_PATTERN):
    """ the matrix block of a z-matrix string
    """
    pattern = matrix_block_pattern(symbol_pattern, key_pattern,
                                   variable_pattern)
    mat_strs = apf.all_captures(pattern, zma_str)
    assert len(mat_strs) == 1
    mat_str, = mat_strs
    return mat_str


def matrix_line_symbol_pattern(symbol_pattern=SYMBOL_PATTERN):
    """ captures the atom symbol from a single line of the matrix block
    """
    return _LSTART + app.capturing(symbol_pattern)


def matrix_line_column_key_entry_patterns(symbol_pattern=SYMBOL_PATTERN,
                                          key_pattern=KEY_PATTERN,
                                          variable_pattern=VARIABLE_PATTERN):
    """ capture each column key and entry from a line of the matrix block
    """
    entry_pattern = _entry_pattern(variable_pattern)
    column_patterns = []
    for col_idx in range(3):
        patterns = ([symbol_pattern] + col_idx * [key_pattern, entry_pattern] +
                    [app.capturing(key_pattern), app.capturing(entry_pattern)])
        pattern = _LSTART + app.LINESPACES.join(patterns)
        column_patterns.append(pattern)
    return tuple(column_patterns)


def matrix_block_pattern(symbol_pattern=SYMBOL_PATTERN,
                         key_pattern=KEY_PATTERN,
                         variable_pattern=VARIABLE_PATTERN):
    """ captures the whole matrix block of a z-matrix
    """
    entry_pattern = _entry_pattern(variable_pattern)

    def _patterns(nentries):
        return [symbol_pattern] + nentries * [key_pattern, entry_pattern]

    line_patterns = [_LSTART + app.LINESPACES.join(_patterns(n)) + _LEND
                     for n in range(0, 4)]
    pattern = app.one_of_these([
        app.NEWLINE.join(line_patterns[:3])
        + app.zero_or_more(app.NEWLINE + line_patterns[3]),
        app.NEWLINE.join(line_patterns[:2]),
        app.NEWLINE.join(line_patterns[:1]),
    ])
    return app.capturing(pattern)


def _entry_pattern(variable_pattern=VARIABLE_PATTERN):
    return app.one_of_these([app.FLOAT, variable_pattern])


# variable block
def variable_values(var_str,
                    variable_pattern=VARIABLE_PATTERN,
                    delim_pattern=DELIM_PATTERN):
    """ the variable values of a variable block, as a dictionary of floats
    """
    pattern = variable_value_pattern(variable_pattern, delim_pattern)
    caps = apf.all_captures(pattern, var_str)
    var_names, var_vals = zip(*caps) if caps else ((), ())
    var_vals = list(map(float, var_vals))
    return dict(zip(var_names, var_vals))


def variable_block(zma_str,
                   variable_pattern=VARIABLE_PATTERN,
                   delim_pattern=DELIM_PATTERN):
    """ the variable block of a z-variable string
    """
    zma_str = apf.remove_empty_lines(zma_str)
    pattern = variable_block_pattern(variable_pattern, delim_pattern)
    var_strs = apf.all_captures(pattern, zma_str)
    assert len(var_strs) == 1
    var_str, = var_strs
    return var_str


def variable_value_pattern(variable_pattern=VARIABLE_PATTERN,
                           delim_pattern=DELIM_PATTERN):
    """ captures the variable name and value from a line of the variable block
    """
    patterns = [app.capturing(variable_pattern), delim_pattern,
                app.capturing(app.FLOAT)]
    pattern = _LSTART + app.LINESPACES.join(patterns) + _LEND
    return pattern


def variable_block_pattern(variable_pattern=VARIABLE_PATTERN,
                           delim_pattern=DELIM_PATTERN):
    """ captures the whole variable block of a z-matrix
    """
    patterns = [variable_pattern, delim_pattern, app.FLOAT]
    line_pattern = _LSTART + app.LINESPACES.join(patterns) + _LEND
    end_pattern = app.one_of_these([app.NEWLINE, app.STRING_END])
    pattern = app.one_or_more(line_pattern + end_pattern)
    return app.capturing(pattern)
