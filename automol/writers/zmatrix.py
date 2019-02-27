""" z-matrix writers
"""


def string(syms, key_mat, name_mat, val_dct, mat_delim=' ', setval_sign='='):
    """ write a z-matrix to a .zmat string
    """
    mat_str = matrix_block_string(syms=syms, key_mat=key_mat,
                                  name_mat=name_mat, delim=mat_delim)
    setval_str = setval_block_string(val_dct=val_dct, setval_sign=setval_sign)
    zma_str = '\n'.join((mat_str, setval_str))
    return zma_str


def matrix_block_string(syms, key_mat, name_mat, delim=' '):
    """ write the .zmat matrix block to a string
    """
    mat_str = ''
    for row_idx, sym in enumerate(syms):
        mat_str += '{:<2s}'.format(sym)
        for col_idx in range(min(row_idx, 3)):
            mat_str += '{}{:>d}'.format(delim, key_mat[row_idx][col_idx])
            mat_str += '{}{:>5s}'.format(delim, name_mat[row_idx][col_idx])
        mat_str += '\n'

    return mat_str


def setval_block_string(val_dct, setval_sign='='):
    """ write the .zmat setval block to a string
    """
    setval_str = ''
    for name, val in val_dct.items():
        setval_str += '{:<5s} {} {:>11.6f}\n'.format(name, setval_sign, val)

    return setval_str
