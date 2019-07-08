""" full z-matrix parsers
"""
import autoparse.find as apf
import autoparse.pattern as app
from autoread.zmatrix.matrix import read as _matrix_read
from autoread.zmatrix.setval import read as _setval_read
from autoread.zmatrix.matrix import block_pattern as _matrix_block_pattern
from autoread.zmatrix.setval import block_pattern as _setval_block_pattern
from autoread import par

KEY_PATTERN = app.UNSIGNED_INTEGER
NAME_PATTERN = app.VARIABLE_NAME
MAT_ENTRY_SEP_PATTERN = app.LINESPACE
SETVAL_START_PATTERN = app.one_or_more(app.NEWLINE)
SETVAL_ENTRY_SEP_PATTERN = app.escape('=')
SETVAL_SEP_PATTERN = app.padded(app.NEWLINE)


def read(string,
         start_ptt=None,
         sym_ptt=par.Pattern.ATOM_SYMBOL,
         key_ptt=KEY_PATTERN,
         name_ptt=NAME_PATTERN,
         val_ptt=par.Pattern.NUMERIC_VALUE,
         mat_entry_start_ptt=None,
         mat_entry_sep_ptt=MAT_ENTRY_SEP_PATTERN,
         mat_entry_end_ptt=None,
         mat_line_start_ptt=None,
         mat_line_end_ptt=None,
         setv_start_ptt=SETVAL_START_PATTERN,
         setv_entry_sep_ptt=SETVAL_ENTRY_SEP_PATTERN,
         setv_entry_start_ptt=None,
         setv_sep_ptt=SETVAL_SEP_PATTERN,
         last=True,
         case=False):
    """ read full z-matrix from a string
    """
    block_ptt_ = block_pattern(sym_ptt=sym_ptt,
                               key_ptt=key_ptt,
                               name_ptt=name_ptt,
                               val_ptt=val_ptt,
                               mat_entry_start_ptt=mat_entry_start_ptt,
                               mat_entry_sep_ptt=mat_entry_sep_ptt,
                               mat_entry_end_ptt=mat_entry_end_ptt,
                               mat_line_start_ptt=mat_line_start_ptt,
                               mat_line_end_ptt=mat_line_end_ptt,
                               setv_start_ptt=setv_start_ptt,
                               setv_entry_sep_ptt=setv_entry_sep_ptt,
                               setv_entry_start_ptt=setv_entry_start_ptt,
                               setv_sep_ptt=setv_sep_ptt,
                               capture_matrix_block=True,
                               capture_setval_block=True)

    if start_ptt is not None:
        block_ptt_ = start_ptt + block_ptt_

    if last:
        cap = apf.last_capture(block_ptt_, string, case=case)
        assert cap is not None
        mat_str, setv_str = cap
    else:
        cap = apf.first_capture(block_ptt_, string, case=case)
        assert cap is not None
        mat_str, setv_str = cap

    syms, key_mat, name_mat = _matrix_read(
        mat_str,
        sym_ptt=sym_ptt,
        key_ptt=key_ptt,
        name_ptt=name_ptt,
        entry_start_ptt=mat_entry_start_ptt,
        entry_sep_ptt=mat_entry_sep_ptt,
        entry_end_ptt=mat_entry_end_ptt,
        line_start_ptt=mat_line_start_ptt,
        line_end_ptt=mat_line_end_ptt)

    if len(syms) == 1:
        val_dct = {}
    else:
        val_dct = _setval_read(
            setv_str,
            name_ptt=name_ptt,
            val_ptt=val_ptt,
            entry_sep_ptt=setv_entry_sep_ptt,
            entry_start_ptt=setv_entry_start_ptt,
            sep_ptt=setv_sep_ptt)

    return syms, key_mat, name_mat, val_dct


def block_pattern(sym_ptt=par.Pattern.ATOM_SYMBOL,
                  key_ptt=KEY_PATTERN,
                  name_ptt=NAME_PATTERN,
                  val_ptt=par.Pattern.NUMERIC_VALUE,
                  mat_entry_start_ptt=None,
                  mat_entry_sep_ptt=MAT_ENTRY_SEP_PATTERN,
                  mat_entry_end_ptt=None,
                  mat_line_start_ptt=None,
                  mat_line_end_ptt=None,
                  setv_start_ptt=SETVAL_START_PATTERN,
                  setv_entry_sep_ptt=SETVAL_ENTRY_SEP_PATTERN,
                  setv_entry_start_ptt=None,
                  setv_sep_ptt=SETVAL_SEP_PATTERN,
                  capture_matrix_block=False,
                  capture_setval_block=False):
    """ full z-matrix pattern
    """
    mat_ptt = _matrix_block_pattern(
        sym_ptt=sym_ptt,
        key_ptt=key_ptt,
        name_ptt=name_ptt,
        entry_start_ptt=mat_entry_start_ptt,
        entry_sep_ptt=mat_entry_sep_ptt,
        entry_end_ptt=mat_entry_end_ptt,
        line_start_ptt=mat_line_start_ptt,
        line_end_ptt=mat_line_end_ptt)
    setv_ptt = app.maybe(_setval_block_pattern(
        name_ptt=name_ptt,
        val_ptt=val_ptt,
        entry_sep_ptt=setv_entry_sep_ptt,
        entry_start_ptt=setv_entry_start_ptt,
        sep_ptt=setv_sep_ptt))
    mat_ptt = app.capturing(mat_ptt) if capture_matrix_block else mat_ptt
    setv_ptt = app.capturing(setv_ptt) if capture_setval_block else setv_ptt
    block_ptt = app.padded(setv_start_ptt).join([mat_ptt, setv_ptt])
    return block_ptt
