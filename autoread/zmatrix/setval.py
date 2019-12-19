""" z-matrix setval block parsers
"""
from autoparse import cast as _cast
import autoparse.pattern as app
import autoparse.find as apf

NAME_PATTERN = app.VARIABLE_NAME
VALUE_PATTERN = app.NUMBER
ENTRY_SEP_PATTERN = app.escape('=')
SEP_PATTERN = app.padded(app.NEWLINE)


def read(string,
         start_ptt=None,
         name_ptt=NAME_PATTERN,
         val_ptt=VALUE_PATTERN,
         entry_sep_ptt=ENTRY_SEP_PATTERN,
         entry_start_ptt=None,
         entry_end_ptt=None,
         sep_ptt=SEP_PATTERN,
         last=True,
         case=False):
    """ read setvalues from a string

    captures variable names and values from the setval block of a z-matrix* and
    returns them as a dictionary

    *works for any series of variable names and values -- they need not be in a
    z-matrix

    :param start_ptt: matches before the start of the setval block
    :type start_ptt: str
    :param name_ptt: matches the variable name in the setval block
    :type name_ptt: str
    :param val_ptt: matches the numeric value in the setval block
    :type name_ptt: str
    :param entry_sep_ptt: matches the separator between a variable name and
        its value, such as the equals sign in 'R1 = 5.00'
    :type entry_sep_ptt: str
    :param entry_start_ptt: matches at the start of a setval entry
    :type entry_start_ptt: str
    :param entry_end_ptt: matches at the end of a setval entry
    :param sep_ptt: matches the separator between setval entries, such as a
        newline or comma
    :param last: capture the last match, instead of the first?
    :type last: bool
    :param case: make the match case-sensitive?
    :type case: bool
    :return: variable names and values
    :rtype: dict
    """
    entry_ptt_ = entry_pattern(
        name_ptt=app.capturing(name_ptt), val_ptt=app.capturing(val_ptt),
        sep_ptt=entry_sep_ptt, start_ptt=entry_start_ptt,
        end_ptt=entry_end_ptt)
    block_ptt_ = app.capturing(block_pattern(
        name_ptt=name_ptt, val_ptt=val_ptt, entry_sep_ptt=entry_sep_ptt,
        entry_start_ptt=entry_start_ptt, entry_end_ptt=entry_end_ptt,
        sep_ptt=sep_ptt))

    block_ptt_ = block_ptt_ if start_ptt is None else start_ptt + block_ptt_

    block_str = (apf.last_capture(block_ptt_, string, case=case) if last else
                 apf.first_capture(block_ptt_, string, case=case))

    caps = apf.all_captures(entry_ptt_, block_str, case=case)
    val_dct = dict(_cast(caps))
    return val_dct


def block_pattern(name_ptt=NAME_PATTERN,
                  val_ptt=VALUE_PATTERN,
                  entry_sep_ptt=ENTRY_SEP_PATTERN,
                  entry_start_ptt=None,
                  entry_end_ptt=None,
                  sep_ptt=SEP_PATTERN):
    """ a single setvalue entry
    """
    entry_ptt = entry_pattern(
        name_ptt=name_ptt,
        val_ptt=val_ptt,
        sep_ptt=entry_sep_ptt,
        start_ptt=entry_start_ptt,
        end_ptt=entry_end_ptt,
    )
    block_ptt = app.series(entry_ptt, app.padded(sep_ptt))
    return block_ptt


def entry_pattern(name_ptt=NAME_PATTERN,
                  val_ptt=VALUE_PATTERN,
                  sep_ptt=ENTRY_SEP_PATTERN,
                  start_ptt=None,
                  end_ptt=None):
    """ a single setvalue entry
    """
    parts = (([] if start_ptt is None else [start_ptt]) +
             [name_ptt] +
             [sep_ptt] +
             [val_ptt] +
             ([] if end_ptt is None else [end_ptt]))

    ptt = app.padded(app.maybe(app.LINESPACES).join(parts))
    return ptt
