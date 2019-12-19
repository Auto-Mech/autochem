""" energy parsers
"""
from autoparse import cast as _cast
import autoparse.find as apf
import autoparse.pattern as app

# VALUE_PATTERN = app.one_of_these([app.FLOAT])
VALUE_PATTERN = app.one_of_these([app.EXPONENTIAL_FLOAT_D, app.FLOAT])
SEP_PATTERN = app.LINESPACES


def read(string,
         start_ptt,
         val_ptt=VALUE_PATTERN,
         last=True,
         case=False):
    """ read energy from a string

    :param start_ptt: matches before the numeric value
    :type start_ptt: str
    :param val_ptt: matches the numeric value
    :type val_ptt: str
    :param last: capture the last match, instead of the first?
    :type last: bool
    :param case: make the match case-sensitive?
    :type case: bool
    :return: energy, converted to float if possible
    """
    ptt_ = pattern(start_ptt=start_ptt, val_ptt=app.capturing(val_ptt))

    ene_str = (apf.last_capture(ptt_, string, case=case) if last else
               apf.first_capture(ptt_, string, case=case))
    ene = _cast(ene_str.replace('D', 'E'))
    return ene


def pattern(start_ptt,
            val_ptt=VALUE_PATTERN):
    """ energy pattern
    """
    return start_ptt + app.lpadded(val_ptt)
