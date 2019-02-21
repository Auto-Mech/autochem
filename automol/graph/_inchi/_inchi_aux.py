""" functions operating on InChI-AuxInfo strings
"""
from autoparse.pattern import escape as _escape
from autoparse.pattern import one_or_more as _one_or_more
from autoparse.pattern import one_of_these as _one_of_these
from autoparse.pattern import named_capturing as _named_capturing
from autoparse.pattern import UNSIGNED_INTEGER as _UNSIGNED_INTEGER
from autoparse.pattern import NONSPACE as _NONSPACE
from autoparse.pattern import STRING_END as _STRING_END
from autoparse.find import first_named_capture as _first_named_capture
from autoparse.find import all_captures as _all_captures

_NONSPACES_NONGREEDY = _one_or_more(_NONSPACE, greedy=False)
_INCHI_SUBLAYER_END = _one_of_these([_escape('/'), _STRING_END])


class PARSE():
    """ InChI-AuxInfo format specifications """
    class NUMBERING():
        """ _ """
        LAYER_KEY = 'all'
        CONTENT_KEY = 'content'

        _START = _escape('/')
        _LAYER = (_escape('N:') +
                  _named_capturing(_NONSPACES_NONGREEDY,
                                   name=CONTENT_KEY))
        _END = _INCHI_SUBLAYER_END

        PATTERN = _START + _named_capturing(_LAYER, name=LAYER_KEY) + _END

        class NUMBER():
            """ _ """
            PATTERN = _UNSIGNED_INTEGER


def sorted_atom_keys(ich_aux):
    """ zero-indexed numbering
    """
    cap_dct = _first_named_capture(PARSE.NUMBERING.PATTERN, ich_aux)
    lyr = cap_dct[PARSE.NUMBERING.CONTENT_KEY]
    atm_keys = tuple(
        map(int, _all_captures(PARSE.NUMBERING.NUMBER.PATTERN, lyr)))
    return atm_keys
