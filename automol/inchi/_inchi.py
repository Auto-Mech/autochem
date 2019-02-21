""" functions operating on InChI strings
"""
import itertools
from string import ascii_lowercase as _ascii_lowercase
from autoparse.pattern import escape as _escape
from autoparse.pattern import named_capturing as _named_capturing
from autoparse.pattern import one_or_more as _one_or_more
from autoparse.pattern import one_of_these as _one_of_these
from autoparse.pattern import not_followed_by as _not_followed_by
from autoparse.pattern import LOWERCASE_LETTER as _LOWERCASE_LETTER
from autoparse.pattern import UNSIGNED_INTEGER as _UNSIGNED_INTEGER
from autoparse.pattern import NONSPACE as _NONSPACE
from autoparse.pattern import STRING_START as _STRING_START
from autoparse.pattern import STRING_END as _STRING_END
from autoparse.find import first_named_capture as _first_named_capture
from autoparse.find import all_captures as _all_captures
from ._rdkit import from_inchi as _rdm_from_inchi
from ._rdkit import to_inchi as _rdm_to_inchi

_NONSPACES_NONGREEDY = _one_or_more(_NONSPACE, greedy=False)
_INCHI_SUBLAYER_END = _one_of_these([_escape('/'), _STRING_END])
_STEREO_UNKNOWN_VAL = 'u'
_STEREO_UNDEFINED_VAL = '?'
_STEREO_MINUS_VAL = '-'
_STEREO_PLUS_VAL = '+'


def _key_layer(key):
    assert key in _ascii_lowercase

    class _KEYxLAYER():
        LAYER_KEY = 'all'
        CONTENT_KEY = 'content'

        _START = _escape('/')
        _LAYER = key + _named_capturing(_NONSPACES_NONGREEDY,
                                        name=CONTENT_KEY)
        _END = _INCHI_SUBLAYER_END

        PATTERN = _START + _named_capturing(_LAYER, name=LAYER_KEY) + _END

    _KEYxLAYER.__name__ = 'KEYxLAYER(\'{:s}\')'.format(key)
    return _KEYxLAYER


class PARSE():
    """ InChI format specifications """

    @classmethod
    def key_layer_(cls, key):
        """ _ """
        return _key_layer(key)

    class PREFIX():
        """ _ """
        LAYER_KEY = 'all'
        CONTENT_KEY = 'content'

        _START = _STRING_START
        _LAYER = (_escape('InChI=') +
                  _named_capturing(_NONSPACES_NONGREEDY,
                                   name=CONTENT_KEY))
        _END = _INCHI_SUBLAYER_END

        PATTERN = _START + _named_capturing(_LAYER, name=LAYER_KEY) + _END

    class FORMULA():
        """ _ """
        LAYER_KEY = 'all'
        CONTENT_KEY = 'content'

        _START = _escape('/')
        _LAYER = (_not_followed_by(_LOWERCASE_LETTER) +
                  _named_capturing(_NONSPACES_NONGREEDY,
                                   name=CONTENT_KEY))
        _END = _INCHI_SUBLAYER_END

        PATTERN = _START + _named_capturing(_LAYER, name=LAYER_KEY) + _END

    class ATOMxSTEREO(_key_layer('t')):
        """ _ """

        class TERM():
            """ _ """
            KEY_KEY = 'key'
            VAL_KEY = 'val'

            PLUS_VAL = _STEREO_PLUS_VAL
            MINUS_VAL = _STEREO_MINUS_VAL
            UNKNOWN_VALS = (_STEREO_UNKNOWN_VAL, _STEREO_UNDEFINED_VAL)
            VALS = (MINUS_VAL, PLUS_VAL) + UNKNOWN_VALS

            _KEY = _UNSIGNED_INTEGER
            _VAL = _one_of_these(list(map(_escape, VALS)))
            PATTERN = (_named_capturing(_KEY, name=KEY_KEY) +
                       _named_capturing(_VAL, name=VAL_KEY))

    class BONDxSTEREO(_key_layer('b')):
        """ _ """

        class TERM():
            """ _ """
            KEY_KEY = 'key'
            VAL_KEY = 'val'

            PLUS_VAL = _STEREO_PLUS_VAL
            MINUS_VAL = _STEREO_MINUS_VAL
            UNKNOWN_VALS = (_STEREO_UNKNOWN_VAL, _STEREO_UNDEFINED_VAL)
            VALS = (MINUS_VAL, PLUS_VAL) + UNKNOWN_VALS

            _KEY = _UNSIGNED_INTEGER + _escape('-') + _UNSIGNED_INTEGER
            _VAL = _one_of_these(list(map(_escape, VALS)))
            PATTERN = (_named_capturing(_KEY, name=KEY_KEY) +
                       _named_capturing(_VAL, name=VAL_KEY))


def recalculate(ich, force_stereo=False):
    """ recalculate InChI string
    """
    _options = '-SUU' if force_stereo else ''
    rdm = _rdm_from_inchi(ich)
    ich = _rdm_to_inchi(rdm, options=_options, with_aux_info=False)
    return ich


def is_closed(ich):
    """ regenerating the InChI string yields the same thing
    """
    return recalculate(ich) == ich


def prefix(ich):
    """ InChI prefix
    """
    cap_dct = _first_named_capture(PARSE.PREFIX.PATTERN, ich)
    assert cap_dct
    pfx = cap_dct[PARSE.PREFIX.LAYER_KEY]
    return pfx


def version(ich):
    """ InChI version
    """
    cap_dct = _first_named_capture(PARSE.PREFIX.PATTERN, ich)
    assert cap_dct
    ver = cap_dct[PARSE.PREFIX.CONTENT_KEY]
    return ver


def formula_layer(ich):
    """ InChI formula
    """
    cap_dct = _first_named_capture(PARSE.FORMULA.PATTERN, ich)
    assert cap_dct
    fml = cap_dct[PARSE.FORMULA.LAYER_KEY]
    return fml


def key_layer(ich, key):
    """ a sublayer from the InChI string, by key
    """
    key_layer_parser = PARSE.key_layer_(key)
    cap_dct = _first_named_capture(key_layer_parser.PATTERN, ich)
    return cap_dct[key_layer_parser.LAYER_KEY] if cap_dct else None


def key_layer_content(ich, key):
    """ a sublayer from the InChI string, by key
    """
    key_layer_parser = PARSE.key_layer_(key)
    cap_dct = _first_named_capture(key_layer_parser.PATTERN, ich)
    return cap_dct[key_layer_parser.CONTENT_KEY] if cap_dct else None


def core_parent(ich):
    """ get the InChI string of the core parent structure
    """
    lyrs = [prefix(ich), formula_layer(ich)]
    for key in ('c', 'h'):
        lyr = key_layer(ich, key=key)
        if lyr is not None:
            lyrs.append(lyr)
    return '/'.join(lyrs)


def atom_stereo_elements(ich):
    """ atom stereo keys and values
    """
    cap_dct = _first_named_capture(PARSE.ATOMxSTEREO.PATTERN, ich)
    ret = ()
    if cap_dct:
        lyr = cap_dct[PARSE.ATOMxSTEREO.LAYER_KEY]
        ret = _all_captures(PARSE.ATOMxSTEREO.TERM.PATTERN, lyr)
    return ret


def bond_stereo_elements(ich):
    """ bond stereo keys and values
    """
    cap_dct = _first_named_capture(PARSE.BONDxSTEREO.PATTERN, ich)
    ret = ()
    if cap_dct:
        lyr = cap_dct[PARSE.BONDxSTEREO.LAYER_KEY]
        ret = _all_captures(PARSE.BONDxSTEREO.TERM.PATTERN, lyr)
    return ret


def known_atom_stereo_elements(ich):
    """ atom stereo keys and values, excluding unknown stereo
    """
    return tuple((key, val) for key, val in atom_stereo_elements(ich)
                 if val not in PARSE.ATOMxSTEREO.TERM.UNKNOWN_VALS)


def known_bond_stereo_elements(ich):
    """ bond stereo keys and values, excluding unknown stereo
    """
    return tuple((key, val) for key, val in bond_stereo_elements(ich)
                 if val not in PARSE.BONDxSTEREO.TERM.UNKNOWN_VALS)


def has_unknown_stereo_elements(ich):
    """ does this InChI string have unknown stereo elements?
    """
    ich_ste = recalculate(ich, force_stereo=True)
    return (atom_stereo_elements(ich_ste) !=
            known_atom_stereo_elements(ich_ste) or
            bond_stereo_elements(ich_ste) !=
            known_bond_stereo_elements(ich_ste))


def compatible_stereoisomers(ich):
    """ expand InChI string to its compatible stereoisomers
    """
    atm_terms_lst = []
    for num, (key, val) in enumerate(
            atom_stereo_elements(recalculate(ich, force_stereo=True))):
        terms = ([key + val] if val not in PARSE.ATOMxSTEREO.TERM.UNKNOWN_VALS
                 else [key + PARSE.ATOMxSTEREO.TERM.MINUS_VAL] if num == 0
                 else [key + PARSE.ATOMxSTEREO.TERM.MINUS_VAL,
                       key + PARSE.ATOMxSTEREO.TERM.PLUS_VAL])
        atm_terms_lst.append(terms)

    bnd_terms_lst = []
    for key, val in bond_stereo_elements(recalculate(ich, force_stereo=True)):
        terms = ([key + val] if val not in PARSE.BONDxSTEREO.TERM.UNKNOWN_VALS
                 else [key + PARSE.BONDxSTEREO.TERM.MINUS_VAL,
                       key + PARSE.BONDxSTEREO.TERM.PLUS_VAL])
        bnd_terms_lst.append(terms)

    ich_cp = core_parent(ich)
    ich_lst = [ich_cp]
    if bnd_terms_lst:
        lyr_lst = ['b' + ','.join(terms)
                   for terms in itertools.product(*bnd_terms_lst)]
        ich_lst = [ich_start + '/' + lyr
                   for ich_start, lyr in itertools.product(ich_lst, lyr_lst)]

    if atm_terms_lst:
        lyr_lst = ['t' + ','.join(terms)
                   for terms in itertools.product(*atm_terms_lst)]
        ich_lst = [ich_start + '/' + lyr
                   for ich_start, lyr in itertools.product(ich_lst, lyr_lst)]

    ich_lst = tuple(map(recalculate, ich_lst))
    return ich_lst
