""" InChIs
"""
import itertools
from string import ascii_lowercase as _ascii_lowercase
import autoparse.pattern as app
import autoparse.find as apf
import automol.convert.inchi

_NONSPACES_NONGREEDY = app.one_or_more(app.NONSPACE, greedy=False)
_INCHI_SUBLAYER_END = app.one_of_these([app.escape('/'), app.STRING_END])
_STEREO_UNKNOWN_VAL = 'u'
_STEREO_UNDEFINED_VAL = '?'
_STEREO_MINUS_VAL = '-'
_STEREO_PLUS_VAL = '+'


def _key_layer_class(key):
    assert key in _ascii_lowercase

    class _KeyLayer():
        LAYER_KEY = 'all'
        CONTENT_KEY = 'content'

        _START = app.escape('/')
        _LAYER = key + app.named_capturing(_NONSPACES_NONGREEDY,
                                           name=CONTENT_KEY)
        _END = _INCHI_SUBLAYER_END

        PATTERN = _START + app.named_capturing(_LAYER, name=LAYER_KEY) + _END

    _KeyLayer.__name__ = 'KeyLayer(\'{:s}\')'.format(key)
    return _KeyLayer


class Parse():
    """ InChI format specifications """

    @classmethod
    def key_layer_(cls, key):
        """ _ """
        return _key_layer_class(key)

    class Prefix():
        """ _ """
        LAYER_KEY = 'all'
        CONTENT_KEY = 'content'

        _START = app.STRING_START
        _LAYER = (app.escape('InChI=') +
                  app.named_capturing(_NONSPACES_NONGREEDY,
                                      name=CONTENT_KEY))
        _END = _INCHI_SUBLAYER_END

        PATTERN = _START + app.named_capturing(_LAYER, name=LAYER_KEY) + _END

    class Formula():
        """ _ """
        LAYER_KEY = 'all'
        CONTENT_KEY = 'content'

        _START = app.escape('/')
        _LAYER = (app.not_followed_by(app.LOWERCASE_LETTER) +
                  app.named_capturing(_NONSPACES_NONGREEDY,
                                      name=CONTENT_KEY))
        _END = _INCHI_SUBLAYER_END

        PATTERN = _START + app.named_capturing(_LAYER, name=LAYER_KEY) + _END

    class AtomStereo(_key_layer_class('t')):
        """ _ """

        class Term():
            """ _ """
            KEY_KEY = 'key'
            VAL_KEY = 'val'

            PLUS_VAL = _STEREO_PLUS_VAL
            MINUS_VAL = _STEREO_MINUS_VAL
            UNKNOWN_VALS = (_STEREO_UNKNOWN_VAL, _STEREO_UNDEFINED_VAL)
            VALS = (MINUS_VAL, PLUS_VAL) + UNKNOWN_VALS

            _KEY = app.UNSIGNED_INTEGER
            _VAL = app.one_of_these(list(map(app.escape, VALS)))
            PATTERN = (app.named_capturing(_KEY, name=KEY_KEY) +
                       app.named_capturing(_VAL, name=VAL_KEY))

    class BondStereo(_key_layer_class('b')):
        """ _ """

        class Term():
            """ _ """
            KEY_KEY = 'key'
            VAL_KEY = 'val'

            PLUS_VAL = _STEREO_PLUS_VAL
            MINUS_VAL = _STEREO_MINUS_VAL
            UNKNOWN_VALS = (_STEREO_UNKNOWN_VAL, _STEREO_UNDEFINED_VAL)
            VALS = (MINUS_VAL, PLUS_VAL) + UNKNOWN_VALS

            _KEY = (app.UNSIGNED_INTEGER + app.escape('-') +
                    app.UNSIGNED_INTEGER)
            _VAL = app.one_of_these(list(map(app.escape, VALS)))
            PATTERN = (app.named_capturing(_KEY, name=KEY_KEY) +
                       app.named_capturing(_VAL, name=VAL_KEY))


def is_closed(ich):
    """ regenerating the InChI string yields the same thing
    """
    return recalculate(ich) == ich


def prefix(ich):
    """ InChI prefix
    """
    cap_dct = apf.first_named_capture(Parse.Prefix.PATTERN, ich)
    assert cap_dct
    pfx = cap_dct[Parse.Prefix.LAYER_KEY]
    return pfx


def version(ich):
    """ InChI version
    """
    cap_dct = apf.first_named_capture(Parse.Prefix.PATTERN, ich)
    assert cap_dct
    ver = cap_dct[Parse.Prefix.CONTENT_KEY]
    return ver


def formula_layer(ich):
    """ InChI formula
    """
    cap_dct = apf.first_named_capture(Parse.Formula.PATTERN, ich)
    assert cap_dct
    fml = cap_dct[Parse.Formula.LAYER_KEY]
    return fml


def key_layer(ich, key):
    """ a sublayer from the InChI string, by key
    """
    key_layer_parser = Parse.key_layer_(key)
    cap_dct = apf.first_named_capture(key_layer_parser.PATTERN, ich)
    return cap_dct[key_layer_parser.LAYER_KEY] if cap_dct else None


def key_layer_content(ich, key):
    """ a sublayer from the InChI string, by key
    """
    key_layer_parser = Parse.key_layer_(key)
    cap_dct = apf.first_named_capture(key_layer_parser.PATTERN, ich)
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
    cap_dct = apf.first_named_capture(Parse.AtomStereo.PATTERN, ich)
    ret = ()
    if cap_dct:
        lyr = cap_dct[Parse.AtomStereo.LAYER_KEY]
        ret = apf.all_captures(Parse.AtomStereo.Term.PATTERN, lyr)
    return ret


def bond_stereo_elements(ich):
    """ bond stereo keys and values
    """
    cap_dct = apf.first_named_capture(Parse.BondStereo.PATTERN, ich)
    ret = ()
    if cap_dct:
        lyr = cap_dct[Parse.BondStereo.LAYER_KEY]
        ret = apf.all_captures(Parse.BondStereo.Term.PATTERN, lyr)
    return ret


def known_atom_stereo_elements(ich):
    """ atom stereo keys and values, excluding unknown stereo
    """
    return tuple((key, val) for key, val in atom_stereo_elements(ich)
                 if val not in Parse.AtomStereo.Term.UNKNOWN_VALS)


def known_bond_stereo_elements(ich):
    """ bond stereo keys and values, excluding unknown stereo
    """
    return tuple((key, val) for key, val in bond_stereo_elements(ich)
                 if val not in Parse.BondStereo.Term.UNKNOWN_VALS)


def has_known_stereo_elements(ich):
    """ does this InChI string have any (known or unknown) stereo elements?
    """
    ich_ste = recalculate(ich, force_stereo=True)
    return (known_atom_stereo_elements(ich_ste) or
            known_bond_stereo_elements(ich_ste))


def has_unknown_stereo_elements(ich):
    """ does this InChI string have unknown stereo elements?
    """
    ich_ste = recalculate(ich, force_stereo=True)
    return (atom_stereo_elements(ich_ste) !=
            known_atom_stereo_elements(ich_ste) or
            bond_stereo_elements(ich_ste) !=
            known_bond_stereo_elements(ich_ste))


def substereomers(ich):
    """ expand InChI string to its compatible stereoisomers
    """
    atm_terms_lst = []
    for num, (key, val) in enumerate(
            atom_stereo_elements(recalculate(ich, force_stereo=True))):
        terms = ([key + val] if val not in Parse.AtomStereo.Term.UNKNOWN_VALS
                 else [key + Parse.AtomStereo.Term.MINUS_VAL] if num == 0
                 else [key + Parse.AtomStereo.Term.MINUS_VAL,
                       key + Parse.AtomStereo.Term.PLUS_VAL])
        atm_terms_lst.append(terms)

    bnd_terms_lst = []
    for key, val in bond_stereo_elements(recalculate(ich, force_stereo=True)):
        terms = ([key + val] if val not in Parse.BondStereo.Term.UNKNOWN_VALS
                 else [key + Parse.BondStereo.Term.MINUS_VAL,
                       key + Parse.BondStereo.Term.PLUS_VAL])
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


# comparisons
def same_connectivity(ich1, ich2):
    """ do these InChI strings have the same connectivity?
    """
    return (formula_layer(ich1) == formula_layer(ich2) and
            key_layer(ich1, 'c') == key_layer(ich2, 'c') and
            key_layer(ich1, 'h') == key_layer(ich2, 'h'))


# conversions
def recalculate(ich, force_stereo=False):
    """ recalculate an inchi string
    """
    return automol.convert.inchi.recalculate(ich, force_stereo=force_stereo)


def geometry(ich):
    """ inchi => geometry
    """
    return automol.convert.inchi.geometry(ich)


def graph(ich):
    """ inchi => graph
    """
    return automol.convert.inchi.graph(ich)


def smiles(ich):
    """ inchi => smiles
    """
    return automol.convert.inchi.smiles(ich)


def inchi_key(ich):
    """ inchi => inchi_key
    """
    return automol.convert.inchi.inchi_key(ich)
