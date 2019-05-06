""" functions operating on InChIKeys
"""
from autoparse.pattern import STRING_START as _STRING_START
from autoparse.pattern import STRING_END as _STRING_END
from autoparse.pattern import UPPERCASE_LETTER as _UPPERCASE_LETTER
from autoparse.pattern import escape as _escape
from autoparse.pattern import named_capturing as _named_capturing
from autoparse.find import has_match as _has_match
from autoparse.find import first_named_capture as _first_named_capture


class PARSE():
    """ InChIKey parser """
    _HASH1 = _UPPERCASE_LETTER * 14
    _HASH2 = _UPPERCASE_LETTER * 8
    _VERSION = _UPPERCASE_LETTER * 2
    _PROT = _UPPERCASE_LETTER

    HASH1_KEY = 'hash1'
    HASH2_KEY = 'hash2'
    VERSION_KEY = 'version'
    PROT_KEY = 'protonation'
    PATTERN = (
        _STRING_START +
        _named_capturing(_HASH1, name=HASH1_KEY) + _escape('-') +
        _named_capturing(_HASH2, name=HASH2_KEY) +
        _named_capturing(_VERSION, name=VERSION_KEY) + _escape('-') +
        _named_capturing(_PROT, name=PROT_KEY) + _STRING_END)


def is_valid(ick):
    """ is this a valid InChIKey?
    """
    assert isinstance(ick, (str, bytes, bytearray))
    return _has_match(PARSE.PATTERN, ick)


def first_hash(ick):
    """ the first hash block, indicating connectivity
    """
    assert is_valid(ick)
    cap_dct = _first_named_capture(PARSE.PATTERN, ick)
    hash1 = cap_dct[PARSE.HASH1_KEY]
    return hash1


def second_hash(ick):
    """ the second hash block, indicating stereochemistry etc.
    """
    assert is_valid(ick)
    cap_dct = _first_named_capture(PARSE.PATTERN, ick)
    hash2 = cap_dct[PARSE.HASH2_KEY]
    return hash2


def version_indicator(ick):
    """ the version indicator following the second hash block
    """
    assert is_valid(ick)
    cap_dct = _first_named_capture(PARSE.PATTERN, ick)
    ver = cap_dct[PARSE.VERSION_KEY]
    return ver


def protonation_indicator(ick):
    """ protonation indicator (last character of the key)
    """
    assert is_valid(ick)
    cap_dct = _first_named_capture(PARSE.PATTERN, ick)
    prot = cap_dct[PARSE.PROT_KEY]
    return prot
