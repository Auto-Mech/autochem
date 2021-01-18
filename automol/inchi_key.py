""" InChIKeys

(could be simplified like I did for automol.inchi, but for now this works)
"""
import autoparse.pattern as app
import autoparse.find as apf


class Parse():
    """ InChIKey parser """
    _HASH1 = app.UPPERCASE_LETTER * 14
    _HASH2 = app.UPPERCASE_LETTER * 8
    _VERSION = app.UPPERCASE_LETTER * 2
    _PROT = app.UPPERCASE_LETTER

    HASH1_KEY = 'hash1'
    HASH2_KEY = 'hash2'
    VERSION_KEY = 'version'
    PROT_KEY = 'protonation'
    PATTERN = (
        app.STRING_START +
        app.named_capturing(_HASH1, name=HASH1_KEY) + app.escape('-') +
        app.named_capturing(_HASH2, name=HASH2_KEY) +
        app.named_capturing(_VERSION, name=VERSION_KEY) + app.escape('-') +
        app.named_capturing(_PROT, name=PROT_KEY) + app.STRING_END)


def is_valid(ick):
    """ Determine if an InChIKey has the proper form.

        :param ick: InChIKey
        :type ick: str
        :rtype: bool
    """
    assert isinstance(ick, (str, bytes, bytearray))
    return apf.has_match(Parse.PATTERN, ick)


def first_hash(ick):
    """ the first hash block, indicating connectivity

        :param ick: InChIKey
        :type ick: str
        :rtype: bool
    """
    assert is_valid(ick)
    cap_dct = apf.first_named_capture(Parse.PATTERN, ick)
    hash1 = cap_dct[Parse.HASH1_KEY]
    return hash1


def second_hash(ick):
    """ the second hash block, indicating stereochemistry etc.

        :param ick: InChIKey
        :type ick: str
        :rtype: bool
    """
    assert is_valid(ick)
    cap_dct = apf.first_named_capture(Parse.PATTERN, ick)
    hash2 = cap_dct[Parse.HASH2_KEY]
    return hash2


def version_indicator(ick):
    """ the version indicator following the second hash block

        :param ick: InChIKey
        :type ick: str
        :rtype: bool
    """
    assert is_valid(ick)
    cap_dct = apf.first_named_capture(Parse.PATTERN, ick)
    ver = cap_dct[Parse.VERSION_KEY]
    return ver


def protonation_indicator(ick):
    """ protonation indicator (last character of the key)

        :param ick: InChIKey
        :type ick: str
        :rtype: bool
    """
    assert is_valid(ick)
    cap_dct = apf.first_named_capture(Parse.PATTERN, ick)
    prot = cap_dct[Parse.PROT_KEY]
    return prot


def second_hash_with_extension(ick):
    """ second hash block with version and protonation indicators

        :param ick: InChIKey
        :type ick: str
        :rtype: bool
    """
    return (second_hash(ick) + version_indicator(ick) + '-' +
            protonation_indicator(ick))
