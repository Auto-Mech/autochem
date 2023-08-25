""" ChIKeys
"""

import pyparsing as pp

DASH = pp.Suppress("-")
HASH1 = pp.Combine(pp.Char(pp.alphas.upper()) * 14)
HASH2 = pp.Combine(pp.Char(pp.alphas.upper()) * 8)
VERSION = pp.Combine(pp.Char(pp.alphas.upper()) * 2)
PROTONATION = pp.Char(pp.alphas.upper())

CHI_KEY = (
    HASH1("hash1")
    + DASH
    + HASH2("hash2")
    + VERSION("version")
    + DASH
    + PROTONATION("protonation")
)


def to_dict(chk: str) -> dict:
    """Split the ChI key into a dictionary of parsed elements

    :param chk: The ChI key
    :type chk: str
    :return: The parsed elements; keys: "hash1", "hash2", "version", "protonation"
    :rtype: dict
    """
    return CHI_KEY.parseString(chk).asDict()


def first_hash(chk):
    """Parse ChIKey for the first hash block, indicating connectivity.

    :param chk: ChIKey
    :type chk: str
    :rtype: str
    """
    chk_dct = to_dict(chk)
    return chk_dct["hash1"]


def second_hash(chk):
    """Parse ChIKey for the second hash block, indicating stereochemistry.

    :param chk: ChIKey
    :type chk: str
    :rtype: str
    """
    chk_dct = to_dict(chk)
    return chk_dct["hash2"]


def version_indicator(chk):
    """Parse ChIKey second-hash block for the ChIKey version indicator.

    :param chk: ChIKey
    :type chk: str
    :rtype: str
    """
    chk_dct = to_dict(chk)
    return chk_dct["version"]


def protonation_indicator(chk):
    """Parse final character of ChIKey for the protonation indicator.

    :param chk: ChIKey
    :type chk: str
    :rtype: str
    """
    chk_dct = to_dict(chk)
    return chk_dct["protonation"]


def second_hash_with_extension(chk):
    """Parse ChIKey second-hash block for version and protonation
    indicators.

    :param chk: ChIKey
    :type chk: str
    :rtype: str
    """
    return second_hash(chk) + version_indicator(chk) + "-" + protonation_indicator(chk)
