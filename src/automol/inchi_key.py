"""Functions operating on InChI and AMChI keys."""

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


def first_hash(chk: str) -> str:
    """Get the first hash block, indicating connectivity.

    :param chk: An InChI or AMChI key
    :return: The first hash block
    """
    return CHI_KEY.parseString(chk).get("hash1")


def second_hash(chk: str) -> str:
    """Parse ChIKey for the second hash block, indicating stereochemistry.

    :param chk: A ChI key
    :return: The second hash block
    """
    return CHI_KEY.parseString(chk).get("hash2")


def version_indicator(chk: str) -> str:
    """Parse ChIKey second-hash block for the ChIKey version indicator.

    :param chk: A ChI key
    :return: ChIKey version indicator
    """
    return CHI_KEY.parseString(chk).get("version")


def protonation_indicator(chk: str) -> str:
    """Parse final character of ChIKey for the protonation indicator.

    :param chk: A ChI key
    :return: Final ChIKey protonation indicator character
    """
    return CHI_KEY.parseString(chk).get("protonation")


def second_hash_with_extension(chk: str) -> str:
    """Parse ChIKey second-hash block for version and protonation
    indicators.

    :param chk: A ChI key
    :return: Version and protonation indicators from second-hash block
    """
    hash2 = CHI_KEY.parseString(chk).get("hash2")
    vers = CHI_KEY.parseString(chk).get("version")
    prot = CHI_KEY.parseString(chk).get("protonation")
    return hash2 + vers + "-" + prot
