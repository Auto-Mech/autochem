"""SMARTS functions."""

from . import rd


# properties
def shape(smarts: str) -> tuple[list[int], list[int]]:
    """Get numbers of reactants and products in SMARTS string.

    :param smarts: SMARTS string
    :return: Reaction shape
    """
    return rd.shape(rd.from_smarts(smarts))


def reactant_count(smarts: str) -> tuple[list[int], list[int]]:
    """Get number of reactants in SMARTS string.

    :param smarts: SMARTS string
    :return: Reaction shape
    """
    return rd.reactant_count(rd.from_smarts(smarts))


def product_count(smarts: str) -> tuple[list[int], list[int]]:
    """Get number of products in SMARTS string.

    :param smarts: SMARTS string
    :return: Reaction shape
    """
    return rd.product_count(rd.from_smarts(smarts))
