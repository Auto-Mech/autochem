"""Determine transition state spin multiplicities."""

import itertools

import numpy
from _collections_abc import Sequence

from ._mult import spin as _spin


def high(rct_mults: Sequence[float], prd_mults: Sequence[float]) -> int:
    """Calculate high-spin multiplicity for a transition state
    from the multiplicities of the reactants and products.

    :param rct_mults: Multiplicites of reactants
    :param prd_mults: Multiplicites of products
    :return: High spin multiplicity
    """
    return min(_high(rct_mults), _high(prd_mults))


def low(rct_mults: Sequence[float], prd_mults: Sequence[float]) -> int:
    """Calculate low-spin multiplicity for a transition state
    from the multiplicities of the reactants and products.

    :param rct_mults: Multiplicites of reactants
    :param prd_mults: Multiplicites of products
    :return: Low spin multiplicity
    """
    return max(_low(rct_mults), _low(prd_mults))


def _high(mults: Sequence[float]) -> int:
    """Obtain the highest spin multiplicity state that can be obtained
    from a set of multiplciities.

    :param mults: Spin multiplicities
    :return: High spin multiplicity
    """
    spns = list(map(_spin, mults))
    hi_spn = sum(spns)
    return int(hi_spn + 1)


def _low(mults: Sequence[float]) -> int:
    """Obtain the lowest spin multiplicity state that can be obtained
    from a set of multiplciities.

    :param mults: spin multiplicities
    :return: Lowest spin multiplicity
    """
    spns = list(map(_spin, mults))
    num = len(spns)
    low_spn = min(
        map(
            abs,
            map(
                sum,
                (
                    numpy.multiply(spns, sgns)
                    for sgns in itertools.product([+1, -1], repeat=num)
                ),
            ),
        )
    )
    return int(low_spn + 1)
