""" Determine transition state spin multiplicities.
"""

import itertools
import numpy
from automol.mult._mult import spin as _spin


def high(rct_mults, prd_mults):
    """ Calculate high-spin multiplicity for a transition state
        from the multiplicities of the reactants and products.

        :param rct_mults: multiplicites of reactants
        :type rct_mults: tuple(float)
        :param prd_mults: multiplicites of products
        :type prd_mults: tuple(float)
        :rtype: int
    """
    return min(_high(rct_mults), _high(prd_mults))


def low(rct_mults, prd_mults):
    """ Calculate low-spin multiplicity for a transition state
        from the multiplicities of the reactants and products.

        :param rct_mults: multiplicites of reactants
        :type rct_mults: tuple(float)
        :param prd_mults: multiplicites of products
        :type prd_mults: tuple(float)
        :rtype: int
    """
    return max(_low(rct_mults), _low(prd_mults))


def _high(mults):
    """ Obtain the highest spin multiplicity state that can be obtained
        from a set of multiplciities.

        :param mults: spin multiplicities
        :type mults: tuple(int)
        :rtype: int
    """
    spns = list(map(_spin, mults))
    hi_spn = sum(spns)
    return int(hi_spn + 1)


def _low(mults):
    """ Obtain the lowest spin multiplicity state that can be obtained
        from a set of multiplciities

        :param mults: spin multiplicities
        :type mults: tuple(int)
        :rtype: int
    """
    spns = list(map(_spin, mults))
    num = len(spns)
    low_spn = min(map(abs, map(sum, (
        numpy.multiply(spns, sgns)
        for sgns in itertools.product([+1, -1], repeat=num)))))
    return int(low_spn + 1)
