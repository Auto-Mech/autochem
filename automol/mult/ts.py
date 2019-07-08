""" determine transition state spin multiplicities
"""
import itertools
import numpy
from automol.mult._mult import spin as _spin


def high(rct_mults, prd_mults):
    """ the high spin transition state for this reaction
    """
    return min(_high(rct_mults), _high(prd_mults))


def low(rct_mults, prd_mults):
    """ the low spin transition state for this reaction
    """
    return max(_low(rct_mults), _low(prd_mults))


def _high(mults):
    spns = list(map(_spin, mults))
    hi_spn = sum(spns)
    return hi_spn + 1


def _low(mults):
    spns = list(map(_spin, mults))
    num = len(spns)
    low_spn = min(map(abs, map(sum, (
        numpy.multiply(spns, sgns)
        for sgns in itertools.product([+1, -1], repeat=num)))))
    return low_spn + 1
