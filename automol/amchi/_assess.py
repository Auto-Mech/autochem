""" Checks that have to be done at the high-level of amchi directory
    since they require other automol imports
"""

import numbers
import automol.graph
from automol.amchi._conv import graph


def is_valid_multiplicity(chi, mul):
    """ is this multiplicity compatible with this amchi string?

    :param chi: ChI string
    :type chi: str
    :param mul: multiplicity
    :type mul: int
    :returns: validity of amchi multiplicity
    :rtype: bool
    """
    assert isinstance(mul, numbers.Integral)
    return mul in automol.graph.possible_spin_multiplicities(
        graph(chi, stereo=False))
