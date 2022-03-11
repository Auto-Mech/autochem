""" Checks that have to be done at the high-level of inchi directory
    since they require other automol imports
"""

import numbers
import automol.graph
import automol.inchi


def is_valid_inchi_multiplicity(ich, mul):
    """ is this multiplicity compatible with this inchi string?

    :param ich: inchi
    :type ich: str
    :param mul: multiplicity
    :type mul: int
    :returns: validity of inchi multiplicity
    :rtype: bool
    """
    assert isinstance(mul, numbers.Integral)
    return mul in automol.graph.possible_spin_multiplicities(
        automol.inchi.graph(ich, stereo=False))
