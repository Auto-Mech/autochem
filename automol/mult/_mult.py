""" Calculate spins (2S) from spin multiplicities.
"""


def spin(mult):
    """ Calculate the spin (2Ms) from using the spin multiplicity
        which is equivalent to the number of unpaired electrons.

        :param mult: multiplicity
        :type mult: int
        :rtype: int
    """
    return mult - 1
