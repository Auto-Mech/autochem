""" Use combining rules to generate energy transfer parameters for A+A
    interactions when what you have are A+B and B+B. This is because A+A
    parameters are often used to represent energy transfer in mechanisms.
"""


def epsilon(ab_eps, bb_eps):
    """ Perform combining rule to get A+A epsilon parameter.
        Output units are whatever those are of the input parameters.

        :param ab_eps: A+B epsilon parameter
        :type ab_eps: float
        :param ab_eps: B+B epsilon parameter
        :type ab_eps: float
        :rtype: float
    """

    if ab_eps is not None and bb_eps is not None:
        aa_eps = ab_eps**2 / bb_eps
    else:
        aa_eps = None

    return aa_eps


def sigma(ab_sig, bb_sig):
    """ Perform combining rule to get A+A sigma parameter.
        Output units are whatever those are of the input parameters.

        :param ab_sig: A+B sigma parameter
        :type ab_sig: float
        :param ab_sig: B+B sigma parameter
        :type ab_sig: float
        :rtype: float
    """

    if ab_sig is not None and bb_sig is not None:
        aa_sig = 2.0 * ab_sig - bb_sig
    else:
        aa_sig = None

    return aa_sig
