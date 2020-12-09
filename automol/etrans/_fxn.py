"""
  Various functions useful for energy transfer calculations
"""

import numpy


# LENNARD JONES
def troe_lj_collision_frequency(eps, sig, red_mass, temp):
    """ Collision Frequency formula from Troe that uses
        Lennard-Jones epsilons and sigma parameters

        :param eps: Target+Bath Lennard-Jones epsilon value (in _)
        :type eps: float
        :param sig: Target+Bath Lennard-Jones sigma value (in _)
        :type sig: float
        :param red_mass: Reduced mass of colliding bodies (in kg)
        :type red_mass: float
        :param temp: Temperature at which the collision occurs
        :type temp: float
        :return zlj: Lennard-Jones collision frequency (in _)
        :rtype: float
    """

    pref1 = 1.0e-14 * numpy.sqrt(
        (8.0 * 1.380603e-23 * temp) / (numpy.pi * red_mass)
    )
    pref2 = 0.7 + 0.52 * (numpy.log(0.69502 * temp / eps) / numpy.log(10))

    zlj = numpy.pi * sig**2 * (pref1 / pref2)

    return zlj


def combine_epsilon(eps1, eps2):
    """ Combine the epsilon values of two species.

        :param eps1: epsilon of species 1
        :type eps1: float
        :param eps2: epsilon of species 2
        :type eps2: float
        :return: eps_comb
        :rtpye: float
    """

    if eps1 is not None and eps2 is not None:
        eps_comb = eps1**2 / eps2
    else:
        eps_comb = None

    return eps_comb


def combine_sigma(sig1, sig2):
    """ Combine the sigma values of two species.

        :param sig1: sigma of species 1
        :type sig1: float
        :param sig2: sigma of species 2
        :type sig2: float
        :return: sig_comb
        :rtpye: float
    """

    if sig1 is not None and sig2 is not None:
        sig_comb = (2.0 * sig1) - sig2
    else:
        sig_comb = None

    return sig_comb
