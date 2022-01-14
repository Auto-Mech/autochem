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
