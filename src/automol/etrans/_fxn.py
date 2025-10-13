"""
  Various functions useful for energy transfer calculations
"""

import numpy
from phydat import phycon


# LENNARD JONES
def troe_lj_collision_frequency(eps, sig, red_mass, temp):
    """ Collision Frequency formula from Troe that uses
        Lennard-Jones epsilons and sigma parameters

        Z = sqrt( (8*kB*T)/(pi*mu) ) * sig^2 * omega
        omega = [ 0.7 + 0.5 log10( (kB*T)/eps ) ]^-1

        :param eps: Target+Bath Lennard-Jones epsilon value (in hart)
        :type eps: float
        :param sig: Target+Bath Lennard-Jones sigma value (in bohr)
        :type sig: float
        :param red_mass: Reduced mass of colliding bodies (in amu)
        :type red_mass: float
        :param temp: Temperature at which the collision occurs
        :type temp: float
        :return zlj: Lennard-Jones collision frequency (in m^3/s)
        :rtype: float
    """

    # Convert LJ params to SI units
    eps *= phycon.EH2J
    sig *= phycon.BOHR2M
    red_mass *= phycon.AMU2KG

    # Create a k_B*T constant since it is in both terms
    kbt = phycon.KB * temp

    # Prefactor term with physical constants and masses
    pref1 = numpy.sqrt((8.0 * kbt) / (numpy.pi * red_mass))

    # Omega integral
    omega_s = (0.7 + 0.52 * (numpy.log10(kbt / eps)))**(-1)

    return pref1 * sig**2 * omega_s
