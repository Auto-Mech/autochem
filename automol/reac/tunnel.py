"""
 New Tunneling Equtions
"""

import math
import numpy
from phydat import phycon


# Tunneling expressions
def transmission_coefficient(enes, valpha, rxn_freq):
    """ Calculate the tunneling transmission coefficients at Es.
    """

    kes = ()
    for ene in enes:
        kes += (_transmission_coefficient(ene, valpha, rxn_freq),)

    return kes


def action(enes, valpha, rxn_freq):
    """ Calculate action at seve
    """

    ses = ()
    for ene in enes:
        ses += (_action(ene, valpha, rxn_freq),)

    return ses


def _transmission_coefficient(ene, valpha, rxn_freq):
    """ Calculate the tunneling transmission coefficient at some energy.

        :param ene: energy to calculate the transmission coefficient
        :param alpha: alpha coefficient in S(E) expansion
        :param rxn_freq: frequency of the reaction mode
    """

    if ene < 0.0:
        denom = math.exp(action(-1.0*ene, valpha, rxn_freq))
    else:
        denom = (-2.0 * numpy.pi * ene) / rxn_freq

    p_e = 1.0 / (1.0 + denom)

    return p_e


def _action(ene, valpha, rxn_freq):
    """ Tunneling action as a second-order expansion of the energy
    """

    rxn_freq *= phycon.WAVEN2EH
    s_e = (
        (2.0 * numpy.pi * ene) / rxn_freq +
        (valpha * ene**2) / rxn_freq**2
    )

    return s_e


# Quantities descrining the surface
def alpha(freqs, cfc_mat, qfc_mat=None):
    """ calculate alpha from expansion
    """

    # Obtain the imaginary frequency and sort other freqs
    ridxs = tuple(idx for idx, freq in enumerate(freqs) if freq < 0.0)
    assert len(ridxs) == 1, ('Freqs should only have one imag')
    ridx = ridxs[0]

    freqs = tuple(abs(freq) * phycon.WAVEN2EH for freq in freqs)
    rfreq = freqs[ridx]

    # Calculate alpha with cubic contributions
    term1 = (5.0 / 3.0) * (
        cfc_mat[ridx, ridx, ridx]**2 / rfreq**3)

    term2 = 0.0
    for i, ifreq in enumerate(freqs):
        vval = cfc_mat[i, ridx, ridx]

        tm1 = (vval**2 / ifreq**3)
        tm2 = 2.0*ifreq**(-2) + (1.0 / (4.0*rfreq**2 + ifreq**2))
        term2 += (tm1 * tm2)

    _alpha = (numpy.pi / 8.0) * (term1 - term2)

    # Calculate quartic contributions to alpha, if available
    if qfc_mat is not None:
        _alpha -= (numpy.pi / 8.0) * (
            qfc_mat[ridx, ridx, ridx, ridx] / rfreq**6)

    return _alpha
