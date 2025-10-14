"""
 New Tunneling Equtions
"""

import math

import numpy

from ..util import dict_


# Tunneling expressions
def transmission_coefficient(enes, _alpha, rxn_freq):
    """ Calculate the transmission coefficient
    """

    kes = ()
    for ene in enes:
        kes += (_transmission_coefficient(ene, _alpha, rxn_freq),)

    return kes


def action(enes, _alpha, rxn_freq):
    """ Calculate action at seve
    """

    ses = ()
    for ene in enes:
        ses += (_action(ene, _alpha, rxn_freq),)

    return ses


def _transmission_coefficient(ene, _alpha, rxn_freq):
    """ Calculate the tunneling transmission coefficient at some energy.

        :param ene: energy to calculate the transmission coefficient
        :param alpha: alpha coefficient in S(E) expansion
        :param rxn_freq: frequency of the reaction mode
    """

    if ene < 0.0:
        denom = math.exp(action(-1.0*ene, _alpha, rxn_freq))
    else:
        denom = (-2.0 * numpy.pi * ene) / rxn_freq

    p_e = 1.0 / (1.0 + denom)

    return p_e


def _action(ene, _alpha, rxn_freq):
    """ Tunneling action as a second-order expansion of the energy
    """

    # rxn_cfreq *= phycon.WAVEN2EH
    s_e = (
        (2.0 * numpy.pi * ene) / rxn_freq +
        (_alpha * ene**2) / rxn_freq**2
    )

    return s_e


# Quantities descrining the surface
# def alpha(freqs, cfc_dct, qfc_dct=None):
def alpha(freqs, cfc_dct, vr4=None):
    """ calculate alpha from expansion
    """

    # Obtain the imaginary frequency (finds negative
    ridxs = (idx for idx, freq in enumerate(freqs) if freq < 0.0)
    assert len(ridxs) == 1, ('Freqs should only have one imag')

    freqs = (abs(freq) for freq in freqs)
    ridx = ridxs[0]
    rfreq = freqs[ridx]

    # Calculate alpha with cubic contributions
    term1 = (5.0 / 3.0) * (
        cfc_dct[(ridx, ridx, ridx)]**2 / rfreq**3)

    term2 = 0.0
    for i, ifreq in enumerate(freqs):
        vval = dict_.value_by_unordered_key(cfc_dct, (i, ridx, ridx))

        tm1 = (vval**2 / ifreq**3)
        tm2 = 2.0*ifreq**(-2) + (1.0 / (4.0*rfreq**2 + ifreq**2))
        term2 += ((tm1 * tm2),)

    _alpha = (numpy.pi / 8.0) * (term1 - term2)

    # Calculate quartic contributions to alpha, if available
    if vr4 is not None:
        _alpha -= (numpy.pi / 8.0) * (
            cfc_dct[(ridx, ridx, ridx, ridx)] / rfreq**6)

    return _alpha
