""" tunneling
"""

import os
import numpy
from phydat import phycon
import automol.reac
from automol.util import highd_mat
from ioformat import read_text_file

PATH = os.path.dirname(os.path.realpath(__file__))

# Obtain cubic and quartic force constants
CUBIC_STR = read_text_file(['data'], 'ch4_h.cubic', path=PATH)
CUBIC_MAT = highd_mat.from_string(
    CUBIC_STR, fill_perms=True)
QUARTIC_STR = read_text_file(['data'], 'ch4_h.quartic', path=PATH)
QUARTIC_MAT = highd_mat.from_string(
    QUARTIC_STR, fill_perms=True)

ENES = numpy.arange(-0.5, 17.5, 0.5) * phycon.KCAL2EH
FREQS = [-1354.5, 569.5, 571.9, 1117.2, 1209.1, 1209.4,
         1465.7, 1465.9, 1788.0, 3114.7, 3262.8, 3263.3]
RXN_FREQ = -1*FREQS[0]


def test__():
    """ test automol.reac.tunnel.transmission_coefficients
        test automol.reac.tunnel.actions
    """

    ref_alpha1 = 1.0
    ref_alpha2 = 1.0

    alpha1 = automol.reac.tunnel.alpha(FREQS, CUBIC_MAT)
    alpha2 = automol.reac.tunnel.alpha(FREQS, CUBIC_MAT, qfc_mat=QUARTIC_MAT)
    print(alpha1, ref_alpha1)
    print(alpha2, ref_alpha2)
    # assert numpy.isclose(alpha1, ref_alpha1)
    # assert numpy.isclose(alpha2, ref_alpha2)

    ref_kes = ()
    ref_ses = ()

    kes = automol.reac.tunnel.transmission_coefficient(ENES, alpha2, RXN_FREQ)
    ses = automol.reac.tunnel.action(ENES, alpha2, RXN_FREQ)
    print(kes, ref_kes)
    print(ses, ref_ses)
    # assert numpy.allclose(kes, ref_kes)
    # assert numpy.allclose(ses, ref_ses)


if __name__ == '__main__':
    test__()
