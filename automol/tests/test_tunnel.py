""" tunneling
"""

import numpy
from phydat import phycon
import automol.reac
import automol.util.highd_mat
from _util import read_file


# Obtain cubic and quartic force constants
CUBIC_STR = read_file(['data'], 'ch4_h.cubic')
CUBIC_MAT = automol.util.highd_mat.from_string(
    CUBIC_STR, fill_perms=True)
QUARTIC_STR = read_file(['data'], 'ch4_h.quartic')
QUARTIC_MAT = automol.util.highd_mat.from_string(
    QUARTIC_STR, fill_perms=True)

ENES = numpy.arange(0.0, 17.0, 0.5) * phycon.KCAL2EH
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
    print(alpha1)
    print(alpha2)
    # assert numpy.isclose(alpha1, ref_alpha1)
    # assert numpy.isclose(alpha2, ref_alpha2)

    ref_kes = ()
    ref_ses = ()

    kes = automol.reac.tunnel.transmission_coefficient(ENES, alpha2, RXN_FREQ)
    ses = automol.reac.tunnel.action(ENES, alpha2, RXN_FREQ)
    # assert numpy.allclose(kes, ref_kes)
    # assert numpy.allclose(ses, ref_ses)
    print(kes)
    print(ses)


if __name__ == '__main__':
    test__()
