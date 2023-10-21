""" tunneling
"""

import os

import numpy
from automol import reac
from automol.util import tensor
from phydat import phycon

PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, "data")


# Obtain cubic and quartic force constants
CUBIC_STR = None
with open(os.path.join(DAT_PATH, "ch4_h.cubic"), mode="r", encoding="utf-8") as fobj:
    CUBIC_STR = fobj.read()

QUARTIC_STR = None
with open(os.path.join(DAT_PATH, "ch4_h.quartic"), mode="r", encoding="utf-8") as fobj:
    QUARTIC_STR = fobj.read()

CUBIC_MAT = tensor.from_string(CUBIC_STR, fill_perms=True)
QUARTIC_MAT = tensor.from_string(QUARTIC_STR, fill_perms=True)

ENES = numpy.arange(-0.5, 17.5, 0.5) * phycon.KCAL2EH
FREQS = [
    -1354.5,
    569.5,
    571.9,
    1117.2,
    1209.1,
    1209.4,
    1465.7,
    1465.9,
    1788.0,
    3114.7,
    3262.8,
    3263.3,
]
RXN_FREQ = -1 * FREQS[0]


def test__():
    """test reac.tunnel.transmission_coefficients
    test reac.tunnel.actions
    """

    ref_alpha1 = 1.0
    ref_alpha2 = 1.0

    alpha1 = reac.tunnel.alpha(FREQS, CUBIC_MAT)
    alpha2 = reac.tunnel.alpha(FREQS, CUBIC_MAT, qfc_mat=QUARTIC_MAT)
    print(alpha1, ref_alpha1)
    print(alpha2, ref_alpha2)
    # assert numpy.isclose(alpha1, ref_alpha1)
    # assert numpy.isclose(alpha2, ref_alpha2)

    ref_kes = ()
    ref_ses = ()

    kes = reac.tunnel.transmission_coefficient(ENES, alpha2, RXN_FREQ)
    ses = reac.tunnel.action(ENES, alpha2, RXN_FREQ)
    print(kes, ref_kes)
    print(ses, ref_ses)
    # assert numpy.allclose(kes, ref_kes)
    # assert numpy.allclose(ses, ref_ses)


if __name__ == "__main__":
    test__()
