""" test reac._pst
"""

import numpy
from automol import reac


KT_PST = 4.0e-10
N_PAR = 6
CN_PAR = 6
MRED = 1.0
TEMP = 300.0


def test__pst():
    """ test reac.pst_kt
        test reac.pst_cn
    """

    pst_kt = reac.pst_kt(N_PAR, MRED, CN_PAR, TEMP)

    pst_cn = reac.pst_cn(KT_PST, N_PAR, MRED, TEMP)
    pst_kt2 = reac.pst_kt(N_PAR, MRED, pst_cn, TEMP)

    assert numpy.isclose(pst_kt, pst_kt2)
