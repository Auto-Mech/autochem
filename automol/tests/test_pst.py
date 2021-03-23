""" PST
"""

import automol.reac


def test__pst():
    """ test pst k(T)  equations
    """

    kt = automol.reac.pst_kt(n_par, mred, cn_par, temp)
    cn = automol.reac.pst_cn(kt_pst, n_pst, mred, temp)
    kt2 = automol.reac.pst_kt(n_par, mred, cn, temp)

    assert numpy.isclose(kt, kt2 )
