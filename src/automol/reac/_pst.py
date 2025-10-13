"""
  Calculate the pre-exponential parameter required to calculate some desired
  Phase-Space Theory rate constant, given an exponential parameter and
  temperature. Using rate expression constant gives functional
  dependence  of k(n, mu, C0, T).
"""

import scipy
import numpy
from phydat import phycon


def pst_kt(n_par, mred, cn_par, temp):
    """ Calculate a rate constant according to Phase-Space Theory.

        :param n_par: exponential parameter
        :param mred: reduced mass ()
        :param cn_par: pre-exponential potential coefficient [in Bohr]
        :param temp: temperature (K)
        :return: k(T)
        :rtype: float

        temp_pst [K] -> [hartree]
        n_par [unitless]
        mred [amu] -> [au]
        cn_par [???] [hartree/bohr^n]
    """

    # k2eh = 0.000003166808534191
    # amu2au = 1.0 / (9.1093837015e-28 * 6.022e23)
    # BOHR2CM = 5.29177e-9

    mred *= phycon.AMU2AU
    temp *= phycon.K2EH

    kt_val = (
        (8.0 * numpy.pi)**(1.0/2.0) *
        ((n_par - 2) / 2)**(2.0/n_par) *
        scipy.special.gamma(1.0 - 2.0/n_par) *
        mred**(-1.0/2.0) *
        cn_par**(2.0/n_par) *
        temp**(1.0/2.0 - 2.0/n_par)
    )
    kt_val *= (phycon.BOHR2CM**3 / phycon.JIFFY)

    return kt_val


def pst_cn(kt_pst, n_pst, mred, t_pst):
    """ Calculate a Cn value to match a k(T) value for PST
  For (n=N, mu=MU, T=T), the target parameter (C0_TGT) needed to obtain
  the target rate constant (k_TGT) can be found via

    C0_TGT**(2/n) = [ k_TGT / k(n=N, mu=MU, C0=1.0, T=300.0) ] * 1.0
    """

    cn_par = numpy.sqrt(kt_pst / pst_kt(n_pst, mred, 1.0, t_pst))**n_pst

    return cn_par
