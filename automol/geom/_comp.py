"""
  Functions used for handling and comparing geometries
"""

import numpy
import automol


def are_torsions_same(geo, geoi, ts_bnds=()):
    """ compare all torsional angle values
    """
    dtol = 0.09
    same_dihed = True
    zma = automol.geom.zmatrix(geo, ts_bnds=ts_bnds)
    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo, ts_bnds=ts_bnds)
    zmai = automol.geom.zmatrix(geoi)
    tors_namesi = automol.geom.zmatrix_torsion_coordinate_names(geoi, ts_bnds=ts_bnds)
    for idx, tors_name in enumerate(tors_names):
        val = automol.zmatrix.values(zma)[tors_name]
        vali = automol.zmatrix.values(zmai)[tors_namesi[idx]]
        valip = vali+2.*numpy.pi
        valim = vali-2.*numpy.pi
        vchk1 = abs(val - vali)
        vchk2 = abs(val - valip)
        vchk3 = abs(val - valim)
        if vchk1 > dtol and vchk2 > dtol and vchk3 > dtol:
            same_dihed = False
    return same_dihed


# Checks
def is_unique(lst, lsti, check_dct):
    """ Compare one of many structure features of a geometry to that of
        a list of geometries to see if it is unique.
    """

    unique = True
    for pair in lsti:
        for key, val in check_dct:
            if key != 'stereo' or val is not None:
                kwargs = {}
            else:
                kwargs = {'chk_arg': val}
            if not CHECK_DCT[key](lst, pair, **kwargs):
                unique = False

    return unique


def _ene(lst, lsti, chk_arg=2e-5):
    """
    """
    ene, enei = lst[0], lsti[0]
    return abs(ene-enei) < chk_arg


def _dist(lst, lsti, chk_arg=3e-1):
    """
    """
    geo, geoi = lst[1], lsti[1]
    return automol.geom.almost_equal_dist_matrix(
        geo, geoi, thresh=chk_arg)


def _tors(lst, lsti, chk_arg=()):
    """
    """
    geo, geoi = lst[1], lsti[1]
    return are_torsions_same(geo, geoi, ts_bnd=chk_arg)


def _stereo(lst, lsti):
    """
    """
    ich, ichi = automol.geom.inchi(lst[1]), automol.geom.inchi(lsti[1])
    return bool(ich == ichi)


def _coloumb(lst, lsti, check_arg=1e-2):
    """
    """
    geo, geoi = lst[1], lsti[1]
    return automol.geom.almost_equal_coulomb_spectrum(
        geo, geoi, rtol=check_arg)


CHECK_DCT = {
    'ene': _ene,
    'dist': _dist,
    'tors': _tors,
    'stereo': _stereo,
    'coloumb': _coloumb
}
