"""
 Build unstable products
"""

import itertools

from phydat import instab_fgrps

from .. import chi as chi_
from .. import geom, graph, zmat
from ._0core import ts_structure
from ._5conv import from_chis, from_zmatrices


# Identify instability products
def instability_product_zmas(zma, stereo=True):
    """Determine if the species has look for functional group attachments that
    could cause molecule instabilities
    """

    ich = geom.chi(zmat.geometry(zma))
    instab_ichs = instability_product_inchis(ich, stereo=stereo)

    if instab_ichs is not None:
        instab_zmas = tuple(chi_.zmatrix(ich) for ich in instab_ichs)
    else:
        instab_zmas = None

    return instab_zmas


def instability_product_inchis(ich, stereo=True):
    """Generate ichs for instable inchi. Also validates if
    the decomposition is valid by ID'd it as a beta-scission
    """

    instab_ichs = None

    gra = graph.explicit(chi_.graph(ich))
    instab_gras = instability_product_graphs(gra, stereo=False)

    if instab_gras:
        instab_ichs = [graph.chi(gra) for gra in instab_gras]
        ste_prd1_ichs = chi_.expand_stereo(instab_ichs[0])
        ste_prd2_ichs = chi_.expand_stereo(instab_ichs[1])
        prd_ichs_lst = itertools.product(ste_prd1_ichs, ste_prd2_ichs)

        for prd_ichs in prd_ichs_lst:
            rxn_objs = from_chis((ich,), prd_ichs, stereo=stereo)
            if rxn_objs is not None:
                instab_ichs = prd_ichs
                break

    return instab_ichs


def instability_product_graphs(gra, stereo=True):
    """Determine if the species has look for functional group attachments that
    could cause molecule instabilities
    """

    # Build graphs for the detection scheme
    rad_grp_dct = graph.radical_group_dct(gra)

    # Check for instability causing functional groups
    prd_gras = ()
    for atm, grps in rad_grp_dct.items():
        if atm in instab_fgrps.DCT:
            fgrps_dct = instab_fgrps.DCT[atm]
            for grp in grps:
                grp_ich = graph.chi(grp, stereo=stereo)
                if grp_ich in fgrps_dct:
                    # If instability found, determine prod of the instability
                    prd_ich = fgrps_dct[grp_ich]
                    prd_geo = chi_.geometry(prd_ich)
                    prd_gra = geom.graph(prd_geo)
                    prd_gras = graph.radical_dissociation_products(gra, prd_gra)
                    break

    return prd_gras


# Build transformation object for instability
def instability_transformation(conn_zma, disconn_zmas):
    """Build the reaction objects for an instability"""

    zrxns = from_zmatrices([conn_zma], disconn_zmas, struc_typ="zmat", stereo=True)
    if zrxns:
        zrxn, *_ = zrxns
        zma = ts_structure(zrxn)
    else:
        zrxn, zma = None, None

    return zrxn, zma
