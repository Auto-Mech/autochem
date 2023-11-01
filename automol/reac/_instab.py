"""
 Build unstable products
"""

import itertools
from phydat import instab_fgrps
import automol.graph
import automol.geom
import automol.chi
import automol.zmat
from automol.graph import radical_dissociation_products
from automol.graph import radical_group_dct
from automol.reac._7conv import with_structures_from_zmatrix


# Identify instability products
def instability_product_zmas(zma, stereo=True):
    """ Determine if the species has look for functional group attachments that
        could cause molecule instabilities
    """

    ich = automol.geom.chi(automol.zmat.geometry(zma))
    instab_ichs = instability_product_inchis(ich, stereo=stereo)

    if instab_ichs is not None:
        instab_zmas = tuple(automol.chi.zmatrix(ich)
                            for ich in instab_ichs)
    else:
        instab_zmas = None

    return instab_zmas


def instability_product_inchis(ich, stereo=True):
    """ Generate ichs for instable inchi. Also validates if
        the decomposition is valid by ID'd it as a beta-scission
    """

    instab_ichs = None

    gra = automol.graph.explicit(automol.chi.graph(ich))
    instab_gras = instability_product_graphs(
        gra, stereo=False)

    if instab_gras:
        instab_ichs = [automol.graph.chi(gra) for gra in instab_gras]
        ste_prd1_ichs = automol.chi.expand_stereo(instab_ichs[0])
        ste_prd2_ichs = automol.chi.expand_stereo(instab_ichs[1])
        prd_ichs_lst = itertools.product(ste_prd1_ichs, ste_prd2_ichs)

        for prd_ichs in prd_ichs_lst:
            rxn_objs = automol.reac.with_structures_from_chi(
                (ich,), prd_ichs, stereo=stereo)
            if rxn_objs is not None:
                instab_ichs = prd_ichs
                break

    return instab_ichs


def instability_product_graphs(gra, stereo=True):
    """ Determine if the species has look for functional group attachments that
        could cause molecule instabilities
    """

    # Build graphs for the detection scheme
    rad_grp_dct = radical_group_dct(gra)

    # Check for instability causing functional groups
    prd_gras = ()
    for atm, grps in rad_grp_dct.items():
        if atm in instab_fgrps.DCT:
            fgrps_dct = instab_fgrps.DCT[atm]
            for grp in grps:
                grp_ich = automol.graph.chi(grp, stereo=stereo)
                if grp_ich in fgrps_dct:
                    # If instability found, determine prod of the instability
                    prd_ich = fgrps_dct[grp_ich]
                    prd_geo = automol.chi.geometry(prd_ich)
                    prd_gra = automol.geom.graph(prd_geo)
                    prd_gras = radical_dissociation_products(gra, prd_gra)
                    break

    return prd_gras


# Build transformation object for instability
def instability_transformation(conn_zma, disconn_zmas):
    """ Build the reaction objects for an instability
    """

    zrxn_objs = with_structures_from_zmatrix(
        [conn_zma], disconn_zmas, zmat=True, stereo=True)
    if zrxn_objs:
        zrxn, zma, _, _ = zrxn_objs[0]
    else:
        zrxn, zma = None, None

    return zrxn, zma
