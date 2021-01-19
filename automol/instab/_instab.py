"""
 Build unstable products
"""

import automol
from phydat import instab_fgrps


def product_zmas(zma):
    """ Determine if the species has look for functional group attachments that
        could cause molecule instabilities
    """

    disconn_zmas = ()
    for gra in product_graphs(automol.zmat.graph(zma)):
        ich = automol.graph.inchi(gra)
        geo_tmp = automol.inchi.geometry(ich)
        zma = automol.geom.zmatrix(geo_tmp)
        disconn_zmas += (zma,)

    return disconn_zmas


def product_graphs(gra):
    """ Determine if the species has look for functional group attachments that
        could cause molecule instabilities
    """

    # Build graphs for the detection scheme
    rad_grp_dct = automol.graph.radical_group_dct(gra)

    # Check for instability causing functional groups
    prd_gras = ()
    for atm, grps in rad_grp_dct.items():
        if atm in instab_fgrps.DCT:
            fgrps, prds = instab_fgrps.DCT[atm]
            for grp in grps:
                grp_ich = automol.graph.inchi(grp)
                if grp_ich in fgrps:
                    # If instability found, determine prod of the instability
                    prd_ich = prds[fgrps.index(grp_ich)]
                    prd_geo = automol.inchi.geometry(prd_ich)
                    prd_gra = automol.geom.graph(prd_geo)
                    prd_gras = automol.graph.radical_dissociation_prods(
                        gra, prd_gra)
                    break

    return prd_gras


# def _instab_info(conn_zma, disconn_zmas):
#     """ Determine keys corresponding to breaking bond for the instabiliity
#     """
#
#     # Get the zma for the connected graph
#     rct_zmas = [conn_zma]
#
#     # Get the zmas used for the identification
#     prd_zmas = disconn_zmas
#
#     # Get the keys
#     ret = beta_scission(rct_zmas, prd_zmas)
#     zma, _, brk_bnd_keys, _, rcts_gra = ret
#
#     return zma, brk_bnd_keys, rcts_gra
