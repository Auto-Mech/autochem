"""
 Build unstable products
"""

from phydat import instab_fgrps
import automol.graph
from automol.reac._util import rxn_objs_from_zmatrix
import automol.geom
import automol.inchi
import automol.zmat
from automol.graph import radical_dissociation_products
from automol.graph import radical_group_dct


def instability_product_zmas(zma):
    """ Determine if the species has look for functional group attachments that
        could cause molecule instabilities
    """

    disconn_zmas = ()
    for gra in instability_product_graphs(automol.zmat.graph(zma)):
        ich = automol.graph.inchi(gra)
        geo_tmp = automol.inchi.geometry(ich)
        zma = automol.geom.zmatrix(geo_tmp)
        disconn_zmas += (zma,)

    return disconn_zmas


def instability_product_graphs(gra):
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
                grp_ich = automol.graph.inchi(grp, stereo=True)
                if grp_ich in fgrps_dct:
                    # If instability found, determine prod of the instability
                    prd_ich = fgrps_dct[grp_ich]
                    prd_geo = automol.inchi.geometry(prd_ich)
                    prd_gra = automol.geom.graph(prd_geo)
                    prd_gras = radical_dissociation_products(gra, prd_gra)
                    break

    return prd_gras


def instability_transformation(conn_zma, disconn_zmas):
    """ Build the reaction objects for an instability
    """

    zrxn_objs = rxn_objs_from_zmatrix(
        [conn_zma], disconn_zmas, indexing='zma', stereo=True)
    if zrxn_objs:
        zrxn, zma, _, _ = zrxn_objs[0]
    else:
        zrxn, zma = None, None

    return zrxn, zma
