"""
 Library to deal unstable species
"""

# init: just check geom
# conf: check ratio of confs
# hr:   check geom along scan,

import automol
import autofile
import elstruct
import mechanalyzer
from lib import filesys
from phydat import instab_fgrps
from automol.zmatrix._unimol_ts import beta_scission


# Check functional groups
def functional_groups_stable(geo, thy_save_fs, mod_thy_info):
    """ look for functional group attachments that could cause
        molecule instabilities
    """

    # Initialize empty set of product graphs
    prd_gras = ()

    # Check for instability causing functional groups
    gra = automol.geom.graph(geo)
    rad_grp_dct = automol.graph.radical_group_dct(gra)
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

    if prd_gras:
        disconn_zmas = []
        for gra in prd_gras:
            ich = automol.graph.inchi(gra)
            geo_tmp = automol.inchi.geometry(ich)
            zma = automol.geom.zmatrix(geo_tmp)
            disconn_zmas.append(zma)
        conn_zma = automol.geom.zmatrix(geo)
        structure.instab.write_instab2(
            conn_zma, disconn_zmas,
            thy_save_fs, mod_thy_info[1:4],
            zma_locs=(0,),
            save_cnf=True)

        stable = False
    else:
        stable = True
    return stable


def _instab_info(conn_zma, disconn_zmas):
    """ Obtain instability info
    """

    # Get the zma for the connected graph
    rct_zmas = [conn_zma]

    # Get the zmas used for the identification
    prd_zmas = disconn_zmas

    # Get the keys
    ret = beta_scission(rct_zmas, prd_zmas)
    zma, _, brk_bnd_keys, _, rcts_gra = ret

    return zma, brk_bnd_keys, rcts_gra


def _disconnected_zmas(disconn_zma):
    """ get graphs
    """

    # Convert to disconnected component graph
    disconn_geo = automol.zmatrix.geometry(disconn_zma)
    disconn_gras = automol.graph.connected_components(
        automol.geom.graph(disconn_geo))

    # Get the zmas
    disconn_zmas = [automol.geom.zmatrix(automol.graph.geometry(gra))
                    for gra in disconn_gras]

    return disconn_zmas
